#!/usr/bin/env python

import sys
import os
import glob

import numpy
import hifive
import h5py
from library import find_binning_expected
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_procs = comm.Get_size()


def main():
    pattern, hdf5_fname, binsize = sys.argv[1:4]
    fnames = glob.glob(pattern)
    if len(fnames) == 0:
        return None
    binsize = int(binsize)
    if rank == 0:
        outfile = h5py.File(hdf5_fname, 'w')
        for i in range(len(fnames)):
            hic = hifive.HiC(fnames[i])
            if i == 0:
                fends = hic.fends['fends'][...]
                chr_indices = hic.fends['chr_indices'][...]
                chromosomes = hic.fends['chromosomes'][...]
                chrom_sizes = hic.fends['chrom_sizes'][...]
                mappings = []
                counts = {}
                for j in range(chr_indices.shape[0] - 1):
                    if chr_indices[j + 1] - chr_indices[j] == 0:
                        mappings.append(None)
                        counts[chromosomes[j]] = None
                        continue
                    start = (fends['mid'][chr_indices[j]] / binsize) * binsize
                    mappings.append((fends['mid'][chr_indices[j]:chr_indices[j + 1]] - start) / binsize)
                    N = mappings[-1][-1] + 1
                    counts[chromosomes[j]] = numpy.zeros(N, dtype=numpy.int64)
                outfile.create_dataset(name='binning_fend_indices', data=hic.binning_fend_indices)
                outfile.create_dataset(name='binning_num_bins', data=hic.binning_num_bins)
                corrections = numpy.zeros((hic.binning_corrections.shape[0], len(fnames)),
                                          dtype=hic.binning_corrections.dtype)
            corrections[:, i] = hic.binning_corrections
            reads = hic.data['cis_data']
            for j in range(chr_indices.shape[0] - 1):
                if mappings[j] is None:
                    continue
                start = hic.data['cis_indices'][chr_indices[j]]
                stop = hic.data['cis_indices'][chr_indices[j + 1]]
                chrom = chromosomes[j]
                counts[chrom] += numpy.bincount(mappings[j][reads[start:stop, 0] - chr_indices[j]],
                                                minlength=counts[chrom].shape[0])
        outfile.create_dataset(name='binning_corrections', data=numpy.median(corrections, axis=1))
        chr2int = {}
        for i, chrom in enumerate(chromosomes):
            chr2int[chrom] = i
        chroms = []
        for i in range(1, 24):
            if str(i) in chromosomes:
                chroms.append(str(i))
        for chrom in ['X', '2L', '2R', '3L', '3R']:
            if chrom in chromosomes:
                chroms.append(chrom)
        lengths = numpy.zeros(len(chroms), dtype=numpy.int32)
        for i, chrom in enumerate(chroms):
            chrint = chr2int[chrom]
            lengths[i] = chrom_sizes[chrint]
            start = (fends['mid'][chr_indices[chrint]] / binsize) * binsize
            stop = ((fends['mid'][chr_indices[chrint + 1] - 1] - 1) / binsize + 1) * binsize
            outfile.attrs['%s.start' % chrom] = start
            outfile.attrs['%s.stop' % chrom] = stop
        outfile.create_dataset(name='chromosomes', data=numpy.array(chroms))
        outfile.create_dataset(name='chrom_sizes', data=lengths)
        outfile.attrs['binsize'] = binsize
        binning_corrections = outfile['binning_corrections'][...]
        binning_num_bins = outfile['binning_num_bins'][...]
        fend_indices = outfile['binning_fend_indices'][...]
        S1, S2, S3 = comm.bcast((binning_corrections.shape, binning_num_bins.shape, fend_indices.shape), root=0)
        chroms = comm.bcast(chroms, root=0)
        chr2int = comm.bcast(chr2int, root=0)
        fends = comm.bcast(fends, root=0)
        chr_indices = comm.bcast(chr_indices, root=0)
    else:
        outfile = None
        S1, S2, S3 = comm.bcast(None, root=0)
        chroms = comm.bcast(None, root=0)
        chr2int = comm.bcast(None, root=0)
        fends = comm.bcast(None, root=0)
        chr_indices = comm.bcast(None, root=0)
        binning_corrections = numpy.zeros(S1, dtype=numpy.float32)
        binning_num_bins = numpy.zeros(S2, dtype=numpy.int32)
        fend_indices = numpy.zeros(S3, dtype=numpy.int32)
        counts = {}
        for chrom in chroms:
            counts[chrom] = None
    if comm is not None:
        comm.Bcast(binning_corrections, root=0)
        comm.Bcast(binning_num_bins, root=0)
        comm.Bcast(fend_indices, root=0)
    for chrom in chroms:
        find_bin_probabilities(chrom, outfile, fends, chr_indices, binsize, chr2int, binning_corrections,
            binning_num_bins, fend_indices, counts[chrom])
    if rank == 0:
        outfile.close()
        print >> sys.stderr, ("\r%s\r") % (" " * 80),

def find_bin_probabilities(chrom, infile, fends, chr_indices, binsize, chr2int, corrections, num_bins, fend_indices,
                           counts):
    if rank == 0:
        if '%s.noise' % (chrom) in infile:
            del infile['%s.noise' % chrom]
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
        print >> sys.stderr, ("\r%s\rFinding chrom %s noise array") % (' '*120, chrom),
    else:
        if comm.recv(source=0, tag=11) == 0:
            return None
    chrint = chr2int[chrom]
    start_fend = chr_indices[chrint]
    stop_fend = chr_indices[chrint + 1]
    if stop_fend - start_fend == 0:
        if rank == 0:
            infile.attrs[chrom] = 'None'
        return None
    mids = fends['mid'][start_fend:stop_fend]
    start = (mids[0] / binsize) * binsize
    stop = ((mids[-1] - 1) / binsize + 1) * binsize
    if stop - start < 1000000:
        if rank == 0:
            infile.attrs[chrom] = 'None'
        return None
    N = (stop - start) / binsize
    mapping = ((mids - start) / binsize).astype(numpy.int32)
    distance_parameters = None
    chrom_mean = 1.0
    expected = numpy.zeros(N * (N + 1) / 2, dtype=numpy.float32)
    if rank == 0:
        indices = list(numpy.triu_indices(N, 0))
        indices[0] = indices[0].astype(numpy.int64)
        indices[1] = indices[1].astype(numpy.int64)
        node_ranges = numpy.round(numpy.linspace(0, indices[0].shape[0], num_procs + 1)).astype(numpy.int32)
        for i in range(1, num_procs):
            comm.send([node_ranges[i], node_ranges[i + 1], indices[0][node_ranges[i]],
                       indices[0][node_ranges[i + 1] - 1], indices[1][node_ranges[i]],
                       indices[1][node_ranges[i + 1] - 1]], dest=i, tag=11)
        start_index = node_ranges[0]
        stop_index = node_ranges[1]
        start1 = indices[0][start_index]
        stop1 = indices[0][stop_index - 1]
        start2 = indices[1][start_index]
        stop2 = indices[1][stop_index - 1]
    else:
        start_index, stop_index, start1, stop1, start2, stop2 = comm.recv(source=0, tag=11)
    startfend1 = numpy.searchsorted(mids, start1 * binsize + start)
    stopfend1 = numpy.searchsorted(mids, (stop1 + 1) * binsize + start)
    startfend2 = numpy.searchsorted(mids, start2 * binsize + start)
    stopfend2 = numpy.searchsorted(mids, (stop2 + 1) * binsize + start)
    find_binning_expected(
        mapping,
        corrections,
        num_bins,
        fend_indices,
        mids,
        distance_parameters,
        expected,
        chrom_mean,
        start_fend,
        startfend1,
        stopfend1,
        startfend2,
        stopfend2,
        start1,
        stop1,
        start2,
        stop2,
        1)
    if rank == 0:
        for i in range(1, num_procs):
            expected[node_ranges[i]:node_ranges[i + 1]] = comm.recv(source=i, tag=11)
        indices = numpy.triu_indices(N, 0)
        expected[numpy.where((counts[indices[0]] == 0) | (counts[indices[1]] == 0))[0]] = 0
        expected = expected.astype(numpy.float64)
        infile.create_dataset(name='%s.noise' % (chrom), data=expected)
        return None
    else:
        comm.send(expected[start_index:stop_index], dest=0, tag=11)
        return None

if __name__ == "__main__":
    main()
