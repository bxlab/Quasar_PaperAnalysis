#!/usr/bin/env python

"""A class for determining HiC data quality."""

import os
import sys
from math import ceil, floor
import argparse

import numpy
import h5py
import hifive
from scipy.optimize import curve_fit
try:
    from mpi4py import MPI
except:
    pass
try:
    from pyx import *
    unit.set(defaultunit="cm")
    text.set(mode="latex")
except:
    pass



class QuasarPseudo(hifive.Quasar):

    """This class performs subsampling and QuASAR transformations for calculating HiC quality.

    .. note::
      This class is also available as hifive.Quasar

    When initialized, this class creates an h5dict in which to store all data associated with this object.
    
    :param filename: The file name of the h5dict to store the QuASAR-transformed data in.
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`Quasar` class object.

    :Attributes: * **file** (*str.*) A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
    """
    def __init__(self, filename, mode='a', silent=False):
        """Create a :class:`Quasar` object."""
        try:
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.num_procs = self.comm.Get_size()
        except:
            self.comm = None
            self.rank = 0
            self.num_procs = 1
        self.file = os.path.abspath(filename)
        self.hic_fname = None
        self.hic_fname2 = None
        self.silent = silent
        self.filetype = 'quasar'
        if mode != "w":
            self.load(mode)
        else:
            if self.rank == 0:
                self.storage = h5py.File(self.file, 'w')
            else:
                self.storage = None
        return None

    def find_transformation(self, hic, hic2, balance, chroms=[], resolutions=[1000000, 200000, 40000],
        coverages=[0, 40000000, 20000000, 10000000, 5000000, 2000000, 1000000], seed=None):
        if self.rank == 0:
            coverages = numpy.array(coverages, dtype=numpy.int64)
            resolutions = numpy.array(resolutions, dtype=numpy.int64)
            resolutions.sort()
            hic_fname = hic.file
            hic_fname2 = hic2.file
            if ((self.hic_fname is not None and self.hic_fname != hic_fname) or
                (self.hic_fname2 is not None and self.hic_fname2 != hic_fname2)):
                for key in self.storage.keys():
                    if key.split('.')[0] in ['valid', 'dist', 'corr']:
                        del self.storage[key]
                    if 'chromosomes' in self.storage:
                        del self.storage['chromosomes']
            self.hic_fname = hic_fname
            self.hic_fname2 = hic_fname2
            self.balance = float(balance)
            if seed is not None:
                RNG = numpy.random.RandomState(seed=seed)
            else:
                RNG = numpy.random.RandomState()

            # load partition information
            if 'binned' in hic.__dict__ and hic.binned is not None:
                temp_mids1 = hic.fends['bins']['mid'][...]
                chr_indices1 = hic.fends['bin_indices'][...]
            else:
                temp_mids1 = hic.fends['fends']['mid'][...]
                chr_indices1 = hic.fends['chr_indices'][...]
            if 'binned' in hic2.__dict__ and hic.binned is not None:
                temp_mids2 = hic2.fends['bins']['mid'][...]
                chr_indices2 = hic2.fends['bin_indices'][...]
            else:
                temp_mids2 = hic2.fends['fends']['mid'][...]
                chr_indices2 = hic2.fends['chr_indices'][...]

            # fill in chromosome list if empty. Otherwise check that all specified chromosomes exist.
            if not isinstance(chroms, list) or len(chroms) == 0:
                chroms = numpy.intersect1d(hic.fends['chromosomes'][...], hic2.fends['chromosomes'][...])
                valid = numpy.ones(chroms.shape[0], dtype=numpy.bool)
                for i in range(chroms.shape[0]):
                    chrint1 = hic.chr2int[chroms[i]]
                    chrint2 = hic2.chr2int[chroms[i]]
                    if (chr_indices1[chrint1 + 1] - chr_indices1[chrint1] == 0 or
                        chr_indices2[chrint2 + 1] - chr_indices2[chrint2] == 0):
                        valid[i] = False
                    elif (hic.data['cis_indices'][chr_indices1[chrint1 + 1]] -
                          hic.data['cis_indices'][chr_indices1[chrint1]] == 0 or
                          hic2.data['cis_indices'][chr_indices2[chrint2 + 1]] -
                          hic2.data['cis_indices'][chr_indices2[chrint2]] == 0):
                        valid[i] = False
                chroms = chroms[valid]

            # Load raw counts
            bounds1 = numpy.zeros((len(chroms), 2), numpy.int64)
            bounds2 = numpy.zeros((len(chroms), 2), numpy.int64)
            for i, chrom in enumerate(chroms):
                chrint1 = hic.chr2int[chrom]
                chrint2 = hic2.chr2int[chrom]
                bounds1[i, 0] = hic.data['cis_indices'][chr_indices1[chrint1]]
                bounds1[i, 1] = hic.data['cis_indices'][chr_indices1[chrint1 + 1]]
                bounds2[i, 0] = hic2.data['cis_indices'][chr_indices2[chrint2]]
                bounds2[i, 1] = hic2.data['cis_indices'][chr_indices2[chrint2 + 1]]
            raw1 = numpy.zeros((numpy.sum(bounds1[:, 1] - bounds1[:, 0]), 3), dtype=numpy.int64)
            indices1 = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            raw2 = numpy.zeros((numpy.sum(bounds2[:, 1] - bounds2[:, 0]), 3), dtype=numpy.int64)
            indices2 = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            mids1 = {}
            mids2 = {}
            starts = numpy.zeros(len(chroms), dtype=numpy.int32)
            for i, chrom in enumerate(chroms):
                chrint1 = hic.chr2int[chrom]
                chrint2 = hic2.chr2int[chrom]
                indices1[i + 1] = indices1[i] + bounds1[i, 1] - bounds1[i, 0]
                indices2[i + 1] = indices2[i] + bounds2[i, 1] - bounds2[i, 0]
                temp = hic.data['cis_data'][bounds1[i, 0]:bounds1[i, 1], :]
                temp[:, :2] -= chr_indices1[chrint1]
                raw1[indices1[i]:indices1[i + 1], :] = temp
                mids1[chrom] = temp_mids1[chr_indices1[chrint1]:chr_indices1[chrint1 + 1]]
                temp = hic2.data['cis_data'][bounds2[i, 0]:bounds2[i, 1], :]
                temp[:, :2] -= chr_indices2[chrint2]
                raw2[indices2[i]:indices2[i + 1], :] = temp
                mids2[chrom] = temp_mids2[chr_indices2[chrint2]:chr_indices2[chrint2 + 1]]
                starts[i] = min(mids1[chrom][0], mids2[chrom][0])

            # only consider coverage levels that are less than or equal to the number of cis reads
            coverages = coverages[numpy.where(min(numpy.sum(raw1[:, 2]), numpy.sum(raw2[:, 2])) >= coverages)]
            coverages.sort()
            if coverages[0] == 0:
                coverages[:-1] = coverages[1:]
                coverages[-1] = 0
            store_coverages = numpy.copy(coverages)
            total_reads1 = numpy.sum(raw1[:, 2])
            total_reads2 = numpy.sum(raw2[:, 2])
            total_reads = (total_reads1 + total_reads2) / 2
            if coverages.shape[0] > 0 and coverages[-1] == 0:
                coverages[-1] = total_reads
            coverages = coverages[::-1]
        else:
            coverages = None
            resolutions = None
        if self.comm is not None:
            coverages = self.comm.bcast(coverages, root=0)
            resolutions = self.comm.bcast(resolutions, root=0)
            chroms = self.comm.bcast(chroms, root=0)
        if coverages.shape[0] == 0:
            return None

        if self.rank == 0:
            # write arguements to h5dict
            if 'chromosomes' in self.storage:
                del self.storage['chromosomes']
            self.storage.create_dataset(name='chromosomes', data=numpy.array(chroms))
            if 'resolutions' in self.storage:
                del self.storage['resolutions']
            self.storage.create_dataset(name='resolutions', data=numpy.array(resolutions))
            if 'coverages' in self.storage:
                del self.storage['coverages']
            self.storage.create_dataset(name='coverages', data=numpy.array(store_coverages))
            if 'starts' in self.storage:
                del self.storage['starts']
            self.storage.create_dataset(name='starts', data=starts)
            self.storage.attrs['total_reads'] = total_reads

            # rebin data to highest resolution for faster processing
            remapped = {}
            new_mids = {}
            new_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                start = (starts[i] / resolutions[0]) * resolutions[0]
                stop = ((max(mids1[chrom][-1], mids1[chrom][-1]) - 1) / resolutions[0] + 1) * resolutions[0]
                N = (stop - start) / resolutions[0]
                mapping = (mids1[chrom] - start) / resolutions[0]
                raw1[indices1[i]:indices1[i + 1], 0] = mapping[raw1[indices1[i]:indices1[i + 1], 0]]
                raw1[indices1[i]:indices1[i + 1], 1] = mapping[raw1[indices1[i]:indices1[i + 1], 1]]
                new_index = numpy.unique(raw1[indices1[i]:indices1[i + 1], 0] * N +
                                         raw1[indices1[i]:indices1[i + 1], 1])
                index = numpy.searchsorted(new_index, raw1[indices1[i]:indices1[i + 1], 0] * N +
                                                      raw1[indices1[i]:indices1[i + 1], 1])
                remapped[chrom] = numpy.zeros((new_index.shape[0], 3), dtype=numpy.int64)
                remapped[chrom][:, 0] = new_index / N
                remapped[chrom][:, 1] = new_index % N
                remapped[chrom][:, 2] = numpy.bincount(index, weights=raw1[indices1[i]:indices1[i + 1], 2])
                new_indices[i + 1] = new_index.shape[0] + new_indices[i]
                new_mids[chrom] = (start + resolutions[0] / 2 + numpy.arange(N) *
                                   resolutions[0]).astype(numpy.int32)
            indices1 = new_indices.astype(numpy.int64)
            mids = new_mids
            raw1 = numpy.zeros((indices1[-1], 3), dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                raw1[indices1[i]:indices1[i + 1], :] = remapped[chrom]
            remapped = {}
            new_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                start = (starts[i] / resolutions[0]) * resolutions[0]
                stop = ((mids[chrom][-1] - 1) / resolutions[0] + 1) * resolutions[0]
                N = (stop - start) / resolutions[0]
                mapping = (mids2[chrom] - start) / resolutions[0]
                raw2[indices2[i]:indices2[i + 1], 0] = mapping[raw2[indices2[i]:indices2[i + 1], 0]]
                raw2[indices2[i]:indices2[i + 1], 1] = mapping[raw2[indices2[i]:indices2[i + 1], 1]]
                new_index = numpy.unique(raw2[indices2[i]:indices2[i + 1], 0] * N +
                                         raw2[indices2[i]:indices2[i + 1], 1])
                index = numpy.searchsorted(new_index, raw2[indices2[i]:indices2[i + 1], 0] * N +
                                                      raw2[indices2[i]:indices2[i + 1], 1])
                remapped[chrom] = numpy.zeros((new_index.shape[0], 3), dtype=numpy.int64)
                remapped[chrom][:, 0] = new_index / N
                remapped[chrom][:, 1] = new_index % N
                remapped[chrom][:, 2] = numpy.bincount(index, weights=raw2[indices2[i]:indices2[i + 1], 2])
                new_indices[i + 1] = new_index.shape[0] + new_indices[i]
            indices2 = new_indices.astype(numpy.int64)
            raw2 = numpy.zeros((indices2[-1], 3), dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                raw2[indices2[i]:indices2[i + 1], :] = remapped[chrom]
            del remapped

        # Transfer mids data
        chrom_ranges = numpy.round(numpy.linspace(0, len(chroms), self.num_procs + 1)).astype(numpy.int32)
        if self.rank == 0:
            for i in range(1, self.num_procs):
                for j in range(chrom_ranges[i], chrom_ranges[i + 1]):
                    chrom = chroms[j]
                    self.comm.send(mids[chrom].shape[0], dest=i, tag=3)
                    self.comm.Send(mids[chrom], dest=i, tag=4)
        else:
            mids = {}
            for i in range(chrom_ranges[self.rank], chrom_ranges[self.rank + 1]):
                chrom = chroms[i]
                N = self.comm.recv(source=0, tag=3)
                mids[chrom] = numpy.zeros(N, dtype=numpy.int32)
                self.comm.Recv(mids[chrom], source=0, tag=4)

        # cycle through coverages
        for c, cov in enumerate(coverages):
            raw_counts = {}
            if self.rank == 0:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rDownsampling to %i coverage") % (' ' * 120, cov),
                if store_coverages[c] == 0:
                    target = int(round(total_reads1 * balance))
                else:
                    target = int(round(cov * balance))
                raw1, indices1 = self._downsample(raw1, indices1, target, RNG)
                raw2, indices2 = self._downsample(raw2, indices2, cov - target, RNG)
                raw, indices = self._combine_counts(raw1, indices1, raw2, indices2)
                for i in range(chrom_ranges[1]):
                    chrom = chroms[i]
                    raw_counts[chrom] = raw[indices[i]:indices[i + 1], :]
                for i in range(1, self.num_procs):
                    for j in range(chrom_ranges[i], chrom_ranges[i + 1]):
                        chrom = chroms[j]
                        self.comm.send(indices[j + 1] - indices[j], dest=i, tag=1)
                        self.comm.Send(raw[indices[j]:indices[j + 1], :], dest=i, tag=2)
            else:
                for i in range(chrom_ranges[self.rank], chrom_ranges[self.rank + 1]):
                    chrom = chroms[i]
                    N = self.comm.recv(source=0, tag=1)
                    raw_counts[chrom] = numpy.zeros((N, 3), dtype=numpy.int64)
                    self.comm.Recv(raw_counts[chrom], source=0, tag=2)

            # cycle through resolutions
            for r, res in enumerate(resolutions):
                dists = {}
                norms = {}
                valids = {}
                corrs = {}

                # For each chromosome, normalize and find distance-corrected matrix
                if self.rank == 0:
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rCoverage %i Resolution %i - Normalizing counts") % (
                            ' ' * 120, cov, res),
                for chrom in raw_counts:
                    norm, dist, valid_rows = self._normalize(chrom, raw_counts[chrom], mids[chrom], res)
                    dists[chrom] = dist
                    norms[chrom] = norm
                    valids[chrom] = valid_rows
                if self.rank == 0:
                    for i in range(1, self.num_procs):
                        self._transfer_dict(dists, 0, i)
                        self._transfer_dict(norms, 0, i)
                        self._transfer_dict(valids, 0, i)
                else:
                    self._transfer_dict(dists, 0, self.rank)
                    self._transfer_dict(norms, 0, self.rank)
                    self._transfer_dict(valids, 0, self.rank)
                    dists = None

                # cycle through chromosomes finding correlation matrices
                for i, chrom in enumerate(chroms):
                    task = True
                    if self.rank == 0:
                        if cov == total_reads:
                            key = '%s.0C.%iR' % (chrom, res)
                        else:
                            key = '%s.%iC.%iR' % (chrom, cov, res)
                        if 'corr.%s' % key in self.storage:
                            task = False
                        if self.comm is not None:
                            task = self.comm.bcast(task, root=0)
                        if task:
                            if not self.silent:
                                print >> sys.stderr, ("\r%s\rCoverage %i Resolution %i - Correlating chrom %s") % (
                                    ' ' * 120, cov, res, chrom),
                            corrs[chrom] = self._find_correlations(norms[chrom], valids[chrom])
                    else:
                        task = self.comm.bcast(task, root=0)
                        if task:
                            self._find_correlations()

                # write resulting matrices to hdf5 file
                if self.rank == 0:
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rCoverage %i Resolution %i - Writing results") % (
                            ' ' * 120, cov, res),
                    for chrom in chroms:
                        if cov == total_reads:
                            key = '%s.0C.%iR' % (chrom, res)
                        else:
                            key = '%s.%iC.%iR' % (chrom, cov, res)
                        if valids[chrom] is None:
                            self.storage.attrs['%s.invalid' % (key)] = True
                        elif chrom in corrs:
                            self.storage.create_dataset(name="valid.%s" % (key),
                                data=valids[chrom])
                            self.storage.create_dataset(name="dist.%s" % (key), data=dists[chrom])
                            self.storage.create_dataset(name="corr.%s" % (key), data=corrs[chrom])
        if self.rank == 0 and not self.silent:
            print >> sys.stderr, ("\r%s\r") % (' ' * 120),
        return None

    def _combine_counts(self, counts1, indices1, counts2, indices2):
        if indices2 is None:
            return counts1, indices1
        elif indices1 is None:
            return counts2, indices2
        new_counts = {}
        indices = numpy.zeros(indices1.shape[0], dtype=indices1.dtype)
        for i in range(indices1.shape[0] - 1):
            if indices1[i + 1] == indices1[i]:
                new_counts[i] = counts2[indices2[i]:indices2[i + 1], :]
            elif indices2[i + 1] == indices2[i]:
                new_counts[i] = counts1[indices1[i]:indices1[i + 1], :]
            else:
                temp1 = counts1[indices1[i]:indices1[i + 1], :]
                temp2 = counts2[indices2[i]:indices2[i + 1], :]
                N = max(numpy.amax(temp1[:, :2]), numpy.amax(temp2[:, :2])) + 1
                index1 = temp1[:, 0] * N + temp1[:, 1]
                index2 = temp2[:, 0] * N + temp2[:, 1]
                index = numpy.unique(numpy.r_[index1, index2])
                new_counts[i] = numpy.zeros((index.shape[0], 3), dtype=numpy.int64)
                new_counts[i][:, 0] = index / N
                new_counts[i][:, 1] = index % N
                new_counts[i][:, 2] = numpy.bincount(numpy.searchsorted(index, index1), weights=temp1[:, 2],
                                                     minlength=index.shape[0])
                new_counts[i][:, 2] += numpy.bincount(numpy.searchsorted(index, index2), weights=temp2[:, 2],
                                                      minlength=index.shape[0]).astype(numpy.int64)
            indices[i + 1] = indices[i] + new_counts[i].shape[0]
        counts = numpy.zeros((indices[-1], 3), dtype=numpy.int64)
        for i in range(indices.shape[0] - 1):
            counts[indices[i]:indices[i + 1], :] = new_counts[i]
        return counts, indices


def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Find a quality score for a HiC dataset"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="hic", type=str, action='store', help="HiFive project")
    parser.add_argument(dest="hic2", type=str, action='store', help="Second HiFive project")
    parser.add_argument(dest="outfile", type=str, action='store', help="Output file")
    parser.add_argument('-r', dest="resolution", type=str, action='store', default='1000000',
        help="Comma-separated list of resolutions to find quality for")
    parser.add_argument('-d', dest="coverage", type=str, action='store', default='0',
        help="Comma-separated list of number of reads to use")
    parser.add_argument('-b', dest="balance", type=float, action='store', default=0.5,
        help="Percentage of reads to draw from first HiC project")
    parser.add_argument('-w', dest="width", type=int, action='store', default=100, help="Number of bins to use")
    parser.add_argument('-s', dest="seed", type=int, action='store', default=1, help="Random seed")
    parser.add_argument('--chroms', dest="chroms", type=str, action='store', default='', help="A Comma-separated list of chromosomes to use. Defaults to Numbered chromosomes up to 22 (fewer if appropriate) and X.")
    return parser

if __name__ == "__main__":
    parser = generate_parser()
    args = parser.parse_args()
    quasar = QuasarPseudo(args.outfile, mode='w')
    hic = hifive.HiC(args.hic)
    hic2 = hifive.HiC(args.hic2)
    resolutions = args.resolution.split(',')
    for i in range(len(resolutions)):
        resolutions[i] = int(resolutions[i])
    coverages = args.coverage.split(',')
    for i in range(len(coverages)):
        coverages[i] = int(coverages[i])
    if args.chroms == '':
        chroms = []
    else:
        chroms = args.chroms.split(',')
    quasar.find_transformation(hic, hic2, args.balance, chroms, resolutions, coverages)
    quasar.save()
    quasar.close()
