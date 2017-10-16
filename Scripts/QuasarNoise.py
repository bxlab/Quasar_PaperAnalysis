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



class QuasarNoise(hifive.Quasar):

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

    def find_transformation(self, hic, noise_fname, chroms=[], resolutions=[1000000, 200000, 40000],
        noises=[0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75], coverage=0, seed=None):
        if self.rank == 0:
            noises = numpy.array(noises, dtype=numpy.float64)
            noises.sort()
            resolutions = numpy.array(resolutions, dtype=numpy.int64)
            resolutions.sort()
            hic_fname = hic.file
            if self.hic_fname is not None and self.hic_fname != hic_fname:
                for key in self.storage.keys():
                    if key.split('.')[0] in ['valid', 'dist', 'corr']:
                        del self.storage[key]
                    if 'chromosomes' in self.storage:
                        del self.storage['chromosomes']
            self.hic_fname = hic_fname
            if seed is not None:
                RNG = numpy.random.RandomState(seed=seed)
            else:
                RNG = numpy.random.RandomState()

            # load partition information
            if 'binned' in hic.__dict__ and hic.binned is not None:
                temp_mids = hic.fends['bins']['mid'][...]
                chr_indices = hic.fends['bin_indices'][...]
            else:
                temp_mids = hic.fends['fends']['mid'][...]
                chr_indices = hic.fends['chr_indices'][...]

            # fill in chromosome list if empty. Otherwise check that all specified chromosomes exist.
            infile = h5py.File(noise_fname, 'r')
            if not isinstance(chroms, list) or len(chroms) == 0:
                chroms = numpy.intersect1d(hic.fends['chromosomes'][...], infile['chromosomes'][...])
                valid = numpy.ones(chroms.shape[0], dtype=numpy.bool)
                for i in range(chroms.shape[0]):
                    chrint = hic.chr2int[chroms[i]]
                    if chr_indices[chrint + 1] - chr_indices[chrint] == 0:
                        valid[i] = False
                    elif (hic.data['cis_indices'][chr_indices[chrint + 1]] -
                         hic.data['cis_indices'][chr_indices[chrint]] == 0):
                        valid[i] = False
                chroms = chroms[valid]

            # Load noise distributions
            print >> sys.stderr, ("\r%s\rLoading noise model") % (' ' * 120),
            noise_dists = {}
            noise_chrom_dist = numpy.zeros(len(chroms), dtype=numpy.float64)
            starts = numpy.zeros(len(chroms), dtype=numpy.int32)
            stops = numpy.zeros(len(chroms), dtype=numpy.int32)
            Ns = numpy.zeros(len(chroms), dtype=numpy.int32)
            for i, chrom in enumerate(chroms):
                M = int((infile['%s.noise' % chrom].shape[0] * 2 + 0.25) ** 0.5 - 0.5)
                indices2 = numpy.triu_indices(M, 0)
                start = infile.attrs['%s.start' % chrom]
                stop = infile.attrs['%s.stop' % chrom]
                starts[i] = (infile.attrs['%s.start' % chrom] / resolutions[0]) * resolutions[0]
                stops[i] = ((infile.attrs['%s.start' % chrom] + M * infile.attrs['binsize'] - 1) /
                           resolutions[0] + 1) * resolutions[0]
                Ns[i] = (stops[i] - starts[i]) / resolutions[0]
                N = Ns[i]
                indices3 = numpy.triu_indices(N, 0)
                mapping = ((numpy.arange(M) * infile.attrs['binsize'] + infile.attrs['%s.start' % chrom] - starts[i]) /
                           resolutions[0]).astype(numpy.int64)
                index = mapping[indices2[0]] * N + mapping[indices2[1]]
                model = infile['%s.noise' % chrom][...].astype(numpy.float64)
                temp = numpy.bincount(index, weights=model, minlength=(N * N)).reshape(N, N)
                model = temp[indices3]
                noise_chrom_dist[i] = numpy.sum(model)
                noise_dists[chrom] = model
                del temp
                del indices2
                del indices3
            noise_chrom_dist = numpy.cumsum(noise_chrom_dist)
            noise_chrom_dist /= noise_chrom_dist[-1]
            for i, chrom in enumerate(chroms):
                noise_dists[chrom] = numpy.cumsum(noise_dists[chrom])
                noise_dists[chrom] /= noise_dists[chrom][-1]
            infile.close()

            # Load raw counts
            bounds = numpy.zeros((len(chroms), 2), numpy.int64)
            for i, chrom in enumerate(chroms):
                chrint = hic.chr2int[chrom]
                bounds[i, 0] = hic.data['cis_indices'][chr_indices[chrint]]
                bounds[i, 1] = hic.data['cis_indices'][chr_indices[chrint + 1]]
            raw = numpy.zeros((numpy.sum(bounds[:, 1] - bounds[:, 0]), 3), dtype=numpy.int64)
            indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            mids = {}
            for i, chrom in enumerate(chroms):
                chrint = hic.chr2int[chrom]
                indices[i + 1] = indices[i] + bounds[i, 1] - bounds[i, 0]
                temp = hic.data['cis_data'][bounds[i, 0]:bounds[i, 1], :]
                temp[:, :2] -= chr_indices[chrint]
                raw[indices[i]:indices[i + 1], :] = temp
                mids[chrom] = temp_mids[chr_indices[chrint]:chr_indices[chrint + 1]]

            total_reads = numpy.sum(raw[:, 2])
            if coverage == 0:
                coverage = total_reads
            current_count = total_reads
            noise_targets = numpy.round(coverage * (1.0 - noises)).astype(numpy.int64)
            noise_counts = numpy.zeros((0, 3), dtype=numpy.int64)
            noise_indices = None
        else:
            noises = None
            resolutions = None
        if self.comm is not None:
            noises = self.comm.bcast(noises, root=0)
            resolutions = self.comm.bcast(resolutions, root=0)
            chroms = self.comm.bcast(chroms, root=0)

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
            self.storage.create_dataset(name='coverages', data=numpy.round(noises * 1000).astype(numpy.int64))
            self.storage.create_dataset(name='noises', data=noises)
            self.storage.attrs['coverage'] = coverage
            if 'starts' in self.storage:
                del self.storage['starts']
            self.storage.create_dataset(name='starts', data=starts)
            self.storage.attrs['total_reads'] = total_reads

            # rebin data to highest resolution for faster processing
            remapped = {}
            new_mids = {}
            new_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                start = starts[i]
                stop = stops[i]
                N = Ns[i]
                mapping = (mids[chrom] - start) / resolutions[0]
                raw[indices[i]:indices[i + 1], 0] = mapping[raw[indices[i]:indices[i + 1], 0]]
                raw[indices[i]:indices[i + 1], 1] = mapping[raw[indices[i]:indices[i + 1], 1]]
                new_index = numpy.unique(raw[indices[i]:indices[i + 1], 0] * N + raw[indices[i]:indices[i + 1], 1])
                index = numpy.searchsorted(new_index, raw[indices[i]:indices[i + 1], 0] * N +
                                                      raw[indices[i]:indices[i + 1], 1])
                remapped[chrom] = numpy.zeros((new_index.shape[0], 3), dtype=numpy.int64)
                remapped[chrom][:, 0] = new_index / N
                remapped[chrom][:, 1] = new_index % N
                remapped[chrom][:, 2] = numpy.bincount(index, weights=raw[indices[i]:indices[i + 1], 2])
                new_indices[i + 1] = new_index.shape[0] + new_indices[i]
                new_mids[chrom] = (start + resolutions[0] / 2 + numpy.arange(N) *
                                   resolutions[0]).astype(numpy.int32)
            indices = new_indices.astype(numpy.int64)
            mids = new_mids
            raw = numpy.zeros((indices[-1], 3), dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                raw[indices[i]:indices[i + 1], :] = remapped[chrom]
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

        # Cycle through desired noise levels
        for h, noise in enumerate(noises):
            raw_counts = {}
            if self.rank == 0:
                print >> sys.stderr, ("\r%s\rDownsampling for noise level %0.1f") % (' ' * 120, noise * 100.),
                if current_count > noise_targets[h]:
                    raw, indices = self._downsample(raw, indices, noise_targets[h], RNG)
                    current_count = noise_targets[h]
                    if noise > 0:
                        noise_counts, noise_indices = self._generate_noise(coverage - current_count, noise_dists,
                            noise_chrom_dist, noise_counts, noise_indices, h, Ns, chroms, RNG)
                        combined_counts, combined_indices = self._combine_counts(raw, indices, noise_counts,
                                                                                 noise_indices)
                    else:
                        combined_counts, combined_indices = raw, indices
                else:
                    combined_counts, combined_indices = raw, indices
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
                        print >> sys.stderr, ("\r%s\rNoise %0.1f Resolution %i - Normalizing counts") % (
                            ' ' * 120, noise * 100.0, res),
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
                        key = '%s.%iC.%iR' % (chrom, int(noise * 1000), res)
                        if 'corr.%s' % key in self.storage:
                            task = False
                        if self.comm is not None:
                            task = self.comm.bcast(task, root=0)
                        if task:
                            if not self.silent:
                                print >> sys.stderr, ("\r%s\rNoise %0.1f Resolution %i - Correlating chrom %s") % (
                                    ' ' * 120, noise * 100.0, res, chrom),
                            corrs[chrom] = self._find_correlations(norms[chrom], valids[chrom])
                    else:
                        task = self.comm.bcast(task, root=0)
                        if task:
                            self._find_correlations()

                # write resulting matrices to hdf5 file
                if self.rank == 0:
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rNoise %0.1f Resolution %i - Writing results") % (
                            ' ' * 120, noise * 100., res),
                    for chrom in chroms:
                        key = '%s.%iC.%iR' % (chrom, int(noise * 1000), res)
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

    def _generate_noise(self, target, dists, chrom_dist, counts, indices, h, Ns, chroms, RNG):
        if indices is not None:
            target -= numpy.sum(counts[:, 2])
        new_counts = self._sample_from_dist(dists, chrom_dist, chroms, target, RNG)
        new_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
        for i, chrom in enumerate(chroms):
            new_indices[i + 1] = new_indices[i] + new_counts[chrom].shape[0]
        all_counts = numpy.zeros((new_indices[-1], 3), dtype=numpy.int64)
        for i, chrom in enumerate(chroms):
            all_counts[new_indices[i]:new_indices[i + 1], :] = new_counts[chrom]
        del new_counts
        if indices is not None:
            all_counts, new_indices = self._combine_counts(counts, indices, all_counts, new_indices)
        return all_counts, new_indices

    def _sample_from_dist(self, dist, chrom_dist, chroms, target, RNG):
        chrom_targets = numpy.bincount(numpy.searchsorted(chrom_dist, RNG.rand(target)), minlength=len(chroms))
        counts = {}
        for i, chrom in enumerate(chroms):
            if chrom_targets[i] > 0:
                N = int((dist[chrom].shape[0] * 2 + 0.25) ** 0.5 - 0.5)
                index = numpy.searchsorted(dist[chrom], RNG.rand(chrom_targets[i]))
                index2 = numpy.unique(index)
                indices = numpy.triu_indices(N, 0)
                counts[chrom] = numpy.zeros((index2.shape[0], 3), dtype=numpy.int64)
                counts[chrom][:, 0] = indices[0][index2]
                counts[chrom][:, 1] = indices[1][index2]
                counts[chrom][:, 2] = numpy.bincount(numpy.searchsorted(index2, index))
            else:
                counts[chrom] = numpy.zeros((0, 3), dtype=numpy.int64)
        return counts

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
    parser.add_argument(dest="model", type=str, action='store', help="Background model file name")
    parser.add_argument(dest="outfile", type=str, action='store', help="Output file")
    parser.add_argument('-r', dest="resolution", type=str, action='store', default='1000000',
        help="Comma-separated list of resolutions to find quality for")
    parser.add_argument('-d', dest="coverage", type=int, action='store', default=0, help="Number of reads to use")
    parser.add_argument('-n', dest="noise", type=str, action='store', default=0,
        help="Comma-separated list of noise levels")
    parser.add_argument('-w', dest="width", type=int, action='store', default=100, help="Number of bins to use")
    parser.add_argument('-s', dest="seed", type=int, action='store', default=1, help="Random seed")
    parser.add_argument('--chroms', dest="chroms", type=str, action='store', default='', help="A Comma-separated list of chromosomes to use. Defaults to Numbered chromosomes up to 22 (fewer if appropriate) and X.")
    return parser

if __name__ == "__main__":
    parser = generate_parser()
    args = parser.parse_args()
    quasar = QuasarNoise(args.outfile, mode='w')
    hic = hifive.HiC(args.hic)
    resolutions = args.resolution.split(',')
    for i in range(len(resolutions)):
        resolutions[i] = int(resolutions[i])
    noises = args.noise.split(',')
    for i in range(len(noises)):
        noises[i] = float(noises[i])
    if args.chroms == '':
        chroms = []
    else:
        chroms = args.chroms.split(',')
    quasar.find_transformation(hic, args.model, chroms, resolutions, noises, args.coverage)
    quasar.save()
    quasar.close()
