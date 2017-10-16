#!/usr/bin/env python

import sys
import glob
import argparse

import h5py
import numpy
from mpi4py import MPI

from find_quasar_replicate_pair_score import correlate_samples, load_data, write_results

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_procs = comm.Get_size()

def main():
    num_procs = comm.Get_size()
    if rank == 0:
        parser = generate_parser()
        args = parser.parse_args()
        args.fnames = glob.glob(args.pattern)
        groups = 1
        while groups * (groups + 1) / 2 <= num_procs:
            groups += 1
        groups -= 1
        num_procs = groups * (groups + 1) / 2
        args.node_ranges = numpy.round(numpy.linspace(0, len(args.fnames), groups + 1)).astype(numpy.int32)
        infile = h5py.File(args.fnames[0], 'r')
        chroms = infile['chromosomes'][...]
        infile.close()
    else:
        args = None
    args = comm.bcast(args, root=0)
    if len(args.fnames) < 2:
        return None
    num_procs = comm.bcast(num_procs, root=0)
    if rank >= num_procs:
        return None
    data = {}
    indices = numpy.triu_indices(args.node_ranges.shape[0] - 1, 0)
    set1 = indices[0][rank]
    set2 = indices[1][rank]
    results = {}
    names1 = []
    names2 = []
    for i in range(args.node_ranges[set1], args.node_ranges[set1 + 1]):
        names1.append('_'.join(args.fnames[i].split('/')[-1].split('.')[:-1]))
    for i in range(args.node_ranges[set2], args.node_ranges[set2 + 1]):
        names2.append('_'.join(args.fnames[i].split('/')[-1].split('.')[:-1]))
    for i in range(args.node_ranges[set1], args.node_ranges[set1 + 1]):
        name1 = '_'.join(args.fnames[i].split('/')[-1].split('.')[:-1])
        data[name1] = load_all_data(args.fnames[i])
        if set1 != set2:
            start = args.node_ranges[set2]
        else:
            start = i + 1
        remaining = names1[(i + 1):] + names2[(start - args.node_ranges[set2]):]
        for j in range(start, args.node_ranges[set2 + 1]):
            if (args.fnames[i].count('Pseudo') + args.fnames[j].count('Pseudo') > 0 and
              args.fnames[i].split('Pseudo')[0].split('Rep')[0] != args.fnames[j].split('Pseudo')[0].split('Rep')[0]):
                remaining = names1[(i + 1):] + names2[(j - args.node_ranges[set2] + 1):]
                continue
            name2 = '_'.join(args.fnames[j].split('/')[-1].split('.')[:-1])
            data[name2] = load_all_data(args.fnames[j])
            keys1 = set(data[name1].keys())
            keys2 = set(data[name2].keys())
            keys = keys1.intersection(keys2)
            for key in keys:
                results[(name1, name2, key[0], key[1])] = correlate_samples(data[name1][key], data[name2][key], key[1])
            remaining = names1[(i + 1):] + names2[(j - args.node_ranges[set2] + 1):]
            if name2 not in remaining:
                del data[name2]
        if name1 not in remaining:
            del data[name1]
    if rank == 0:
        for i in range(1, num_procs):
            results.update(comm.recv(source=i, tag=i))
        write_results(results, chroms, args)
    else:
        comm.send(results, dest=0, tag=rank)

def load_all_data(fname):
    infile = h5py.File(fname, 'r')
    starts = infile['starts'][...]
    chroms = infile['chromosomes'][...]
    keys = []
    for key in infile.keys():
        temp = key.split('.')
        if len(temp) == 4:
            if temp[2] == '0C':
                keys.append(tuple([int(temp[2][:-1]), int(temp[3][:-1])]))
    keys = set(keys)
    data = {}
    for key in keys:
        data[key] = load_data(infile, chroms, key[0], key[1])
    infile.close()
    return data

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Find a quality score for a HiC dataset"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="pattern", type=str, action='store', help="Second quasar file pattern")
    parser.add_argument(dest="output", type=str, action='store', help="Results file name")
    return parser

if  __name__ == "__main__":
    main()
