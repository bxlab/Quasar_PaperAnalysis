#!/usr/bin/env python

import sys
import glob
import argparse

import h5py
import numpy

def main():
    parser = generate_parser()
    args = parser.parse_args()
    infile = h5py.File(args.input, 'r')
    if 'chromosomes' not in infile:
        infile.close()
        return None
    chroms = infile['chromosomes'][...]
    keys = []
    for key in infile.keys():
        temp = key.split('.')
        if len(temp) == 4:
            keys.append(tuple([int(temp[2][:-1]), int(temp[3][:-1])]))
    keys = set(keys)
    results = {}
    for key in keys:
        data = load_data(infile, chroms, key[0], key[1])
        if data is None:
            continue
        results[key] = find_results(data)
    infile.close()
    write_results(results, chroms, args)

def load_data(infile, chroms, cov, res):
    data = {}
    for i, chrom in enumerate(chroms):
        if chrom not in chroms:
            continue
        if '%s.%iC.%iR.invalid' % (chrom, cov, res) in infile.attrs:
            continue
        dist = infile['dist.%s.%iC.%iR' % (chrom, cov, res)][...]
        valid_rows = infile['valid.%s.%iC.%iR' % (chrom, cov, res)][...]
        corr = infile['corr.%s.%iC.%iR' % (chrom, cov, res)][...]
        valid = numpy.zeros(corr.shape, dtype=numpy.bool)
        N, M = corr.shape
        valid = numpy.zeros((N, M), dtype=numpy.int32)
        for i in range(min(N - 1, M)):
            P = N - i - 1
            valid[:P, i] = valid_rows[(i + 1):] * valid_rows[:P]
        valid[numpy.where((numpy.abs(dist) == numpy.inf) | (numpy.abs(corr) == numpy.inf))] = 0
        trans = numpy.zeros((N, M), dtype=numpy.float64)
        where = numpy.where(valid)
        trans[where] = corr[where] * dist[where]
        data[chrom] = [corr, dist, trans, valid]
    if len(data) == 0:
        return None
    return data

def find_results(data):
    chroms = data.keys()
    chroms.sort()
    results = numpy.zeros(len(chroms) + 1, dtype=numpy.float64)
    temp = numpy.zeros(4, dtype=numpy.float64)
    for j, chrom in enumerate(chroms):
        where = numpy.where(data[chrom][3])
        if where[0].shape[0] == 0:
            continue
        dist = numpy.sum(data[chrom][1][where])
        corr = numpy.sum(data[chrom][0][where])
        corrdist = numpy.sum(data[chrom][2][where])
        N = where[0].shape[0]
        if N == 0:
            continue
        results[j] = corrdist / dist - corr / N
        temp += [corrdist, dist, corr, N]
    results[-1] = temp[0] / temp[1] - temp[2] / temp[3]
    return results

def write_results(results, chroms, args):
    keys = results.keys()
    keys.sort()
    chroms.sort()
    output = open(args.output, 'w')
    temp = "Coverage\tResolution\tAll"
    for chrom in chroms:
        temp += "\t%s" % chrom
    print >> output, temp
    for key in keys:
        temp = [str(key[0]), str(key[1]), str(results[key][-1])]
        for i in range(results[key].shape[0] - 1):
            temp.append(str(results[key][i]))
        print >> output, '\t'.join(temp)
    output.close()

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Find a quality score for a HiC dataset"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="input", type=str, action='store', help="Quasar file name")
    parser.add_argument(dest="output", type=str, action='store', help="Results file name")
    return parser

if  __name__ == "__main__":
    main()
