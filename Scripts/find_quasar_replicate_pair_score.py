#!/usr/bin/env python

import sys
import glob
import argparse

import h5py
import numpy

def main():
    parser = generate_parser()
    args = parser.parse_args()
    infile1 = h5py.File(args.input1, 'r')
    infile2 = h5py.File(args.input2, 'r')
    if 'chromosomes' not in infile1 or 'chromosomes' not in infile2:
        infile1.close()
        infile2.close()
        return None
    chroms = numpy.intersect1d(infile1['chromosomes'][...], infile2['chromosomes'][...])
    keys1 = []
    for key in infile1.keys():
        temp = key.split('.')
        if len(temp) == 4:
            keys1.append(tuple([int(temp[2][:-1]), int(temp[3][:-1])]))
    keys1 = set(keys1)
    keys2 = []
    for key in infile2.keys():
        temp = key.split('.')
        if len(temp) == 4:
            keys2.append(tuple([int(temp[2][:-1]), int(temp[3][:-1])]))
    keys2 = set(keys2)
    keys = keys1.intersection(keys2)
    results = {}
    for key in keys:
        data1 = load_data(infile1, chroms, key[0], key[1])
        data2 = load_data(infile2, chroms, key[0], key[1])
        if data1 is None or data2 is None:
            continue
        name1 = '.'.join(args.input1.split('/')[-1].split('.')[:-1])
        name2 = '.'.join(args.input2.split('/')[-1].split('.')[:-1])
        results[(name1, name2, key[0], key[1])] = correlate_samples(data1, data2, key[1])
    write_results(results, chroms, args)
    infile1.close()

def load_data(infile, chroms, cov, res):
    starts = infile['starts'][...]
    data = {}
    for i, chrom in enumerate(chroms):
        if 'valid.%s.%iC.%iR' % (chrom, cov, res) not in infile:
            data[chrom] = None
            continue
        start = (starts[i] / res) * res
        dist = (infile['dist.%s.%iC.%iR' % (chrom, cov, res)][...] + 1) ** 0.5
        valid_rows = infile['valid.%s.%iC.%iR' % (chrom, cov, res)][...]
        corr = infile['corr.%s.%iC.%iR' % (chrom, cov, res)][...]
        valid = numpy.zeros(corr.shape, dtype=numpy.bool)
        N, M = corr.shape
        valid = numpy.zeros((N, M), dtype=numpy.int32)
        for i in range(min(N - 1, M)):
            P = N - i - 1
            valid[:P, i] = valid_rows[(i + 1):] * valid_rows[:P]
        valid[numpy.where((numpy.abs(dist) == numpy.inf) | (numpy.abs(corr) == numpy.inf))] = False
        trans = numpy.zeros(corr.shape, dtype=numpy.float64)
        where = numpy.where(valid)
        trans[where] = corr[where] * dist[where]
        data[chrom] = [start, trans, valid]
    return data

def correlate_samples(data1, data2, res):
    chroms = data1.keys()
    chroms.sort()
    results = numpy.zeros(len(chroms) + 1, dtype=numpy.float64)
    temp = numpy.zeros(6, dtype=numpy.float64)
    for j, chrom in enumerate(chroms):
        if chrom not in data2 or data1[chrom] is None or data2[chrom] is None:
            continue
        start1 = data1[chrom][0]
        start2 = data2[chrom][0]
        if start2 > start1:
            start1 = (start2 - start1) / res
            start2 = 0
        elif start2 > start1:
            start2 = (start1 - start2) / res
            start1 = 0
        else:
            start1 = 0
            start2 = 0
        stop1 = data1[chrom][1].shape[0] - start1
        stop2 = data2[chrom][1].shape[0] - start2
        if stop2 > stop1:
            stop2 = stop1 + start2
            stop1 = start1 + stop1
        else:
            stop1 = stop2 + start1
            stop2 = start2 + stop2
        stop3 = min(data1[chrom][1].shape[1], data2[chrom][1].shape[1])
        valid1 = data1[chrom][2][start1:stop1, :stop3]
        valid2 = data2[chrom][2][start2:stop2, :stop3]
        valid = numpy.where(valid1 & valid2)
        if valid[0].shape[0] == 0:
            continue
        trans1 = data1[chrom][1][start1:stop1, :stop3][valid]
        trans2 = data2[chrom][1][start2:stop2, :stop3][valid]
        X = numpy.sum(trans1)
        Y = numpy.sum(trans2)
        X2 = numpy.sum(trans1 ** 2.0)
        Y2 = numpy.sum(trans2 ** 2.0)
        XY = numpy.sum(trans1 * trans2)
        N = valid[0].shape[0]
        if N == 0:
            continue
        temp += [X, Y, X2, Y2, XY, N]
        Xmu = X / N
        Ymu = Y / N
        X2mu = X2 / N
        Y2mu = Y2 / N
        XYmu = XY / N
        if Xmu ** 2.0 > X2mu or Ymu ** 2.0 > Y2mu:
            continue
        Xstd = (X2mu - Xmu ** 2.0) ** 0.5
        Ystd = (Y2mu - Ymu ** 2.0) ** 0.5
        if Xstd == 0 or Ystd == 0:
            continue
        results[j] = (XYmu - Xmu * Ymu) / (Xstd * Ystd)
    Xmu = temp[0] / temp[5]
    Ymu = temp[1] / temp[5]
    X2mu = temp[2] / temp[5]
    Y2mu = temp[3] / temp[5]
    XYmu = temp[4] / temp[5]
    Xstd = (X2mu - Xmu ** 2.0) ** 0.5
    Ystd = (Y2mu - Ymu ** 2.0) ** 0.5
    results[-1] = (XYmu - Xmu * Ymu) / (Xstd * Ystd)
    return results

def write_results(results, chroms, args):
    chroms.sort()
    output = open(args.output, 'w')
    temp = "Sample1\tSample2\tCoverage\tResolution\tAll"
    for chrom in chroms:
        temp += "\t%s" % chrom
    print >> output, temp
    keys = results.keys()
    keys.sort()
    for key in keys:
        if results[key][-1] == 0.0:
            continue
        temp = [key[0], key[1], str(key[2]), str(key[3]), str(results[key][-1])]
        for i in range(results[key].shape[0] - 1):
            temp.append(str(results[key][i]))
        print >> output, '\t'.join(temp)
    output.close()

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Find a quality score for a HiC dataset"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="input1", type=str, action='store', help="First quasar file name")
    parser.add_argument(dest="input2", type=str, action='store', help="Second quasar file name")
    parser.add_argument(dest="output", type=str, action='store', help="Results file name")
    return parser

if  __name__ == "__main__":
    main()
