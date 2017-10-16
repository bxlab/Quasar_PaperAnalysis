#!/usr/bin/env python

import sys
import os
from math import ceil

import numpy
import h5py
from pyx import *
from PIL import Image
from scipy.optimize import curve_fit

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\usepackage{amsmath}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")
text.preamble(r"\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}")

page_width = 17.0
page_height = 22.5
text_size = -3

A_width = (page_width - 0.5) * 0.5
B_width = A_width

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_3.pdf" % base_dir
    c = canvas.canvas()
    c.insert(plot_A(A_width))
    c.insert(plot_B(B_width), [trafo.translate(A_width + 0.4, 0)])
    c.writePDFfile(out_fname)

def plot_A(width):
    hoffset = 0.7
    voffset = 0.55
    spacer = 0.25
    pwidth = (width - hoffset - 2 * spacer) / 3.0
    data_fname = "%s/Results/quasar_noise_replicate_results.txt" % base_dir
    data = load_data(data_fname)
    data['param'] /= 1000.0
    c = canvas.canvas()
    c.text(0, pwidth * 3 + spacer * 2 + voffset + 0.3, 'A', [text.halign.left, text.valign.top, text.size(0)])
    resolutions = {}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        where = numpy.where(data['genome'] == gen)[0]
        resolutions[gen] = numpy.unique(data['resolution'][where])
    genomes = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        where = numpy.where(data['genome'] == gen)[0]
        minX = 0.0
        maxX = 0.75
        for j in range(resolutions[gen].shape[0]):
            minY = 0.0
            maxY = 1.05
            where = numpy.where((data['genome'] == gen) & (data['resolution'] == resolutions[gen][j]))[0]
            p = canvas.canvas()
            for k in numpy.linspace(0, 60, 4):
                label = "%i" % k
                X = (k / 100.0 - minX) / (maxX - minX) * pwidth
                if j == 2:
                    p.stroke(path.line(X, 0, X, -0.05))
                    p.text(X, -0.1, label, [text.halign.center, text.valign.top, text.size(text_size)])
                p.stroke(path.line(X, 0, X, pwidth), [color.gray(0.9)])
            for k in numpy.linspace(0, 1, 6):
                Y = (k - minY) / (maxY - minY) * pwidth
                if i == 0:
                    p.stroke(path.line(0, Y, -0.05, Y))
                    p.text(-0.1, Y, "%0.1f" % k, [text.halign.right, text.valign.middle, text.size(text_size)])
                p.stroke(path.line(0, Y, pwidth, Y), [color.gray(0.9)])
            p.insert(plot_A_res(data[where], minX, maxX, minY, maxY, pwidth, resolutions[gen][j]))
            if resolutions[gen][j] >= 1000000:
                label = '%iMb' % (resolutions[gen][j] / 1000000)
            else:
                label = '%iKb' % (resolutions[gen][j] / 1000)
            p.text(0.1, 0.1, label, [text.valign.bottom, text.halign.left, text.size(text_size)])
            if j == 0:
                p.text(pwidth * 0.5, pwidth + 0.3, genomes[gen],
                    [text.halign.center, text.valign.top, text.size(text_size)])
            c.insert(p, [trafo.translate(hoffset + (pwidth + spacer) * i, voffset + (pwidth + spacer) * (2 - j))])
    c.text(0, voffset + 1.5 * pwidth + spacer, 'Replicate Score',
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + 1.5 * pwidth + spacer, 0, r'\% Noise Reads',
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    return c

def plot_B(width):
    hoffset = 0.7
    voffset = 0.55
    spacer = 0.25
    pwidth = (width - hoffset - 2 * spacer) / 3.0
    data_fname = "%s/Results/quasar_coverage_replicate_results.txt" % base_dir
    data = load_data(data_fname)
    data = data[numpy.where(data['param'] > 0)]
    data['param'] = numpy.log10(data['param'])
    c = canvas.canvas()
    c.text(0, pwidth * 3 + spacer * 2 + voffset + 0.3, 'B', [text.halign.left, text.valign.top, text.size(0)])
    resolutions = {}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        where = numpy.where(data['genome'] == gen)[0]
        resolutions[gen] = numpy.unique(data['resolution'][where])
    genomes = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        where = numpy.where(data['genome'] == gen)[0]
        minX = numpy.amin(data['param'][where])
        maxX = numpy.amax(data['param'][where])
        for j in range(resolutions[gen].shape[0]):
            minY = 0.0
            maxY = 1.05
            where = numpy.where((data['genome'] == gen) & (data['resolution'] == resolutions[gen][j]))[0]
            p = canvas.canvas()
            coverages = numpy.unique(data['param'][where])
            start = 0
            if gen == 'dm6':
                start = 1
            stop = ceil((coverages[-1] - coverages[0]) / numpy.log10(4))
            cov_labels = numpy.linspace(coverages[start] - numpy.log10(4),
                                        coverages[start] + numpy.log10(4) * stop, stop + 2)
            for k in cov_labels[1:-1]:
                cov = numpy.round(10 ** k)
                if cov >= 40000000:
                    continue
                if cov >= 1000000:
                    label = "%i" % (cov / 1000000)
                else:
                    label = "%0.2f" % (cov / 1000000.)
                X = (k - minX) / (maxX - minX) * pwidth
                if j == 2:
                    p.stroke(path.line(X, 0, X, -0.085))
                    p.text(X, -0.1, label, [text.halign.center, text.valign.top, text.size(text_size)])
                    span = cov_labels[1] - cov_labels[0]
                p.stroke(path.line(X, 0, X, pwidth), [color.gray(0.9)])
            if j == 2:
                for k in range(cov_labels.shape[0] - 1):
                    for m in numpy.linspace(10 ** cov_labels[k], 10 ** cov_labels[k + 1], 9)[1:-1]:
                        X1 = (numpy.log10(m) - minX) / (maxX - minX) * pwidth
                        if X1 >= 0 and X1 < pwidth:
                            p.stroke(path.line(X1, 0, X1, -0.05))
            for k in numpy.linspace(0, 1, 6):
                Y = (k - minY) / (maxY - minY) * pwidth
                if i == 0:
                    p.stroke(path.line(0, Y, -0.05, Y))
                    p.text(-0.1, Y, "%0.1f" % k, [text.halign.right, text.valign.middle, text.size(text_size)])
                p.stroke(path.line(0, Y, pwidth, Y), [color.gray(0.9)])
            p.insert(plot_A_res(data[where], minX, maxX, minY, maxY, pwidth))
            if resolutions[gen][j] >= 1000000:
                label = '%iMb' % (resolutions[gen][j] / 1000000)
            else:
                label = '%iKb' % (resolutions[gen][j] / 1000)
            p.text(pwidth - 0.1, 0.1, label, [text.valign.bottom, text.halign.right, text.size(text_size)])
            if j == 0:
                p.text(pwidth * 0.5, pwidth + 0.3, genomes[gen],
                    [text.halign.center, text.valign.top, text.size(text_size)])
            c.insert(p, [trafo.translate(hoffset + (pwidth + spacer) * i, voffset + (pwidth + spacer) * (2 - j))])
    c.text(0, voffset + 1.5 * pwidth + spacer, 'Replicate Score',
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + 1.5 * pwidth + spacer, 0, 'Millions of Reads',
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    return c

def plot_A_res(data, minX, maxX, minY, maxY, width, temp=False):
    c = canvas.canvas()
    samples = numpy.unique(data['sample1'])
    ramp = Ramp(0.0, len(samples) - 1)
    for h, sample in enumerate(samples):
        Xs = data['param'][numpy.where(data['sample1'] == sample)]
        noi = numpy.copy(Xs)
        Ys = data['score'][numpy.where(data['sample1'] == sample)]
        order = numpy.argsort(Xs)
        Xs = Xs[order]
        Ys = Ys[order]
        if numpy.amax(Ys[1:]) > Ys[0] and temp:
            if Ys[-1] / Ys[0] < 1.0:
                print sample, temp, numpy.amax(Ys/Ys[0]), noi[numpy.where((Ys[:-1]/Ys[0] > 1.0) & (Ys[1:]/Ys[0] < 1.0))[0][0] + 1]
            else:
                print sample, temp, numpy.amax(Ys/Ys[0]), noi[-1]
        Xs = (Xs - minX) / (maxX - minX) * width
        Ys = (Ys - minY) / (maxY - minY) * width
        lpath = path.path(path.moveto(Xs[0], Ys[0]))
        for i in range(1, Xs.shape[0]):
            lpath.append(path.lineto(Xs[i], Ys[i]))
        c.stroke(lpath, [ramp.get_rgb_color(h)])
    c.stroke(path.rect(0, 0, width, width))
    return c

def load_data(fname):
    data = []
    infile = open(fname)
    line = infile.readline()
    line = infile.readline()
    while line:
        temp = line.split('\t')
        data.append((temp[0], temp[1], temp[2], int(temp[3]), int(temp[4]), float(temp[5])))
        line = infile.readline()
    data = numpy.array(data, dtype=numpy.dtype([('sample1', 'S80'), ('sample2', 'S80'), ('genome', 'S4'),
        ('param', numpy.float64), ('resolution', numpy.int32), ('score', numpy.float32)]))
    return data

class Ramp():
    def __init__(self, minval, maxval):
        self.min = float(minval)
        self.max = float(maxval)
        self.breaks = numpy.linspace(minval, maxval, 8)
        self.span = self.breaks[1] - self.breaks[0]
        self.colors = numpy.array([
            [0, 0, 0],
            [26, 26, 89],
            [77, 38, 166],
            [153, 51, 128],
            [255, 64, 38],
            [230, 140, 0],
            [230, 191, 26],
            [230, 230, 128]], dtype=numpy.uint32)
        self.base = 256 ** 3 * 255
        self.scale = numpy.array([1, 256, 256 ** 2], dtype=numpy.uint32)

    def get_color(self, val):
        index = numpy.searchsorted(self.breaks[1:-1], val)
        frac = numpy.maximum(0.0, numpy.minimum(1.0, (val - self.breaks[index]) / self.span))
        cvals = (self.colors[index, :] * (1.0 - frac) +
                 self.colors[index + 1, :] * frac) / 255.0
        return color.rgb(r=cvals[0], g=cvals[1], b=cvals[2])

    def get_rgb_color(self, val):
        index = numpy.searchsorted(self.breaks[1:-1], val)
        frac = max(0.0, min(1.0, (val - self.breaks[index]) / self.span))
        cvals = (self.colors[index] * (1.0 - frac) + self.colors[index + 1] * frac) / 255.0
        return color.rgb(r=cvals[0], g=cvals[1], b=cvals[2])


if __name__ == "__main__":
    main()
