#!/usr/bin/env python

import sys
import os
from math import ceil, floor

import numpy
import h5py
from pyx import *
from PIL import Image
from scipy.stats import linregress

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\usepackage{amsmath}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")
text.preamble(r"\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}")

page_width = 16.9
page_height = 22.5
text_size = -3

spacer = 0.4
A_width = (page_width - spacer) / 2.0
pwidth = (A_width - 1.45) / 3.0
A_height = pwidth * 3 + 1.85

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_S2.pdf" % base_dir
    c = canvas.canvas()
    data_fname = "%s/Results/quasar_coverage_quality_results.txt"%  base_dir
    data = load_data(data_fname)
    ideal_fname = "%s/Results/quality_coverage_params.txt" % base_dir
    ideal = load_ideal(ideal_fname)
    sample_fname = "%s/Data/samples.txt" % base_dir
    sample_data = load_sample_data(sample_fname)
    letters = ['A', 'B', 'C', 'D']
    for i, var in enumerate(['insert_size', 'same_fragment', 'failed_cut']):
        c.insert(plot_subfigure(A_width, var, letters[i], data, ideal, sample_data),
            [trafo.translate((A_width + spacer) * (i % 2), (A_height + spacer) * (1 - i / 2))])
    c.writePDFfile(out_fname)

def plot_subfigure(width, var, letter, data, ideal, sample_data):
    hoffset = 0.55
    hoffset2 = 0.3
    spacer = 0.3
    spacer2 = 0.45
    voffset = 0.6
    voffset2 = 0.35
    pwidth = (width - hoffset - hoffset2 - spacer * 2) / 3.0
    pheight = pwidth
    height = pheight * 3 + voffset + voffset2 + spacer2 * 2
    Y_data = {}
    if var == 'ratio':
        for key in sample_data:
            Y_data[key] = sample_data[key]['valid_cis_reads'] / float(sample_data[key]['valid_trans_reads'] + sample_data[key]['valid_cis_reads'])
    else:
        for key in sample_data:
            Y_data[key] = float(sample_data[key][var]) / float(sample_data[key]['total_reads'])
    c = canvas.canvas()
    genomes = {'hg38': 'Human, 1Mb', 'mm10': 'Mouse, 1Mb', 'dm6': 'Fruit Fly, 50Kb'}
    c.text(0, height, letter, [text.halign.left, text.valign.top, text.size(0)])
    where = numpy.where((((data['genome'] == 'dm6') & (data['resolution'] == 100000)) |
                         ((data['genome'] != 'dm6') & (data['resolution'] == 1000000))) &
                         (data['param'] == 0))[0]
    c.insert(plot_panel(width, Y_data, data[where]), [trafo.translate(0, voffset + 2 * (pheight + spacer2))])
    where = numpy.where(((data['genome'] == 'dm6') & (data['resolution'] == 100000) & (data['param'] == 1000000)) |
                        ((data['genome'] != 'dm6') & (data['resolution'] == 1000000) & (data['param'] == 10000000)))[0]
    c.insert(plot_panel(width, Y_data, data[where]), [trafo.translate(0, voffset + (pheight + spacer2))])
    c.insert(plot_panel(width, Y_data, ideal), [trafo.translate(0, voffset)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        c.text(hoffset + pwidth * (i + 0.5) + spacer * i, height - 0.025, genomes[gen],
            [text.halign.center, text.valign.top, text.size(text_size)])
    c.text(hoffset + pwidth * 1.5 + spacer, 0, 'Quality Score',
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    labels = {
        'ratio': 'Valid Cis / Valid Trans Read Ratio',
        'insert_size': r'\% Insert Too Large Reads',
        'same_fragment': r'\% Circular Fragment Reads',
        'failed_cut': r'\% Putative Failed Restriction Cut Reads',
    }
    c.text(0, voffset + pheight * 1.5 + spacer2, labels[var],
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    for i, label in enumerate(['Max Score', 'Uniform Coverage', 'Predicted Max']):
        c.text(width, voffset + pheight * (2.5 - i) + spacer2 * (2 - i), label,
            [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(-90)])
    return c

def plot_panel(width, Y_data, data):
    hoffset = 0.55
    hoffset2 = 0.3
    spacer = 0.3
    spacer2 = 0.45
    voffset = 0.6
    voffset2 = 0.35
    pwidth = (width - hoffset - hoffset2 - spacer * 2) / 3.0
    pheight = pwidth
    c = canvas.canvas()
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        samples = numpy.unique(data['sample'][numpy.where(data['genome'] == gen)])
        ramp = Ramp(0.0, len(samples) - 1.0)
        c1 = canvas.canvas([canvas.clip(path.rect(0, 0, pwidth, pheight))])
        c2 = canvas.canvas()
        c2.insert(c1)
        Xs = []
        Ys = []
        for k, sample in enumerate(samples):
            Xs.append(data['score'][numpy.where(data['sample'] == sample)][0])
            Ys.append(Y_data[sample])
        Xs = numpy.array(Xs)
        Ys = numpy.array(Ys)
        minX = numpy.amin(Xs)
        maxX = numpy.amax(Xs)
        minY = numpy.amin(Ys)
        maxY = numpy.amax(Ys)
        spanX = maxX - minX
        spanY = maxY - minY
        minX -= spanX * 0.05
        maxX += spanX * 0.05
        minY -= spanY * 0.05
        maxY += spanY * 0.2
        spanX = maxX - minX
        spanY = maxY - minY
        Xs = (Xs - minX) / (maxX - minX) * pwidth
        Ys = (Ys - minY) / (maxY - minY) * pheight
        slope, intercept, r, pval = linregress(Xs, Ys)[:4]
        c1.stroke(path.line(0, intercept, pwidth, pwidth * slope + intercept))
        c2.text(0.05, pheight - 0.05, r"r=%0.2f" % (r),
            [text.halign.left, text.valign.top, text.size(text_size)])
        c2.text(0.05, pheight - 0.3, r"p=%0.1e" % (pval),
            [text.halign.left, text.valign.top, text.size(text_size)])
        for k in range(Xs.shape[0]):
            c1.fill(path.circle(Xs[k], Ys[k], 0.04), [ramp.get_rgb_color(k)])
        c2.stroke(path.rect(0, 0, pwidth, pheight))
        scale = 10 ** -ceil(numpy.log10(maxX) - 2)
        val1 = ceil((minX + 0.15 * spanX) * scale) / scale
        X = (val1 - minX) / spanX * pwidth
        c2.stroke(path.line(X, 0, X, -0.05))
        c2.text(X, -0.1, "%0.1e" % val1, [text.valign.top, text.halign.center, text.size(text_size)])
        val2 = floor((maxX - 0.15 * spanX) * scale) / scale
        X = (val2 - minX) / spanX * pwidth
        c2.stroke(path.line(X, 0, X, -0.05))
        c2.text(X, -0.1, "%0.1e" % val2, [text.valign.top, text.halign.center, text.size(text_size)])
        if i == 0:
            scale = 10 ** -ceil(numpy.log10(maxY) - 2)
            val1 = ceil((minY + 0.1 * spanY) * scale) / scale
            Y = (val1 - minY) / spanY * pheight
            c2.stroke(path.line(0, Y, -0.05, Y))
            label = "%0.1e" % val1
            label = label.replace('e-0', 'e-').replace('e+0', 'e')
            c2.text(-0.1, Y, label,
                [text.valign.bottom, text.halign.center, text.size(text_size), trafo.rotate(90)])
            val2 = floor((maxY - 0.1 * spanY) * scale) / scale
            label = "%0.1e" % val2
            label = label.replace('e-0', 'e-').replace('e+0', 'e')
            Y = (val2 - minY) / spanY * pheight
            c2.stroke(path.line(0, Y, -0.05, Y))
            c2.text(-0.1, Y, label,
                [text.valign.bottom, text.halign.center, text.size(text_size), trafo.rotate(90)])
        c.insert(c2, [trafo.translate(hoffset + i * (pwidth + spacer), 0)])
    return c

    c.text(hoffset + pwidth * 1.5 + hoffset2, 0, 'Quality Score',
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    labels = {
        'ratio': 'Valid Cis / Valid Trans Read Ratio',
        'insert_size': r'\% Insert Too Large Reads',
        'same_fragment': r'\% Circular Fragment Reads',
        'failed_cut': r'\% Putative Failed Restriction Cut Reads',
    }
    c.text(0, voffset + pheight * 1.5 + spacer, labels[var],
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    Xs = data[numpy.where(data)]
    return c

def load_sample_data(fname):
    data = {}
    labels = None
    for line in open(fname):
        temp = line.rstrip('\n').split('\t')
        if labels is None:
            labels = temp[1:]
            continue
        temp1 = {}
        for i, val in enumerate(temp[1:]):
            try:
                temp1[labels[i]] = int(val)
            except:
                temp1[labels[i]] = val
        data[temp[0]] = temp1
    return data

def load_data(fname):
    data = []
    infile = open(fname)
    line = infile.readline()
    line = infile.readline()
    while line:
        temp = line.split('\t')
        data.append((temp[0], temp[1], float(temp[2]), int(temp[3]), float(temp[4])))
        line = infile.readline()
    data = numpy.array(data, dtype=numpy.dtype([('sample', 'S80'), ('genome', 'S4'), ('param', numpy.float64),
                       ('resolution', numpy.float32), ('score', numpy.float32)]))
    return data

def load_ideal(fname):
    data = []
    infile = open(fname)
    line = infile.readline()
    line = infile.readline()
    while line:
        temp = line.split('\t')
        data.append((temp[0], temp[1], float(temp[2])))
        line = infile.readline()
    data = numpy.array(data, dtype=numpy.dtype([('sample', 'S80'), ('genome', 'S80'), ('score', numpy.float64)]))
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
