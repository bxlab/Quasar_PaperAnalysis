#!/usr/bin/env python

import sys
import os
from math import ceil, floor

import numpy
import h5py
from pyx import *
from PIL import Image
from scipy.stats import linregress, spearmanr

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\usepackage{amsmath}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")
text.preamble(r"\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}")

page_width = 16.95
page_height = 22.5
text_size = -3

A_width = (page_width - 0.7) * 0.5

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_S1.pdf" % base_dir
    c = canvas.canvas()
    c.insert(plot_A(A_width))
    c.insert(plot_B(A_width), [trafo.translate(A_width + 0.6, 0)])
    c.writePDFfile(out_fname)

def plot_A(width):
    hoffset = 0.6
    voffset = 0.6
    voffset2 = 0.3
    spacer = 0.7
    pwidth = (width - hoffset - spacer * 2) / 3.0
    pheight = pwidth
    height = pheight * 3 + spacer * 2 + voffset + voffset2
    data_fname = "%s/Results/quasar_coverage_quality_results.txt" % base_dir
    data = load_data(data_fname)
    c = canvas.canvas()
    genomes = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    c.text(0, height, 'A', [text.halign.left, text.valign.top, text.size(0)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        samples = numpy.unique(data['sample'][numpy.where(data['genome'] == gen)])
        if gen == 'dm6':
            reses1 = [100000, 100000, 20000]
            reses2 = [20000, 4000, 4000]
            mincov = 1000000
        else:
            reses1 = [1000000, 1000000, 200000]
            reses2 = [200000, 40000, 40000]
            mincov = 10000000
        for j in range(len(reses1)):
            res1 = reses1[j]
            res2 = reses2[j]
            c1 = canvas.canvas([canvas.clip(path.rect(0, 0, pwidth, pheight))])
            c2 = canvas.canvas()
            c2.insert(c1)
            Xs = []
            Ys = []
            vsamples = []
            for k, sample in enumerate(samples):
                if numpy.where((data['sample'] == sample) & (data['param'] == mincov))[0].shape[0] == 0:
                    continue
                Xs.append(data['score'][numpy.where((data['resolution'] == res1) & (data['sample'] == sample) &
                    (data['param'] == 0))][0])
                Ys.append(data['score'][numpy.where((data['resolution'] == res2) & (data['sample'] == sample) &
                    (data['param'] == 0))][0])
                vsamples.append(sample)
            valid = numpy.zeros(data.shape[0], dtype=numpy.bool)
            for sample in vsamples:
                valid[numpy.where((data['sample'] == sample) & (data['param'] == 0))] = True
            temp = data[valid]
            minX = numpy.amin(temp['score'][numpy.where((temp['resolution'] == res1) & (temp['genome'] == gen))])
            maxX = numpy.amax(temp['score'][numpy.where((temp['resolution'] == res1) & (temp['genome'] == gen))])
            minY = numpy.amin(temp['score'][numpy.where((temp['resolution'] == res2) & (temp['genome'] == gen))])
            maxY = numpy.amax(temp['score'][numpy.where((temp['resolution'] == res2) & (temp['genome'] == gen))])
            spanX = maxX - minX
            spanY = maxY - minY
            minX -= spanX * 0.05
            maxX += spanX * 0.05
            minY -= spanY * 0.05
            maxY += spanY * 0.05
            spanX = maxX - minX
            spanY = maxY - minY
            ramp = Ramp(0.0, len(vsamples) - 1.0)
            Xs = numpy.array(Xs)
            Ys = numpy.array(Ys)
            Xs = (Xs - minX) / (maxX - minX) * pwidth
            Ys = (Ys - minY) / (maxY - minY) * pheight
            r, pval = spearmanr(Xs, Ys)
            c2.text(0.05, pheight - 0.05, r"r=%0.2f" % (r),
                [text.halign.left, text.valign.top, text.size(text_size)])
            c2.text(0.05, pheight - 0.3, r"p=%0.1e" % (pval),
                [text.halign.left, text.valign.top, text.size(text_size)])
            for k in range(Xs.shape[0]):
                c1.fill(path.circle(Xs[k], Ys[k], 0.04), [ramp.get_rgb_color(k)])
            if res1 >= 1000000:
                label1 = '%iMb' % (res1 / 1000000)
            else:
                label1 = '%iKb' % (res1 / 1000)
            if res2 >= 1000000:
                label2 = '%iMb' % (res2 / 1000000)
            else:
                label2 = '%iKb' % (res2 / 1000)
            c2.text(pwidth * 0.5, -0.5, "%s Qualities" % label1,
                [text.halign.center, text.valign.bottom, text.size(text_size)])
            c2.text(-0.5, pheight * 0.5, "%s Qualities" % label2,
                [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
            c2.stroke(path.rect(0, 0, pwidth, pheight))
            if j == 0:
                c2.text(pwidth * 0.5, pheight + 0.3, genomes[gen],
                    [text.halign.center, text.valign.top, text.size(text_size)])
            for k in [minX + 0.15 * spanX, maxX - 0.15 * spanX]:
                f = floor(numpy.log10(k)) - 1
                val = round(k / 10 ** f) * 10 ** f
                label = "%0.1e" % val
                label = label.replace("e+0", "e").replace("e-0", "e-")
                X = (val - minX) / (maxX - minX) * pwidth
                c2.stroke(path.line(X, 0, X, -0.05))
                c2.text(X, -0.1, label, [text.halign.center, text.valign.top, text.size(text_size)])
            for k in [minY + 0.15 * spanY, maxY - 0.15 * spanY]:
                f = floor(numpy.log10(k)) - 1
                val = round(k / 10 ** f) * 10 ** f
                label = "%0.1e" % val
                label = label.replace("e+0", "e").replace("e-0", "e-")
                Y = (val - minY) / (maxY - minY) * pwidth
                c2.stroke(path.line(0, Y, -0.05, Y))
                c2.text(-0.1, Y, label,
                    [text.halign.center, text.valign.bottom, text.size(text_size), trafo.rotate(90)])
            c.insert(c2, [trafo.translate(hoffset + i * (pwidth + spacer), voffset + (2 - j) * (pheight + spacer))])
    Xs = data[numpy.where(data)]
    return c

def plot_B(width):
    hoffset = 0.6
    voffset = 0.6
    voffset2 = 0.3
    spacer = 0.7
    pwidth = (width - hoffset - spacer * 2) / 3.0
    pheight = pwidth
    height = pheight * 3 + spacer * 2 + voffset + voffset2
    data_fname = "%s/Results/quasar_coverage_quality_results.txt" % base_dir
    data = load_data(data_fname)
    data = data[numpy.where(((data['genome'] != 'dm6') & (data['param'] == 10000000)) |
                            ((data['genome'] == 'dm6') & (data['param'] == 1000000)))]
    c = canvas.canvas()
    genomes = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    c.text(0, height, 'B', [text.halign.left, text.valign.top, text.size(0)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        samples = numpy.unique(data['sample'][numpy.where(data['genome'] == gen)])
        ramp = Ramp(0.0, len(samples) - 1.0)
        if gen == 'dm6':
            reses1 = [100000, 100000, 20000]
            reses2 = [20000, 4000, 4000]
        else:
            reses1 = [1000000, 1000000, 200000]
            reses2 = [200000, 40000, 40000]
        for j in range(len(reses1)):
            res1 = reses1[j]
            res2 = reses2[j]
            c1 = canvas.canvas([canvas.clip(path.rect(0, 0, pwidth, pheight))])
            c2 = canvas.canvas()
            c2.insert(c1)
            minX = numpy.amin(data['score'][numpy.where((data['resolution'] == res1) & (data['genome'] == gen))])
            maxX = numpy.amax(data['score'][numpy.where((data['resolution'] == res1) & (data['genome'] == gen))])
            minY = numpy.amin(data['score'][numpy.where((data['resolution'] == res2) & (data['genome'] == gen))])
            maxY = numpy.amax(data['score'][numpy.where((data['resolution'] == res2) & (data['genome'] == gen))])
            spanX = maxX - minX
            spanY = maxY - minY
            minX -= spanX * 0.05
            maxX += spanX * 0.05
            minY -= spanY * 0.05
            maxY += spanY * 0.05
            spanX = maxX - minX
            spanY = maxY - minY
            Xs = []
            Ys = []
            for k, sample in enumerate(samples):
                Xs.append(data['score'][numpy.where((data['resolution'] == res1) & (data['sample'] == sample))][0])
                Ys.append(data['score'][numpy.where((data['resolution'] == res2) & (data['sample'] == sample))][0])
            Xs = numpy.array(Xs)
            Ys = numpy.array(Ys)
            Xs = (Xs - minX) / (maxX - minX) * pwidth
            Ys = (Ys - minY) / (maxY - minY) * pheight
            r, pval = spearmanr(Xs, Ys)
            c2.text(0.05, pheight - 0.05, r"r=%0.2f" % (r),
                [text.halign.left, text.valign.top, text.size(text_size)])
            c2.text(0.05, pheight - 0.3, r"p=%0.1e" % (pval),
                [text.halign.left, text.valign.top, text.size(text_size)])
            for k in range(Xs.shape[0]):
                c1.fill(path.circle(Xs[k], Ys[k], 0.04), [ramp.get_rgb_color(k)])
            if res1 >= 1000000:
                label1 = '%iMb' % (res1 / 1000000)
            else:
                label1 = '%iKb' % (res1 / 1000)
            if res2 >= 1000000:
                label2 = '%iMb' % (res2 / 1000000)
            else:
                label2 = '%iKb' % (res2 / 1000)
            c2.text(pwidth * 0.5, -0.5, "%s Qualities" % label1,
                [text.halign.center, text.valign.bottom, text.size(text_size)])
            c2.text(-0.5, pheight * 0.5, "%s Qualities" % label2,
                [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
            c2.stroke(path.rect(0, 0, pwidth, pheight))
            if j == 0:
                c2.text(pwidth * 0.5, pheight + 0.3, genomes[gen],
                    [text.halign.center, text.valign.top, text.size(text_size)])
            for k in [minX + 0.15 * spanX, maxX - 0.15 * spanX]:
                f = floor(numpy.log10(k)) - 1
                val = round(k / 10 ** f) * 10 ** f
                label = "%0.1e" % val
                label = label.replace("e+0", "e").replace("e-0", "e-")
                X = (val - minX) / (maxX - minX) * pwidth
                c2.stroke(path.line(X, 0, X, -0.05))
                c2.text(X, -0.1, label, [text.halign.center, text.valign.top, text.size(text_size)])
            for k in [minY + 0.15 * spanY, maxY - 0.15 * spanY]:
                f = floor(numpy.log10(k)) - 1
                val = round(k / 10 ** f) * 10 ** f
                label = "%0.1e" % val
                label = label.replace("e+0", "e").replace("e-0", "e-")
                Y = (val - minY) / (maxY - minY) * pwidth
                c2.stroke(path.line(0, Y, -0.05, Y))
                c2.text(-0.1, Y, label,
                    [text.halign.center, text.valign.bottom, text.size(text_size), trafo.rotate(90)])
            c.insert(c2, [trafo.translate(hoffset + i * (pwidth + spacer), voffset + (2 - j) * (pheight + spacer))])
    Xs = data[numpy.where(data)]
    return c

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
