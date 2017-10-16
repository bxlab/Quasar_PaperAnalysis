#!/usr/bin/env python

import sys
import os
from math import ceil, floor

import numpy
import h5py
from pyx import *
from scipy.stats import linregress
from scipy.optimize import curve_fit

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\usepackage{amsmath}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")
text.preamble(r"\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}")

page_width = 16.9
page_height = 22.4
text_size = -3
hoffset = 3.2
voffset = 0.75
voffset2 = 0.25
pwidth = (page_width - hoffset) / 4.0
pheight = page_height - voffset - voffset2 * 3

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_S3.pdf" % base_dir
    c = canvas.canvas()
    data_fname = "%s/Results/quasar_replicate_results.txt" % base_dir
    data = load_rdata(data_fname)
    resolutions = {}
    sizes = numpy.zeros(4, dtype=numpy.float64)
    all_samples = {}
    genomes = ['hg38', 'mm10', 'dm6']
    total = 0
    for i, gen in enumerate(genomes):
        all_samples[gen] = []
        sizes[i + 1] = sizes[i]
        where = numpy.where(data['genome'] == gen)
        resolutions[gen] = numpy.unique(data['resolution'][numpy.where(data['genome'] == gen)])
        samples = numpy.unique(numpy.r_[data['sample1'][where], data['sample2'][where]])
        for samp in samples:
            if samp.count('Pseudo') == 0:
                sizes[i + 1] += 1
                all_samples[gen].append(samp)
    bounds = numpy.zeros((4, 2), dtype=numpy.float64)
    bounds[:, 0] = 0
    bounds[:, 1] = 1.0
    spans = bounds[:, 1] - bounds[:, 0]
    bounds[:, 0] -= spans * 0.05
    bounds[:, 1] += spans * 0.05
    total = sizes[-1]
    heights = (sizes[1:] - sizes[:-1]) / sizes[-1] * pheight + voffset2
    span = pheight / total
    pos = page_height - span * 0.2 - voffset2
    pos2 = page_height
    gen_label = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    for i, gen in enumerate(genomes):
        for j in range(len(all_samples[gen])):
            temp = all_samples[gen][j].split('_')
            c.text(hoffset - 0.01, pos, ' '.join([temp[0], temp[2], temp[4]]),
                [text.halign.right, text.valign.top, text.size(text_size - 1)])
            pos -= span
        pos -= voffset2
        pos2 -= heights[i]
        c.text(0, pos2 + (heights[i] - voffset2) * 0.5, gen_label[gen],
            [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
        for j in range(resolutions[gen].shape[0]):
            where = numpy.where((data['genome'] == gen) & (data['resolution'] == resolutions[gen][j]))[0]
            c.insert(plot_samples(data[where], heights[i], resolutions[gen][j], all_samples[gen], bounds[j, :],
                ((i == 2) | ((i == 1) & (j == 3)))), [trafo.translate(hoffset + pwidth * j, pos2)])
    c.insert(plot_key(heights[-1] * 0.6), [trafo.translate(hoffset + pwidth * 3 + 0.5, 0.5)])
    c.writePDFfile(out_fname)

def plot_samples(data, height, res, samples, bounds, Xlab):
    pheight = height - voffset2
    step = pheight / float(len(samples))
    indices = {}
    for i, sample in enumerate(samples):
        indices[sample] = pheight - (i + 0.5) * step
    c = canvas.canvas()
    maxX = bounds[1]
    minX = bounds[0]
    span = maxX - minX
    for i in range(int(ceil(minX * 10)), int(floor(maxX * 10) + 1), 2):
        val = i / 10.0
        X = (val - minX) / span * pwidth
        if Xlab:
            c.stroke(path.line(X, 0, X, -0.05))
            c.text(X, -0.1, "%0.1f" % val,
                [text.halign.left, text.valign.middle, text.size(text_size), trafo.rotate(-90)])
        c.stroke(path.line(X, 0, X, pheight), [color.gray(0.9)])
    if Xlab:
        c.text(pwidth * 0.5, -voffset, "Replicate Score",
            [text.halign.center, text.valign.bottom, text.size(text_size)])
    ramp = Ramp(0.0, 4.0)
    strokes = [
        path.circle(0, 0, 0.08),
        path.path(path.moveto(-0.1, -0.0867), path.lineto(0, 0.0867), path.lineto(0.1, -0.0867), path.closepath()),
        path.rect(-0.075, -0.075, 0.15, 0.15),
        path.path(path.moveto(-0.1, 0), path.lineto(0, 0.1), path.lineto(0.1, 0),
                  path.lineto(0, -0.1), path.closepath()),
    ]
    points = [[], [], [], []]

    for i in range(data.shape[0]):
        name1 = data['sample1'][i]
        name2 = data['sample2'][i]
        t1 = name1.split('_')[2].split('H1-')[-1].split('-Dil')[0].split('-IS')[0]
        t2 = name2.split('_')[2].split('H1-')[-1].split('-Dil')[0].split('-IS')[0]
        val = (data['score'][i] - minX) / span * pwidth
        if name1.count('Pseudo') + name2.count('Pseudo') > 0:
            jitter = numpy.random.normal(0, step * 0.1)
            index = 1
        elif name1.split('Rep')[0] == name2.split('Rep')[0]:
            jitter = 0
            index = 3
        elif t1 == t2:
            jitter = numpy.random.normal(0, step * 0.1)
            index = 2
        else:
            jitter = numpy.random.normal(0, step * 0.1)
            index = 0
        if name1 in indices:
            Y = indices[name1]
            points[index].append([val, Y + jitter])
        if name2 in indices:
            Y = indices[name2]
            points[index].append([val, Y + jitter])
    for i in range(len(points)):
        for j in range(len(points[i])):
            if i == 0:
                c.stroke(strokes[i],
                    [trafo.scale(0.6), trafo.translate(points[i][j][0], points[i][j][1]), ramp.get_rgb_color(i)])
            else:
                c.stroke(strokes[i], [trafo.scale(0.6), style.linewidth.thin, deco.filled([ramp.get_rgb_color(i)]),
                    trafo.translate(points[i][j][0], points[i][j][1])])
    c.stroke(path.rect(0, 0, pwidth, pheight))
    if res >= 1000:
        label = "%iMb" % (res / 1000)
    else:
        label = "%iKb" % (res)
    c.text(pwidth * 0.5, height - 0.02, label, [text.halign.center, text.valign.top, text.size(text_size)])
    for i in range(1, len(samples)):
        c.stroke(path.line(0, step * i, pwidth, step * i), [style.linestyle.dotted])
    return c

def plot_key(height):
    ramp = Ramp(0.0, 4.0)
    strokes = [
        path.circle(0, 0, 0.08),
        path.path(path.moveto(-0.1, -0.0867), path.lineto(0, 0.0867), path.lineto(0.1, -0.0867), path.closepath()),
        path.rect(-0.075, -0.075, 0.15, 0.15),
        path.path(path.moveto(-0.1, 0), path.lineto(0, 0.1), path.lineto(0.1, 0),
                  path.lineto(0, -0.1), path.closepath()),
    ]
    c = canvas.canvas()
    step = height / 5.0
    for i, label in enumerate(['Unrelated', 'Pseudo Replicate', 'Same Tissue Type', 'Biological Replicate']):
        if i == 0:
            c.stroke(strokes[i], [ramp.get_rgb_color(i), trafo.translate(0, step * (i + 0.5))])
        else:
            c.stroke(strokes[i], [style.linewidth.thin, deco.filled([ramp.get_rgb_color(i)]),
                trafo.translate(0, step * (i + 0.5))])
        c.text(0.2, step * (i + 0.5) + 0.08, label, [text.halign.left, text.valign.top, text.size(text_size)])
    return c

def load_rdata(fname):
    data = []
    infile = open(fname)
    line = infile.readline()
    line = infile.readline()
    while line:
        temp = line.rstrip('\n').split('\t')
        data.append((temp[0], temp[1], temp[2], int(temp[3]), int(temp[4]) / 1000, float(temp[5])))
        line = infile.readline()
    data = numpy.array(data, dtype=numpy.dtype([('sample1', 'S50'), ('sample2', 'S50'), ('genome', 'S4'),
        ('coverage', numpy.int32), ('resolution', numpy.int32), ('score', numpy.float64)]))
    for i in range(data.shape[0]):
        if data['sample1'][i] > data['sample2'][i]:
            data['sample1'][i], data['sample2'][i] = data['sample2'][i], data['sample1'][i]
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
