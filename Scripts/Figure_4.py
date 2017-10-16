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

A_width = (page_width - 0.5) * 0.5 + 0.3
A_height = page_height
B_width = page_width - 0.5 - A_width

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_4.pdf" % base_dir
    c = canvas.canvas()
    data_fname = "%s/Results/quasar_replicate_results.txt" % base_dir
    data = load_rdata(data_fname)
    sizes = []
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        sizes.append(0)
        where = numpy.where(data['genome'] == gen)
        samples = numpy.unique(numpy.r_[data['sample1'][where], data['sample2'][where]])
        for samp in samples:
            if samp.count('Pseudo') == 0:
                sizes[-1] += 1
    B_height = A_height * sizes[1] / float(sizes[0])
    C_height = A_height * sizes[2] / float(sizes[0])
    key_height = 2.0
    D_height = A_height - B_height - C_height - key_height - 1.0
    c.insert(plot_A(A_width, A_height, data))
    c.insert(plot_B(B_width, B_height, data), [trafo.translate(A_width + 0.5, 0)])
    c.insert(plot_C(B_width, C_height, data), [trafo.translate(A_width + 0.5, -B_height - 0.5)])
    c.insert(plot_D(B_width, D_height, data), [trafo.translate(A_width + 0.5, -A_height)])
    c.insert(plot_key(key_height), [trafo.translate(A_width + 3.0, -A_height + D_height + 0.5)])
    c.writePDFfile(out_fname)

def plot_A(width, height, all_data):
    data = all_data[numpy.where(all_data['genome'] == 'hg38')]
    pheight = height - 0.3
    hoffset = 3.2
    c = plot_samples(data, hoffset, width, pheight)
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(0, -0.3)])
    c1.text(hoffset + (width - hoffset) * 0.5, 0, 'Human', [text.halign.center, text.valign.top, text.size(text_size)])
    c1.text(0, 0, 'A', [text.valign.top, text.halign.left, text.size(0)])
    return c1

def plot_B(width, height, all_data):
    data = all_data[numpy.where(all_data['genome'] == 'mm10')]
    pheight = height - 0.3
    hoffset = 2.6
    c = plot_samples(data, hoffset, width, pheight)
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(0, -0.3)])
    c1.text(hoffset + (width - hoffset) * 0.5, 0, 'Mouse', [text.halign.center, text.valign.top, text.size(text_size)])
    c1.text(0, 0, 'B', [text.valign.top, text.halign.left, text.size(0)])
    return c1

def plot_C(width, height, all_data):
    data = all_data[numpy.where(all_data['genome'] == 'dm6')]
    pheight = height - 0.3
    hoffset = 2.6
    c = plot_samples(data, hoffset, width, pheight)
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(0, -0.3)])
    c1.text(hoffset + (width - hoffset) * 0.5, 0, 'Fruit Fly',
        [text.halign.center, text.valign.top, text.size(text_size)])
    c1.text(0, 0, 'C', [text.valign.top, text.halign.left, text.size(0)])
    return c1

def plot_D(width, height, rdata):
    hoffset = 0.8
    voffset = 0.4
    pwidth = width - hoffset
    pheight = height - voffset
    qdata_fname = "%s/Results/quasar_quality_results.txt" % base_dir
    qdata = load_qdata(qdata_fname)
    Xs = [[], [], []]
    Ys = [[], [], []]
    reses = {'hg38': 1000000, 'mm10': 1000000, 'dm6': 100000}
    indices = {'hg38': 0, 'mm10': 1, 'dm6': 2}
    ramp = Ramp(0.0, 4.0)
    colors = {'hg38': ramp.get_rgb_color(1.0), 'mm10': ramp.get_rgb_color(2.0), 'dm6': ramp.get_rgb_color(3.0)}
    for sample in numpy.unique(qdata['sample']):
        if sample.count('Rep2') > 0:
            continue
        s2 = sample.replace('Rep1', 'Rep2')
        qwhere = numpy.where(qdata['sample'] == sample)[0]
        gen = qdata['genome'][qwhere[0]]
        res = reses[gen]
        rwhere = numpy.where((((rdata['sample1'] == sample) & (rdata['sample2'] == s2)) | 
                             ((rdata['sample2'] == sample) & (rdata['sample1'] == s2))) &
                             (rdata['coverage'] == 0) & (rdata['resolution'] == 0))[0]
        if rwhere.shape[0] == 0:
            continue
        index = indices[gen]
        score1 = qdata['score'][numpy.where((qdata['sample'] == sample) & (qdata['param'] == 0) &
                                            (qdata['resolution'] == res))][0]
        score2 = qdata['score'][numpy.where((qdata['sample'] == s2) & (qdata['param'] == 0) &
                                            (qdata['resolution'] == res))][0]
        Xs[index].append(min(score1,score2))
        Ys[index].append(rdata['score'][rwhere][0])
    all_Xs = Xs[0] + Xs[1] + Xs[2]
    all_Ys = Ys[0] + Ys[1] + Ys[2]
    all_Xs = numpy.array(all_Xs)
    all_Ys = numpy.array(all_Ys)
    minX = 0.0#numpy.amin(all_Xs)
    maxX = numpy.amax(all_Xs)
    minY = numpy.amin(all_Ys)
    maxY = numpy.amax(all_Ys)
    spanX = maxX - minX
    spanY = maxY - minY
    #minX -= spanX * 0.05
    maxX += spanX * 0.05
    minY -= spanY * 0.05
    maxY += spanY * 0.05
    spanX = maxX - minX
    spanY = maxY - minY
    c = canvas.canvas()
    c1 = canvas.canvas([canvas.clip(path.rect(0, 0, pwidth, pheight))])
    def f(x, a, b):
        return a * x + b
    for i in range(len(Xs)):
        Xs[i] = numpy.array(Xs[i])
        Ys[i] = numpy.array(Ys[i])
        Xs[i] = (Xs[i] - minX) / spanX * pwidth
        Ys[i] = (Ys[i] - minY) / spanY * pheight
        for j in range(Xs[i].shape[0]):
            c1.stroke(path.circle(Xs[i][j], Ys[i][j], 0.05),
                [ramp.get_rgb_color(i + 1.0)])
    c.insert(c1, [trafo.translate(hoffset, voffset)])
    c.stroke(path.rect(hoffset, voffset, pwidth, pheight))
    c.text(0.1, voffset + pheight * 0.5, "Replicate Score",
        [trafo.rotate(90), text.valign.top, text.halign.center, text.size(text_size)])
    c.text(hoffset + pwidth * 0.5, 0, "Min Quality Score",
        [text.valign.bottom, text.halign.center, text.size(text_size)])
    c.text(0, height, 'D', [text.halign.left, text.valign.top, text.size(0)])
    val = ceil(minX * 10.0) / 10.0
    X = (val - minX) / spanX * pwidth + hoffset
    c.stroke(path.line(X, voffset, X, voffset - 0.05))
    c.text(X, voffset - 0.1, "%0.1f" % val, [text.halign.center, text.valign.top, text.size(text_size)])
    val = floor(maxX * 10.0) / 10.0
    X = (val - minX) / spanX * pwidth + hoffset
    c.stroke(path.line(X, voffset, X, voffset - 0.05))
    c.text(X, voffset - 0.1, "%0.1f" % val, [text.halign.center, text.valign.top, text.size(text_size)])
    val = ceil(minY * 10.0) / 10.0
    Y = (val - minY) / spanY * pheight + voffset
    c.stroke(path.line(hoffset, Y, hoffset - 0.05, Y))
    c.text(hoffset - 0.1, Y, "%0.1f" % val, [text.halign.right, text.valign.middle, text.size(text_size)])
    val = floor(maxY * 10.0) / 10.0
    Y = (val - minY) / spanY * pheight + voffset
    c.stroke(path.line(hoffset, Y, hoffset - 0.05, Y))
    c.text(hoffset - 0.1, Y, "%0.1f" % val, [text.halign.right, text.valign.middle, text.size(text_size)])
    labels = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        label = labels[gen]
        pcolor = ramp.get_rgb_color(i + 1.0)
        X = hoffset + pwidth * 0.75
        Y = voffset + 0.3 * (3 - i)
        c.stroke(path.circle(X, Y, 0.05), [pcolor])
        c.text(X + 0.15, Y + 0.11, label, [text.halign.left, text.valign.top, text.size(text_size)])
    return c

def plot_samples(data, hoffset, width, height):
    pwidth = width - hoffset
    voffset = 0.6
    pheight = height - voffset
    samples = numpy.unique(numpy.r_[data['sample1'], data['sample2']])
    valid = numpy.zeros(samples.shape[0], dtype=bool)
    for i, samp in enumerate(samples):
        if samp.count('Pseudo') == 0:
            valid[i] = True
    samples = samples[valid]
    step = pheight / float(samples.shape[0])
    indices = {}
    for i, sample in enumerate(samples):
        indices[sample] = -(i + 0.5) * step
    c = canvas.canvas()
    maxX = numpy.amax(data['score'])
    minX = numpy.amin(data['score'])
    span = maxX - minX
    minX -= span * 0.05
    maxX += span * 0.05
    span = maxX - minX
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
            points[index].append([hoffset + val, Y + jitter])
        if name2 in indices:
            Y = indices[name2]
            points[index].append([hoffset + val, Y + jitter])
    scale = 10.0
    for i in range(int(ceil(minX * scale)), int(floor(maxX * scale)) + 1):
        X = (i / scale - minX) / span * pwidth + hoffset
        c.stroke(path.line(X, -pheight, X, -pheight - 0.05))
        c.text(X, -pheight - 0.1, "%0.1f" % (i / scale), [text.halign.center, text.valign.top, text.size(text_size)])
        c.stroke(path.line(X, -pheight, X, 0), [color.gray(0.9)])
    for i in range(len(points)):
        for j in range(len(points[i])):
            if i == 0:
                c.stroke(strokes[i],
                    [trafo.scale(0.6), trafo.translate(points[i][j][0], points[i][j][1]), ramp.get_rgb_color(i)])
            else:
                c.stroke(strokes[i], [trafo.scale(0.6), style.linewidth.thin, deco.filled([ramp.get_rgb_color(i)]),
                    trafo.translate(points[i][j][0], points[i][j][1])])
    for sample in samples:
        temp = sample.split('_')
        if len(temp) > 4:
            name = '%s %s %s' % (temp[0], temp[2], temp[4])
        else:
            name = '%s %s %s' % (temp[0], temp[2], temp[3])
        c.text(hoffset - 0.05, indices[sample], name,
            [text.halign.right, text.valign.middle, text.size(text_size - 1)])
    c.stroke(path.rect(hoffset, 0, pwidth, -pheight))
    for i in range(1, samples.shape[0]):
        c.stroke(path.line(hoffset, -step * i, hoffset + pwidth, -step * i), [style.linestyle.dotted])
    c.text(hoffset + pwidth * 0.5, -height, "Replicate Score",
        [text.halign.center, text.valign.bottom, text.size(text_size)])
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
    where = numpy.where(data['score'] == 1)
    remove = []
    for i in where[0]:
        if data['sample1'][i] not in remove and data['sample2'][i] not in remove:
            remove.append(data['sample1'][i])
    for sample in remove:
        mask = (data['sample1'] != sample)
        mask2 = (data['sample2'] != sample)
        data = data[numpy.where((data['sample1'] != sample) & (data['sample2'] != sample))]
    temp = []
    for i in range(data.shape[0]):
        temp.append((data['sample1'][i], data['sample2'][i], data['coverage'][i], data['genome'][i]))
    temp = list(set(temp))
    new_data = numpy.zeros(len(temp), dtype=numpy.dtype([('sample1', 'S50'), ('sample2', 'S50'), ('genome', 'S4'),
        ('coverage', numpy.int32), ('resolution', numpy.int32), ('score', numpy.float64)]))
    for i in range(new_data.shape[0]):
        new_data['sample1'][i] = temp[i][0]
        new_data['sample2'][i] = temp[i][1]
        new_data['coverage'][i] = temp[i][2]
        new_data['genome'][i] = temp[i][3]
        new_data['score'][i] = numpy.amax(data['score'][numpy.where(
            (data['sample1'] == temp[i][0]) & (data['sample2'] == temp[i][1]) &
            (data['coverage'] == temp[i][2]))])
    return new_data

def load_qdata(fname):
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
