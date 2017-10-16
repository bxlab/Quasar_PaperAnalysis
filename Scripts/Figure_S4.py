#!/usr/bin/env python

import sys
import os
from math import ceil, floor

import numpy
import h5py
from pyx import *
from scipy.stats import linregress
from scipy.optimize import curve_fit
import hifive
from PIL import Image

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{subscript}")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\usepackage{amsmath}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")
text.preamble(r"\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}")

page_width = 16.8
page_height = 22.4
A_width = (page_width - 0.8) / 4.0
A_height = 19 / 3.0
C_width = A_width
B_width = A_width * 2
C_height = C_width - 0.2
B_height = C_height
A_height = C_height
text_size = -3

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_S4.pdf" % base_dir
    c = canvas.canvas()
    data_fname = "%s/Results/quasar_replicate_results.txt" % base_dir
    data = load_rdata(data_fname)
    qdata_fname = "%s/Results/quasar_quality_results.txt" % base_dir
    qdata = load_qdata(qdata_fname)
    c = plot_all(data, qdata)
    c.writePDFfile(out_fname)

def plot_all(rdata, qdata):
    hoffset = 0.5
    voffset = 0.5
    pwidth = (A_width - 0.4 - hoffset) / 3.0
    pheight = pwidth
    height = pheight + voffset + 0.2
    resolutions = {}
    sizes = numpy.zeros(4, dtype=numpy.float64)
    all_samples = {}
    genomes = ['hg38', 'mm10', 'dm6']
    total = 0
    for i, gen in enumerate(genomes):
        all_samples[gen] = []
        sizes[i + 1] = sizes[i]
        where = numpy.where(rdata['genome'] == gen)
        resolutions[gen] = numpy.unique(rdata['resolution'][numpy.where(rdata['genome'] == gen)])
        samples = numpy.unique(numpy.r_[rdata['sample1'][where], rdata['sample2'][where]])
        for samp in samples:
            if samp.count('Pseudo') == 0:
                sizes[i + 1] += 1
                all_samples[gen].append(samp)
    rep_data = []
    for i, gen in enumerate(genomes):
        for j in range(resolutions[gen].shape[0]):
            where = numpy.where((rdata['genome'] == gen) & (rdata['resolution'] == resolutions[gen][j]))[0]
            where2 = numpy.where((qdata['genome'] == gen) & (qdata['resolution'] == resolutions[gen][j]))[0]
            rep_data += compile_data(rdata[where], qdata[where2], resolutions[gen][j], gen, all_samples[gen])
    data = numpy.array(rep_data, dtype=numpy.dtype([('sample', 'S80'), ('genome', 'S10'), ('resolution', numpy.int32), ('rscore', numpy.float64), ('qscore1', numpy.float64), ('qscore2', numpy.float64), ('pass', numpy.bool)]))
    data = data[numpy.argsort(data['sample'])]
    results = {}
    all_results = []
    for i, cutoff0 in enumerate(numpy.linspace(0.75, 0.99, 25)):
        rcutoff = cutoff0
        old_rcutoff = rcutoff + 1
        pos = 0
        while rcutoff != old_rcutoff:
            old_rcutoff = rcutoff
            qcutoff, qscore = find_qcutoff(data, rcutoff)
            rcutoff, rscore = find_rcutoff(data, qcutoff)
        results[(rcutoff, qcutoff)] = (qscore + rscore, qscore, rscore)
        all_results.append((cutoff0, qscore + rscore, qcutoff, rcutoff))
    results2 = []
    for key, score in results.iteritems():
        results2.append((key[0], key[1], score[0], score[1], score[2]))
    results2 = numpy.array(results2)
    all_results = numpy.array(all_results)
    where = numpy.where(results2[:, 2] == numpy.amin(results2[:, 2]))[0][0]
    best_cutoffs = numpy.copy(results2[where, :2])
    results2[where, :] = numpy.inf
    where = numpy.where(results2[:, 2] == numpy.amin(results2[:, 2]))[0][0]
    good_cutoffs = numpy.copy(results2[where, :2])
    c = canvas.canvas()
    c.insert(plot_A(all_results, A_width, A_height), [trafo.translate(0, 0)])
    c.insert(plot_B(data, good_cutoffs, B_width, B_height), [trafo.translate(A_width + 0.4, 0)])
    c.insert(plot_C(data, good_cutoffs, C_width, C_height), [trafo.translate(A_width + B_width + 0.8, 0)])
    return c

def plot_A(data, width, height):
    hoffset = 0.6
    hoffset2 = 0.85
    voffset = 0.4
    pwidth = width - hoffset - hoffset2
    pheight = height - voffset
    c = canvas.canvas()
    c1 = canvas.canvas()
    c.insert(c1, [trafo.translate(hoffset, voffset)])
    data = data[numpy.argsort(data[:, 0]), :]
    Xs = numpy.round(data[:, 0].astype(numpy.float64) * 100.0) / 100.0
    minX = Xs[0]
    maxX = Xs[-1]
    Y1s = data[:, 2]
    Y2s = data[:, 3]
    Y3s = data[:, 1]
    minY1 = min(numpy.amin(Y1s), numpy.amin(Y2s))
    maxY1 = max(numpy.amax(Y1s), numpy.amax(Y2s))
    minY3 = numpy.amin(Y3s)
    maxY3 = numpy.amax(Y3s)
    span1 = maxY1 - minY1
    minY1 = max(0, minY1 - 0.05 * span1)
    maxY1 += 0.05 * span1
    span3 = maxY3 - minY3
    minY3 -= 0.05 * span3
    maxY3 += 0.05 * span3
    Xs = (Xs - minX) / (maxX - minX) * pwidth
    Y1s = (Y1s - minY1) / (maxY1 - minY1) * pheight
    Y2s = (Y2s - minY1) / (maxY1 - minY1) * pheight
    Y3s = (Y3s - minY3) / (maxY3 - minY3) * pheight
    ramp = Ramp(0., 3.)
    for i, Ys in enumerate([Y1s, Y2s, Y3s]):
        lpath = path.path(path.moveto(Xs[0], Ys[0]))
        for j in range(1, Xs.shape[0]):
            lpath.append(path.lineto(Xs[j], Ys[j]))
        c1.stroke(lpath, [ramp.get_rgb_color(i)])
    c1.stroke(path.rect(0, 0, pwidth, pheight))
    c.text(0, height, 'A', [text.halign.left, text.valign.top, text.size(0)])
    for i in [minX, maxX]:
        val = i
        X = (val - minX) / (maxX - minX) * pwidth
        c1.stroke(path.line(X, 0, X, -0.05))
        c1.text(X, -0.1, "%0.2f" % val, [text.valign.top, text.halign.center, text.size(text_size)])
    for i in [ceil(minY1 * 100) / 100, 0.85]:
        val = i
        Y = (val - minY1) / (maxY1 - minY1) * pheight
        c1.stroke(path.line(0, Y, -0.05, Y))
        c1.text(-0.1, Y, "%0.2f" % val, [text.valign.middle, text.halign.right, text.size(text_size)])
    for i in [ceil(minY3 * 100) / 100, floor(maxY3 * 100) / 100]:
        val = i
        Y = (val - minY3) / (maxY3 - minY3) * pheight
        c1.stroke(path.line(pwidth, Y, pwidth + 0.05, Y))
        c1.text(pwidth + 0.1, Y, "%0.2f" % val, [text.valign.middle, text.halign.left, text.size(text_size)])
    c1.text(pwidth * 0.5, -voffset, r"Initial Cutoff\textsubscript{r}",
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    c1.text(-hoffset, pheight * 0.5, r"Optimized Cutoff",
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c1.text(pwidth + hoffset2, pheight * 0.5, r"Gini Impurity Score",
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(-90)])
    c2 = canvas.canvas()
    for i, label in enumerate([r'cutoff\textsubscript{q}', r'cutoff\textsubscript{r}', ['combined', 'score']]):
        if i == 2:
            Y = -0.75
        else:
            Y = -0.3 * i
        c2.stroke(path.line(0, Y, 0.3, Y), [ramp.get_rgb_color(i)])
        if i == 2:
            c2.text(0.35, Y + 0.25, label[0], [text.halign.left, text.valign.top, text.size(text_size)])
            c2.text(0.35, Y - 0.03, label[1], [text.halign.left, text.valign.top, text.size(text_size)])
        else:
            c2.text(0.35, Y + 0.12, label, [text.halign.left, text.valign.top, text.size(text_size)])
    c1.insert(c2, [trafo.translate(0.1, pheight * 0.75)])
    return c

def plot_B(data, cutoffs, width, height):
    hoffset = 0.4
    voffset = 0.5
    spacer = 0.5
    pwidth = (width - spacer - hoffset) / 2.0
    pheight = height - voffset - 0.3
    c = canvas.canvas()
    labels1 = data['qscore1'] >= cutoffs[1]
    labels2 = data['qscore2'] >= cutoffs[1]
    where = numpy.where(labels1 == labels2)[0]
    c1, minX, maxX = plot_histograms(data['rscore'][where], labels1[where], cutoffs[0], pwidth, pheight)
    val = ceil(minX * 10.0) / 10.0
    X = (val - minX) / (maxX - minX) * pwidth
    c1.stroke(path.line(X, 0, X, -0.05))
    c1.text(X, -0.1, "%0.1f" % val, [text.valign.top, text.halign.center, text.size(text_size)])
    val = floor(maxX * 10.0) / 10.0
    X = (val - minX) / (maxX - minX) * pwidth
    c1.stroke(path.line(X, 0, X, -0.05))
    c1.text(X, -0.1, "%0.1f" % val, [text.valign.top, text.halign.center, text.size(text_size)])
    c.insert(c1, [trafo.translate(hoffset, voffset)])
    labels = numpy.r_[data['rscore'] > cutoffs[0], data['rscore'] > cutoffs[0]]
    c2, minX, maxX = plot_histograms(numpy.r_[data['qscore1'], data['qscore2']], labels, cutoffs[1], pwidth, pheight)
    val = ceil(minX * 10.0) / 10.0
    X = (val - minX) / (maxX - minX) * pwidth
    c2.stroke(path.line(X, 0, X, -0.05))
    c2.text(X, -0.1, "%0.1f" % val, [text.valign.top, text.halign.center, text.size(text_size)])
    val = floor(maxX * 10.0) / 10.0
    X = (val - minX) / (maxX - minX) * pwidth
    c2.stroke(path.line(X, 0, X, -0.05))
    c2.text(X, -0.1, "%0.1f" % val, [text.valign.top, text.halign.center, text.size(text_size)])
    c.insert(c2, [trafo.translate(hoffset + pwidth + spacer, voffset)])
    c1.text(pwidth * 0.5, pheight + 0.25, "Replicate",
        [text.halign.center, text.valign.top, text.size(text_size)])
    c2.text(pwidth * 0.5, pheight + 0.25, "Quality",
        [text.halign.center, text.valign.top, text.size(text_size)])
    c.text(hoffset + pwidth + spacer * 0.5, 0, "Score",
         [text.halign.center, text.valign.bottom, text.size(text_size)])
    c.text(0, voffset + pheight * 0.5, 'Density',
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(0, height + 0.05, 'B', [text.halign.left, text.valign.top, text.size(0)])
    c1.stroke(path.line(0.9, pheight - 0.4, 1.2, pheight - 0.4))
    c1.text(1.3, pheight - 0.28, r'< cutoff\textsubscript{q}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    c1.stroke(path.line(0.9, pheight - 0.7, 1.2, pheight - 0.7), [color.rgb.red])
    c1.text(1.3, pheight - 0.58, r'> cutoff\textsubscript{q}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    c2.stroke(path.line(0.9, pheight - 0.4, 1.2, pheight - 0.4))
    c2.text(1.3, pheight - 0.28, r'< cutoff\textsubscript{r}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    c2.stroke(path.line(0.9, pheight - 0.7, 1.2, pheight - 0.7), [color.rgb.red])
    c2.text(1.3, pheight - 0.58, r'> cutoff\textsubscript{r}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    return c

def plot_C(data, cutoffs, width, height):
    hoffset = 0.8
    voffset = 0.6
    pwidth = width - hoffset
    pheight = height - voffset
    Xs = []
    Ys = []
    indices = []
    samples = numpy.unique(data['sample'])
    for i, sample in enumerate(samples):
        where = numpy.where(data['sample'] == sample)[0]
        reses = numpy.log10(data['resolution'][where])
        X = find_intercept(data['rscore'][where], reses, cutoffs[0])
        if X is None:
            continue
        Y1 = find_intercept(data['qscore1'][where], reses, cutoffs[1])
        if Y1 is not None:
            Xs.append(X)
            Ys.append(Y1)
            indices.append(i)
        Y2 = find_intercept(data['qscore2'][where], reses, cutoffs[1])
        if Y2 is not None:
            Xs.append(X)
            Ys.append(Y2)
            indices.append(i)
        #print sample, Y1, Y2
    Xs = numpy.array(Xs)
    Ys = numpy.array(Ys)
    minX = min(numpy.amin(Xs), numpy.amin(Ys))
    maxX = max(numpy.amax(Xs), numpy.amax(Ys))
    span = maxX - minX
    minX -= 0.05 * span
    maxX += 0.05 * span
    span = maxX - minX
    c = canvas.canvas()
    c1 = canvas.canvas([canvas.clip(path.rect(0, 0, pwidth, pheight))])
    c.insert(c1, [trafo.translate(hoffset, voffset)])
    ramp = Ramp(0.0, len(samples) - 1.0)
    A, B, R, P = linregress(Xs, Ys)[:4]
    c1.stroke(path.line(0, (minX * A + B - minX) / span * pheight, pwidth, (maxX * A + B - minX) / span * pheight))
    for i in range(Xs.shape[0]):
        X = (Xs[i] - minX) / span * pwidth
        Y = (Ys[i] - minX) / span * pheight
        c1.stroke(path.circle(X, Y, 0.04), [ramp.get_rgb_color(indices[i])])
    c.stroke(path.rect(hoffset, voffset, pwidth, pheight))
    for i in range(int(ceil(minX)), int(floor(maxX)) + 1):
        X = (i - minX) / span * pwidth + hoffset
        Y = (i - minX) / span * pheight + voffset
        c.stroke(path.line(X, voffset, X, voffset - 0.085))
        c.stroke(path.line(hoffset, Y, hoffset - 0.085, Y))
        c.text(X, voffset - 0.1, '1e%i' % i, [text.valign.top, text.halign.center, text.size(text_size)])
        c.text(hoffset - 0.1, Y, '1e%i' % i, [text.valign.middle, text.halign.right, text.size(text_size)])
    for i in range(int(ceil(minX)) - 1, int(floor(maxX)) + 2):
        for j in numpy.linspace(10 ** i, 10 ** (i + 1), 9)[1:-1]:
            X = (numpy.log10(j) - minX) / span * pwidth
            Y = (numpy.log10(j) - minX) / span * pheight
            if X > 0 and X < pwidth:
                c.stroke(path.line(X + hoffset, voffset, X + hoffset, voffset - 0.05))
            if Y > 0 and Y < pheight:
                c.stroke(path.line(hoffset, Y + voffset, hoffset - 0.05, Y + voffset))
    c.text(0, voffset + pheight * 0.5, r"Resolution\textsubscript{q}",
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + pwidth * 0.5, 0.25, r"Resolution\textsubscript{r}",
        [text.halign.center, text.valign.top, text.size(text_size)])
    c.text(0, height + 0.05, 'C', [text.halign.left, text.valign.top, text.size(0)])
    c1.text(0.05, pheight - 0.05, "R=%0.2f" % R, [text.halign.left, text.valign.top, text.size(text_size)])
    c1.text(0.05, pheight - 0.3, "p=%0.1e" % P, [text.halign.left, text.valign.top, text.size(text_size)])
    return c

def compile_data(rdata, qdata, res, gen, samples):
    results = []
    for i in range(rdata.shape[0]):
        name1 = rdata['sample1'][i]
        name2 = rdata['sample2'][i]
        t1 = name1.split('_')[2].split('H1-')[-1].split('-Dil')[0].split('-IS')[0]
        t2 = name2.split('_')[2].split('H1-')[-1].split('-Dil')[0].split('-IS')[0]
        if name1.split('Rep')[0] == name2.split('Rep')[0]:
            samp = name1.split('_Rep')[0]
            w1 = numpy.where(qdata['sample'] == '%s_Rep1' % (samp))[0][0]
            w2 = numpy.where(qdata['sample'] == '%s_Rep2' % (samp))[0][0]
            results.append((samp, gen, res, rdata['score'][i], qdata['score'][w1], qdata['score'][w2], True))
    return results

def plot_histograms(scores, labels, cutoff, width, height):
    minX = 0.0#numpy.amin(scores)
    maxX = numpy.amax(scores)
    where = numpy.where(labels)[0]
    hist1, bins = numpy.histogram(scores[where], range=(minX, maxX), density=True, bins=15)
    mids = (bins[:-1] + bins[1:]) / 2.0
    where = numpy.where(labels == False)[0]
    hist2 = numpy.histogram(scores[where], range=(minX, maxX), density=True, bins=15)[0]
    minY = 0.0
    maxY = max(numpy.amax(hist1), numpy.amax(hist2))
    maxY += (maxY - minY) * 0.05
    c = canvas.canvas()
    c1 = canvas.canvas()
    c.insert(c1)
    step = width / hist1.shape[0]
    step1 = step / 3.0
    Xs = (mids - bins[0]) / (bins[-1] - bins[0]) * width
    Ys1 = (hist1 - minY) / (maxY - minY) * height
    Ys2 = (hist2 - minY) / (maxY - minY) * height
    lpath1 = path.path(path.moveto(Xs[0], Ys1[0]))
    lpath2 = path.path(path.moveto(Xs[0], Ys2[0]))
    for i in range(1, hist1.shape[0]):
        lpath1.append(path.lineto(Xs[i], Ys1[i]))
        lpath2.append(path.lineto(Xs[i], Ys2[i]))
    c1.stroke(lpath1, [color.rgb.red])
    c1.stroke(lpath2, [color.rgb.black])
    X = (cutoff - minX) / (maxX - minX) * width
    c1.stroke(path.line(X, 0, X, height), [style.linestyle.dotted])
    c1.stroke(path.rect(0, 0, width, height))
    for h, i in enumerate([ceil(minY), floor(maxY)]):
        pos = (i - minY) / (maxY - minY) * height
        c.stroke(path.line(0, pos, -0.05, pos))
        c.text(-0.1, pos, '%i' % int(i), [text.halign.right, text.valign.middle, text.size(text_size)])
    return c, bins[0], bins[-1]

def find_qcutoff(data, rcutoff):
    labels = data['rscore'] >= rcutoff
    return find_cutoff(numpy.r_[data['qscore1'], data['qscore2']], numpy.r_[labels, labels])

def find_rcutoff(data, qcutoff):
    labels1 = data['qscore1'] >= qcutoff
    labels2 = data['qscore2'] >= qcutoff
    #where = numpy.where(labels1 == labels2)[0]
    #return find_cutoff(data['rscore'][where], labels1[where])
    return find_cutoff(data['rscore'], (labels1 & labels2))

def find_cutoff(data, labels):
    temp =  numpy.copy(data)
    temp.sort()
    mids = (temp[1:] + temp[:-1]) / 2.0
    scores = numpy.zeros(mids.shape[0], dtype=numpy.float64)
    for i in range(mids.shape[0]):
        where = numpy.where(data <= mids[i])[0]
        p0 = numpy.where(labels[where])[0].shape[0] / float(where.shape[0])
        p1 = 1.0 - p0
        scores[i] += p0 * p1 * where.shape[0]
        where = numpy.where(data > mids[i])[0]
        p0 = numpy.where(labels[where])[0].shape[0] / float(where.shape[0])
        p1 = 1.0 - p0
        scores[i] += p0 * p1 * where.shape[0]
    scores /= mids.shape[0]
    return mids[numpy.where(scores == numpy.amin(scores))[0][0]], numpy.amin(scores)

def find_intercept(scores, reses, cutoff):
    order = numpy.argsort(reses)
    Xs = numpy.copy(reses)[order]
    Ys = numpy.copy(scores)[order]
    pos = 0
    while pos < Xs.shape[0] - 1 and cutoff > Ys[pos + 1]:
        pos += 1
    if pos == Xs.shape[0] - 1 or cutoff < Ys[pos]:
        return None
    A = (Ys[pos + 1] - Ys[pos]) / (Xs[pos + 1] - Xs[pos])
    B = Ys[pos] - Xs[pos] * A
    return (cutoff - B) / A

def plot_hm(data, maxscore):
    ramp = Ramp(0, 1.0)
    score = data[:, :, 1] / maxscore
    where = numpy.where(data[:, :, 0] > 0)
    hm = numpy.zeros(score.shape, dtype=numpy.uint32)
    hm.fill(int('FF888888', 16))
    hm[where] = ramp.get_color(score[where])
    return Image.frombuffer('RGBA', (hm.shape[1], hm.shape[0]), hm, 'raw', 'RGBA', 0, 1)

def load_rdata(fname):
    data = []
    infile = open(fname)
    line = infile.readline()
    line = infile.readline()
    while line:
        temp = line.rstrip('\n').split('\t')
        data.append((temp[0], temp[1], temp[2], int(temp[3]), int(temp[4]), float(temp[5])))
        line = infile.readline()
    data = numpy.array(data, dtype=numpy.dtype([('sample1', 'S50'), ('sample2', 'S50'), ('genome', 'S4'),
        ('coverage', numpy.int32), ('resolution', numpy.int32), ('score', numpy.float64)]))
    for i in range(data.shape[0]):
        if data['sample1'][i] > data['sample2'][i]:
            data['sample1'][i], data['sample2'][i] = data['sample2'][i], data['sample1'][i]
    return data

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
        self.scale = numpy.array([1, 256, 256 ** 2], dtype=numpy.uint32).reshape(1, -1)

    def get_color(self, val):
        if isinstance(val, int):
            index = numpy.searchsorted(self.breaks[1:-1], val)
            frac = max(0.0, min(1.0, (val - self.breaks[index]) / self.span))
            cvals = numpy.round(self.colors[index] * (1.0 - frac) + self.colors[index + 1] * frac).astype(numpy.uint32)
            return self.base + numpy.sum(cvals * self.scale)
        else:
            index = numpy.searchsorted(self.breaks[1:-1], val)
            frac = numpy.maximum(0.0, numpy.minimum(1.0, (val - self.breaks[index]) / self.span)).reshape(-1, 1)
            cvals = numpy.round(self.colors[index, :].reshape(-1, 3) * (1.0 - frac) +
                                self.colors[index + 1, :].reshape(-1, 3) * frac).astype(numpy.uint32)
            return self.base + numpy.sum(cvals * self.scale.reshape(1, -1), axis=1)

    def get_rgb_color(self, val):
        index = numpy.searchsorted(self.breaks[1:-1], val)
        frac = max(0.0, min(1.0, (val - self.breaks[index]) / self.span))
        cvals = (self.colors[index] * (1.0 - frac) + self.colors[index + 1] * frac) / 255.0
        return color.rgb(r=cvals[0], g=cvals[1], b=cvals[2])

if __name__ == "__main__":
    main()
