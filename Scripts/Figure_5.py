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
A_width = (page_width - 0.3) * 0.5
A_height = 19 / 3.0
C_width = (A_width - 0.4) / 3.0
B_width = A_width - C_width - 0.4
C_height = C_width - 0.2
B_height = C_height
D_width = A_width
D_height = A_height + B_height + 0.25
text_size = -3
hoffset = 3.2
voffset = 0.5
voffset2 = 0.25

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_5.pdf" % base_dir
    c = canvas.canvas()
    data_fname = "%s/Results/quasar_replicate_results.txt" % base_dir
    data = load_rdata(data_fname)
    qdata_fname = "%s/Results/quasar_quality_results.txt" % base_dir
    qdata = load_qdata(qdata_fname)
    c = canvas.canvas()
    c.insert(plot_A(data, qdata))
    c.writePDFfile(out_fname)

def plot_A(rdata, qdata):
    hoffset = 0.85
    voffset = 1.05
    voffset2 = 0.3
    hspacer = 0.3
    vspacer = 0.2
    pwidth = (A_width - hspacer * 2 - hoffset) / 3.0
    pheight = (A_height - voffset2 - voffset - vspacer) / 2.0
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
    for i, cutoff0 in enumerate(numpy.linspace(0.75, 0.99, 50)):
        rcutoff = cutoff0
        old_rcutoff = rcutoff + 1
        pos = 0
        while rcutoff != old_rcutoff:
            old_rcutoff = rcutoff
            qcutoff, qscore = find_qcutoff(data, rcutoff)
            rcutoff, rscore = find_rcutoff(data, qcutoff)
        results[(rcutoff, qcutoff)] = (qscore + rscore, qscore, rscore)
    results2 = []
    for key, score in results.iteritems():
        results2.append((key[0], key[1], score[0], score[1], score[2]))
    results2 = numpy.array(results2)
    where = numpy.where(results2[:, 2] == numpy.amin(results2[:, 2]))[0][0]
    best_cutoffs = numpy.copy(results2[where, :2])
    results2[where, :] = numpy.inf
    where = numpy.where(results2[:, 2] == numpy.amin(results2[:, 2]))[0][0]
    good_cutoffs = numpy.copy(results2[where, :2])
    print best_cutoffs
    print good_cutoffs
    c = canvas.canvas()
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        c1, c2 = plot_A_gen(data, gen, pwidth, pheight, best_cutoffs, good_cutoffs)
        c.insert(c1, [trafo.translate(hoffset + (pwidth + hspacer) * i, pheight + vspacer + voffset)])
        c.insert(c2, [trafo.translate(hoffset + (pwidth + hspacer) * i, voffset)])
    c.text(0, voffset + vspacer + pheight * 1.5, 'Replicate Score',
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(0, voffset + pheight * 0.5, 'Quality Score',
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + pwidth * 1.5 + hspacer, 0, 'Resolution',
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    c.text(0, A_height, 'A', [text.halign.left, text.valign.top, text.size(0)])
    c1 = canvas.canvas()
    c1.insert(c, [trafo.translate(0, B_height + 0.25)])
    c1.insert(plot_B(data, best_cutoffs, B_width, B_height))
    c1.insert(plot_C(data, best_cutoffs, C_width, C_height), [trafo.translate(B_width + 0.4, 0)])
    c1.insert(plot_D(D_width, D_height), [trafo.translate(A_width + 0.4, 0)])
    return c1

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
    c1.stroke(path.line(0.5, pheight - 0.4, 0.8, pheight - 0.4))
    c1.text(0.9, pheight - 0.28, r'< cutoff\textsubscript{q}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    c1.stroke(path.line(0.5, pheight - 0.7, 0.8, pheight - 0.7), [color.rgb.red])
    c1.text(0.9, pheight - 0.58, r'> cutoff\textsubscript{q}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    c2.stroke(path.line(0.5, pheight - 0.4, 0.8, pheight - 0.4))
    c2.text(0.9, pheight - 0.28, r'< cutoff\textsubscript{r}',
        [text.halign.left, text.valign.top, text.size(text_size)])
    c2.stroke(path.line(0.5, pheight - 0.7, 0.8, pheight - 0.7), [color.rgb.red])
    c2.text(0.9, pheight - 0.58, r'> cutoff\textsubscript{r}',
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
    ramp = Ramp(0.0, (len(samples) - 1.0) * 1.15)
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

def plot_D(width, height):
    fnames = [
        '%s/Data/GSE52457_hg38_H1-MesenchymalStem_HindIII_Rep2_bin.hcp' % base_dir,
        '%s/Data/GSE52457_hg38_H1-NPC_HindIII_Rep2_bin.hcp' % base_dir,
        '%s/Data/Encode_hg38_LNCaP-clone-FGC_HindIII_Rep2_bin.hcp' % base_dir,
    ]
    hic = []
    for fname in fnames:
        hic.append(hifive.HiC(fname, silent=True))
    hspacer = 0.1
    vspacer = 0.1
    hoffset = 0.5
    voffset = 0.55
    voffset2 = 0.5
    pwidth = (width - hoffset - 2 * hspacer) / 3.0
    pheight = (height - voffset - voffset2 - 2 * vspacer) / 3.0
    pwidth = min(pwidth, pheight)
    pheight = pwidth
    hspacer = (width - hoffset - pwidth * 3) / 2.0
    vspacer = (height - voffset - voffset2 - pheight * 3) / 2.0
    c = canvas.canvas()
    reses = [10000, 100000, 1000000]
    ares = [4.01543392754, 4.99605762231, 5.9148410993]
    mid = 75000000
    key_array = numpy.zeros((1000, 5), dtype=numpy.uint32)
    key_array.shape = (key_array.shape[1], key_array.shape[0])
    ramp = Ramp(0., 999.)
    key_array[:, :] = ramp.get_color(numpy.arange(1000)).reshape(1, -1)
    key_bm = Image.frombuffer('RGBA', (key_array.shape[1], key_array.shape[0]), key_array, 'raw', 'RGBA', 0, 1)
    key = canvas.canvas()
    key.insert(bitmap.bitmap(0, 0, key_bm, width=(pwidth - 0.7), height=0.2))
    key.stroke(path.rect(0, 0, pwidth - 0.7, 0.2))
    for i, res1 in enumerate(reses):
        data = []
        span = res1 * 200
        start = max(0, mid - span/2)
        stop = min(150000000, mid + span/2)
        maxscore = 0.0
        for j, res2 in enumerate(reses):
            data.append(hic[j].cis_heatmap(chrom='1', start=start, stop=stop, binsize=res1,
                                           datatype='raw', arraytype='full'))
            where = numpy.where(numpy.amin(data[-1], axis=2) == 0)
            data[-1][where[0], where[1], 1] = 0
            where = numpy.where(numpy.amin(data[-1], axis=2) > 0)
            data[-1][where[0], where[1], 1] = numpy.log10(data[-1][where[0], where[1], 0])
            maxscore = max(maxscore, numpy.amax(data[-1][:, :, 1]))
        maxlabel = int(ceil(10**maxscore))
        n = len(str(maxlabel)) - 1
        maxlabel = int(ceil(maxlabel / 10.0 ** n) * 10 ** n)
        maxscore = numpy.log10(maxlabel)
        for j, res2 in enumerate(reses):
            c.insert(bitmap.bitmap(0, 0, plot_hm(data[j], maxscore), width=pwidth),
                [trafo.translate(hoffset + (pwidth + hspacer) * i, voffset + (pheight + vspacer) * (2 - j))])
            if i == 0:
                eres = str(int(10 ** ares[j]))
                label = []
                while len(eres) > 3:
                    label = [eres[-3:]] + label
                    eres = eres[:-3]
                label = [eres] + label
                c.text(0.25, voffset + pwidth * (2.5 - j) + (2 - j) * vspacer, '%s bp' % ','.join(label),
                [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
        if res1 >= 1000000:
            label = "%iMb" % (res1 / 1000000)
        else:
            label = "%iKb" % (res1 / 1000)
        c.text(hoffset + pwidth * (i + 0.5) + hspacer * i, height - 0.25, label,
            [text.halign.center, text.valign.top, text.size(text_size)])
        c.insert(key, [trafo.translate(hoffset + (pwidth + hspacer) * i + 0.2, 0.3)])
        c.text(hoffset + (pwidth + hspacer) * i + 0.15, 0.4, '1',
            [text.halign.right, text.valign.middle, text.size(text_size)])
        label = '%0.0e' % maxlabel
        c.text(hoffset + (pwidth + hspacer) * i + pwidth - 0.45, 0.4, label.replace('e+0','e'),
            [text.halign.left, text.valign.middle, text.size(text_size)])
    c.text(0, height, 'D', [text.halign.left, text.valign.top, text.size(0)])
    c.text(hoffset + pwidth * 1.5 + hspacer, height, "Resolution",
        [text.halign.center, text.valign.top, text.size(text_size)])
    c.text(0, voffset + 1.5 * pheight + vspacer, "Sample Resolution at Quality Cutoff",
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + pwidth * 1.5 + hspacer, 0, "Reads per bin",
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    c.fill(path.rect(hoffset + 0.3, 0.05, 0.2, 0.2), [color.gray(0.5)])
    c.text(hoffset + 0.55, 0.05, 'No reads', [text.halign.left, text.valign.bottom, text.size(text_size)])
    return c

def plot_A_gen(data, gen, width, height, best, good):
    temp = numpy.copy(data[numpy.where(data['genome'] == gen)])
    resolutions = numpy.unique(temp['resolution'])
    minX = numpy.amin(numpy.log10(temp['resolution']))
    maxX = numpy.amax(numpy.log10(temp['resolution']))
    samples = numpy.unique(temp['sample'])
    ramp = Ramp(0.0, (samples.shape[0] - 1.0) * 1.15)
    minY = 0.0
    maxY = 1.09
    c0 = canvas.canvas()
    c = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])
    for i in numpy.linspace(minY, floor(maxY), 6):
        pos = (i - minY) / (maxY - minY) * height
        if gen == 'hg38':
            c0.stroke(path.line(0, pos, -0.05, pos))
            c0.text(-0.1, pos, '%0.1f' % i, [text.halign.right, text.valign.middle, text.size(text_size)])
        c0.stroke(path.line(0, pos, width, pos), [color.gray(0.9)])
    for i, sample in enumerate(samples):
        where = numpy.where(temp['sample'] == sample)[0]
        Xs = numpy.log10(temp['resolution'][where])
        Ys = temp['rscore'][where]
        order = numpy.argsort(Xs)
        Xs = (Xs[order] - minX) / (maxX - minX) * width
        Ys = (Ys[order] - minY) / (maxY - minY) * height
        lpath = path.path(path.moveto(Xs[0], Ys[0]))
        for j in range(1, Xs.shape[0]):
            lpath.append(path.lineto(Xs[j], Ys[j]))
        c.stroke(lpath, [ramp.get_rgb_color(i)])
    Y = (best[0] - minY) / (maxY - minY) * height
    c.stroke(path.line(0, Y, width, Y), [style.linestyle.dashed])
    Y = (good[0] - minY) / (maxY - minY) * height
    c.stroke(path.line(0, Y, width, Y), [style.linestyle.dotted])
    gen_label = {'hg38': 'Human', 'mm10': 'Mouse', 'dm6': 'Fruit Fly'}
    c0.text(width * 0.5, height + 0.275, gen_label[gen], [text.valign.top, text.halign.center, text.size(text_size)])
    c2 = canvas.canvas()
    c1 = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])
    minY = 0.0
    maxY = max(numpy.amax(data['qscore1']), numpy.amax(data['qscore2']))
    for i, sample in enumerate(samples):
        where = numpy.where(temp['sample'] == sample)[0]
        Xs = numpy.log10(temp['resolution'][where])
        Y1s = temp['qscore1'][where]
        Y2s = temp['qscore2'][where]
        order = numpy.argsort(Xs)
        Xs = (Xs[order] - minX) / (maxX - minX) * width
        Y1s = (Y1s[order] - minY) / (maxY - minY) * height
        Y2s = (Y2s[order] - minY) / (maxY - minY) * height
        lpath1 = path.path(path.moveto(Xs[0], Y1s[0]))
        lpath2 = path.path(path.moveto(Xs[0], Y2s[0]))
        for j in range(1, Xs.shape[0]):
            lpath1.append(path.lineto(Xs[j], Y1s[j]))
            lpath2.append(path.lineto(Xs[j], Y2s[j]))
        c1.stroke(lpath1, [ramp.get_rgb_color(i)])
        c1.stroke(lpath2, [ramp.get_rgb_color(i)])
    Y = (best[1] - minY) / (maxY - minY) * height
    c1.stroke(path.line(0, Y, width, Y), [style.linestyle.dashed])
    Y = (good[1] - minY) / (maxY - minY) * height
    c1.stroke(path.line(0, Y, width, Y), [style.linestyle.dotted])
    for i in numpy.linspace(ceil(minY * 10.0) / 10.0, floor(maxY * 10.0) / 10.0, 6):
        pos = (i - minY) / (maxY - minY) * height
        if gen == 'hg38':
            c2.stroke(path.line(0, pos, -0.05, pos))
            c2.text(-0.1, pos, '%0.2f' % i, [text.halign.right, text.valign.middle, text.size(text_size)])
        c2.stroke(path.line(0, pos, width, pos), [color.gray(0.9)])
    for i, res in enumerate(resolutions):
        if res >= 1000000:
            label = '%iMb' % (res / 1000000)
        else:
            label = '%iKb' % (res / 1000)
        if i == 0:
            halign = text.halign.left
            shift = -0.075
        elif i == resolutions.shape[0] - 1:
            halign = text.halign.right
            shift = 0.075
        else:
            halign = text.halign.center
            shift = 0
        pos = (numpy.log10(res) - minX) / (maxX - minX) * width
        c2.stroke(path.line(pos, 0, pos, -0.05))
        #c2.text(pos + shift, -0.1, label, [halign, text.valign.top, text.size(text_size)])
        c2.text(pos, -0.1, label, [halign.right, text.valign.middle, text.size(text_size), trafo.rotate(90)])
        c0.stroke(path.line(pos, 0, pos, height), [color.gray(0.9)])
        c2.stroke(path.line(pos, 0, pos, height), [color.gray(0.9)])
    c0.insert(c)
    c2.insert(c1)
    c0.stroke(path.rect(0, 0, width, height))
    c2.stroke(path.rect(0, 0, width, height))
    return c0, c2

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
    minX = 0.0
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
        self.colors = numpy.array([
            [0, 0, 0],
            [26, 26, 89],
            [77, 38, 166],
            [153, 51, 128],
            [255, 64, 38],
            [230, 140, 0],
            [230, 191, 26],
            [230, 230, 128],
            [255, 255, 255]], dtype=numpy.uint32)
        self.breaks = numpy.linspace(minval, maxval, self.colors.shape[0])
        self.span = self.breaks[1] - self.breaks[0]
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
