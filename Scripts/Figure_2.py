#!/usr/bin/env python

import sys
import os
from math import ceil, floor

import numpy
import h5py
from pyx import *
from PIL import Image
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
page_height = 22.5
text_size = -3

A_width = (page_width - 0.4) * 0.5
A_height = A_width + 0.05
C_width = A_width
C_height = A_height
D_width = A_width * 1.3
D_height = (A_height + C_height) * 0.5
B_width = page_width - D_width - 0.4
B_height = A_height + C_height - D_height

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_2.pdf" % base_dir
    output = open('%s/Results/quality_coverage_params.txt' % base_dir, 'w')
    output.close()
    c = canvas.canvas()
    c.insert(plot_A(A_width, A_height))
    c.insert(plot_B(B_width, B_height), [trafo.translate(0, -A_height - 0.4)])
    c.insert(plot_C(C_width, C_height), [trafo.translate(A_width + 0.4, 0)])
    c.insert(plot_D(D_width, D_height), [trafo.translate(B_width + 0.4, -A_height - 0.4)])
    c.writePDFfile(out_fname)

def plot_A(width, height):
    hoffset = 0.8
    voffset = 0.55
    voffset2 = 0.3
    spacer = 0.2
    pwidth = (width - hoffset - 2 * spacer) / 3.0
    pwidth2 = pwidth * 0.3
    pheight = (height - voffset - voffset2 - 2 * spacer) / 3.0
    pheight2 = pheight * 0.3
    data_fname = "%s/Results/quasar_noise_quality_results.txt" % base_dir
    data = load_data(data_fname)
    data['param'] /= 10.0
    c = canvas.canvas()
    c.text(0, height, 'A', [text.halign.left, text.valign.top, text.size(0)])
    for i, gen in enumerate(['Human', 'Mouse', 'Fruit Fly']):
        X0 = hoffset + i * (pwidth + spacer)
        for j in range(0, 75, 20):
            X = j / 75.0 * pwidth + X0
            c.stroke(path.line(X, voffset, X, voffset - 0.05))
            c.text(X, voffset - 0.1, j, [text.halign.center, text.valign.top, text.size(text_size)])
            for k in range(3):
                Y = voffset + (pheight + spacer) * k
                c.stroke(path.line(X, Y, X, Y + pheight), [color.gray(0.9)])
        c.text(hoffset + pwidth * (i + 0.5) + spacer * i, height - 0.05, gen,
            [text.valign.top, text.halign.center, text.size(text_size)])
    for i in range(3):
        p = plot_A_res(data, pwidth, pwidth2, pheight, pheight2, spacer, i)
        c.insert(p, [trafo.translate(hoffset, voffset + (pheight + spacer) * (2 - i))])
    c.text(hoffset + spacer + pwidth * 1.5, 0, r"\% Noise Reads",
        [text.valign.bottom, text.halign.center, text.size(text_size)])
    c.text(0, voffset + pheight * 1.5 + spacer, r"\% of Raw Quality Score",
        [text.valign.top, text.halign.center, text.size(text_size), trafo.rotate(90)])
    return c

def plot_A_res(data, pwidth, pwidth2, pheight, pheight2, spacer, h):
    c = canvas.canvas()
    reses = {}
    for gen in ['hg38', 'mm10', 'dm6']:
        reses[gen] = numpy.unique(data['resolution'][numpy.where(data['genome'] == gen)])
    samples = numpy.unique(data['sample'])
    for sample in samples:
        gen = data['genome'][numpy.where(data['sample'] == sample)[0][0]]
        where = numpy.where((data['sample'] == sample) & (data['resolution'] == reses[gen][h]))
        where2 = numpy.where((data['sample'] == sample) & (data['resolution'] == reses[gen][h]) &
                             (data['param'] == 0))
        data['score'][where] /= data['score'][where2][0]
    max_X = 75.0
    max_X2 = [2.0, 2.0, 2.0][2 - h]
    max_score = 1.05#numpy.amax(data['score'])
    min_score = 0
    span = max_score - min_score
    #max_score += span * 0.02
    span = max_score - min_score
    max_score2 = 1.005#max_score
    span2 = max_X2 / max_X * span
    min_score2 = max_score2 - span2
    scale = 5.0
    for i in range(int(ceil(scale * min_score)), int(scale * max_score) + 1):
        val = i / scale
        Y = (val - min_score) / span * pheight
        c.stroke(path.line(0, Y, -0.05, Y))
        c.text(-0.05, Y, "%i" % int(val * 100), [text.halign.right, text.valign.middle, text.size(text_size)])
        for j in range(3):
            X = (pwidth + spacer) * j
            c.stroke(path.line(X, Y, X + pwidth, Y), [color.gray(0.9)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        res = reses[gen][h]
        p = plot_lines(data[numpy.where((data['genome'] == gen) & (data['resolution'] == res))],
            None, res, min_score, max_score, pwidth, pheight)
        c.insert(p, [trafo.translate((pwidth + spacer) * i, 0)])
        X = 0.3
        Y = 0.4
        where = numpy.where((data['genome'] == gen) & (data['resolution'] == res) &
            (data['param'] <= max_X2))
        p.fill(path.rect(X, Y, pwidth2, pheight2), [color.rgb.white])
        p1 = canvas.canvas()
        Y2 = (1 - min_score2) / span2 * pheight2
        p1.stroke(path.line(0, Y2, pwidth2, Y2), [color.gray(0.9)])
        p1.stroke(path.line(0, Y2, -0.05, Y2))
        p1.insert(plot_lines(data[where], None, None, min_score2, max_score2, pwidth2, pheight2))
        p.insert(p1, [trafo.translate(X, Y)])
        p.stroke(path.rect(0, (min_score2 - min_score) / span * pheight, max_X2 / max_X * pwidth,
            span2 / span * pheight))
        p.stroke(path.line(pwidth * (max_X2 / max_X), (min_score2 - min_score + span2) / span * pheight,
            X + pwidth2, Y + pheight2))
        p.stroke(path.line(0, (min_score2 - min_score) / span * pheight, X, Y))
        p.stroke(path.line(X, Y, X, Y - 0.05))
        p.text(X, Y - 0.1, "0", [text.halign.center, text.valign.top, text.size(text_size)])
        p.stroke(path.line(X + pwidth2, Y, X + pheight2, Y - 0.05))
        p.text(X + pwidth2, Y - 0.1, "%i" % max_X2, [text.halign.center, text.valign.top, text.size(text_size)])
        p.stroke(path.rect(0, 0, pwidth, pheight))
        p.stroke(path.rect(X, Y, pwidth2, pheight2))
        if res >= 1000000:
            label = '%iMb' % (res / 1000000)
        else:
            label = '%iKb' % (res / 1000)
        p.text(pwidth - 0.1, 0.1, label,
            [text.halign.right, text.valign.bottom, text.size(text_size)])
    return c

def plot_B(width, height):
    spacer = 0.2
    hoffset = 0.85
    hoffset2 = 0.25
    voffset = 1.1
    c = canvas.canvas()
    pwidth = width - hoffset - hoffset2
    pheight = (height - voffset - 2 * spacer) / 3.0
    data_fname = "%s/Results/quasar_heterogeneity_results.txt" % base_dir
    data = load_het_data(data_fname)
    samples1 = numpy.unique(data['sample1'])
    samples2 = []
    for sample in samples1:
        samples2.append(data['sample2'][numpy.where(data['sample1'] == sample)[0][0]])
    ramp = Ramp(0.0, len(samples1) + 1)
    for i, res in enumerate([1000000, 200000, 40000]):
        c1 = canvas.canvas()
        where = numpy.where(data['resolution'] == res)[0]
        minY = numpy.amin(data['score'][where])
        maxY = numpy.amax(data['score'][where])
        spanY = maxY - minY
        minY -= 0.05 * spanY
        maxY += 0.05 * spanY
        spanY = maxY - minY
        for j in numpy.linspace(0, 1, 6):
            X = j * pwidth
            if i == 0:
                c1.stroke(path.line(X, 0, X, -0.1))
                c1.text(X, -0.15, "%i:%i" % (int(j * 100), 100 - int(j * 100)),
                    [text.valign.middle, text.halign.right, text.size(text_size), trafo.rotate(90)])
            c1.stroke(path.line(X, 0, X, pheight), [color.gray(0.9)])
        scale = 200.0
        while int(maxY * scale) - int(ceil(minY * scale)) > 5:
            scale /= 2
        for j in range(int(ceil(minY * scale)), int(maxY * scale) + 1):
            val = j / scale
            Y = (val - minY) / spanY * pheight
            c1.stroke(path.line(0, Y, -0.05, Y))
            c1.text(-0.05, Y, "%0.2f" % val, [text.halign.right, text.valign.middle, text.size(text_size)])
            c1.stroke(path.line(0, Y, pwidth, Y), [color.gray(0.9)])
        for j, sample in enumerate(samples1):
            where = numpy.where((data['resolution'] == res) & (data['sample1'] == sample))[0]
            Xs = data['balance'][where]
            Ys = data['score'][where]
            order = numpy.argsort(Xs)
            Xs = Xs[order] * pwidth
            Ys = (Ys[order] - minY) / spanY * pheight
            lpath = path.path(path.moveto(Xs[0], Ys[0]))
            for k in range(1, Xs.shape[0]):
                lpath.append(path.lineto(Xs[k], Ys[k]))
            c1.stroke(lpath, [ramp.get_rgb_color(j)])
        c1.stroke(path.rect(0, 0, pwidth, pheight))
        if res >= 1000000:
            label = '%iMb' % (res / 1000000)
        else:
            label = '%iKb' % (res / 1000)
        c1.text(pwidth + hoffset2, pheight * 0.5, label,
            [text.halign.center, text.valign.top, trafo.rotate(-90), text.size(text_size)])
        c.insert(c1, [trafo.translate(hoffset, voffset + (pheight + spacer) * i)])
    c.text(0, pheight * 1.5 + spacer + voffset, r"Quality Score",
        [text.valign.top, text.halign.center, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + pwidth * 0.5 + spacer, 0, r"Sample Percent Contributions",
        [text.valign.bottom, text.halign.center, text.size(text_size)])
    c.text(0, height, 'B', [text.halign.left, text.valign.top, text.size(0)])
    return c

def plot_C(width, height):
    hoffset = 0.8
    voffset = 0.55
    voffset2 = 0.3
    spacer = 0.2
    pwidth = (width - hoffset - 2 * spacer) / 3.0
    pheight = (height - voffset - voffset2 - 2 * spacer) / 3.0
    data_fname = "%s/Results/quasar_coverage_quality_results.txt" % base_dir
    sample_fname = "%s/Data/samples.txt" % base_dir
    data = load_data(data_fname)
    sample_data = load_sample_data(sample_fname)
    where = numpy.where(data['param'] == 0)[0]
    for i in where:
        data['param'][i] = sample_data[data['sample'][i]]['valid_cis_reads']
    data = data[numpy.where((data['genome'] != 'hg38') | (data['param'] >= 1000000))]
    data['param'] = numpy.log2(data['param'])
    c = canvas.canvas()
    reses = {}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        reses[gen] = numpy.unique(data['resolution'][numpy.where(data['genome'] == gen)])
        min_X = numpy.amin(data['param'][numpy.where(data['genome'] == gen)])
        max_X = numpy.amax(data['param'][numpy.where(data['genome'] == gen)])
        X0 = hoffset + (pwidth + spacer) * i 
        for j in range(max(-2, int(min_X - numpy.log2(1000000))), int(max_X - numpy.log2(1000000)) + 1, 2):
            if j < 0:
                label = "%0.2f" % (2 ** j)
            else:
                label = "%i" % (2 ** j)
            val = j + numpy.log2(1000000)
            X = X0 + (val - min_X) / (max_X - min_X) * pwidth
            c.text(X, voffset - 0.1, label, [text.halign.center, text.valign.top, text.size(text_size)])
            c.stroke(path.line(X, voffset, X, voffset - 0.085))
            for k in range(3):
                Y = voffset + (pheight + spacer) * k
                c.stroke(path.line(X, Y, X, Y + pheight), [color.gray(0.9)])
        for j in range(max(-2, int(min_X - numpy.log2(1000000))) - 2, int(max_X - numpy.log2(1000000)) + 3, 2):
            Xs = numpy.linspace(2 ** j  * 1000000, 2 ** (j + 2) * 1000000, 9)[1:-1]
            for x in Xs:
                X = (numpy.log2(x) - min_X) / (max_X - min_X) * pwidth
                if X > 0 and X < pwidth:
                    c.stroke(path.line(X + X0, voffset, X + X0, voffset - 0.05))
    c.text(0, height, 'C', [text.halign.left, text.valign.top, text.size(0)])
    for i in range(3):
        p = plot_C_res(data, pwidth, pheight, spacer, i + 1)
        c.insert(p, [trafo.translate(hoffset, voffset + (pheight + spacer) * (2 - i))])
    for i, gen in enumerate(['Human', 'Mouse', 'Fruit Fly']):
        X0 = hoffset + i * (pwidth + spacer)
        c.text(hoffset + pwidth * (i + 0.5) + spacer * i, height - 0.025, gen,
            [text.valign.top, text.halign.center, text.size(text_size)])
    c.text(hoffset + spacer + pwidth * 1.5, 0, r"Millions of Reads",
        [text.valign.bottom, text.halign.center, text.size(text_size)])
    c.text(0, voffset + pheight * 1.5 + spacer, r"Quality Score",
        [text.valign.top, text.halign.center, text.size(text_size), trafo.rotate(90)])
    return c

def plot_C_res(data, pwidth, pheight, spacer, h):
    reses = {}
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        reses[gen] = numpy.unique(data['resolution'][numpy.where(data['genome'] == gen)])
    c = canvas.canvas()
    min_score = 1.0
    max_score = 0.05
    max_X = {'hg38': -numpy.inf, 'mm10': -numpy.inf, 'dm6': -numpy.inf}
    min_X = {'hg38': numpy.inf, 'mm10': numpy.inf, 'dm6': numpy.inf}
    samples = numpy.unique(data['sample'])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        res = reses[gen][h - (gen == 'dm6')]
        where = numpy.where((data['genome'] == gen) & (data['resolution'] == res))
        min_score = min(min_score, numpy.amin(data['score'][where]))
        max_score = max(max_score, numpy.amax(data['score'][where]))
        min_X[gen] = min(min_X[gen], numpy.amin(data['param'][where]))
        max_X[gen] = max(max_X[gen], numpy.amax(data['param'][where]))
    span = max_score - min_score
    min_score -= span * 0.05
    max_score += span * 0.05
    span = max_score - min_score
    scale = 50.0
    if int(max_score * scale) - int(ceil(min_score * scale)) > 5:
        scale /= 2
    for i in range(int(ceil(min_score * scale)), int(max_score * scale) + 1):
        val = i / scale
        Y = (val - min_score) / span * pheight
        c.stroke(path.line(0, Y, -0.05, Y))
        c.text(-0.05, Y, "%0.2f" % val, [text.halign.right, text.valign.middle, text.size(text_size)])
        for j in range(3):
            X = (pwidth + spacer) * j
            c.stroke(path.line(X, Y, X + pwidth, Y), [color.gray(0.9)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        if gen == 'dm6':
            scale = 1000000.
        else:
            scale = 1000000.
        res = reses[gen][h - (gen == 'dm6')]
        p0 = plot_estimated_lines(data[numpy.where((data['genome'] == gen) & (data['resolution'] == res))],
            min_score, max_score, pwidth, pheight, ((gen == 'dm6') & (res == 100000)) | ((gen != 'dm6') &
            (res == 1000000)), gen=gen)
        c.insert(p0, [trafo.translate((pwidth + spacer) * i, 0)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        if gen == 'dm6':
            scale = 1000000.
        else:
            scale = 1000000.
        res = reses[gen][h - (gen == 'dm6')]
        p = plot_lines(data[numpy.where((data['genome'] == gen) & (data['resolution'] == res))],
            None, None, min_score, max_score, pwidth, pheight)
        c.insert(p, [trafo.translate((pwidth + spacer) * i, 0)])
        p.stroke(path.rect(0, 0, pwidth, pheight))
        if res >= 1000000:
            label = '%iMb' % (res / 1000000)
        else:
            label = '%iKb' % (res / 1000)
        p.text(0.1, pheight - 0.1, label,
            [text.halign.left, text.valign.top, text.size(text_size)])
    return c

def plot_D(width, height):   
    data_fname = "%s/Results/quasar_coverage_quality_results.txt" % base_dir
    data = load_data(data_fname)
    ideal_fname = "%s/Results/quality_coverage_params.txt" % base_dir
    ideal = load_ideal(ideal_fname)
    sample_fname = "%s/Data/samples.txt" % base_dir
    sample_data = load_sample_data(sample_fname)
    hoffset = 0.55
    hoffset2 = 0.25
    spacer = 0.3
    spacer2 = 0.45
    voffset = 0.6
    voffset2 = 0.3
    pwidth = (width - hoffset - hoffset2 - spacer * 2) / 3.0
    pheight = (height - voffset - voffset2 - spacer2 * 2) / 3.0
    Y_data = {}
    for key in sample_data:
        Y_data[key] = sample_data[key]['valid_cis_reads'] / float(sample_data[key]['valid_trans_reads'] + sample_data[key]['valid_cis_reads'])
    c = canvas.canvas()
    genomes = {'hg38': 'Human, 1Mb', 'mm10': 'Mouse, 1Mb', 'dm6': 'Fruit Fly, 50Kb'}
    where = numpy.where((((data['genome'] == 'dm6') & (data['resolution'] == 100000)) |
                         ((data['genome'] != 'dm6') & (data['resolution'] == 1000000))) &
                         (data['param'] == 0))[0]
    c.insert(plot_panel(width, pheight, Y_data, data[where]), [trafo.translate(0, voffset + 2 * (pheight + spacer2))])
    where = numpy.where(((data['genome'] == 'dm6') & (data['resolution'] == 100000) & (data['param'] == 1000000)) |
                        ((data['genome'] != 'dm6') & (data['resolution'] == 1000000) & (data['param'] == 10000000)))[0]
    c.insert(plot_panel(width, pheight, Y_data, data[where]), [trafo.translate(0, voffset + (pheight + spacer2))])
    c.insert(plot_panel(width, pheight, Y_data, ideal), [trafo.translate(0, voffset)])
    for i, gen in enumerate(['hg38', 'mm10', 'dm6']):
        c.text(hoffset + pwidth * (i + 0.5) + spacer * i, height - 0.025, genomes[gen],
            [text.halign.center, text.valign.top, text.size(text_size)])
    c.text(hoffset + pwidth * 1.5 + spacer, 0, 'Quality Score',
        [text.halign.center, text.valign.bottom, text.size(text_size)])
    c.text(0, voffset + pheight * 1.5 + spacer2, 'Valid Cis / Total Valid Read Ratio',
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    for i, label in enumerate(['Max Score', 'Uniform Coverage', 'Predicted Max']):
        c.text(width, voffset + pheight * (2.5 - i) + spacer2 * (2 - i), label,
            [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(-90)])
    c.text(0, height, 'D', [text.halign.left, text.valign.top, text.size(0)])
    return c

def plot_panel(width, pheight, Y_data, data):
    hoffset = 0.55
    hoffset2 = 0.25
    spacer = 0.3
    spacer2 = 0.45
    voffset = 0.6
    voffset2 = 0.3
    pwidth = (width - hoffset - hoffset2 - spacer * 2) / 3.0
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
        maxY += spanY * 0.1
        spanX = maxX - minX
        spanY = maxY - minY
        Xs = (Xs - minX) / (maxX - minX) * pwidth
        Ys = (Ys - minY) / (maxY - minY) * pheight
        slope, intercept, r, pval = linregress(Xs, Ys)[:4]
        c1.stroke(path.line(0, intercept, pwidth, pwidth * slope + intercept))
        c2.text(pwidth - 0.05, 0.65, r"r=%0.2f" % (r),
            [text.halign.right, text.valign.top, text.size(text_size)])
        c2.text(pwidth - 0.05, 0.35, r"p=%0.1e" % (pval),
            [text.halign.right, text.valign.top, text.size(text_size)])
        for k in range(Xs.shape[0]):
            c1.fill(path.circle(Xs[k], Ys[k], 0.04), [ramp.get_rgb_color(k)])
        c2.stroke(path.rect(0, 0, pwidth, pheight))
        scale = 10 ** -ceil(numpy.log10(maxX) - 2)
        val1 = ceil((minX + 0.12 * spanX) * scale) / scale
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
            label = "%0.1e" % val1
            label = label.replace('e-0', 'e-').replace('e+0', 'e')
            c2.stroke(path.line(0, Y, -0.05, Y))
            c2.text(-0.1, Y, label,
                [text.valign.bottom, text.halign.center, text.size(text_size), trafo.rotate(90)])
            val2 = floor((maxY - 0.1 * spanY) * scale) / scale
            Y = (val2 - minY) / spanY * pheight
            label = "%0.1e" % val2
            label = label.replace('e-0', 'e-').replace('e+0', 'e')
            c2.stroke(path.line(0, Y, -0.05, Y))
            c2.text(-0.1, Y, label,
                [text.valign.bottom, text.halign.center, text.size(text_size), trafo.rotate(90)])
        c.insert(c2, [trafo.translate(hoffset + i * (pwidth + spacer), 0)])
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

def load_het_data(fname):
    data = []
    infile = open(fname)
    line = infile.readline()
    line = infile.readline()
    while line:
        temp = line.split('\t')
        data.append((temp[0], temp[1], float(temp[2]), int(temp[3]), int(temp[4]), float(temp[5])))
        line = infile.readline()
    data = numpy.array(data, dtype=numpy.dtype([('sample1', 'S80'), ('sample2', 'S80'), ('balance', numpy.float64),
                       ('param', numpy.float64), ('resolution', numpy.float64), ('score', numpy.float64)]))
    return data

def plot_lines(data, gen, res, minY, maxY, width, height):
    scores = numpy.copy(data)
    samples = numpy.unique(scores['sample'])
    spanY = maxY - minY
    print minY, maxY
    minX = numpy.amin(data['param'])
    maxX = numpy.amax(data['param'])
    spanX = maxX - minX
    ramp = Ramp(0, len(samples) - 1)
    c = canvas.canvas()
    c1 = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])
    for h, sample in enumerate(samples):
        temp = scores[numpy.where(scores['sample'] == sample)]
        if temp.shape[0] == 0:
            continue
        Xs = temp['param']
        Ys = temp['score']
        order = numpy.argsort(Xs)
        Xs = Xs[order]
        Ys = Ys[order]
        try:
            if numpy.amax(Ys[1:]) > Ys[0] and Xs[1] > Xs[0] and Ys[0] == 1 and res is not None:
                print sample, res, list(Ys), Xs[numpy.where((Ys[:-1] > 1.0) & (Ys[1:] < 1.0))[0][0] + 1]
        except:
            pass
        Xs = (Xs - minX) / spanX * width
        Ys = (Ys - minY) / spanY * height
        lpath = path.path(path.moveto(Xs[0], Ys[0]))
        for i in range(1, Xs.shape[0]):
            lpath.append(path.lineto(Xs[i], Ys[i]))
        c1.stroke(lpath, [ramp.get_rgb_color(h)])
    c.insert(c1)
    if gen is not None:
        if res >= 1000000:
            label = "%i Mb" % (res / 1000000)
        else:
            label = "%i Kb" % (res / 1000)
        genome = {'hg38': "Human", 'mm10': 'Mouse', 'dm6': 'Fruit fly'}
        c.text(width * 0.5, height + 0.25, "%s, %s" % (genome[gen], label),
            [text.halign.center, text.valign.top, text.size(text_size)])
    return c

def plot_estimated_lines(data, minY, maxY, width, height, write_params=False, gen=None):
    if write_params:
        output = open('%s/Results/quality_coverage_params.txt' % base_dir, 'a')
    scores = numpy.copy(data)
    samples = numpy.unique(scores['sample'])
    spanY = maxY - minY
    minX = numpy.amin(data['param'])
    maxX = numpy.amax(data['param'])
    spanX = maxX - minX
    ramp = Ramp(0, len(samples) - 1)
    c = canvas.canvas()
    c1 = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])

    def f(x, x0, k, L):
        return L / (1 + numpy.exp(-k * (x - x0)))

    for h, sample in enumerate(samples):
        temp = scores[numpy.where(scores['sample'] == sample)]
        if temp.shape[0] == 0:
            continue
        Xs = temp['param']
        Ys = temp['score']
        order = numpy.argsort(Xs)
        Xs = Xs[order]
        Ys = Ys[order]
        if (Xs[-1] < numpy.log2(8000000) and gen != 'dm6') or (Xs[-1] < numpy.log2(1000000) and gen == 'dm6'):
            #print sample
            continue
        try:
            params = curve_fit(f, Xs, Ys, p0=(Xs[-1], 0.5, Ys[-1] * 2), maxfev=(5000*Xs.shape[0]),
                 bounds=((-numpy.inf, -numpy.inf, 0), (numpy.inf, numpy.inf, 2)))[0]
            Xs = numpy.linspace(minX, maxX, 20)
            Ys = f(Xs, *params)
            Xs = (Xs - minX) / spanX * width
            Ys = (Ys - minY) / spanY * height
            lpath = path.path(path.moveto(Xs[0], Ys[0]))
            for i in range(1, Xs.shape[0]):
                lpath.append(path.lineto(Xs[i], Ys[i]))
            c1.stroke(lpath, [color.gray(0.7)])
            if write_params and params[2] < 2.0:
                print >> output, "%s\t%s\t%f\t%f\t%f" % (sample, gen, params[2], params[0], params[1])
        except:
            #print sample
            pass
    if write_params:
        output.close()
    c.insert(c1)
    return c

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
