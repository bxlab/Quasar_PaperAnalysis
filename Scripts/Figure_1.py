#!/usr/bin/env python

import sys
import os
from math import ceil

import numpy
import h5py
from pyx import *
from PIL import Image

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{uarial}")
text.preamble(r"\usepackage{amsmath}")
text.preamble(r"\usepackage{dsfont}")
text.preamble(r"\usepackage{graphicx}")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
text.preamble(r"\renewcommand*\encodingdefault{T1}")
text.preamble(r"\newcommand{\textoverline}[1]{$\overline{\mbox{#1}}$}")

page_width = 17.3
page_height = 22.4
text_size = -3

space = 0.3
A_width = 2.2
hm_width = (page_width - space * 2 - 0.4 - A_width - 0.6) / 4.4
B_width = hm_width * 3 + 0.2 + 0.3
B_height = hm_width * 3 + 0.2 + 0.25
A_height = B_height
C_width = page_width - A_width - B_width - space * 2 - 0.2
C_height = hm_width * 2 + 0.2 + 0.25
D_width = C_width - 1.
D_height = A_height - C_height - 0.5

base_dir = "."

def main():
    out_fname = "%s/Figures/Figure_1.pdf" % base_dir
    c = canvas.canvas()
    c.insert(plot_A(A_width, A_height))
    c.insert(plot_B(B_width, B_height), [trafo.translate(A_width + space - 0.2, 0)])
    c.insert(plot_C(C_width), [trafo.translate(A_width + B_width + space * 2 - 0.2, 0)])
    c.insert(plot_key(), [trafo.translate(page_width - C_width, -A_height)])
    c.insert(plot_D(D_width, D_height), [trafo.translate(page_width - D_width - 0.4, -A_height)])
    c.writePDFfile(out_fname)

def plot_A(awidth, height):
    width = 4.0
    c = canvas.canvas([canvas.clip(path.rect(0, 0, awidth, -height))])
    c.text(0, 0, 'A', [text.halign.left, text.valign.top, text.size(0)])
    pheight = 2.85
    p = numpy.array([
        [width * 0.2, pheight * 0.3],
        [width * 0.45, pheight * 0.1],
        [width * 0.45, pheight * 0.35],
        [width * 0.8 - 0.1, pheight * 0.8 - 0.1]
    ])
    arrow = path.path(path.moveto(0, -0.2), path.lineto(0.1, -0.3), path.lineto(0, 0), path.lineto(-0.1, -0.3), path.closepath())
    c1 = canvas.canvas()
    c1.stroke(path.path(path.moveto(width * 0.25, pheight * 0.5),
        path.curveto(width * 0.25 - 0.1, pheight * 0.5 - 0.4, p[0, 0], p[0, 1] + 0.4, p[0, 0], p[0, 1]),
        path.curveto(p[0, 0], p[0, 1] - 0.6, p[1, 0] - 0.6, p[1, 1], p[1, 0], p[1, 1]),
        path.curveto(p[1, 0] + 0.6, p[1, 1], p[2, 0] + 0.6, p[2, 1] - 0.2, p[2, 0], p[2, 1]),
        path.curveto(p[2, 0] - 1.0, p[2, 1] + 0.2, p[3, 0] - 0.6, p[3, 1] + 0.2, p[3, 0], p[3, 1]),
        path.curveto(p[3, 0] + 0.1, p[3, 1] - 0.4, p[3,0] + 0.1, p[3,1] - 0.2, p[3,0], p[3,1] - 0.5)),
        [style.linewidth.THIck, color.rgb(r=(255/255.), g=(64/255.), b=(38/255.))])
    c1.stroke(path.line(p[0, 0], p[0, 1], p[1, 0], p[1, 1]), [color.gray(0.0), style.linewidth.THICK])
    c1.stroke(path.line(p[0, 0], p[0, 1], p[2, 0], p[2, 1]), [color.gray(0.2), style.linewidth.THICk])
    c1.stroke(path.line(p[1, 0], p[1, 1], p[2, 0], p[2, 1]), [color.gray(0.3), style.linewidth.THICk])
    c1.stroke(path.line(p[0, 0], p[0, 1], p[3, 0], p[3, 1]), [color.gray(0.9), style.linewidth.Thick])
    c1.stroke(path.line(p[1, 0], p[1, 1], p[3, 0], p[3, 1]), [color.gray(0.8), style.linewidth.thick])
    c1.stroke(path.line(p[2, 0], p[2, 1], p[3, 0], p[3, 1]), [color.gray(0.75), style.linewidth.Thick])
    c1.stroke(path.circle(p[0, 0], p[0, 1], 0.1), [deco.filled([color.rgb.white])])
    c1.stroke(path.circle(p[1, 0], p[1, 1], 0.1), [deco.filled([color.rgb.white])])
    c1.stroke(path.circle(p[2, 0], p[2, 1], 0.1), [deco.filled([color.rgb.white])])
    c1.stroke(path.circle(p[3, 0], p[3, 1], 0.1), [deco.filled([color.rgb.white])])
    c2 = canvas.canvas()
    c2.insert(c1, [trafo.rotate(30), trafo.translate(0.1, -3.8)])
    c2.text(awidth * 0.5, 0, 'Consistent', [text.valign.top, text.halign.center, text.size(text_size)])
    c2.fill(arrow, [trafo.rotate(90), trafo.translate(1.6, -1.7)])
    c.insert(c2, [trafo.translate(0, -0.1)])
    c1 = canvas.canvas()
    c1.stroke(path.path(path.moveto(width * 0.25, pheight * 0.5),
        path.curveto(width * 0.25 - 0.1, pheight * 0.5 - 0.4, p[0, 0], p[0, 1] + 0.4, p[0, 0], p[0, 1]),
        path.curveto(p[0, 0], p[0, 1] - 0.6, p[1, 0] - 0.6, p[1, 1], p[1, 0], p[1, 1]),
        path.curveto(p[1, 0] + 0.6, p[1, 1], p[2, 0] + 0.6, p[2, 1] - 0.2, p[2, 0], p[2, 1]),
        path.curveto(p[2, 0] - 1.0, p[2, 1] + 0.2, p[3, 0] - 0.6, p[3, 1] + 0.2, p[3, 0], p[3, 1]),
        path.curveto(p[3, 0] + 0.1, p[3, 1] - 0.4, p[3,0] + 0.1, p[3,1] - 0.2, p[3,0], p[3,1] - 0.5)),
        [style.linewidth.THIck, color.rgb(r=(255/255.), g=(64/255.), b=(38/255.))])
    c1.stroke(path.line(p[1, 0], p[1, 1], p[3, 0], p[3, 1]), [color.gray(0.8), style.linewidth.THICK])
    c1.stroke(path.line(p[0, 0], p[0, 1], p[1, 0], p[1, 1]), [color.gray(0.0), style.linewidth.THICK])
    c1.stroke(path.line(p[0, 0], p[0, 1], p[2, 0], p[2, 1]), [color.gray(0.2), style.linewidth.THICk])
    c1.stroke(path.line(p[1, 0], p[1, 1], p[2, 0], p[2, 1]), [color.gray(0.3), style.linewidth.THICk])
    c1.stroke(path.line(p[0, 0], p[0, 1], p[3, 0], p[3, 1]), [color.gray(0.9), style.linewidth.Thick])
    c1.stroke(path.line(p[2, 0], p[2, 1], p[3, 0], p[3, 1]), [color.gray(0.75), style.linewidth.Thick])
    c1.stroke(path.circle(p[0, 0], p[0, 1], 0.1), [deco.filled([color.rgb.white])])
    c1.stroke(path.circle(p[1, 0], p[1, 1], 0.1), [deco.filled([color.rgb.white])])
    c1.stroke(path.circle(p[2, 0], p[2, 1], 0.1), [deco.filled([color.rgb.white])])
    c1.stroke(path.circle(p[3, 0], p[3, 1], 0.1), [deco.filled([color.rgb.white])])
    c2 = canvas.canvas()
    c2.insert(c1, [trafo.rotate(30), trafo.translate(0.1, -3.8)])
    c2.text(awidth * 0.5, 0, 'Inconsistent', [text.valign.top, text.halign.center, text.size(text_size)])
    c2.fill(arrow, [trafo.rotate(90), trafo.translate(1.6, -1.7)])
    c.insert(c2, [trafo.translate(0, -height * 1.0 / 3.0 - 0.55)])
    key = canvas.canvas()
    c1 = canvas.canvas()
    step = width / 20.0
    c1.stroke(path.circle(0, 1.5, 0.1))
    c1.text(0.3, 1.5, "Hi-C Bin", [text.halign.left, text.valign.middle, text.size(text_size)])
    c1.stroke(path.line(-0.15, 1.2, 0.15, 1.2), [style.linewidth.THIck, color.rgb.red])
    c1.text(0.3, 1.2, "Chromatin", [text.halign.left, text.valign.middle, text.size(text_size)])
    c1.stroke(path.line(-0.15, 0.775, 0.15, 0.775), [style.linewidth.THIck])
    c1.text(0.3, 0.9, "Correlation /", [text.halign.left, text.valign.middle, text.size(text_size)])
    c1.text(0.3, 0.65, "Enrichment", [text.halign.left, text.valign.middle, text.size(text_size)])
    key.insert(c1, [trafo.translate(awidth * 0.5 - step * 4, -0.35)])

    key.stroke(path.line(awidth * 0.5 - step * 3.07, 0.05, awidth * 0.5 - step * 3.07, -0.1), [style.linewidth.normal])
    key.stroke(path.line(awidth * 0.5 - step * 2.07, 0.05, awidth * 0.5 - step * 2.07, -0.1), [style.linewidth.thick])
    key.stroke(path.line(awidth * 0.5 - step * 1.07, 0.05, awidth * 0.5 - step * 1.07, -0.1), [style.linewidth.Thick])
    key.stroke(path.line(awidth * 0.5 - step * 0.07, 0.05, awidth * 0.5 - step * 0.07, -0.1), [style.linewidth.THick])
    key.stroke(path.line(awidth * 0.5 + step * 0.93, 0.05, awidth * 0.5 + step * 0.93, -0.1), [style.linewidth.THIck])
    key.stroke(path.line(awidth * 0.5 + step * 1.93, 0.05, awidth * 0.5 + step * 1.93, -0.1), [style.linewidth.THICk])
    key.stroke(path.line(awidth * 0.5 + step * 2.93, 0.05, awidth * 0.5 + step * 2.93, -0.1), [style.linewidth.THICK])
    key.text(awidth * 0.5 - step * 3.05 - 0.1, -0.065, '-',
        [text.halign.right, text.valign.middle, text.size(text_size)])
    key.text(awidth * 0.5 + step * 3.05 + 0.1, -0.025, '+',
        [text.halign.left, text.valign.middle, text.size(text_size)])
    key.text(awidth * 0.5, -0.15, 'Enrichment', [text.halign.center, text.valign.top, text.size(text_size)])
    strength = numpy.zeros((256, 25), dtype=numpy.uint32)
    strength.shape = (25, 256)
    strength += 256 ** 3 * 255 + (256 ** 2 + 257) * numpy.arange(256)[::-1].astype(numpy.uint32).reshape(1, -1)
    pilImage = Image.frombuffer('RGBA', (strength.shape[1], strength.shape[0]), strength, 'raw', 'RGBA', 0, 1)
    key.insert(bitmap.bitmap(0, 0, pilImage, width=(step * 6.1)), [trafo.translate(awidth * 0.5 - step * 3.05, -0.6)])
    key.stroke(path.rect(awidth * 0.5 - step * 3.05, -0.6, step * 6.1, step * 6.1 / 256. * 25.))
    key.text(awidth * 0.5 - step * 3.05 - 0.1, -0.6 + step * 6.1 / 256. * 12.5 - 0.04, '-', [text.halign.right, text.valign.middle, text.size(text_size)])
    key.text(awidth * 0.5 + step * 3.05 + 0.1, -0.6 + step * 6.1 / 256. * 12.5, '+', [text.halign.left, text.valign.middle, text.size(text_size)])
    key.text(awidth * 0.5, -0.65, 'Correlation', [text.halign.center, text.valign.top, text.size(text_size)])
    c2 = canvas.canvas()
    c2.insert(key, [trafo.translate(0, 0.9)])
    c.insert(c2, [trafo.translate(0, -height)])
    return c

def plot_B(width, height):
    spacer = 0.1
    hoffset = 0.3
    voffset = 0.3
    hm_width = (width - 2 * spacer - hoffset) / 3.0
    cstart = 55000000
    cstop = cstart + 21000000
    c = canvas.canvas()
    c.text(0, 0, 'B', [text.halign.left, text.valign.top, text.size(0)])
    in_fname = "%s/Results/Encode_hg38_A549_HindIII_Rep1.quasar" % base_dir
    infile = h5py.File(in_fname, 'r')
    start = infile['starts'][numpy.where(infile['chromosomes'][...] == '1')[0][0]]
    for i, res in enumerate([1000000, 200000, 40000]):
        M = infile['dist.1.0C.%iR' % res].shape[1]
        raw = compact_to_full(infile['dist.1.0C.%iR' % res], start, res, cstart, cstop)
        corr = compact_to_full(infile['corr.1.0C.%iR' % res], start, res, cstart, cstop)
        valid = row_to_full(infile['valid.1.0C.%iR' % res], corr, M, start, res, cstart, cstop)
        raw *= valid
        raw = (raw + 1) ** 0.5
        corr *= valid
        trans = raw * corr
        c1 = canvas.canvas()
        c1.insert(plot_hm(trans, valid, hm_width),
            [trafo.translate(0, 0)])
        c1.insert(plot_hm(corr, valid, hm_width),
            [trafo.translate(0, hm_width + spacer)])
        c1.insert(plot_hm(raw, valid, hm_width),
            [trafo.translate(0, (hm_width + spacer) * 2)])
        if res >= 1000000:
            label = '%iMb' % (res / 1000000)
        else:
            label = '%iKb' % (res / 1000)
        c1.text(hm_width * 0.5, 3 * hm_width + 2 * spacer + voffset - 0.025, label,
            [text.halign.center, text.valign.top, text.size(text_size)])
        if i == 2:
            for j, label in enumerate(['matrix S', 'matrix C', 'matrix T']):
                c1.text(0.2, 0.2 + (hm_width + spacer) * (2 - j), label, [text.halign.left, text.valign.bottom, text.size(text_size), color.rgb.white])
        c.insert(c1, [trafo.translate(hoffset + (hm_width + spacer) * i, -voffset - hm_width * 3 - spacer * 2)])
    infile.close()
    for i, label in enumerate(['Raw', 'Correlation', 'Transformed']):
        c.text(0.025, -voffset - (i + 0.5) * hm_width - spacer * i, label,
            [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    return c

def plot_C(width):
    spacer = 0.1
    hoffset = 0.3
    voffset = 0.3
    hm_width = (width - spacer - hoffset) / 2.0
    in_fnames = [
        "%s/Results/GSE35156_mm10_ES_HindIII_Rep1.quasar" % base_dir,
        "%s/Results/GSE60494_mm10_ES_MboI_Rep2.quasar" % base_dir,
    ]
    o_cstart = 75000000
    o_cstop = o_cstart + 25000000
    c = canvas.canvas()
    c.text(0, 0, 'C', [text.halign.left, text.valign.top, text.size(0)])
    for h, fname in enumerate(in_fnames):
        infile = h5py.File(fname, 'r')
        start = infile['starts'][numpy.where(infile['chromosomes'][...] == '1')[0][0]]
        cstart = o_cstart
        cstop = o_cstop
        for i, res in enumerate([200000, 40000, 10000]):
            M = infile['dist.1.0C.%iR' % res].shape[1]
            raw = compact_to_full(infile['dist.1.0C.%iR' % res], start, res, cstart, cstop)
            corr = compact_to_full(infile['corr.1.0C.%iR' % res], start, res, cstart, cstop)
            valid = row_to_full(infile['valid.1.0C.%iR' % res], corr, M, start, res, cstart, cstop)
            raw *= valid
            raw = (raw + 1) ** 0.5
            corr *= valid
            trans = raw * corr
            c.insert(plot_hm(trans, valid, hm_width),
                [trafo.translate(hoffset + (hm_width + spacer) * h,
                 -voffset - hm_width * (i + 1) - spacer * i)])
            span = cstop - cstart
            cstop = cstart + span / 4
            if h == 0:
                if res < 1000000:
                    label = "%iKb" % (res / 1000)
                else:
                    label = "%iMb" % (res / 1000000)
                c.text(0.025, -voffset - hm_width * (i + 0.5) - spacer * i, label,
                    [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
        c.text(hoffset + (h + 0.5) * hm_width + spacer * h, -0.025, "Sample %i" % (h + 1),
            [text.halign.center, text.valign.top, text.size(text_size)])
    return c

def plot_D(width, height):
    hoffset = 0.9
    voffset = 0.5
    voffset2 = 0.25
    pwidth = width - hoffset
    pheight = height -voffset - voffset2
    quality_fname = "%s/Results/quasar_quality_results.txt" % base_dir
    data = load_quality_data(quality_fname)
    c = canvas.canvas()
    c1 = canvas.canvas()
    min_score = 0
    max_score = numpy.amax(data['score'])
    span = max_score - min_score
    max_score += span * 0.05
    span = max_score - min_score
    base = numpy.amin(data['resolution'])
    min_X = 0.0
    max_X = numpy.log(numpy.amax(data['resolution']) / base)
    spanX = max_X - min_X
    samples = numpy.unique(data['sample'])
    ramp = Ramp(0, 1.5)
    for i in [10000, 40000, 200000, 1000000]:
        X = (numpy.log(i / base) - min_X) / spanX * pwidth
        c1.stroke(path.line(X, 0, X, -0.05))
        if i >= 1000000:
            label = "%iMb" % (i / 1000000)
            c1.text(X, -0.1, label, [text.halign.right, text.valign.top, text.size(text_size)])
        else:
            label = "%iKb" % (i / 1000)
            c1.text(X, -0.1, label, [text.halign.center, text.valign.top, text.size(text_size)])
    scale = 20.0
    for i in range(int(ceil(min_score * scale)), int(max_score * scale) + 1):
        Y = (i / scale - min_score) / span * pheight
        c1.stroke(path.line(0, Y, -0.05, Y))
        c1.text(-0.1, Y, "%0.2f" % (i / scale), [text.halign.right, text.valign.middle, text.size(text_size)])
    for i, sample in enumerate(samples):
        where = numpy.where(data['sample'] == sample)[0]
        Xs = numpy.log(data['resolution'][where] / base)
        Ys = data['score'][where]
        order = numpy.argsort(Xs)
        Xs = Xs[order]
        Ys = Ys[order]
        lpath = path.path(path.moveto((Xs[0] - min_X) / spanX * pwidth, (Ys[0] - min_score) / span * pheight))
        for j in range(1, Xs.shape[0]):
            X = (Xs[j] - min_X) / spanX * pwidth
            Y = (Ys[j] - min_score) / span * pheight
            lpath.append(path.lineto(X, Y))
        c1.stroke(lpath, [ramp.get_rgb_color(i)])
    c1.stroke(path.rect(0, 0, pwidth, pheight))
    c.insert(c1, [trafo.translate(hoffset, voffset)])
    c.text(hoffset + pwidth * 0.5, 0, "Resolution", [text.halign.center, text.valign.bottom, text.size(text_size)])
    c.text(0, voffset + pheight * 0.5, "Quality Score",
        [text.halign.center, text.valign.top, text.size(text_size), trafo.rotate(90)])
    c.text(hoffset + 0.05, height, "Sample 1",
        [text.valign.top, text.halign.left, text.size(text_size), ramp.get_rgb_color(0)])
    c.text(hoffset + pwidth - 0.05, height, "Sample 2",
        [text.valign.top, text.halign.right, text.size(text_size), ramp.get_rgb_color(1)])
    return c

def plot_key():
    key_val = numpy.zeros((500, 100), dtype=numpy.float64)
    key_val += numpy.arange(key_val.shape[0]).reshape(-1, 1)
    key = canvas.canvas()
    kw = 1.5
    kh = kw * 100. / 500.0
    key.insert(plot_hm(key_val, numpy.ones(key_val.shape), kw, False, True),
        [trafo.rotate(90), trafo.translate(0.5 * kh, 0)])
    key.stroke(path.rect(-kh * 0.5, 0, kh, kw))
    key.text(0, -0.05, "Low",
        [text.halign.center, text.valign.top, text.size(text_size)])
    key.text(0, kw + 0.3, "High",
        [text.halign.center, text.valign.top, text.size(text_size)])
    key.fill(path.rect(-0.3, -0.7, 0.2, 0.2),
        [color.gray(0.5)])
    key.text(0, -0.35, 'No', [text.halign.left, text.valign.top, text.size(text_size)])
    key.text(0, -0.6, 'Data', [text.halign.left, text.valign.top, text.size(text_size)])
    c = canvas.canvas()
    c.insert(key, [trafo.translate(0.1, 0.8)])
    return c

def plot_hm_from_infile(infile, key, width, full=False):
    data = infile[key][...]
    valid = data > 0
    data[valid] = numpy.log(data[valid])
    return plot_hm(data, valid, width, rescale=False, full=full)

def plot_multires_hms(infile, key, width): 
    data = infile[key][...]
    plots = []
    hm_width = width * 0.7
    step = (width - hm_width) / 2.0
    for i in range(3):
        N = data.shape[0]
        L = 5 ** i
        M = N / L
        data1 = numpy.zeros((M, M), dtype=numpy.float64)
        for j in range(L):
            for k in range(L):
                data1 += data[j::L, k::L]
        valid = data1 > 0
        data1[valid] = numpy.log(data1[valid])
        plots.append(plot_hm(data1, valid, hm_width, rescale=False, full=True))
        plots[-1].stroke(path.rect(0, 0, hm_width, hm_width))
    c = canvas.canvas()
    for i in range(len(plots))[::-1]:
        c.insert(plots[i], [trafo.translate((2 - i) * step, i * step)])
    return c    

def compact_to_full(cdata, cstart, res, start, stop):
    starti = (start - cstart) / res
    stopi = (stop - cstart) / res
    N = stopi - starti
    fdata = numpy.zeros((N, N), dtype=numpy.float64)
    for i in range(min(N, cdata.shape[1])):
        temp1 = numpy.arange(N - i)
        temp2 = numpy.arange(i, N)
        fdata[temp1, temp2] = cdata[starti:(stopi - i), i]
    indices = numpy.triu_indices(N, 1)
    fdata[indices[1], indices[0]] = fdata[indices]
    return fdata

def row_to_full(cdata, corr, M, cstart, res, start, stop):
    starti = (start - cstart) / res
    stopi = (stop - cstart) / res
    N = stopi - starti
    fdata = numpy.zeros((N, N), dtype=numpy.float64)
    for i in range(min(N, M)):
        temp1 = numpy.arange(N - i)
        temp2 = numpy.arange(i, N)
        fdata[temp1, temp2] = cdata[starti:(stopi - i)] * cdata[(starti + i):(stopi)]
    indices = numpy.triu_indices(N, 1)
    fdata[indices[1], indices[0]] = fdata[indices]
    #fdata[numpy.where(corr == numpy.inf)] = 0
    return fdata

def plot_hm(data, valid, width, rescale=True, full=False):
    where = numpy.where(valid)
    scores = data[where]
    scores.sort()
    if full:
        ramp = Ramp(scores[0], scores[-1])
    else:
        ramp = Ramp(scores[int(0.01 * scores.shape[0])], scores[int(0.99 * scores.shape[0])])
    hm = numpy.zeros(data.shape, dtype=numpy.uint32)
    hm.shape = (hm.shape[1], hm.shape[0])
    hm.fill(int('FF888888', 16))
    hm[where[1], where[0]] = ramp.get_color(data[where])
    if rescale:
        scale = 200.0 / hm.shape[0]
        if scale > 1:
            M = int(numpy.ceil(scale))
            N = hm.shape[0] * M
            new_hm = numpy.zeros((N, N), dtype=numpy.uint32)
            for i in range(M):
                for j in range(M):
                    new_hm[i::M, j::M] = hm
            hm = new_hm
    pilImage = Image.frombuffer('RGBA', (hm.shape[1], hm.shape[0]), hm, 'raw', 'RGBA', 0, 1)
    c = canvas.canvas()
    c.insert(bitmap.bitmap(0, 0, pilImage, width=width))
    return c

def load_quality_data(fname):
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
    data = data[numpy.where(((data['sample'] == 'GSE35156_mm10_ES_HindIII_Rep1') |
                             (data['sample'] == 'GSE60494_mm10_ES_MboI_Rep2')) &
                             (data['param'] == 0))]
    return data

class Ramp():
    def __init__(self, min, max):
        self.min = min
        self.max = max
        self.breaks = numpy.linspace(min, max, 9)
        self.span = self.breaks[1] - self.breaks[0]
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
        self.base = 256 ** 3 * 255
        self.scale = numpy.array([1, 256, 256 ** 2], dtype=numpy.uint32)

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
