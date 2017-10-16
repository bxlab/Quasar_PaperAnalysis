#!/usr/bin/env python

import sys
import numpy

base_dir = "."

def main():
    out_fname = "%s/Figures/table_S1.tex" % base_dir
    data_fname = "%s/Data/samples.txt" % base_dir
    data = load_sample_data(data_fname)
    output = open(out_fname, 'w')
    print >> output, r"\begin{landscape}"
    print >> output, r"\begin{ThreePartTable}"
    print >> output, r"\small"
    #print >> output, r"\begin{ltabulary}{ p{5.2cm}p{0.9cm}p{0.9cm}p{0.9cm}p{1.2cm}p{1.2cm}p{1.2cm}p{1.2cm}p{1.2cm}p{1.2cm} }"
    print >> output, r"\begin{ltabulary}{ p{5.4cm}p{1.2cm}p{1.0cm}p{1.6cm}p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm}p{1.5cm} }"
    print >> output, r"\caption{\textbf{Datasets and associated read statistics used for analysis.}}\\"
    print >> output, r" \toprule"
    print >> output, r" Sample Name & Source & Genome & Restriction Enzyme & Total Reads & Cis Reads & Trans Reads & Insert Too Large & Circular Fragments & Failed RE Cut \\"
    print >> output, r" \toprule"
    print >> output, r" \midrule"
    print >> output, r" \endfirsthead"
    print >> output, r" \multicolumn{10}{l}{\footnotesize\itshape\tablename~\thetable: continued from previous page} \\"
    print >> output, r" \toprule"
    print >> output, r" Sample Name & Source & Genome & Restriction Enzyme & Total Reads & Cis Reads & Trans Reads & Insert Too Large & Circular Fragments & Failed RE Cut \\"
    print >> output, r" \toprule"
    print >> output, r" \midrule"
    print >> output, r" \endhead"
    print >> output, r" \midrule"
    print >> output, r" \multicolumn{10}{r}{\footnotesize\itshape\tablename~\thetable: continued on following page} \\"
    print >> output, r" \endfoot"
    print >> output, r" \bottomrule"
    print >> output, r" \endlastfoot"
    lines = []
    for key in data:
        temp = key.split('_')
        name = "%s %s %s" % (temp[0], temp[2], temp[4])
        genome = temp[1]
        re = temp[3]
        if temp[0][:3] == 'GSE':
            source = r"{\href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s}{GEO}}" % data[key]['GEO_ID']
        else:
            source = r"{\href{https://www.encodeproject.org/biosamples/%s/}{Encode}}" % data[key]['GEO_ID']
        lines.append([name, source, genome, re, str2num(data[key]['total_reads']),
            str2num(data[key]['valid_cis_reads']), str2num(data[key]['valid_trans_reads']),
            str2num(data[key]['insert_size']), str2num(data[key]['same_fragment']),
            str2num(data[key]['failed_cut'])])
    temp = numpy.array(lines)
    order = numpy.lexsort((temp[:, 0], temp[:, 2]))
    for i in order:
        print >> output, r" %s \\" % (" & ".join(lines[i]))
        print >> output, r" \hline"
    print >> output, r"\end{ltabulary}"
    print >> output, r"\end{ThreePartTable}"
    print >> output, r"\end{landscape}"
    output.close()

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
            temp1[labels[i]] = val
        data[temp[0]] = temp1
    return data

def str2num(s):
    n = []
    while len(s) > 3:
        n = [s[-3:]] + n
        s = s[:-3]
    if len(s) > 0:
        n = [s] + n
    return ','.join(n)


if __name__ == "__main__":
    main()
