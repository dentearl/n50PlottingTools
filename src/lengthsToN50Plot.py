#!/usr/bin/env python
"""
lengthsToN50Plot.py
12 January 2012
dent earl, dearl(a)soe ucsc edu

This script takes an arbitrary number of files which consist of lengths,
one per line, and it then produces a figure showing
the cumulative plot of the N statistics for the files.

"""
##############################
# Copyright (C) 2009-2011 by
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedict.paten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC).
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
import matplotlib
matplotlib.use('Agg')
#####
# the param pdf.fonttype allows for text to be editable in Illustrator.
# Use either Output Type 3 (Type3) or Type 42 (TrueType)
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
import numpy
from optparse import OptionParser
import os
import sys
import re


class LengthObj:
    """ Used to hold the length length information from one file
    """
    def __init__(self, name='', lengths=None, preSorted=False):
        if lengths is None:
            lengths = []
        if not isinstance(lengths, numpy.ndarray):
            raise RuntimeError('Input `lengths\' must be of type '
                               'numpy.ndarray, not %s.' % lengths.__class__)
        self.name = name
        self.lengths = lengths
        self.n = numpy.alen(lengths)
        if not preSorted:
            self.lengths = numpy.sort(self.lengths)
            self.lengths = self.lengths[::-1] # reverse
        # xData is paired one-to-one with lengths
        self.xData = numpy.cumsum(self.lengths, dtype=numpy.float64)

def initOptions(parser):
    parser.add_option(
        '--genomeLength', dest='genomeLength', type='float',
        help=('Total length of the genome. default=Max of sum of '
              'input file lengths.'))
    parser.add_option('--title', dest='title',
                      type='string', default='N-Statistics',
                      help='Title of the plot. default=%default')
    parser.add_option('--linear', dest='linear', default=False,
                      action='store_true',
                      help='Puts y axis into linear scale. default=%default')
    parser.add_option(
        '--n50Line', dest='n50Line', default=False, action='store_true',
        help=('Adds straight lines from the y-axis to the curves. '
              'default=%default'))
    parser.add_option(
        '--xlabel', dest='xlabel', type='string',
        default='Cumulative length proportional to genome length',
        help='Label on the x-axis. default=%default')
    parser.add_option('--reportN50Values', dest='reportN50Values',
                      default=False, action='store_true',
                      help='prints n50 values to stdout. default=%default')
    parser.add_option(
        '--preSorted', dest='preSorted', default=False, action='store_true',
        help=('Switches off the sort step, this can help with enormous numbers '
              'of reads. default=%default'))
    parser.add_option(
        '--dpi', dest='dpi', default=300, type='int',
        help=('Dots per inch of the output, if --outFormat is all or png. '
              'default=%default'))
    parser.add_option('--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps]. default=%default')
    parser.add_option(
        '--out', dest='out', default='myPlot', type='string',
        help=('path/filename where figure will be created. '
              'No extension needed. default=%default'))

def checkOptions(options, args, parser):
    options.log = True
    if options.linear:
        options.log = False
    if len(args) is 0:
        parser.error('Specify input files on command line.\n')
    for a in args:
        if not os.path.exists(a):
            parser.error('File %s does not exist.\n' % a)
    if options.dpi < 72:
        parser.error('--dpi %d less than screen res, 72. Must be >= 72.'
                     % options.dpi)
    if options.outFormat not in ('pdf', 'png', 'eps', 'all'):
        parser.error('Unrecognized output format: %s. Choose one from: '
                     'pdf png eps all.' % options.outFormat)
    if (options.out.endswith('.png') or options.out.endswith('.pdf') or
        options.out.endswith('.eps')):
        options.out = options.out[:-4]

def initImage(width, height, options):
    """
    initImage takes a width and height and returns
    both a fig and pdf object. options must contain outFormat,
    and dpi
    """
    pdf = None
    if options.outFormat == 'pdf' or options.outFormat == 'all':
        pdf = pltBack.PdfPages(options.out + '.pdf')
    fig = plt.figure(figsize=(width, height), dpi=options.dpi, facecolor='w')
    return (fig, pdf)

def writeImage(fig, pdf, options):
    """
    writeImage assumes options contains outFormat and dpi.
    """
    if options.outFormat == 'pdf':
        fig.savefig(pdf, format = 'pdf')
        pdf.close()
    elif options.outFormat == 'png':
        fig.savefig(options.out + '.png', format='png', dpi=options.dpi)
    elif options.outFormat == 'all':
        fig.savefig(pdf, format='pdf')
        pdf.close()
        fig.savefig(options.out + '.png', format='png', dpi=options.dpi)
        fig.savefig(options.out + '.eps', format='eps')
    elif options.outFormat == 'eps':
        fig.savefig(options.out + '.eps', format='eps')

def readFile(filename):
    with open(filename) as f:
        for filelen, l in enumerate(f, 1):
            pass
    d = numpy.zeros(filelen, dtype=numpy.uint64)
    with open(filename) as f:
        for i, line in enumerate(f):
            d[i] = int(line.strip())
    return d

def establishAxis(fig, options):
    """ create one axes per chromosome
    """
    options.axLeft = 0.1
    options.axWidth = 0.85
    options.axBottom = 0.15
    options.axHeight = 0.75
    ax = fig.add_axes([options.axLeft, options.axBottom,
                       options.axWidth, options.axHeight])
    return ax

def drawData(data, ax, options):
    colorList = ['#1f77b4', # dark blue
                 '#aec7e8', # light blue
                 '#ff7f0e', # bright orange
                 '#ffbb78', # light orange
                 '#4B4C5E', # dark slate gray
                 '#9edae5', # light blue
                 '#7F80AB', # purple-ish slate blue
                 '#c7c7c7', # light gray
                 '#9467bd', # dark purple
                 '#c5b0d5'  # light purple
                 ]
    labels = map(lambda x: x.name, data)
    lineStyleList = ['-', '--', '-.', ':']
    lineStyleIndex = -1
    ax.set_title(options.title)
    # create the N50 line
    globalMin = min(map(lambda x: x.lengths[-1], data))
    if options.n50Line:
        color50 = (0.4, 0.4, 0.4)
        for d in data:
            # horizontal lines
            ax.add_line(lines.Line2D(xdata=[0.0, 0.5],
                                     ydata=[nValue(d, 0.5),
                                            nValue(d, 0.5)],
                                     color=color50,
                                     linewidth=0.75,
                                     linestyle=':'))
    if options.reportN50Values:
        for d in data:
            print '%9s  #lens: %10d' % (d.name, d.n)
            print '%9s cumlen: %10d' % (d.name, numpy.sum(d.lengths))
            for n in [10, 50, 90, 95]:
                print '%9s %6s: %10d' % (d.name, 'n%d' % n, nValue(d, n / 100.0))

    for i, d in enumerate(data, 0):
        j = i % len(colorList)
        if j == 0 : lineStyleIndex = (lineStyleIndex + 1) % len(lineStyleList)
        ax.plot(d.xData, d.lengths,
                color = colorList[j],
                linestyle = lineStyleList[lineStyleIndex],
                label = labels[i])
    for loc, spine in ax.spines.iteritems():
        if loc in ['left','bottom']:
            spine.set_position(('outward', 10)) # outward by 10 points
        elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine
        else:
            raise ValueError('unknown spine location: %s' % loc)

    if options.log:
        ax.set_yscale('log')
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(1, 10)))
    plt.ylabel('Length')

    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.xaxis.set_ticklabels([0, '', '', '', '', 0.5, '', '', '', '', 1.0])
    # turn off ticks where there is no spine
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.xlabel(options.xlabel)
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(handles, labels)
    plt.setp(leg.get_texts(), fontsize='x-small') # legend fontsize
    leg._drawFrame = False

def nValue(a, x):
    """ given a dict populated as the processData returned dicts, `a'
    and a a float x in the range (0, 1.0) nValue returns the Nx value
    """
    if not isinstance(a, LengthObj):
        raise RuntimeError('Type of `a\' must be LengthObj, not %s'
                           % x.__class__)
    if not isinstance(x, float):
        raise RuntimeError('Type of `x\' must be float, not %s' % x.__class__)
    if not (0.0 < x < 1.0):
        raise RuntimeError('Value of `x\' must be in (0.0, 1.0), not %f' % x)
    n = numpy.sum(a.xData > x)
    if n == 0:
        return 0
    return a.lengths[- n]

def processData(lengthObjList, options):
    """ processData() takes the list of LengthObjs and sorts
    their lengths, if necessary,
    """
    if options.genomeLength is None:
        options.genomeLength = max(map(lambda x: x.xData[-1], lengthObjList))
    for l in lengthObjList:
        l.xData = numpy.divide(l.xData, float(options.genomeLength))

def main():
    usage = ('usage: %prog [options] lengths1.txt lengths2.txt lengths3.txt '
             '...\n\n'
             '%prog takes as many lengths files as you offer, some options\n'
             'and then produces a figure showing all N-values.')
    parser = OptionParser(usage=usage)
    initOptions(parser)

    options, args = parser.parse_args()
    checkOptions(options, args, parser)

    lengthObjList = []
    for a in args:
        # each input length file is transformed into a LengthObj
        lengthObjList.append(LengthObj(a, readFile(a), options.preSorted))

    processData(lengthObjList, options)

    fig, pdf = initImage(8.0, 5.0, options)
    ax = establishAxis(fig, options)
    drawData(lengthObjList, ax, options)
    writeImage(fig, pdf, options)

if __name__ == '__main__':
    main()
