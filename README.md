# N50 Plotting Tools
Scripts to create N50 plots from files containing lengths of sequences

## Author
[Dent Earl](https://github.com/dentearl/)

## Dependencies
* python >= 2.6
* numpy
* matplotlib

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Usage
lengthsToN50Plot.py --title=TITLE lengths1.txt lengths2.txt lengths3.txt ...

lengthsToN50Plot.py takes as many lengths files as you offer, the size of the genome
(--size), a title (--title) and then produces an N50 style figure.

Options:
  -h, --help            show this help message and exit
  --genomeLength=GENOMELENGTH
                        Total length of the genome.
  --title=TITLE         Title of the plot.
  --log                 Puts y axis into log scale. default=False
  --n50Line             Adds straight lines from the y-axis to the curves. default=False
  --xlabel=XLABEL       Label on the x-axis. default=Cumulative length proportional to genome length
  --reportN50Values     prints n50 values to stdout. default=False
  --dpi=DPI             Dots per inch of the output, if --outFormat is all or png. default=300
  --outFormat=OUTFORMAT
                        output format [pdf|png|all|eps]. default=pdf
  --out=OUT             filename where figure will be created. No extension needed. default=myPlot

## Example
<code>
$ bin/fastaLengthSummarizer.py < genomeContigs.fa  > contigLengths.txt
$ bin/fastaLengthSummarizer.py < genomeScaffolds.fa  > scaffoldLengths.txt
$ bin/lengthsToN50Plot.py --title "My genome's N-stats" --n50Line --log contigLengths.txt scaffoldLengths.txt
</code>
