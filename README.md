# N50 Plotting Tools
Scripts to create N50 plots from files containing lengths of sequences

## Author
[Dent Earl](https://github.com/dentearl/)

## Dependencies
* python >= 2.6
* [matplotlib](http://matplotlib.sourceforge.net/)

## Installation
1. Download the package.
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

## Usage
<code>$ lengthsToN50Plot.py [options] lengths1.txt lengths2.txt lengths3.txt ... </code>

lengthsToN50Plot.py takes as many lengths files as you offer and then produces an N50 style figure.


Options:

* <code>-h, --help</code>            show this help message and exit
* <code>--genomeLength=GENOMELENGTH</code>Total length of the genome. default=Max of Sum of input file lengths.
* <code>--title=TITLE</code>         Title of the plot.
* <code>--log</code>                 Puts y axis into log scale. default=False
* <code>--n50Line</code>             Adds straight lines from the y-axis to the curves. default=False
* <code>--xlabel=XLABEL</code>       Label on the x-axis. default=Cumulative length proportional to genome length
* <code>--reportN50Values</code>     prints n50 values to stdout. default=False
* <code>--dpi=DPI</code>             Dots per inch of the output, if --outFormat is all or png. default=300
* <code>--outFormat=OUTFORMAT</code> output format [pdf|png|all|eps]. default=pdf
* <code>--out=OUT</code>             path/filename where figure will be created. No extension needed. default=myPlot


## Example
<code>$ bin/fastaLengthSummarizer.py < genomeContigs.fa  > contigLengths.txt</code>
<code>$ bin/fastaLengthSummarizer.py < genomeScaffolds.fa  > scaffoldLengths.txt</code>
<code>$ bin/lengthsToN50Plot.py --title "My genome's N-stats" --n50Line --log contigLengths.txt scaffoldLengths.txt</code>

