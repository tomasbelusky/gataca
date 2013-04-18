gataca
======
Detection of genome variations by combination of read pair and split read methods. Depth of coverage is also used for increase of accuracy. Variations are set into clusters for better sensitivity.
WORK IN PROGRESS.

manual
------
    Usage: gataca.py [OPTIONS] <sample.bam> <reference.fasta>

    Options:
      -h, --help                          show this help message and exit

      Input/output:
        -r STR, --region=STR              specify region (chr:from-to) of your interest, default: whole genome
        -o STR, --output=STR              name of output VCF file, default: standard output

      Reads:
        -p STR, --policy=STR              set how reads were sequenced (fr, rf) [fr]
        -q INT, --min_quality=INT         minimal mapping Phred quality score of read [30]
        -l INT, --min_length=INT          minimal length of split part [10]

      Depth of coverage:
        -w INT, --window_size=INT         size of window for getting coverage [100]
        -c STR, --coverage=STR            interval (min,max) of accepted coverage in windows, default: estimate from reads
        -a FLOAT, --coverage_core=FLOAT   core of windows from which minimal and maximal allowed coverage will be estimated [0.1]
        -u INT, --min_coverage_count=INT  minimal number of windows in core [1000]

      Insert size:
        -i STR, --insert_size=STR         interval (min,max) of accepted size between reads, default: estimate from reads
        -n INT, --insert_reads=INT        number of reads from which insert size will be estimated [50000]
        -e FLOAT, --insert_core=FLOAT     core of reads_num from which minimal and maximal insert size will be estimated [0.1]
        -m INT, --min_insert_count=INT    minimal number of reads in core [1000]

      Variations:
        -v FLOAT, --min_confidence=FLOAT   minimal confidence about variation [0.3]

requirements
------------

 * Python libraries
   * [pysam](http://code.google.com/p/pysam/) - manipulating with reads and reference genome
   * [bx-python](https://bitbucket.org/james_taylor/bx-python/wiki/Home) - fast implementation of finding overlapped intervals in tree
 * Tools
   * [bwa](http://bio-bwa.sourceforge.net) - remapping of clipped sequences in reads
