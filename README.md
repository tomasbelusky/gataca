gataca
======
Detection of genome variations by combination of read pair and split read methods. Depth of coverage is also used for increase of accuracy. Variations are set into clusters for better sensitivity.
WORK IN PROGRESS.

requirements
------------

* Python libraries
    * [pysam](http://code.google.com/p/pysam/) - manipulating with reads and reference genome
    * [bx-python](https://bitbucket.org/james_taylor/bx-python/wiki/Home) - fast implementation of finding overlapped intervals in tree
* Tools
    * [bwa](http://bio-bwa.sourceforge.net) - remapping of clipped sequences in reads
