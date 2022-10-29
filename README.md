# A Tn-Seq toolbox (TnBox)

Transposon sequencing requires the creation of a transposon insertion library, which will contain a group of mutants that collectively have transposon insertions in all non-essential genes. The library is grown under the condition that is of interest. Mutants with transposons inserted in genes required for growth under the test condition will diminish in frequency from the population. To identify genes being lost, sequences encompassing the transposon ends are amplified by PCR and sequenced to determine the location and abundance of each insertion mutation.

## Installation

**Requirements**

- [Python] [python]
- [Burrows Wheeler Aligner][bwa]
- [Samtools][samtools]

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[python]: https://www.python.org/
[conda]: https://www.anaconda.com/

**Tutorial**

The easiest way to set up TnBox is to download the depedencies using the python package manager [Anaconda] [conda]. 
