# A Tn-Seq toolbox (TnBox)

Transposon sequencing requires the creation of a transposon insertion library, which will contain a group of mutants that collectively have transposon insertions in all non-essential genes. The library is grown under the condition that is of interest. Mutants with transposons inserted in genes required for growth under the test condition will diminish in frequency from the population. To identify genes being lost, sequences encompassing the transposon ends are amplified by PCR and sequenced to determine the location and abundance of each insertion mutation.

## Requirements

**Installation**

- [Burrows Wheeler Aligner][bwa]
- [Samtools][samtools]
- [Qualimap][quali]
- R ([tidyverse][tidy], [data.table][d.table], [patchwork][patch])

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[tidy]: https://www.tidyverse.org/
[d.table]: https://github.com/Rdatatable/data.table
[patch]: https://github.com/thomasp85/patchwork
[quali]: http://qualimap.bioinfo.cipf.es/

**Input Data**
