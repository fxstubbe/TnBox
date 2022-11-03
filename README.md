# TnBox
## _A Tn-seq Toolbox_


TnBox is a user friendly toolbox for the analysis of high throughput mutant libraries (Tn-seq). 
It provides a set of tools allowing the user to parse and align the sequencing reads on a prokaryotic genome. It then implements two algorithms (Rslide & TnIF) for the study of gene essentiality.  


## Installation
The easiest way to set up TnBox is to download the depedencies using the python package manager [Anaconda][conda]. Each Anaconda release comes with its version of python. Please, download python > 3
##### Dependencies

- [Burrow Wheeler Aligner][bwa] 
- [Samtools][samtools]
- [Git][git]

```sh
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c anaconda git
```
##### Python packages

Similarly, using conda to manage python
- [Numpy][python]
- [Pandas][bwa] 
- [Matplotlib][bwa] 
- [Searborn][samtools]

```sh
conda install -c anaconda numpy
conda install -c anaconda pandas
conda install -c anaconda matplotlib
conda install -c anaconda searborn
```



## Analysis


## License

MIT

**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[python]: https://www.python.org/
[conda]: https://www.anaconda.com/
[git]: https://git-scm.com/

[pandas]: https://pandas.pydata.org/
[numpy]: https://numpy.org/
[plt]: https://matplotlib.org/
[sns]: https://seaborn.pydata.org/



