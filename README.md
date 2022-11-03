# TnBox
## _A Tn-seq Toolbox_


TnBox is a user friendly toolbox for the analysis of high throughput mutant libraries (Tn-seq). 
It provides a set of tools allowing the user to parse and align the sequencing reads on a prokaryotic genome. It then implements two algorithms (Rslide & TnIF) for the study of gene essentiality.  

For support, questions or requests, please contact: francois-xavier.stubbe@unamur.be

## Installation & quick start

TnBox comes as a simple python program that works on both macOS and Linux systems. However, a few dependencies need to be installed and set in your PATH. Similarly, a couple python packages need to be downloaded. There are many ways to install these but the eaisest way might be to use the python package manager [Anaconda][conda]. Since Anaconda comes along a version of python, it is not necessary to have python pre-installed on your machine.

#### Install dependencies with conda on linux/macOS systems

- [Burrow Wheeler Aligner][bwa] 
- [Samtools][samtools]
- [Bedtools][bedtools]
- [Git][git]

Paste the following commands in a terminal

```sh
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c anaconda git
```
#### Install python packages with conda on linux/macOS systems

- [Numpy][python]
- [Pandas][bwa] 
- [Matplotlib][bwa] 
- [Searborn][samtools]

Paste the following commands in a terminal

```sh
conda install -c anaconda numpy
conda install -c anaconda pandas
conda install -c anaconda matplotlib
conda install -c anaconda searborn
```

#### Setting up and testing TnBox

##### Clone TnBox

To install TnBox, clone this repository using git : 
```sh
git clone https://github.com/fxstubbe/TnBox
```
Alternatively, you can also download TnBox using the green button at page's top.

##### Setting up TnBox

If you correctly installed all the required repositories and packages, TnBox is ready to use as it is. To start the program, you need to cd into TnBox folder and TnBox.py.

Assuming you've installed TnBox on the desktop; paste the following commands :

```sh
cd Desktop/TnBox/
python TnBox.py
```

On the first use, TnBox will automatically clone [BBMAP][bbmap] into the TnBox folder. This steps is crucial and if not executed properly will lead to failure. If an interface opens up : congrats ! TnBox is now woriking on your computer

## Analysis

TnBox interface is divided into 4 main panels : 
- 1. Process libraries
- 2. Get transposons insertion sites
- 3. Indexing
- 4. Explore

Each panel is designed to fulfill a step of the analysis pipeline

### 1. Process libraries

The first step into any high throughput sequencing project is to map the sequencing reads (.fastq.gz) onto a a genome of interest. 



## License

MIT

**Free Software, Hell Yeah!**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

[bwa]: https://sourceforge.net/projects/bio-bwa/files/
[samtools]: http://www.htslib.org/
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[python]: https://www.python.org/
[conda]: https://www.anaconda.com/
[git]: https://git-scm.com/

[pandas]: https://pandas.pydata.org/
[numpy]: https://numpy.org/
[plt]: https://matplotlib.org/
[sns]: https://seaborn.pydata.org/

[bbmap]: https://github.com/BioInfoTools/BBMap

