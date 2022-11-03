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

![](https://github.com/fxstubbe/TnBox/blob/main/Images/TnBox.png)

### 1. Process libraries

The first step int any high throughput sequencing project is to map the sequencing reads (.fastq.gz) onto a a genome of interest. To do so, reads must first be trimmed of adapter and quality trimmed. Then, the processed reads are aliged onto the genome using an aligner. TnBox uses [BBMAP][bbmap] for the filtering and [bwa][bwa] for the alignment. 

#### Transposon fishing

Depending on the sequencing platform, your sequencing read will start right after the inserted transposon or will include the transposon. If the latter is true, tick the Parse miniTn5 tickbox. This will allow TnBox to map only transposon containing reads. Skipping this 

#### Select or Add a reference

The first thing to do is import a genome of interest. Genomes can be downloaded from NCBI. For example, you can download Brucella abortus genome [here][abortus]. A few common references are provided with TnBox. Next time you'll use TnBox, the references you added will still be available.


#### Add sequencing files (.fastq.gz)

Add your files. Your sequecing reads must be with the file extension .fastq.gz otherwise TnBox will fail. Tn-seq are usually sequence in single-end but if your data has been sequenced i paired-end, then only include the forward reads (R1) as this is where the transposon lies. Similarly, sometimes your sequencing reads will be split over mutiple file. You can simply concaten them together as follow : 

```sh
cat file_1 file_2 > concatenated_file
```

#### Start the anaysis

Click on start. That's it.

This process is time consuming. It largely depends on the amount of reads, the genome size ... but also your machine performance. Be patient, make yourself a cofee and come back later. TnBox provides a visual cue of where it is in the process. When the library name is orange, it's under process. When it turns green, TnBox is done with this library.

### 2. Get transposons insertion sites 

#### Choose or Add a reference (.gff)
Now that the reads have been aligned to the genome of interest, it's time extract where the transposons are. The first step is to select the appropirate reference or the add your reference of choice in gff format. For example, you can dowload Brucella abortus reference [here][abortus]. 

#### Define metrics for the algorithms

**Trim end** 

Since many genes can tolerate insertions in their 5' or 3', it is wise to remove those extremities while looking for insertion sites. By default, TnBox trims 10% of the transcript length on both 5' and 3' ends. 
This parameter is used for both the Rslide and TnIF algorithms (see below)

**R window**, **Slide**

While using the RSlide algorithm, a sliding window of size n (R window) is moved with an increment p (slide) over the genome. By default, the R window is set to be 100 as described in [Potemberg et al.,][] [JF et al.,][]. 

#### Select files and algorithm

When processing the libraries, TnBox generated 2 types of files :
i. TA files, where only 5' end of each mapped read has been counted. This is equivalent to the exact insertion site. This is the recommended method when looking for essential genes.
ii. TnIF files, where the whole mapped read has been counted. Even if this method is less sensitive to detect essential genes, it performs better when using the TnIF algorithm which has for goal to detect essentiality variations across several conditions as described in [Potemberg et al., ][] 


### 3. Indexing

### 4. Explore


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

[abortus]:https://www.ncbi.nlm.nih.gov/genome/520
