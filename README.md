[![Build Status](https://travis-ci.org/walaj/svaba.svg?branch=master)](https://travis-ci.org/walaj/svaba)

## *SvAbA* - Structural variation and indel detection by assembly

[Pronounced Sah-Bah or Svah-bah. This project was formerly "Snowman"]

**License:** [GNU GPLv3][license] 

Table of contents
=================

  * [Installation](#gh-md-toc)
  * [Description](#description)
  * [Output file description](#output-file-description)
  * [Filtering and refiltering](#filtering-and-refiltering)
  * [Recipes and examples](#recipes-and-examples)
    * [Whole genome somatic SV and indel detection](#whole-genome-somatic-sv-and-indel-detection)
    * [Whole genome germline SV and indel detection](#whole-genome-germline-sv-and-indel-detection)
    * [Targeted (exome) detection](#targeted-detection)
    * [Targeted local assembly](#targeted-local-assembly)
    * [Assemble all reads](#assemble-all-reads)
    * [Runtime snapshot](#snapshot-of-where-svaba-run-is-currently-operating)
    * [Debug a local assembly and produce assembly graph](#debug-a-local-assembly-and-produce-the-assembly-graph)
    * [View all of the ASCII alignments](#view-all-of-the-ascii-alignments)
    * [View a particular contig with read-to-contig alignments](#view-a-particular-contig-with-read-to-contig-alignments)
    * [Make a function to sort and index the contigs](#make-a-function-to-sort-and-index-contigs)
  * [Attributions](#attributions)


Installation
------------
We recommend compiling with GCC-4.8 or greater. We have successfully compiled on RHEL6, CentOS with GCC-4.8, 4.9 and 5.1, and MacOSX with Clang (Apple LLVM version 7.0.2 (clang-700.1.81))

```
git clone --recursive https://github.com/walaj/svaba
cd svaba
./configure
make
make install

## QUICK START (eg run tumor / normal on Chr22, with 4 cores)
bin/svaba -t tumor.bam -n normal.bam -k 22 -G ref.fa -a test_id -p -4

## get help
svaba --help
svaba run --help
```

SvAbA uses the [SeqLib][seqlib] API for BAM access, BWA-MEM alignments, interval trees and operations,
and several other auxillary operations.

Description
-----------

SvAbA is a method for detecting structural variants in sequencing data using genome-wide local assembly. Under the hood, 
SvAbA uses a custom implementation of SGA (String Graph Assembler) by Jared Simpson, and BWA-MEM by Heng Li. Contigs are assembled
for every 25kb window (with some small overlap) for every region in the genome. The default is to use only clipped, discordant, 
unmapped and indel reads, although this can be customized to any set of reads at the command line using [VariantBam][vbam] rules. 
These contigs are then immediately aligned to the reference with BWA-MEM and parsed to identify variants. Sequencing reads are likewise 
realigned to the contigs with BWA-MEM, and variants are scored 

SvAbA is currently configured to provide indel and rearrangement calls (and anything "in between"). It can jointly call any number of BAM/CRAM/SAM files,
and has built-in support for case-control experiments (e.g. tumor/normal, or trios or quads). In case/control mode, 
any number of cases and controls (but min of 1 case) can be input, and 
will jointly assemble all sequences together. If both a case and control are present, variants are output separately in "somatic" and "germline" VCFs. 
If only a single BAM is present (input with the ``-t`` flag), a single SV and a single indel
VCF will be emitted.

A BWA-MEM index reference genome must also be supplied with ``-G``.

<img src="https://github.com/broadinstitute/SvAbASV/blob/master/gitfig_schematic.png"
width=800/>

Output file description
-----------------------

##### ``*.bps.txt.gz``
Raw, unfiltered variants. This file is parsed at the end to produce the VCF files. With the bps.txt.gz,
one can define a new set of filteirng criteria (depending on sensitivity/specificity needs) using ``svaba refilter``. 

##### ``*.contigs.bam``
All assembly contigs as aligned to the reference with BWA-MEM. Note that this is an unsorted file. To view in IGV,
it must be first sorted and indexed (e.g. ``samtools sort -m 8G id.contigs.bam id.sort && samtools index id.sort.bam``)

##### ``*.discordants.txt.gz``
Information on all clusters of discordant reads identified with 2+ reads. 

##### ``*.log``
Log file giving run-time information, including CPU and Wall time (and how it was partitioned among the tasks), number of 
reads retrieved and contigs assembled for each region.

##### ``*.alignments.txt.gz``
An ASCII plot of variant-supporting contigs and the BWA-MEM alignment of reads to the contigs. This file is incredibly
useful for debugging and visually inspecting the exact information SvAbA saw when it performed the variant-calling. This file
is typically quite large. The recommended usage is to identify the contig name of your variant of interest first from the VCF file 
(SCTG=contig_name). Then do ``gunzip -c id.alignment.txt.gz | grep contig_name > plot.txt``. It is highly recommended that you 
view in a text editor with line truncation turned OFF, so as to not jumble the alignments.

<img src="https://github.com/broadinstitute/SvAbASV/blob/master/gitfig_ascii.png"
width=800/>

##### ``*.bad_mate_regions.bed``
A BED file of regions that were suspected of having poor alignment quality. When encountered, these regions are excluded from future
mate-lookups. This prevents excessive mate-lookups to centromeres, etc.

##### ``*.vcf``
VCF of rearrangements and indels parsed from bps.txt.gz and with a somatic_score == 1 (somatic) or 0 (germline) and quality == PASS. *NOTE* that 
the cutoff for rearrangement vs indel is taken from BWA-MEM, whether it produces a single gapped-alignment 
or two separate alignments. This is an arbitrary cutoff, just as there is no clear consensus distinction between what 
constitutes an "indel" and a "structural variant". The unfiltered VCF files include non-PASS variants. 

Filtering and Refiltering
-----------------------

SvAbA performs a series of log-likelihood calculations for each variant. The purpose is to first classify a variant as real vs artifact, 
and then to determine if the variant is somatic or germline. These log-likelihoods are output in the VCF and bps.txt.gz file and described here:
* ``LOD (LO)`` - Log of the odds that variant is real vs artifact. For indels, the likelihood of an artifact read is proportional to the length of local repeats (repeating units up to 5 long per unit)
* ``LR`` - Log of the odds that the variant has allelic fraction (AF) of 0 or >=0.5. This is used for somatic vs germline classification
* ``SL`` - Scaled LOD. LOD scores is heuristically scaled as: (min(Mapping quality #1, Mapping quality #2) - 2 * NM) / 60 * LOD

SvAbA can refilter the bps.txt.gz file to produce new VCFs with different stringency cutoffs. To run, the following are required:
* ``-b`` - a BAM from the original run, which is used just for its header
* ``-i`` - input bps.txt.gz file

Examples and recipes
--------------------

#### Whole genome somatic sv and indel detection 
```
wget "https://data.broadinstitute.org/snowman/dbsnp_indel.vcf" ## get a DBSNP known indel file
DBSNP=dbsnp_indel.vcf
CORES=8 ## set any number of cores
REF=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta
## -a is any string you like, which gives the run a unique ID
svaba run -t $TUM_BAM -n $NORM_BAM -p $CORES -D $DBSNP -a somatic_run -G $REF
```

#### Whole genome germline sv and indel detection
```
## Set -I to not do mate-region lookup if mates are mapped to different chromosome.
##   This is appropriate for germline-analysis, where we don't have a built-in control
##   to against mapping artifacts, and we don't want to get bogged down with mate-pair
##   lookups.
## Set -L to 6 which means that 6 or more mate reads must be clustered to 
##   trigger a mate lookup. This also reduces spurious lookups as above, and is more 
##   appropriate the expected ALT counts found in a germline sample 
##   (as opposed to impure, subclonal events in cancer that may have few discordant reads).
svaba run -t $GERMLINE_BAM -p $CORES -L 6 -I -a germline_run -G $REF
```

#### Targeted detection
```
## eg targets.bed is a set of exome capture regions
svaba run -t $BAM -k targets.bed -a exome_cap -G $REF
```

#### Targeted local assembly
```
## -k can be a chromosome, a samtools/IGV style string 
##     (e.g. 1:1,000,000-2,000,000), or a BED file
k=chr17:7,541,145-7,621,399
svaba run -t $TUM_BAM -n $NORM_BAM -p $CORES -k $k  -a TP53 -G $REF
```

#### Assemble all reads
```
## default behavior is just assemble clipped/discordant/unmapped/gapped reads
## This can be overridden with -r all flag
svaba run -t $BAM -r all -g $REF
```

#### Snapshot of where svaba run is currently operating
```
tail somatic_run.log
```

#### Debug a local assembly and produce the assembly graph
```
k=chr17:7,541,145-7,621,399
svaba run -t $BAM -a local_test -k $k --write-asqg

## plot the graph
$GIT/SvAbASV/R/snow-asqg.R
```

#### View all of the ASCII alignments 
```
## Make a read-only and no-line-wrapping version of emacs.
## Very useful for *.alignments.txt.gz files
function ev { 
  emacs $1 --eval '(setq buffer-read-only t)' -nw --eval '(setq truncate-lines t)';
  }
ev somatic_run.alignments.txt.gz 
```

#### View a particular contig with read to contig alignments
```
gunzip -c somatic_run.alignments.txt.gz | grep c_1_123456789_123476789 > c_1_123456789_123476789.alignments.txt
ev c_1_123456789_123476789.alignments.txt
```

#### Make a function to sort and index contigs
```
function sai() {
  if [[ -f $1.contigs.bam ]]; then
     samtools sort -m 4G $1.contigs.bam -o $1.contigs.sort.bam
     mv $1.contigs.sort.bam $1.contigs.bam
     samtools index $1.contigs.bam
  fi
}
## for example, for somatic_run.contigs.bam:
sai somatic_run
```


Attributions
============

SvAbA is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org) --  Rameen Berkoukhim lab -- Dana Farber Cancer Institute, Boston, MA. 

This project was developed in collaboration with the Cancer Genome Analysis team at the Broad Institute. Particular thanks to:
* Cheng-Zhong Zhang (Matthew Meyerson Lab)
* Marcin Imielinski (http://vivo.med.cornell.edu/display/cwid-mai9037)

Additional thanks to Jared Simpson for SGA, Heng Li for htslib and BWA, and for the other developers whose  
code contributed to [SeqLib](https://github.com/walaj/SeqLib).

[vbam]: https://github.com/walaj/VariantBam

[license]: https://github.com/walaj/svaba/blob/master/LICENSE

[seqlib]: https://github.com/walaj/SeqLib
