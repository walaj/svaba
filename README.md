[![Build Status](https://travis-ci.org/broadinstitute/SnowmanSV.svg?branch=master)](https://travis-ci.org/broadinstitute/SnowmanSV)

Snowman - Structural Variation Detection by Rolling Local Assembly
==================================================================

**License:** [GNU GPLv3][license] 

Installation
------------
We recommend compiling with GCC-4.8 or greater. We have successfully compiled on RHEL6, CentOS with GCC-4.8, 4.9 and 5.1, and MacOSX with Clang (Apple LLVM version 7.0.2 (clang-700.1.81))

```
git clone --recursive https://github.com/broadinstitute/SnowmanSV
cd SnowmanSV
./configure
make

############### QUICK START ############### 
# run tumor / normal on Chr22, with 4 cores
SnowmanSV/src/Snowman/snowman -t tumor.bam -n normal.bam -k 22 -G ref.fa -a test_id -p -4

## get help
snowman --help
snowman run --help
```

Description
-----------

Snowman is a method for detecting structural variants in sequencing data using genome-wide local assembly. Under the hood, 
Snowman uses a custom implementation of SGA (String Graph Assembler) by Jared Simpson, and BWA-MEM by Heng Li. Contigs are assembled
for every 10kb window (with some small overlap) for every region in the genome. The default is to use only clipped, discordant, 
unmapped and indel reads, although this can be customized to any set of reads at the command line using [VariantBam][vbam] rules. 
These contigs are then immediately aligned to the reference with BWA-MEM and parsed to identify variants. Sequencing reads are likewise 
realigned to the contigs with BWA-MEM, and variants are scored 

Scope and Inputs
----------------

Snowman is currently configured to provide indel and rearrangement calls (and anything "in between"). It is setup to joint calling of any number of BAM/CRAM/SAM files,
and with built-in support for case-control experiments (e.g. tumor/normal, or trios or quads). In case/control mode, 
any number of cases and controls (but min of 1 case) can be input, and asseembly
will jointly assemble them all. If both a case and control are present, variants are output separately in "somatic" and "germline" VCFs. 
If only a single BAM is present, input with the ``-t`` flag (case). 
In this case, the results will contain all calls, with no germline/somatic designation.

A BWA-MEM index reference genome must also be supplied with ``-G``.

Output file description
-----------------------

##### ``*.bps.txt.gz``
Raw, unfiltered variants. This file is parsed at the end to produce the VCF files. With the bps.txt.gz,
one can define a new set of filteirng criteria (depending on sensitivity/specificity needs) using ``snowman refilter``. 

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
useful for debugging and visually inspecting the exact information Snowman saw when it performed the variant-calling. This file
is typically quite large. The recommended usage is to identify the contig name of your variant of interest first from the VCF file 
(SCTG=contig_name). Then do ``gunzip -c id.alignment.txt.gz | grep contig_name > plot.txt``. It is highly recommended that you 
view in a text editor with line truncation turned OFF, so as to not jumble the alignments.

##### ``*.bad_mate_regions.bed``
A BED file of regions that were suspected of having poor alignment quality. When encountered, these regions are excluded from future
mate-lookups. This prevents excessive mate-lookups to centromeres, etc.

##### ``*.vcf``
VCF of rearrangements and indels parsed from bps.txt.gz and with a somatic_score == 1 (somatic) or 0 (germline) and quality == PASS. *NOTE* that 
the cutoff for rearrangement vs indel is taken from BWA-MEM, whether it produces a single gapped-alignment 
or two separate alignments. This is an arbitrary cutoff, just as there is no clear consensus distinction between what 
constitutes an "indel" and a "structural variant". The unfiltered VCF files include non-pass variants. 

Filtering / Refiltering
-----------------------

Snowman performs a series of log-likelihood calculations for each variant. The purpose is to first classify a variant as real vs artifact, 
and then to determine if the variant is somatic or germline. These log-likelihoods are output in the VCF and bps.txt.gz file and described here:
* ``LOD (LO)`` - Log of the odds that variant is real vs artifact. For indels, the likelihood of an artifact read is proportional to the length of local repeats (repeating units up to 5 long per unit)
* ``LR`` - Log of the odds that the variant has allelic fraction (AF) of 0 or >=0.5. This is used for somatic vs germline classification
* ``SL`` - Scaled LOD. LOD scores is heuristically scaled as: (min(Mapping quality #1, Mapping quality #2) - 2 * NM) / 60 * LOD

Snowman can refilter the bps.txt.gz file to produce new VCFs with different stringency cutoffs. To run, the following are required:
* ``-b`` - a BAM from the original run, which is used just for its header
* ``-i`` - input bps.txt.gz file


[vbam]: https://github.com/jwalabroad/VariantBam

[license]: https://github.com/broadinstitute/SnowmanSV/blob/master/LICENSE
