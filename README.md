[![Build Status](https://magnum.travis-ci.com/broadinstitute/SnowmanSV.svg?token=QTnp48gNXtKQKRDpquf3&branch=master)](https://magnum.travis-ci.com/broadinstitute/SnowmanSV)

Snowman - Structural Variation Detection by Rolling Local Assembly
==================================================================

**License:** [GNU GPLv3][license] 

Installation
------------
We recommend compiling with GCC-4.8 or greater. We have successfully compiled on RHEL6, CentOS with GCC-4.8, 4.9 and 5.1, and MacOSX with Clang (Apple LLVM version 7.0.2 (clang-700.1.81))

```
### if on Broad Institute servers, add GCC-4.9
reuse -q GCC-4.9

############## DOWNLOAD BOOST (if not installed) ###############
wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
tar -xvzf boost_1_61_0.tar.gz
## we only user header-only libraries, so no compiling of Boost is needed

############### DOWNLOAD VARIANT BAM ############### 
git clone --recursive https://github.com/broadinstitute/SnowmanSVgit
cd SnowmanSV

############### COMPILE AND INSTALL ###############
./configure --with-boost=<path_to_boost>  ## e.g. ~/boost_1_61_0
make

############### QUICK START ############### 
# run tumor / normal on Chr22, with 4 cores
SnowmanSV/src/Snowman/snowman -t tumor.bam -n normal.bam -k 22 -G ref.fa -a test_id -p -4

## get help
snowman --help
```

Description
-----------

Snowman is a method for detecting structural variants in sequencing data using genome-wide local assembly. Under the hood, 
Snowman uses a custom implementation of SGA (String Graph Assembler) by Jared Simpson, and BWA-MEM by Heng Li. Contigs are assembled
for every 10kb window (with some small overlap) for every region in the genome. The default is to use only clipped, discordant, 
unmapped and indel reads, although this can be customized to any set of reads at the command line using [VariantBam][vbam] rules. 
These contigs are then immediately aligned to the reference with BWA-MEM and parsed to identify variants. Sequencing reads are likewise 
realigned to the contigs with BWA-MEM, and variants are scored 

Scope
-----

Snowman is currently configured to provide indel and rearrangement calls (and anything "in between"). It has been most widely tested
as a somatic variant caller, and outputs separate VCFs for somatic and germline. If only a single BAM is present, input with the ``-t`` flag. 
In this case, the results will contain all calls, with no germline/somatic designation.

Required Inputs
---------------

Any number of BAM/SAM/CRAM files can be supplied at once. Snowman uses random access of the BAMs to obtain pair-mate reads,
and so requires the files to be indexed (``samtools index``). Tumor BAMs are input with ``-t`` and normal with ``-n``. The order
does not matter. At least one "tumor" BAM is required, although this could be just a single germline sample, or paired with a set of parents input
with ``-n`` to search for de novo alterations. All assemblies are done
jointly, combining reads across all of the BAM files. The source of the variant support reads is then tracked during realignment of reads to 
assembly contigs. A BWA indexed reference genome must be supplied as well (``-G``). 

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

##### ``*.cigarmap.txt.gz``
Information on clusters of indel alignments as identified from the original BAMs.

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

Auxillary Tools
---------------

# ``snowman benchmark``

Snowman ships with tools for testing the assemblies and for generating in-silico tumor genomes for testing.

##### ``--sim-breaks-power``
Simulates contigs containing structural variants.

```
### Simulate SVs and indels from the reference.
-R <num_rearrangements> -X <num_indels> --add-scrambled-inserts

```


##### ``--split-bam``
Simple tool for splitting a BAM file randomly into smaller fragments. This is useful for generating a sub-sampled BAM, 
and ensuring that each of the subsampled BAMs have different reads in them. ``split-bam`` will use the read name (and a seed 
provied with ``-s``) to generate the random numbers, which ensures that read-pairs are always sent to the same BAM.

```
### split BAM into 25% and 50% pieces at certain regions
snowman benchmark --split-bam -b <in.bam> -f 0.2,0.5 -k <regions.bed>
### split the entire BAM
snowman benchmark --split-bam -b <in.bam> -f 0.2,0.5 
```

##### --realign-test
Test the abiltiy of BWA-MEM to realign contigs of different lengths, with errors. Useful for gauging how often true indels will be missed 
during the BWA-MEM realignment phase.
```
#### 
snowman benchmark --realign-test
```

##### --realign-sv-test
Test the abiltiy of BWA-MEM to realign SV contigs of different lengths, with errors. Useful for gauging how often true rearrangements will be missed 
during the BWA-MEM realignment phase.
```
#### 
snowman benchmark --realign-sv-test
```

[vbam]: https://github.com/jwalabroad/VariantBam

[license]: https://github.com/broadinstitute/variant-bam/blob/master/LICENSE
