bin_PROGRAMS = snowmanutils

snowmanutils_CPPFLAGS = \
	-I$(top_srcdir)/src/SGA/Util \
	-I$(top_srcdir)/src/SGA/Bigraph \
	-I$(top_srcdir)/src/SGA/SuffixTools \
	-I$(top_srcdir)/src/SGA/StringGraph \
	-I$(top_srcdir)/src/SGA/Algorithm \
	-I$(top_srcdir)/src/SGA/SQG \
	-I$(top_srcdir)/src/SGA/SGA \
	-I$(top_srcdir)/src/Snowman \
	-I$(top_srcdir)/SeqLib \
	-I$(top_srcdir)/SeqLib/htslib

snowmanutils_LDADD = \
	$(top_builddir)/src/SGA/SGA/libsga.a \
	$(top_builddir)/src/SGA/StringGraph/libstringgraph.a \
	$(top_builddir)/src/SGA/Algorithm/libalgorithm.a \
	$(top_builddir)/src/SGA/SuffixTools/libsuffixtools.a \
	$(top_builddir)/src/SGA/Bigraph/libbigraph.a \
	$(top_builddir)/src/SGA/Util/libutil.a \
	$(top_builddir)/src/SGA/SQG/libsqg.a \
	$(top_builddir)/SeqLib/src/libseqlib.a \
	$(top_builddir)/SeqLib/bwa/libbwa.a \
	$(top_builddir)/SeqLib/htslib/libhts.a \
	$(top_builddir)/SeqLib/fermi-lite/libfml.a

snowmanutils_LDFLAGS = -pthread -std=c++11

snowmanutils_SOURCES = snowutils.cpp splitcounter.cpp ../Snowman/BreakPoint2.cpp ../Snowman/AlignedContig2.cpp \
		../Snowman/DiscordantCluster.cpp ../Snowman/DBSnpFilter.cpp ../Snowman/PONFilter.cpp \
		../Snowman/SnowmanUtils.cpp ../Snowman/SnowmanAssemblerEngine.cpp ../Snowman/vcf.cpp ../Snowman/DiscordantRealigner.cpp \
		../Snowman/SnowmanOverlapAlgorithm.cpp ../Snowman/SnowmanASQG.cpp \
		../Snowman/SnowmanAssemble.cpp \
		../Snowman/KmerFilter.cpp ../Snowman/SnowmanBamWalker.cpp  \
		SeqFrag.cpp SimGenome.cpp BamSplitter.cpp SimTrainerWalker.cpp \
		Fractions.cpp ../Snowman/STCoverage.cpp ../Snowman/Histogram.cpp ../Snowman/BamStats.cpp \
		PowerLawSim.cpp AssemblyBamWalker.cpp assembly2vcf.cpp
