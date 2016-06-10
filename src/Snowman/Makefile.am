bin_PROGRAMS = snowman

snowman_CPPFLAGS = \
	-I$(top_srcdir)/src/SGA/Util \
	-I$(top_srcdir)/src/SGA/Bigraph \
	-I$(top_srcdir)/src/SGA/SuffixTools \
	-I$(top_srcdir)/src/SGA/StringGraph \
##	-I$(top_srcdir)/src/SGA/Concurrency \
	-I$(top_srcdir)/src/SGA/Algorithm \
	-I$(top_srcdir)/src/SGA/SQG \
##	-I$(top_srcdir)/src/SGA/Scaffold \
##	-I$(top_srcdir)/src/SGA/GraphDiff \
	-I$(top_srcdir)/src/SGA/SGA \
	-I$(top_srcdir)/SnowTools/src \
	-I$(top_srcdir)/SnowTools \
	-I$(top_srcdir)/SnowTools/htslib \
	-I$(top_srcdir)/SnowTools/multifast-v1.4.2 \
	-I$(top_srcdir)/SnowTools/blat/inc \
	-I$(top_srcdir)/SnowTools/blat/jkOwnLib -Wno-write-strings

snowman_LDADD = \
	$(top_builddir)/src/SGA/SGA/libsga.a \
	$(top_builddir)/SnowTools/bwa/libbwa.a \
	$(top_builddir)/SnowTools/htslib/libhts.a \
	$(top_builddir)/SnowTools/src/libsnowtools.a \
##	$(top_builddir)/src/SGA/Scaffold/libscaffold.a \
##	$(top_builddir)/src/SGA/GraphDiff/libgraphdiff.a \
	$(top_builddir)/src/SGA/StringGraph/libstringgraph.a \
##	$(top_builddir)/src/SGA/Concurrency/libconcurrency.a \
	$(top_builddir)/src/SGA/Algorithm/libalgorithm.a \
	$(top_builddir)/src/SGA/SuffixTools/libsuffixtools.a \
	$(top_builddir)/src/SGA/Bigraph/libbigraph.a \
	$(top_builddir)/src/SGA/Util/libutil.a \
	$(top_builddir)/src/SGA/SQG/libsqg.a \
	$(top_builddir)/SnowTools/multifast-v1.4.2/ahocorasick/build/libahocorasick.a \
	$(top_builddir)/SnowTools/blat/lib/jkOwnLib.a $(top_builddir)/SnowTools/blat/lib/jkweb.a 
##        @boost_lib@/libboost_regex.a @boost_lib@/libboost_filesystem.a @boost_lib@/libboost_system.a


snowman_LDFLAGS = -pthread -std=c++11

snowman_SOURCES = run_snowman2.cpp splitcounter.cpp BreakPoint2.cpp AlignedContig2.cpp DiscordantCluster.cpp DBSnpFilter.cpp PONFilter.cpp \
		SnowmanUtils.cpp snowman.cpp workqueue.h SnowmanAssemblerEngine.cpp vcf.cpp SnowmanRead.cpp \
		SnowmanOverlapAlgorithm.cpp SnowmanASQG.cpp \
		SnowmanAssemble.h SnowmanAssemble.cpp genpon.cpp \
		KmerFilter.cpp SnowmanBamWalker.cpp benchmark.cpp PonWalker.cpp ReadSim.cpp \
		SeqFrag.cpp SimGenome.cpp BamSplitter.cpp SimTrainerWalker.cpp refilter.cpp LearnBamParams.cpp \
		PowerLawSim.cpp AssemblyBamWalker.cpp assembly2vcf.cpp

