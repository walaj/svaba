bin_PROGRAMS = snowman

snowman_CPPFLAGS = \
	-I$(top_srcdir)/SGA/Util \
	-I$(top_srcdir)/SGA/Bigraph \
	-I$(top_srcdir)/SGA/SuffixTools \
	-I$(top_srcdir)/SGA/StringGraph \
	-I$(top_srcdir)/SGA/Concurrency \
	-I$(top_srcdir)/SGA/Algorithm \
	-I$(top_srcdir)/SGA/SQG \
	-I$(top_srcdir)/SGA/Scaffold \
	-I$(top_srcdir)/SGA/GraphDiff \
	-I$(top_srcdir)/SGA/SGA \
	-I$(top_srcdir)/SGA/Thirdparty \
	-I$(top_srcdir)/SnowUtils

snowman_LDADD = \
	$(top_builddir)/SGA/SGA/libsga.a \
	$(top_builddir)/SGA/Scaffold/libscaffold.a \
	$(top_builddir)/SGA/GraphDiff/libgraphdiff.a \
	$(top_builddir)/SGA/StringGraph/libstringgraph.a \
	$(top_builddir)/SGA/Concurrency/libconcurrency.a \
	$(top_builddir)/SGA/Algorithm/libalgorithm.a \
	$(top_builddir)/SGA/SuffixTools/libsuffixtools.a \
	$(top_builddir)/SGA/Bigraph/libbigraph.a \
	$(top_builddir)/SGA/Util/libutil.a \
	$(top_builddir)/SGA/SQG/libsqg.a \
	$(top_builddir)/SGA/Thirdparty/libthirdparty.a \
	$(top_builddir)/SnowUtils/libsnowutils.a 


snowman_LDFLAGS = -pthread -std=c++11

snowman_SOURCES = snowman.cpp workqueue.h grabReads.h grabReads.cpp \
                gather.h gather.cpp AlignedContig.h AlignedContig.cpp \
		kseq.h VCFRecord.h vcf.h vcf.cpp \
		nozzle.h nozzle.cpp zfstream.h zfstream.cpp bamprocess.h bamprocess.cpp \
		SnowmanOverlapAlgorithm.h SnowmanOverlapAlgorithm.cpp SnowmanASQG.h SnowmanASQG.cpp \
		SnowmanAssemble.h SnowmanAssemble.cpp 
