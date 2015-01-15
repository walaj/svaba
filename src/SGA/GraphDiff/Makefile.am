noinst_LIBRARIES = libgraphdiff.a

libgraphdiff_a_CPPFLAGS = \
	-I$(top_srcdir)/SGA/Bigraph \
	-I$(top_srcdir)/SGA/Util \
	-I$(top_srcdir)/SGA/StringGraph \
	-I$(top_srcdir)/SGA/SQG \
	-I$(top_srcdir)/SGA/SuffixTools \
	-I$(top_srcdir)/SGA/Thirdparty \
	-I$(top_srcdir)/SGA/Algorithm \
    -I$(top_srcdir)/SGA/Concurrency

libgraphdiff_a_SOURCES = \
        DindelUtil.h DindelUtil.cpp \
        DindelRealignWindow.h DindelRealignWindow.cpp \
        DindelHMM.h DindelHMM.cpp \
        OverlapHaplotypeBuilder.h OverlapHaplotypeBuilder.cpp \
        StringHaplotypeBuilder.h StringHaplotypeBuilder.cpp \
        DeBruijnHaplotypeBuilder.h DeBruijnHaplotypeBuilder.cpp \
        PairedDeBruijnHaplotypeBuilder.h PairedDeBruijnHaplotypeBuilder.cpp \
        HapgenUtil.h HapgenUtil.cpp \
        VCFTester.h VCFTester.cpp \
        GraphCompare.h GraphCompare.cpp \
        OverlapExtractorWithCorrection.h OverlapExtractorWithCorrection.cpp \
        DBGPathGuide.h DBGPathGuide.cpp \
        IOverlapExtractor.h
