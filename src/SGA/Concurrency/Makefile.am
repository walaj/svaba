noinst_LIBRARIES = libconcurrency.a

libconcurrency_a_CPPFLAGS = \
	-I$(top_srcdir)/SGA/Bigraph \
	-I$(top_srcdir)/SGA/Util \
	-I$(top_srcdir)/SGA/SuffixTools \
	-I$(top_srcdir)/SGA/Thirdparty \
	-I$(top_srcdir)/SGA/Algorithm

libconcurrency_a_SOURCES = \
        OverlapProcess.h OverlapProcess.cpp \
        RmdupProcess.h RmdupProcess.cpp \
        SequenceProcessFramework.h \
        SequenceWorkItem.h \
        ThreadWorker.h \
		MkqsThread.h
