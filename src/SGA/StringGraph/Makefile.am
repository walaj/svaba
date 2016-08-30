noinst_LIBRARIES = libstringgraph.a

libstringgraph_a_CPPFLAGS = \
	-I$(top_srcdir)/src/SGA/Bigraph \
	-I$(top_srcdir)/src/SGA/Util \
	-I$(top_srcdir)/src/SGA/Thirdparty \
	-I$(top_srcdir)/src/SGA/Algorithm \
	-I$(top_srcdir)/src/SGA/SQG \
	-I$(top_srcdir)/SeqLib

libstringgraph_a_SOURCES = \
        SGUtil.cpp SGUtil.h \
        SGAlgorithms.cpp SGAlgorithms.h \
        SGVisitors.h SGVisitors.cpp \
        CompleteOverlapSet.h CompleteOverlapSet.cpp \
        RemovalAlgorithm.h RemovalAlgorithm.cpp \
		SGSearch.h SGSearch.cpp \
		GraphSearchTree.h \
		SGWalk.h SGWalk.cpp

