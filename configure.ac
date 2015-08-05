#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowman, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/Snowman/snowman.cpp])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE([enable])
##m4_include([m4/m4_ax_openmp.m4])

# Checks for programs.
AC_PROG_CXX ## test for cpp compiler
AC_PROG_CC ## test for C compiler
AC_PROG_RANLIB ## required if libraries are built in package

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )

AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))
if test "$with_boost" -a -d "$with_boost"; then
    boost_include="-I$with_boost"
else
    boostr=/xchip/gistic/Jeremiah/boost_1_58_0
    boost_include="-I$boostr"
fi

### CHECK FOR SAMTOOLS
AC_ARG_WITH(samtools, AS_HELP_STRING([--with-samtools=PATH],
	[specify directory containing the Samtools library]))
if test "$with_samtools" -a -d "$with_samtools"; then
    samtools_ldflags="-L$with_samtools"
    samtools_include="-Iwith_samtools"
    samtools_libs="-lbam"
else
    samtools_ldflags="-L /xchip/gistic/Jeremiah/GIT/samtools"
    samtools_include="-I /xchip/gistic/Jeremiah/GIT/samtools"
    samtools_libs="-lbam"
fi  

### CHECK FOR SNOW TOOLS
#AC_ARG_WITH(snowtools, AS_HELP_STRING([--with-snowtools=PATH],
#	[specify directory containing the SnowTools library]))
##if test "$with_snowtools" -a -d "$with_snowtools"; then
#    snowtools_ldflags="-L$with_snowtools"
#    snowtools_include="-I$with_snowtools -I$with_snowtools/multifast-v1.4.2 -I$with_snowtools/htslib"
#    snowtools_libs="-lsnowtools"
#else
#    broad_st="/xchip/gistic/Jeremiah/GIT/SnowTools/src"
#    snowtools_ldflags="-L$broad_st -Wl,-rpath,$broad_st/htslib -L$broad_st/htslib -L$broad_st/bwa"
#    snowtools_include="-I$broad_st -I$broad_st/multifast-v1.4.2 -I$broad_st/htslib"
#    snowtools_libs="-lhts -lbwa -lsnowtools"
#fi  

AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $samtools_include $boost_include")
AC_SUBST(LDFLAGS, "$samtools_ldflags $LDFLAGS")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS $samtools_libs")

# Check for hash table headers
AC_CHECK_HEADERS([ \
		   google/sparse_hash_set google/sparse_hash_map \
		   unordered_map tr1/unordered_map ext/hash_map \
		   unordered_set tr1/unordered_set ext/hash_set \
])

# Make sure the bamtools headers can be found
#AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

# Make sure Boost can be found
##AC_CHECK_HEADERS([boost/algorithm/string.hpp],,[AC_MSG_ERROR([The Boost library must be installed (). You can specify its path with the --with-boost=PATH option])])

# Make sure the seqan headers can be found
##AC_CHECK_HEADERS([seqan/basic.h],,[AC_MSG_WARN([The seqan library must be installed (http://trac.seqan.de) for full functionality. You can specify its path with the --with-seqan=PATH option])])

# Make sure the HTSlib headers can be found
#AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTSlib library must be installed. You can specify its path with the --with-htslib=PATH option])])

## Make sure the snowtools headers can be found
#AC_CHECK_HEADERS([SnowTools/SnowTools.h],,[AC_MSG_ERROR([The SnowTools library must be installed. You can specify its path with the --with-snowtools=PATH option])])

## Make sure the samtools headers can be found
##AC_CHECK_HEADERS([sam.h],,[AC_MSG_ERROR([The Samtools library must be installed. You can specify its path with the --with-samtools=PATH option])])

AC_CONFIG_FILES([Makefile
                SnowTools/src/Makefile
		src/SGA/Util/Makefile
		src/SGA/SQG/Makefile
		src/SGA/Bigraph/Makefile
		src/SGA/Algorithm/Makefile
		src/SGA/StringGraph/Makefile
		src/SGA/SuffixTools/Makefile
                src/SGA/SGA/Makefile
                src/Snowman/Makefile])

AC_OUTPUT
