#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowman, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Snowman/snowman.cpp])
AC_CONFIG_HEADER([config.h])
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

# Check for the bamtools library pathz
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here
AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

if test "$with_bamtools" -a -d "$with_bamtools"; then
    bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib -Wl,-rpath,$with_bamtools/lib/bamtools -L$with_bamtools/lib/bamtools"
    bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
else
    bamt=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21
    bamtools_ldflags="-Wl,-rpath,$bamt/lib -L$bamt/lib -Wl,-rpath,$bamt/lib/bamtools -L$bamt/lib/bamtools"
    bamtools_include="-I$bamt/include -I$bamt/include/bamtools"
fi

##AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
##	[specify directory containing the Boost library ()]))

##if test "$with_boost" -a -d "$with_boost"; then
##    boost_ldflags="-Wl,-rpath,$with_boost/lib -L$with_boost/lib"
##    boost_include="-I$with_boost/include"
##else
##    bamt=/broad/software/free/Linux/redhat_5_x86_64/pkgs/pezmaster31_bamtools-6708a21
##    boostr=/broad/software/free/Linux/redhat_5_x86_64/pkgs/boost_1.54.0
##    bamtools_ldflags="-Wl,-rpath,$boostr/lib -L$boostr/lib"
##    bamtools_include="-I$boostr/include"
##fi

# Check for the seqan library paths
AC_ARG_WITH(seqan, AS_HELP_STRING([--with-seqan=PATH],
	[specify directory containing the seqan library (http://trac.seqan.de)]))

if test "$with_seqan" -a -d "$with_seqan"; then
    seqan_include="-I$with_seqan/include"
##else 
##    seqan_include="-I/xchip/gistic/Jeremiah/seqan-trunk/core/include"
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
AC_ARG_WITH(snowtools, AS_HELP_STRING([--with-snowtools=PATH],
	[specify directory containing the SnowTools library]))
if test "$with_snowtools" -a -d "$with_snowtools"; then
    snowtools_ldflags="-L$with_snowtools"
    snowtools_include="-Iwith_snowtools"
    snowtools_libs="-lsnowtools"
else
    snowtools_ldflags="-L /xchip/gistic/Jeremiah/GIT/isva/SnowTools/src"
    snowtools_include="-I /xchip/gistic/Jeremiah/GIT/isva/SnowTools/src"
    snowtools_libs="-lsnowtools"
fi  

# Check for the hts library paths
AC_ARG_WITH(htslib, AS_HELP_STRING([--with-htslib=PATH],
	[specify directory containing the Heng Li's HTSlib)]))

if test "$with_htslib" -a -d "$with_htslib"; then
    htslib_include="-I$with_htslib"
    htslib_ldflags="-Wl,-rpath,$with_htslib -L$with_htslib"
else
    htslibr=/xchip/gistic/Jeremiah/htslib-1.1	
    htslib_include="-I$htslibr"
    htslib_ldflags="-Wl,-rpath,$htslibr -L$htslibr"
fi

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11")
AC_SUBST(CXXFLAGS, "-O3")
AC_SUBST(CFLAGS, "-O3")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $snowtools_include $htslib_include $bamtools_include $samtools_include $seqan_include")
AC_SUBST(LDFLAGS, "$snowtools_ldflags $bamtools_ldflags $htslib_ldflags $samtools_ldflags $LDFLAGS")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS $snowtools_libs -lbamtools -lhts $samtools_libs")


# Check for hash table headers
AC_CHECK_HEADERS([ \
		   google/sparse_hash_set google/sparse_hash_map \
		   unordered_map tr1/unordered_map ext/hash_map \
		   unordered_set tr1/unordered_set ext/hash_set \
])

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

# Make sure Boost can be found
##AC_CHECK_HEADERS([boost/algorithm/string.hpp],,[AC_MSG_ERROR([The Boost library must be installed (). You can specify its path with the --with-boost=PATH option])])

# Make sure the seqan headers can be found
AC_CHECK_HEADERS([seqan/basic.h],,[AC_MSG_WARN([The seqan library must be installed (http://trac.seqan.de) for full functionality. You can specify its path with the --with-seqan=PATH option])])

# Make sure the HTSlib headers can be found
AC_CHECK_HEADERS([htslib/hts.h],,[AC_MSG_ERROR([The HTSlib library must be installed. You can specify its path with the --with-htslib=PATH option])])

## Make sure the snowtools headers can be found
AC_CHECK_HEADERS([SnowTools/SnowTools.h],,[AC_MSG_ERROR([The SnowTools library must be installed. You can specify its path with the --with-snowtools=PATH option])])

## Make sure the samtools headers can be found
##AC_CHECK_HEADERS([sam.h],,[AC_MSG_ERROR([The Samtools library must be installed. You can specify its path with the --with-samtools=PATH option])])

AC_CONFIG_FILES([Makefile
		SGA/Thirdparty/Makefile
		SGA/Util/Makefile
		SGA/SQG/Makefile
		SGA/Bigraph/Makefile
		SGA/Algorithm/Makefile
		SGA/StringGraph/Makefile
		SGA/Concurrency/Makefile
		SGA/SuffixTools/Makefile
                SGA/GraphDiff/Makefile
                SGA/Scaffold/Makefile
                SGA/SGA/Makefile
                SnowUtils/Makefile
		BWA/Makefile
		BedTools/Makefile
                Snowman/Makefile])

AC_OUTPUT
