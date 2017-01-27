# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(svaba, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/SVaBA/svaba.cpp])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE([disable])
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

##AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
##	[specify directory containing the boost library)]))
##if test "$with_boost" -a -d "$with_boost"; then
##    boost_include="-I$with_boost"
##    AC_CHECK_FILE("$with_boost/stage/lib", STAGE_PATH=1, STAGE_PATH=0)
##    if test ${STAGE_PATH} = 1; then
##        boost_lib="$with_boost/stage/lib"  
##    else
##        boost_lib="$with_boost/lib"        
##    fi
##    boost_lib_flag="-L/$boost_lib"
##    AC_SUBST(boost_lib)
##fi

AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")
##AC_SUBST(CXXFLAGS, "-O3")
##AC_SUBST(CFLAGS, "-O3")
##AC_SUBST(CPPFLAGS, "$CPPFLAGS $boost_include")
AC_SUBST(LDFLAGS, "$LDFLAGS -pthread -std=c++11")

# We always need to specify to link to certain libraries
##AC_SUBST(LIBS, "$LIBS") # libboost_regex.a libboost_filesystem.a libboost_system.a")

# Check for hash table headers
AC_CHECK_HEADERS([ \
		   google/sparse_hash_set google/sparse_hash_map \
		   unordered_map tr1/unordered_map ext/hash_map \
		   unordered_set tr1/unordered_set ext/hash_set \
])

# Make sure the boost headers can be found
##AC_CHECK_HEADERS([boost/regex.hpp],[AC_SUBST(CPPFLAGS, "$CPPFLAGS -DHAVE_BOOST")],[AC_MSG_ERROR([The Boost library must be installed for SVaBA. Specify its path with the --with-boost=PATH option])])

AC_CONFIG_FILES([Makefile
                SeqLib/src/Makefile
		src/SGA/Util/Makefile
		src/SGA/SQG/Makefile
		src/SGA/Bigraph/Makefile
		src/SGA/Algorithm/Makefile
		src/SGA/StringGraph/Makefile
		src/SGA/SuffixTools/Makefile
                src/SGA/SGA/Makefile
                src/SVaBA/Makefile])

AC_OUTPUT
