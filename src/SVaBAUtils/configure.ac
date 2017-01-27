# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowman, 1.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([snowmanutils.cpp])
AC_CONFIG_HEADER([config.h])
AM_MAINTAINER_MODE([disable])

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

AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-g -Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0")

# Check for hash table headers
AC_CHECK_HEADERS([ \
		   google/sparse_hash_set google/sparse_hash_map \
		   unordered_map tr1/unordered_map ext/hash_map \
		   unordered_set tr1/unordered_set ext/hash_set \
])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
