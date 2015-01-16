# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(snowtools, 0.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign) ## foreign tells each automake Makefile.am not to worry about GNU standards
AC_CONFIG_SRCDIR([snowtools.cpp]) ## check that we are looking at right dir
AC_CONFIG_HEADER([config.h])
##m4_include([m4/m4_ax_openmp.m4])

AC_PROG_CXX ## test for cpp compiler
AC_PROG_CC ## test for C compiler
AC_PROG_RANLIB ## required if libraries are built in package

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for the bamtools library pathz
# Bamtools has two different path formats
# If it is built in place the files will be in {path}/lib/ and {path}/include/
# If it is installed (with make install) it will be in {path}/include/bamtools
# and {path}/lib/bamtools. Both variants are handled here
AC_ARG_WITH(bamtools, AS_HELP_STRING([--with-bamtools=PATH],
	[specify directory containing the bamtools library (http://github.com/pezmaster31/bamtools)]))

if test "$with_bamtools" -a -d "$with_bamtools"; then
    bamtools_ldflags="-Wl,-rpath,$with_bamtools/lib -L$with_bamtools/lib"
    bamtools_include="-I$with_bamtools/include -I$with_bamtools/include/bamtools"
fi

## check for snowman library path
AC_ARG_WITH(snowutils, AS_HELP_STRING([--with-snowutils=PATH],
	[specify directory containing the snowutils library]))

if test "$with_snowutils" -a -d "$with_snowutils"; then
    snowutils_ldflags="-Wl,-rpath,$with_snowutils -L$with_snowutils"
    snowutils_include="-I$with_snowutils"
fi

AC_MSG_WARN([$snowutils_ldflags])
AC_MSG_WARN([$snowutils_include])

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
AC_SUBST(CPPFLAGS, "$CPPFLAGS $bamtools_include $snowutils_include")
AC_SUBST(LDFLAGS, "$bamtools_ldflags $snowutils_ldflags $LDFLAGS")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS -lbamtools")
# Add a link to snowutils library
AC_SUBST(LIBS, "$LIBS -lsnowutils")

# Make sure the bamtools headers can be found
AC_CHECK_HEADERS([api/BamReader.h],,[AC_MSG_ERROR([The bamtools library must be installed (http://github.com/pezmaster31/bamtools). You can specify its path with the --with-bamtools=PATH option])])

AC_MSG_WARN([$CPPFLAGS])
AC_MSG_WARN([$LDFLAGS])
# Make sure the snowutil headers can be found
AC_CHECK_HEADERS([SnowCommon.h],,[AC_MSG_ERROR([The snowutils library must be installed. You can specify its path with the --with-snowutils=PATH option])])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
