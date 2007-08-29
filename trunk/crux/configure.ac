dnl 
dnl <Copyright = jasone>
dnl <License>
dnl 
dnl Version: Crux <Version = crux>
dnl 
dnl Process this file with autoconf to produce a configure script.
AC_INIT([Makefile.in])

dnl 
dnl Version requirements for external software we depend on.
dnl 
PYTHON_VERSION_MIN=2.3.3

srcroot=$srcdir
AC_SUBST([srcroot])
abs_srcroot=`cd "$srcdir"; pwd`
AC_SUBST([abs_srcroot])

objroot=.
AC_SUBST([objroot])
abs_objroot=`pwd`
AC_SUBST([abs_objroot])

dnl Munge install path variables.
if test "x$prefix" = "xNONE" ; then
  prefix="/usr/local"
fi
if test "x$exec_prefix" = "xNONE" ; then
  exec_prefix=$prefix
fi
PREFIX=$prefix
AC_SUBST([PREFIX])
BINDIR=`eval echo $bindir`
BINDIR=`eval echo $BINDIR`
AC_SUBST([BINDIR])
DATADIR=`eval echo $datadir`
DATADIR=`eval echo $DATADIR`
AC_SUBST([DATADIR])
INCLUDEDIR=`eval echo $includedir`
INCLUDEDIR=`eval echo $INCLUDEDIR`
AC_SUBST([INCLUDEDIR])
MANDIR=`eval echo $mandir`
MANDIR=`eval echo $MANDIR`
AC_SUBST([MANDIR])

cfgoutputs="Makefile setup.py bin/crux pkg/crux/__init__.py"
cfghdrs="config.stamp include/_cruxmodule_defs.h"
bins=
mods=

dnl Prevent autoconf from molesting CFLAGS.
CFLAGS=$CFLAGS
AC_PROG_CC
AC_PROG_CPP

dnl Platform-specific settings.
AC_CANONICAL_HOST
case "${host_cpu}" in
  i[[3456]]86)
	AC_DEFINE([CxmCpuIa32])
	extra_compile_args="'-msse2'"
	;;
  x86_64)
	AC_DEFINE([CxmCpuAmd64])
	extra_compile_args="'-msse2'"
	;;
  powerpc)
	AC_DEFINE([CxmCpuPpc])
	extra_compile_args="'-faltivec'"
	;;
  *)
	extra_compile_args=""
	;;
esac
AC_SUBST([extra_compile_args])

dnl Disable rules that do automatic regeneration of configure output by default.
AC_ARG_ENABLE([autogen], [  --enable-autogen        Automatically regenerate configure output],
if test "x$enable_autogen" = "xno" ; then
  enable_autogen="0"
else
  enable_autogen="1"
fi
,
enable_autogen="0"
)
AC_SUBST([enable_autogen])

AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], , [$PATH])
AC_PATH_PROG([LD], [ld], , [$PATH])
AC_PATH_PROG([TEST], [test], , [$PATH])
AC_PATH_PROG([ECHO], [echo], , [$PATH])
AC_PATH_PROG([PERL], [perl], , [$PATH])
AC_PATH_PROG([ACLOCAL], [aclocal], , [$PATH])
AC_PATH_PROG([AUTOCONF], [autoconf], , [$PATH])
AC_PATH_PROG([AUTOHEADER], [autoheader], , [$PATH])

dnl Do not compile with debugging by default.
AC_ARG_ENABLE([debug], [  --enable-debug          Build debugging code],
[if test "x$enable_debug" = "xno" ; then
  enable_debug="0"
else
  enable_debug="1"
fi
],
[enable_debug="0"]
)
if test "x$enable_debug" = "x1" ; then
  python_opts=""
  AC_DEFINE([CxmDebug])
  AC_DEFINE([CxmAssertions])
else
  python_opts=" -O"
fi
AC_SUBST([enable_debug])

dnl Do not silence irrelevant compiler warnings by default, since enabling this
dnl option incurs a slight performance penalty.
AC_ARG_ENABLE([cc-silence], [  --enable-cc-silence     Silence irrelevant compiler warnings],
[if test "x$enable_cc_silence" = "xno" ; then
  enable_cc_silence="0"
else
  enable_cc_silence="1"
fi
],
[enable_cc_silence="0"]
)
if test "x$enable_cc_silence" = "x1" ; then
  AC_DEFINE([CxmCcSilence])
fi
AC_SUBST([enable_cc_silence])

dnl 
dnl Dependency configuration.
dnl 
dnl Get path to python.
AC_ARG_WITH([python], [  --with-python           python program],
if test "x$with_python" = "xno" ; then
  PYTHON=
else
  PYTHON="$with_python"
fi
,
  PYTHON=
)
if test "x$PYTHON" != "x" ; then
   if test ! -x "${PYTHON}" ; then
     AC_MSG_ERROR([Unusable or missing python: ${PYTHON}])
   fi
else
  AC_PATH_PROG([PYTHON], [python], , [$PATH])
  if test "x$PYTHON" = "x" ; then
    AC_MSG_ERROR([Cannot build without python])
  fi
fi

dnl Make sure the version of python found is sufficient.
AC_MSG_CHECKING([for python >= $PYTHON_VERSION_MIN])
PYTHON_FOUND=`2>&1 $PYTHON -V | awk '{print $2}'`
if test "x$PYTHON_FOUND" != "x<Version" ; then
  PYTHON_MAJOR=`echo ${PYTHON_FOUND} | tr . " " | awk '{print $1}'`
  PYTHON_MINOR=`echo ${PYTHON_FOUND} | tr . " " | awk '{print $2}' | tr a-z " " |awk '{print $1}'`
  PYTHON_BRANCH=`echo ${PYTHON_FOUND} | tr . " " | awk '{print $3}' | tr a-z " " |awk '{print $1}'`
  if test "x$PYTHON_BRANCH" = "x" ; then
    PYTHON_BRANCH=0
  fi
  PYTHON_MAJOR_MIN=`echo $PYTHON_VERSION_MIN | tr . " " | awk '{print $1}'`
  PYTHON_MINOR_MIN=`echo $PYTHON_VERSION_MIN | tr . " " | awk '{print $2}'`
  PYTHON_BRANCH_MIN=`echo $PYTHON_VERSION_MIN | tr . " " | awk '{print $3}'`
  if test $PYTHON_MAJOR -gt $PYTHON_MAJOR_MIN \
       -o $PYTHON_MAJOR -eq $PYTHON_MAJOR_MIN \
       -a $PYTHON_MINOR -gt $PYTHON_MINOR_MIN \
       -o $PYTHON_MAJOR -eq $PYTHON_MAJOR_MIN \
       -a $PYTHON_MINOR -eq $PYTHON_MINOR_MIN \
       -a $PYTHON_BRANCH -ge $PYTHON_BRANCH_MIN ; then
    AC_MSG_RESULT([$PYTHON_MAJOR.$PYTHON_MINOR.$PYTHON_BRANCH])
    have_python="1"
  else
    AC_MSG_RESULT([no])
    have_python="0"
  fi
else
  AC_MSG_RESULT([unknown development version])
  have_python="1"
fi
if test "x$have_python" = "x0" ; then
  AC_MSG_ERROR([Cannot build without python])
fi
AC_SUBST([python_opts])
AC_SUBST([PYTHON])

AC_SUBST([bins])
AC_SUBST([mods])

dnl 
dnl crux configuration.
dnl 
if echo "<Version = crux>" | grep Version >/dev/null ; then
  crux_version=0.0.0
else
  crux_version="<Version = crux>"
fi
AC_SUBST([crux_version])

dnl Packages.
crux_packages=""
for p in crux ; do
  crux_packages="${crux_packages} '${p}',"
done
crux_packages="[[${crux_packages}]]"
AC_SUBST([crux_packages])

dnl C sources.
crux_sources=""
for s in `find ${srcroot}/src -type f |grep "\.c$"` ; do
  crux_sources="${crux_sources} '${s}',"
done
crux_sources="[[${crux_sources}]]"
AC_SUBST([crux_sources])

dnl C headers.
crux_headers=""
for h in `find ${srcroot}/include -type f |grep "\.h$"` \
         $objroot/_cruxmodule_defs.h.in ${srcroot}/src/asprintf_c ; do
  crux_headers="${crux_headers} ${h}"
done
AC_SUBST([crux_headers])

dnl Python scripts.
crux_scripts=""
for s in crux ; do
  crux_scripts="${crux_scripts} '${objroot}/bin/${s}',"
done
crux_scripts="[[${crux_scripts}]]"
AC_SUBST([crux_scripts])

dnl Crux tests.
crux_tests=""
for t in `find ${srcroot}/test -type f |grep "\.crux$"`; do
  crux_tests="${crux_tests} ${t}"
done
AC_SUBST([crux_tests])

dnl Check for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_HEADER_STDBOOL

dnl Check for library functions.
AC_CHECK_FUNCS([asprintf])

dnl Process .in files.
AC_SUBST([cfghdrs])
AC_CONFIG_HEADER([$cfghdrs])

dnl crux is handled separately, since it must be executable.
AC_CONFIG_FILES([$cfgoutputs])
AC_CONFIG_FILES([crux], [chmod +x crux])
cfgoutputs="$cfgoutputs crux"
AC_SUBST([cfgoutputs])
AC_OUTPUT

dnl 
dnl Print out the results of configuration.
dnl 
AC_MSG_RESULT([===============================================================================])
AC_MSG_RESULT([Crux version  : $crux_version])
AC_MSG_RESULT([])
AC_MSG_RESULT([PYTHON        : ${PYTHON}])
AC_MSG_RESULT([CC            : ${CC}])
dnl AC_MSG_RESULT([CPPFLAGS      : ${CPPFLAGS}])
dnl AC_MSG_RESULT([CFLAGS        : ${CFLAGS}])
dnl AC_MSG_RESULT([LDFLAGS       : ${LDFLAGS}])
dnl AC_MSG_RESULT([LIBS          : ${LIBS}])
AC_MSG_RESULT([])
AC_MSG_RESULT([PREFIX        : ${PREFIX}])
AC_MSG_RESULT([BINDIR        : ${BINDIR}])
AC_MSG_RESULT([DATADIR       : ${DATADIR}])
AC_MSG_RESULT([INCLUDEDIR    : ${INCLUDEDIR}])
AC_MSG_RESULT([MANDIR        : ${MANDIR}])
AC_MSG_RESULT([])
AC_MSG_RESULT([srcroot       : ${srcroot}])
AC_MSG_RESULT([abs_srcroot   : ${abs_srcroot}])
AC_MSG_RESULT([objroot       : ${objroot}])
AC_MSG_RESULT([abs_objroot   : ${abs_objroot}])
AC_MSG_RESULT([])
AC_MSG_RESULT([autogen       : ${enable_autogen}])
AC_MSG_RESULT([debug         : ${enable_debug}])
AC_MSG_RESULT([cc-silence    : ${enable_cc_silence}])
AC_MSG_RESULT([===============================================================================])
