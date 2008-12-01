dnl Process this file with autoconf to produce a configure script.
AC_INIT([Makefile.in])

dnl 
dnl Version requirements for external software we depend on.
dnl 
PYTHON_VERSION_MIN=2.6.0
CYTHON_VERSION_MIN=0.10.0

srcroot=$srcdir
if test "x${srcroot}" = "x." ; then
  srcroot=""
else
  srcroot="${srcroot}/"
fi
AC_SUBST([srcroot])
abs_srcroot="`cd \"${srcdir}\"; pwd`/"
AC_SUBST([abs_srcroot])

objroot=""
AC_SUBST([objroot])
abs_objroot="`pwd`/"
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
MANDIR=`eval echo $mandir`
MANDIR=`eval echo $MANDIR`
AC_SUBST([MANDIR])

cfgoutputs="Makefile"
cfghdrs=""

dnl If CFLAGS isn't defined and using gcc, set CFLAGS to something reasonable.
dnl Otherwise, just prevent autoconf from molesting CFLAGS.
CFLAGS=$CFLAGS
AC_PROG_CC
if test "x$CFLAGS" = "x" ; then
  no_CFLAGS="yes"
fi
if test "x$no_CFLAGS" = "xyes" -a "x$GCC" = "xyes" ; then
 CFLAGS="-std=gnu99 -Wall -pipe -g3"
fi
dnl Append EXTRA_CFLAGS to CFLAGS, if defined.
if test "x$EXTRA_CFLAGS" != "x" ; then
  CFLAGS="$CFLAGS $EXTRA_CFLAGS"
fi
AC_PROG_CPP

AC_CANONICAL_HOST
dnl CPU-specific settings.
case "${host_cpu}" in
  i[[3456]]86)
	AC_DEFINE([CxmCpuIa32])
	CFLAGS="${CFLAGS} -msse2"
	;;
  x86_64)
	AC_DEFINE([CxmCpuAmd64])
	CFLAGS="${CFLAGS} -msse2"
	;;
  powerpc)
	AC_DEFINE([CxmCpuPpc])
	CFLAGS="${CFLAGS} -faltivec"
	;;
  *)
	;;
esac

dnl Platform-specific settings.  abi and RPATH can probably be determined
dnl programmatically, but doing so is error-prone, which makes it generally
dnl not worth the trouble.
dnl 
dnl Define cpp macros in CPPFLAGS, rather than doing AC_DEFINE(macro), since the
dnl definitions need to be seen before any headers are included, which is a pain
dnl to make happen otherwise.
case "${host}" in
  *-*-darwin*)
	CFLAGS="$CFLAGS -fno-common -no-cpp-precomp"
	abi="macho"
	RPATH=""
	;;
  *-*-freebsd*)
	CFLAGS="$CFLAGS"
	abi="elf"
	RPATH="-Wl,-rpath,"
	;;
  *-*-linux*)
	CFLAGS="$CFLAGS"
	abi="elf"
	RPATH="-Wl,-rpath,"
	;;
  *-*-netbsd*)
	AC_MSG_CHECKING([ABI])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#ifdef __ELF__
/* ELF */
#else
#error aout
#endif
]])],
                          [CFLAGS="$CFLAGS"; abi="elf"],
                          [abi="aout"])
	AC_MSG_RESULT([$abi])
	RPATH="-Wl,-rpath,"
	;;
  *-*-solaris2*)
	CFLAGS="$CFLAGS"
	abi="elf"
	RPATH="-Wl,-R,"
	dnl Solaris needs this for sigwait().
	CPPFLAGS="$CPPFLAGS -D_POSIX_PTHREAD_SEMANTICS"
	LIBS="$LIBS -lposix4 -lsocket -lnsl"
	;;
  *)
	AC_MSG_RESULT([Unsupported operating system: ${host}])
	abi="elf"
	RPATH="-Wl,-rpath,"
	;;
esac
AC_SUBST([abi])
AC_SUBST([RPATH])

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

dnl Disable Cython code generation by default.
AC_ARG_ENABLE([cython], [  --enable-cython         Compile .pyx/.pxd sources to .c outputs],
if test "x$enable_cython" = "xyes" ; then
  enable_cython="1"
else
  enable_cython="0"
fi
,
enable_cython="0"
)
AC_SUBST([enable_cython])

AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], , [$PATH])
AC_PATH_PROG([LD], [ld], , [$PATH])
AC_PATH_PROG([TEST], [test], , [$PATH])
AC_PATH_PROG([ECHO], [echo], , [$PATH])
AC_PATH_PROG([PERL], [perl], , [$PATH])
AC_PATH_PROG([AUTOCONF], [autoconf], , [$PATH])

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

dnl Only optimize if not debugging.
if test "x$enable_debug" = "x0" -a "x$no_CFLAGS" = "xyes" ; then
  dnl Make sure that an optimization flag was not specified in EXTRA_CFLAGS.
  if test "x$GCC" = "xyes" ; then
    echo "$EXTRA_CFLAGS" | grep "\-O" >/dev/null || CFLAGS="$CFLAGS -O2"
  else
    echo "$EXTRA_CFLAGS" | grep "\-O" >/dev/null || CFLAGS="$CFLAGS -O"
  fi
fi

dnl Support inlines by default.
AC_ARG_ENABLE([inlines],
              [  --disable-inlines       Do not use inline functions],
[if test "x$enable_inlines" = "xno" ; then
  enable_inlines="0"
else
  enable_inlines="1"
fi
],
[enable_inlines="1"]
)
dnl inlines are incompatible with debug.
if test "x$enable_debug" = "x1" ; then
  enable_inlines="0"
fi
if test "x$enable_inlines" = "x1" ; then
  AC_DEFINE([CxmUseInlines])
fi
AC_SUBST([enable_inlines])

dnl ============================================================================
dnl Make sure the version of python found is sufficient.

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

AC_MSG_CHECKING([for python-config])
python_config="${PYTHON}-config"
if test -x "${python_config}" ; then
  AC_MSG_RESULT([${python_config}])
else
  AC_MSG_ERROR([Cannot build without python-config])
fi
python_config_includes=`${python_config} --includes`
AC_SUBST([python_config_includes])
python_config_cflags=`${python_config} --cflags`
AC_SUBST([python_config_cflags])
python_config_ldflags=`${python_config} --ldflags`
AC_SUBST([python_config_ldflags])

AC_MSG_CHECKING([python site-packages path])
1>conftest cat << EOF
import distutils.sysconfig
print distutils.sysconfig.get_python_lib()
EOF
site_packages=`${PYTHON} conftest`
rm -f conftest
AC_MSG_RESULT([${site_packages}])
AC_SUBST([site_packages])

dnl ============================================================================
dnl Make sure the version of cython found is sufficient.
if test "x${enable_cython}" = "x1" ; then
  dnl Get path to cython.
  AC_ARG_WITH([cython], [  --with-cython           cython program],
  if test "x$with_cython" = "xno" ; then
    CYTHON=
  else
    CYTHON="$with_cython"
  fi
  ,
    CYTHON=
  )
  if test "x$CYTHON" != "x" ; then
     if test ! -x "${CYTHON}" ; then
       AC_MSG_ERROR([Unusable or missing cython: ${CYTHON}])
     fi
  else
    AC_PATH_PROG([CYTHON], [cython], , [$PATH])
    if test "x$CYTHON" = "x" ; then
      AC_MSG_ERROR([Cannot build without cython])
    fi
  fi

  AC_MSG_CHECKING([for cython >= $CYTHON_VERSION_MIN])
  CYTHON_FOUND=`2>&1 $CYTHON -V | awk '{print $3}'`
  if test "x$CYTHON_FOUND" != "x<Version" ; then
    CYTHON_MAJOR=`echo ${CYTHON_FOUND} | tr . " " | awk '{print $1}'`
    CYTHON_MINOR=`echo ${CYTHON_FOUND} | tr . " " | awk '{print $2}' | tr a-z " " |awk '{print $1}'`
    CYTHON_BRANCH=`echo ${CYTHON_FOUND} | tr . " " | awk '{print $3}' | tr a-z " " |awk '{print $1}'`
    if test "x$CYTHON_BRANCH" = "x" ; then
      CYTHON_BRANCH=0
    fi
    CYTHON_MAJOR_MIN=`echo $CYTHON_VERSION_MIN | tr . " " | awk '{print $1}'`
    CYTHON_MINOR_MIN=`echo $CYTHON_VERSION_MIN | tr . " " | awk '{print $2}'`
    CYTHON_BRANCH_MIN=`echo $CYTHON_VERSION_MIN | tr . " " | awk '{print $3}'`
    if test $CYTHON_MAJOR -gt $CYTHON_MAJOR_MIN \
         -o $CYTHON_MAJOR -eq $CYTHON_MAJOR_MIN \
         -a $CYTHON_MINOR -gt $CYTHON_MINOR_MIN \
         -o $CYTHON_MAJOR -eq $CYTHON_MAJOR_MIN \
         -a $CYTHON_MINOR -eq $CYTHON_MINOR_MIN \
         -a $CYTHON_BRANCH -ge $CYTHON_BRANCH_MIN ; then
      AC_MSG_RESULT([$CYTHON_MAJOR.$CYTHON_MINOR.$CYTHON_BRANCH])
      have_cython="1"
    else
      AC_MSG_RESULT([no])
      have_cython="0"
    fi
  else
    AC_MSG_RESULT([unknown development version])
    have_cython="1"
  fi
  if test "x$have_cython" = "x0" ; then
    AC_MSG_ERROR([Cannot build without cython])
  fi
  AC_SUBST([CYTHON])
fi

dnl ============================================================================
dnl Crux configuration.
dnl 
if echo "<Version = crux>" | grep Version >/dev/null ; then
  crux_version=0.0.0
else
  crux_version="<Version = crux>"
fi
AC_SUBST([crux_version])

dnl ============================================================================
dnl Generate file lists.

dnl C headers.
chdrs=""
for h in `find ${srcroot}lib -type f |grep "\.h$" |sort` ; do
  chdrs="${chdrs} ${h}"
done
for h in `find ${srcroot}lib -type f |grep "\.h\.in$" |sort` ; do
  chdr="`dirname ${h}`/`basename ${h} .in`"
  chdrs="${chdrs} ${chdr}"
  cfghdrs="${cfghdrs} ${chdr}"
done
AC_SUBST([chdrs])

dnl C sources.
csrcs=""
for s in `find ${srcroot}lib -type f |grep "Cx[[A-Za-z0-9]]*\.c$" |sort` ; do
  csrcs="${csrcs} ${s}"
done
AC_SUBST([csrcs])

dnl Cython headers.
pxds=""
for d in `find ${srcroot}pkg -type f |grep "\.pxd$" |sort` ; do
  pxds="${pxds} ${d}"
done
AC_SUBST([pxds])

dnl Cython module/package sources.
pyxs=""
pkgs=""
for s in `find ${srcroot}pkg -type f |grep "\.pyx$" |sort` ; do
  if test "`basename ${s} .pyx`" != "__init__" ; then
    pyxs="${pyxs} ${s}"
  else
    pkgs="${pkgs} ${s}"
  fi
done
AC_SUBST([pyxs])
AC_SUBST([pkgs])

dnl Python modules.
pys=""
for s in `find ${srcroot}pkg -type f |grep "\.py$" |sort` ; do
  pys="${pys} ${s}"
done
AC_SUBST([pys])

dnl Python scripts.
bins=""
for b in `find ${srcroot}bin -type f |sort` ; do
  bins="{$bins} ${b}"
done
AC_SUBST([bins])

dnl Tests.
tests=""
for t in `find ${srcroot}test -type f |grep "\.crux$" |sort`; do
  tests="${tests} ${t}"
done
AC_SUBST([tests])

dnl ============================================================================
dnl Check for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_HEADER_STDBOOL

dnl Check for library functions.

dnl Process .in files.
AC_SUBST([cfghdrs])
AC_CONFIG_HEADER([$cfghdrs cfghdrs.stamp])

AC_CONFIG_FILES([$cfgoutputs cfgoutputs.stamp])
AC_DEFUN([CX_CONFIG_BINS], [
AC_CONFIG_FILES([$1], [chmod +x $1])
cfgoutputs="${cfgoutputs} $1"
])
CX_CONFIG_BINS([cfg])
CX_CONFIG_BINS([cydep])
AC_SUBST([cfgoutputs])
AC_OUTPUT

dnl ============================================================================
dnl Print out the results of configuration.
AC_MSG_RESULT([===============================================================================])
AC_MSG_RESULT([Crux version  : $crux_version])
AC_MSG_RESULT([])
AC_MSG_RESULT([PYTHON        : ${PYTHON}])
AC_MSG_RESULT([CYTHON        : ${CYTHON}])
AC_MSG_RESULT([CC            : ${CC}])
AC_MSG_RESULT([CPPFLAGS      : ${CPPFLAGS}])
AC_MSG_RESULT([CFLAGS        : ${CFLAGS}])
AC_MSG_RESULT([LDFLAGS       : ${LDFLAGS}])
AC_MSG_RESULT([LIBS          : ${LIBS}])
AC_MSG_RESULT([])
AC_MSG_RESULT([PREFIX        : ${PREFIX}])
AC_MSG_RESULT([BINDIR        : ${BINDIR}])
AC_MSG_RESULT([DATADIR       : ${DATADIR}])
AC_MSG_RESULT([MANDIR        : ${MANDIR}])
AC_MSG_RESULT([])
AC_MSG_RESULT([srcroot       : ${srcroot}])
AC_MSG_RESULT([abs_srcroot   : ${abs_srcroot}])
AC_MSG_RESULT([objroot       : ${objroot}])
AC_MSG_RESULT([abs_objroot   : ${abs_objroot}])
AC_MSG_RESULT([])
AC_MSG_RESULT([autogen       : ${enable_autogen}])
AC_MSG_RESULT([cython        : ${enable_cython}])
AC_MSG_RESULT([debug         : ${enable_debug}])
AC_MSG_RESULT([cc-silence    : ${enable_cc_silence}])
AC_MSG_RESULT([inlines       : ${enable_inlines}])
AC_MSG_RESULT([===============================================================================])
