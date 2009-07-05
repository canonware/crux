dnl Process this file with autoconf to produce a configure script.
AC_INIT([Makefile.in])

dnl 
dnl Version requirements for external software we depend on.
dnl 
PYTHON_VERSION_MIN=2.5.0
CYTHON_VERSION_MIN=0.11.0
FLEX_VERSION_MIN=2.5.35
MPI4PY_VERSION_MIN=1.0.0

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
dnl Ideally we would tell gcc to generate code that is optimized for the CPU it
dnl is compiled on, but -march=native is only available in the more recent
dnl versions of gcc, and autoconf does not provide an easy way to test what C
dnl compiler flags are compatible.
dnl 
dnl if test "x$GCC" = "xyes" ; then
dnl   CFLAGS="${CFLAGS} -march=native"
dnl fi
AC_PROG_CPP

dnl sizeof_ptr is needed by the build system.
AC_CHECK_SIZEOF([void *])
sizeof_ptr="${ac_cv_sizeof_void_p}"
AC_SUBST([sizeof_ptr])

AC_CANONICAL_HOST
dnl CPU-specific settings.
FPUCFLAGS=""
MTCPPFLAGS="-DMEXP=19937"
MTCFLAGS="-fno-strict-aliasing"
case "${host_cpu}" in
  i[[3456]]86)
	AC_DEFINE([CxmCpuIa32])
	CFLAGS="${CFLAGS} -msse2"
	FPUCFLAGS="-mfpmath=sse"
	MTCPPFLAGS="${MTCPPFLAGS} -DHAVE_SSE2=1"
	;;
  x86_64)
	AC_DEFINE([CxmCpuAmd64])
	CFLAGS="${CFLAGS} -msse2"
	FPUCFLAGS="-mfpmath=sse"
	MTCPPFLAGS="${MTCPPFLAGS} -DHAVE_SSE2=1"
	;;
  powerpc)
	AC_DEFINE([CxmCpuPpc])
	CFLAGS="${CFLAGS} -faltivec"
	MTCPPFLAGS="${MTCPPFLAGS} -DHAVE_ALTIVEC=1"
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
	AC_DEFINE([CxmOsDarwin])
	CFLAGS="$CFLAGS -fno-common -no-cpp-precomp"
	abi="macho"
	RPATH=""
	;;
  *-*-freebsd*)
	AC_DEFINE([CxmOsBSD])
	CFLAGS="$CFLAGS"
	abi="elf"
	RPATH="-Wl,-rpath,"
	;;
  *-*-linux*)
	AC_DEFINE([CxmOsLinux])
	CFLAGS="$CFLAGS"
	abi="elf"
	RPATH="-Wl,-rpath,"
	;;
  *-*-netbsd*)
	AC_DEFINE([CxmOsBSD])
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
	AC_DEFINE([CxmOsSolaris])
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

dnl Support optional additions to rpath.
AC_ARG_WITH([rpath],
  [AS_HELP_STRING([--with-rpath=<rpath>], [colon-separated rpath (ELF systems only)])],
if test "x$with_rpath" = "xno" ; then
  RPATH_EXTRA=
else
  RPATH_EXTRA="`echo $with_rpath | tr \":\" \" \"`"
fi,
  RPATH_EXTRA=
)
AC_SUBST([RPATH_EXTRA])

dnl Disable rules that do automatic regeneration of configure output by default.
AC_ARG_ENABLE([autogen],
  [AS_HELP_STRING[--enable-autogen], [Automatically regenerate configure output])],
if test "x$enable_autogen" = "xno" ; then
  enable_autogen="0"
else
  enable_autogen="1"
fi
,
enable_autogen="0"
)
AC_SUBST([enable_autogen])

dnl Disable flex code generation by default.
AC_ARG_ENABLE([flex],
  [AS_HELP_STRING([--enable-flex], [Compile .l sources to C outputs])],
if test "x$enable_flex" = "xyes" ; then
  enable_flex="1"
else
  enable_flex="0"
  fi
,
enable_flex="0"
)
AC_SUBST([enable_flex])

AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], , [$PATH])
AC_PATH_PROG([LD], [ld], , [$PATH])
AC_PATH_PROG([TEST], [test], , [$PATH])
AC_PATH_PROG([grep], [grep], , [$PATH])
AC_PATH_PROG([AUTOCONF], [autoconf], , [$PATH])

dnl Do not silence irrelevant compiler warnings by default, since enabling this
dnl option incurs a performance penalty.
AC_ARG_ENABLE([cc-silence],
  [AS_HELP_STRING([--enable-cc-silence], [Silence irrelevant compiler warnings])],
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

dnl Do not compile with debugging by default.
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug], [Build debugging code])],
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
  dnl Disabling assertions causes Cython-generated C sources to contain unused
  dnl labels, which causes compiler warnings.  Therefore, leave assertions
  dnl enabled if --enable-cc-silence was specified.
  if test "x$enable_cc_silence" = "x0" ; then
    CPPFLAGS="${CPPFLAGS} -DPYREX_WITHOUT_ASSERTIONS"
  fi
fi
AC_SUBST([enable_debug])

dnl Only optimize if not debugging.
if test "x$enable_debug" = "x0" -a "x$no_CFLAGS" = "xyes" ; then
  dnl Make sure that an optimization flag was not specified in EXTRA_CFLAGS.
  if test "x$GCC" = "xyes" ; then
    echo "$EXTRA_CFLAGS" | grep "\-O" >/dev/null || CFLAGS="$CFLAGS -O3 -funroll-loops -fomit-frame-pointer"
    MTCFLAGS="${MTCFLAGS} --param max-inline-insns-single=1800 --param inline-unit-growth=500 --param large-function-growth=900"
  else
    echo "$EXTRA_CFLAGS" | grep "\-O" >/dev/null || CFLAGS="$CFLAGS -O"
  fi
fi
AC_SUBST([MTCPPFLAGS])
AC_SUBST([MTCFLAGS])

dnl Support inlines by default.
AC_ARG_ENABLE([inlines],
  [AS_HELP_STRING([--disable-inlines], [Do not use inline functions])],
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
dnl Configure pthreads.

AC_CHECK_HEADERS([pthread.h], , [AC_MSG_ERROR([pthread.h is missing])])
AC_CHECK_LIB([pthread], [pthread_create], [LIBS="$LIBS -lpthread"],
  [AC_CHECK_LIB([c_r], [pthread_create], 
    [LIBS="$LIBS -pthread"],
    [AC_MSG_ERROR([libpthread appears non-functional])])])

CPPFLAGS="$CPPFLAGS -D_REENTRANT"

dnl ============================================================================
dnl Make sure the version of python found is sufficient.

dnl Get path to python.
AC_ARG_WITH([python], [AS_HELP_STRING([--with-python=<path>], [python program])],
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

dnl CX_VERSION_CHECK([prog], [version], [version_min], [result])
AC_DEFUN([CX_VERSION_CHECK], [
CX_VERSION_PROG=$1
CX_VERSION=$2
CX_VERSION_MIN=$3
AC_MSG_CHECKING([for $CX_VERSION_PROG >= $CX_VERSION_MIN])
CX_VERSION_MAJOR=`echo ${CX_VERSION} | tr . " " | awk '{print [$]1}'`
CX_VERSION_MINOR=`echo ${CX_VERSION} | tr . " " | awk '{print [$]2}' | tr a-z " " |awk '{print [$]1}'`
CX_VERSION_BRANCH=`echo ${CX_VERSION} | tr . " " | awk '{print [$]3}' | tr a-z " " |awk '{print [$]1}'`
if test "x$CX_VERSION_BRANCH" = "x" ; then
  CX_VERSION_BRANCH=0
fi
CX_VERSION_MAJOR_MIN=`echo $CX_VERSION_MIN | tr . " " | awk '{print [$]1}'`
CX_VERSION_MINOR_MIN=`echo $CX_VERSION_MIN | tr . " " | awk '{print [$]2}'`
CX_VERSION_BRANCH_MIN=`echo $CX_VERSION_MIN | tr . " " | awk '{print [$]3}'`
if test $CX_VERSION_MAJOR -gt $CX_VERSION_MAJOR_MIN \
     -o $CX_VERSION_MAJOR -eq $CX_VERSION_MAJOR_MIN \
     -a $CX_VERSION_MINOR -gt $CX_VERSION_MINOR_MIN \
     -o $CX_VERSION_MAJOR -eq $CX_VERSION_MAJOR_MIN \
     -a $CX_VERSION_MINOR -eq $CX_VERSION_MINOR_MIN \
     -a $CX_VERSION_BRANCH -ge $CX_VERSION_BRANCH_MIN ; then
  AC_MSG_RESULT([$CX_VERSION_MAJOR.$CX_VERSION_MINOR.$CX_VERSION_BRANCH])
  $4="1"
else
  AC_MSG_RESULT([no ($CX_VERSION_MAJOR.$CX_VERSION_MINOR.$CX_VERSION_BRANCH)])
  $4="0"
fi
])

PYTHON_FOUND=`2>&1 $PYTHON -V | awk '{print $2}'`
CX_VERSION_CHECK([python], [$PYTHON_FOUND], [$PYTHON_VERSION_MIN],
  [have_python])
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
site_packages="`${PYTHON} conftest`"
rm -f conftest
AC_MSG_RESULT([${site_packages}])
AC_SUBST([site_packages])

dnl ============================================================================
dnl Make sure the version of cython found is sufficient.

dnl Get path to cython.
AC_ARG_WITH([cython], [AS_HELP_STRING([--with-cython=<path>], [cython program])],
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
  if test -x "${objroot}contrib/cython/bin/cython" ; then
    dnl Use the Cython that is shipped with Crux releases.
    CYTHON="${abs_objroot}cython"
  else
    dnl Try to find an installed version of Cython.
    AC_PATH_PROG([CYTHON], [cython], , [$PATH])
    if test "x$CYTHON" = "x" ; then
      AC_MSG_ERROR([Cannot build without cython])
    fi
  fi
fi

dnl Avoid testing the version if using contrib/cython, since the wrapper script
dnl hasn't been created yet.
if test "x${CYTHON}" != "x${abs_objroot}cython" ; then
  CYTHON_FOUND=`2>&1 $CYTHON -V | awk '{print $3}'`
  CX_VERSION_CHECK([cython], [$CYTHON_FOUND], [$CYTHON_VERSION_MIN],
    [have_cython])
  if test "x$have_cython" = "x0" ; then
    AC_MSG_ERROR([Incompatible cython])
  fi
fi
AC_SUBST([CYTHON])

dnl ============================================================================
dnl Use mpi4py by default, if it can be found.
AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--disable-mpi], [Do not use MPI])],
[if test "x$enable_mpi" = "xno" ; then
  enable_mpi="0"
else
  enable_mpi="1"
fi
],
[enable_mpi="1"]
)
if test "x$enable_mpi" = "x1" ; then
  AC_CHECK_HEADERS([mpi.h], , [enable_mpi=0])
fi
mpi4py_includes=""
if test "x$enable_mpi" = "x1" ; then
  dnl Check whether mpi4py is installed.
  AC_MSG_CHECKING([for mpi4py])
  [
  1>conftest cat << EOF
try:
    import mpi4py
    print "1"
except:
    print "0"
EOF
]
  enable_mpi="`${PYTHON} conftest`"
  rm -f conftest
  if test "x$enable_mpi" = "x1" ; then
    AC_MSG_RESULT([yes])
    dnl Get the mpi4py version.
    [
    1>conftest cat << EOF
import mpi4py
print mpi4py.__version__
EOF
]
    MPI4PY_FOUND="`${PYTHON} conftest`"
    rm -f conftest
    CX_VERSION_CHECK([mpi4py], [$MPI4PY_FOUND], [$MPI4PY_VERSION_MIN],
      [enable_mpi])
  else
    AC_MSG_RESULT([no])
  fi

  if test "x$enable_mpi" = "x1" ; then
    dnl Get the include path.
    AC_MSG_CHECKING([mpi4py include path])
    [
    1>conftest cat << EOF
import mpi4py
print mpi4py.get_include()
EOF
]
    mpi4py_includes="`${PYTHON} conftest`"
    rm -f conftest
    AC_MSG_RESULT([$mpi4py_includes])
  fi
fi
AC_SUBST([enable_mpi])
AC_SUBST([mpi4py_includes])

dnl ============================================================================
dnl Make sure the version of flex found is sufficient.
if test "x${enable_flex}" = "x1" ; then
  dnl Get path to flex.
  AC_ARG_WITH([flex], [AS_HELP_STRING([--with-flex=<path>], [flex program])],
  if test "x$with_flex" = "xno" ; then
    FLEX=
  else
    FLEX="$with_flex"
  fi
  ,
    FLEX=
  )
  if test "x$FLEX" != "x" ; then
     if test ! -x "${FLEX}" ; then
       AC_MSG_ERROR([Unusable or missing flex: ${FLEX}])
     fi
  else
    AC_PATH_PROG([FLEX], [flex], , [$PATH])
    if test "x$FLEX" = "x" ; then
      AC_MSG_ERROR([Cannot build without flex])
    fi
  fi

  FLEX_FOUND=`2>&1 $FLEX --version | awk '{print $2}'`
  CX_VERSION_CHECK([flex], [$FLEX_FOUND], [$FLEX_VERSION_MIN],
    [have_flex])
  if test "x$have_flex" = "x0" ; then
    AC_MSG_ERROR([Cannot build without flex])
  fi
  AC_SUBST([FLEX])
fi

dnl ============================================================================
dnl Crux configuration.
dnl 
crux_version=`cat ${srcroot}VERSION`
AC_SUBST([crux_version])

dnl ============================================================================
dnl Generate file lists.

dnl C headers.
for h in `cd ./${srcroot}; find lib -type f |grep "\.h\.in$" |sort` ; do
  chdr="`dirname ${h}`/`basename ${h} .in`"
  cfghdrs="${cfghdrs} ${chdr}"
done

dnl Flex sources.  csrcs uses lsrcs to avoid multiple listing, so lsrcs must be
dnl generated first.
lsrcs=""
for f in `find ${srcroot}lib -type f |grep "/Cx[[A-Za-z0-9]]*\.l$" |sort` ; do
  lsrcs="${lsrcs} ${f}"
done
AC_SUBST([lsrcs])

dnl C sources.
csrcs=""
for s in `find ${srcroot}lib -type f |grep "/[[A-Za-z0-9]]*\.c$" |sort` ; do
  is_fsrc="0"
  for f in ${lsrcs} ; do
    c="`dirname ${f}`/`basename ${f} .l`.c"
    if test "x${s}" = "x${c}" ; then
      is_fsrc="1"
    fi
  done
  if test "x$is_fsrc" = "x0" ; then
    csrcs="${csrcs} ${s}"
  fi
done
AC_SUBST([csrcs])

dnl Fortran sources.
fsrcs=""
for s in `find ${srcroot}lib -type f |grep "/[[A-Za-z0-9]]*\.f$" |sort` ; do
  fsrcs="${fsrcs} ${s}"
done
AC_SUBST([fsrcs])

dnl Cython headers that aren't associated with .pyx files.
libCx_pxds=""
for d in `find ${srcroot}lib -type f |grep "/[[A-Za-z0-9_]]*\.pxd$" |sort` ; do
  libCx_pxds="${libCx_pxds} ${d}"
done
AC_SUBST([libCx_pxds])

dnl Cython module/package sources.
pyxs=""
pkgs=""
for s in `find ${srcroot}pkg -type f |grep "/[[A-Za-z0-9_]]*\.pyx$" |sort` ; do
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
for s in `find ${srcroot}pkg -type f |grep "/[[A-Za-z0-9_]]*\.py$" |sort` ; do
  pys="${pys} ${s}"
done
AC_SUBST([pys])

dnl Executables.
bins=""
for b in `find ${srcroot}bin -type f |grep "/[[A-Za-z0-9]]*$" |sort` ; do
  bins="${bins} ${b}"
done
AC_SUBST([bins])

dnl Tests.
tests=""
for t in `cd ./${srcroot}; find test -type f |grep "/[[A-Za-z0-9_]]*\.cx$" |sort`; do
  tests="${tests} ${t}"
done
AC_SUBST([tests])

dnl ============================================================================
dnl Check for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_HEADER_STDBOOL

AC_CHECK_LIB([m], [log], [], [AC_MSG_ERROR([Cannot build without libm])])

AC_CHECK_FUNC([posix_memalign],
  [have_posix_memalign=1; AC_DEFINE([CxmHavePosixMemalign])],
  [have_posix_memalign=0])
AC_SUBST([have_posix_memalign])

AC_MSG_CHECKING([for mallopt])
AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <malloc.h>
]], [[
    mallopt(M_MMAP_THRESHOLD, 128*1024);

    return 0;
]])],
              AC_MSG_RESULT([yes])
              AC_DEFINE([CxmHaveMallopt]),
              AC_MSG_RESULT([no]))

dnl ============================================================================
dnl Use system ATLAS/LAPACK by default.
AC_ARG_ENABLE([sys_linalg],
  [AS_HELP_STRING([--disable-sys-linalg],
   [Do not use system ATLAS/LAPACK libraries])],
[if test "x$enable_sys_linalg" = "xno" ; then
  enable_sys_linalg="0"
else
  enable_sys_linalg="1"
fi
],
[enable_sys_linalg="1"]
)
did_fc="0"
if test "x$enable_sys_linalg" = "x1" ; then
  TLIBS="${LIBS}"

  LIBS="-framework vecLib $LIBS"
  AC_CHECK_FUNC([dgetrf],
    [have_vecLib="yes"],
    [have_vecLib="no"; LIBS="$TLIBS"])
  if test "x$have_vecLib" = "xno" ; then
    if test "x$did_fc" = "x0" ; then
      did_fc="1"
      AC_PROG_FC([gfortran])
      AC_FC_LIBRARY_LDFLAGS
    fi

    if test "x$enable_sys_linalg" = "x1" ; then
      AC_CHECK_LIB([atlas], [ATL_dscal], [], [])
    fi

    if test "x$enable_sys_linalg" = "x1" ; then
      TLIBS="${LIBS}"
      LIBS="${LIBS} ${FCLIBS}"
      AC_F77_FUNC([dscal])
      AC_CHECK_LIB([f77blas], $dscal, [],
        [LIBS="${TLIBS}"; enable_sys_linalg="0"])
      if test "x$enable_sys_linalg" = "x0" ; then
        AC_CHECK_LIB([blas], $dscal, [enable_sys_linalg="1"])
      fi
    fi

    if test "x$enable_sys_linalg" = "x1" ; then
      TLIBS="${LIBS}"
      LIBS="${LIBS} ${FCLIBS}"
      AC_F77_FUNC([dgetrf])
      AC_CHECK_LIB([lapack], $dgetrf, [],
        [LIBS="${TLIBS}"; enable_sys_linalg="0"])
      if test "x$enable_sys_linalg" = "x0" ; then
        AC_CHECK_LIB([lapack_atlas], [clapack_dgetrf],
          [enable_sys_linalg="1"; LIBS="-llapack_atlas -llapack $LIBS"])
      fi
    fi
  fi

  if test "x$enable_sys_linalg" = "x0" ; then
    dnl Restore LIBS.
    LIBS="${TLIBS}"
  fi
fi

if test "x$enable_sys_linalg" = "x0" ; then
  if test "x$did_fc" = "x0" ; then
    did_fc="1"
    AC_PROG_FC([gfortran])
    AC_FC_LIBRARY_LDFLAGS
  fi
  dnl Check whether gfortran is available.
  if test "x$FC" = "x" ; then
    AC_MSG_ERROR([Fortran compiler required to build BLAS and LAPACK routines])
  fi
  dnl Set gfortran-specific flags.
dnl See notes re: -march=native for CFLAGS above.
dnl   echo "${FC}" | grep "gfortran">/dev/null \
dnl     && FCFLAGS="${FCFLAGS} -march=native"
  dnl Enable FPU-related [F]CFLAGS.  It is dangerous to do this if we are using
  dnl system libraries (ATLAS/LAPACK) that may have been compiled to use a
  dnl different FPU.
  CFLAGS="${CFLAGS} ${FPUCFLAGS}"
  FCFLAGS="${FCFLAGS} ${FPUCFLAGS}"
fi
AC_SUBST([enable_sys_linalg])

dnl Process .in files.
AC_SUBST([cfghdrs])
AC_CONFIG_HEADERS([$cfghdrs cfghdrs.stamp])

AC_CONFIG_FILES([$cfgoutputs cfgoutputs.stamp])
AC_DEFUN([CX_CONFIG_BINS], [
AC_CONFIG_FILES([$1], [chmod +x $1])
cfgoutputs="${cfgoutputs} $1"
])
CX_CONFIG_BINS([cfg])
CX_CONFIG_BINS([cydep])
CX_CONFIG_BINS([cython])
AC_SUBST([cfgoutputs])
AC_OUTPUT

dnl ============================================================================
dnl Print out the results of configuration.
AC_MSG_RESULT([===============================================================================])
AC_MSG_RESULT([Crux version  : $crux_version])
AC_MSG_RESULT([])
AC_MSG_RESULT([PYTHON        : ${PYTHON}])
AC_MSG_RESULT([CYTHON        : ${CYTHON}])
AC_MSG_RESULT([FLEX          : ${FLEX}])
AC_MSG_RESULT([CC            : ${CC}])
AC_MSG_RESULT([CPPFLAGS      : ${CPPFLAGS}])
AC_MSG_RESULT([MTCPPFLAGS    : ${MTCPPFLAGS}])
AC_MSG_RESULT([CFLAGS        : ${CFLAGS}])
AC_MSG_RESULT([MTCFLAGS      : ${MTCFLAGS}])
AC_MSG_RESULT([FC            : ${FC}])
AC_MSG_RESULT([FCFLAGS       : ${FCFLAGS}])
AC_MSG_RESULT([LDFLAGS       : ${LDFLAGS}])
AC_MSG_RESULT([LIBS          : ${LIBS}])
AC_MSG_RESULT([FCLIBS        : ${FCLIBS}])
AC_MSG_RESULT([RPATH_EXTRA   : ${RPATH_EXTRA}])
AC_MSG_RESULT([])
AC_MSG_RESULT([site_packages : ${site_packages}])
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
AC_MSG_RESULT([flex          : ${enable_flex}])
AC_MSG_RESULT([debug         : ${enable_debug}])
AC_MSG_RESULT([cc-silence    : ${enable_cc_silence}])
AC_MSG_RESULT([inlines       : ${enable_inlines}])
AC_MSG_RESULT([sys-linalg    : ${enable_sys_linalg}])
AC_MSG_RESULT([mpi           : ${enable_mpi}])
AC_MSG_RESULT([===============================================================================])
