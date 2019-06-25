AC_INIT([phylogen], [1.0], [dmulder@suse.com])
AM_INIT_AUTOMAKE()
LT_INIT
AC_PROG_CXX
AM_PATH_PYTHON([3.6])
AX_PKG_SWIG([], [], AC_MSG_ERROR([SWIG is required to build.]))
AX_SWIG_PYTHON
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
