AC_INIT([hybrid_coal],[beta],[sha.joe.zhu@gmail.com])
#AC_INIT([hybrid_sim],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[joe.zhu@pg.canterbury.ac.nz])
dnl AC_CONFIG_SRCDIR([simu.c])

AM_INIT_AUTOMAKE
AC_PROG_INSTALL

AC_CANONICAL_HOST

# Checks for programs.
#AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for libraries.
#AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_LANG(C++) 
#AC_CHECK_HEADERS([iostream string sstream fstream vector iomanip valarray math.h stdio.h ctime])
AC_CHECK_HEADERS([iostream string sstream fstream vector iomanip valarray math.h stdio.h ctime gmp.h])
# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([pow])
#AC_CHECK_FUNCS([factorial])
AC_HEADER_STDBOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

