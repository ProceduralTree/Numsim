/* adio/include/romioconf.h.  Generated from romioconf.h.in by configure.  */
/* adio/include/romioconf.h.in.  Generated from configure.ac by autoheader.  */

/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */
#ifndef ROMIOCONF_H_INCLUDED
#define ROMIOCONF_H_INCLUDED
#include <mplconfig.h>

#include "romioconf-undefs.h"


/* BGQ platform */
/* #undef BGQPLATFORM */

/* Define to 1 if you have the <aio.h> header file. */
/* #undef HAVE_AIO_H */

/* Define to 1 if you have the <aio-lite.h> header file. */
/* #undef HAVE_AIO_LITE_H */

/* Define to 1 if you have the <daos_types.h> header file. */
/* #undef HAVE_DAOS_TYPES_H */

/* Define to 1 if you have the declaration of `MPI_COMBINER_HINDEXED_BLOCK',
   and to 0 if you don't. */
#define HAVE_DECL_MPI_COMBINER_HINDEXED_BLOCK 1

/* Define to 1 if you have the declaration of `pwrite', and to 0 if you don't.
   */
#define HAVE_DECL_PWRITE 1

/* Define to 1 if you have the <dirent.h> header file. */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `fsync' function. */
#define HAVE_FSYNC 1

/* Define to 1 if you have the `ftruncate' function. */
#define HAVE_FTRUNCATE 1

/* Define if GNU __attribute__ is supported */
#define HAVE_GCC_ATTRIBUTE 1

/* Define to 1 if you have the <gpfs_fcntl.h> header file. */
/* #undef HAVE_GPFS_FCNTL_H */

/* Define to 1 if you have the <gpfs.h> header file. */
/* #undef HAVE_GPFS_H */

/* Define to 1 if you have the <gurt/hash.h> header file. */
/* #undef HAVE_GURT_HASH_H */

/* Define to 1 if you have the <ime_native.h> header file. */
/* #undef HAVE_IME_NATIVE_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if int smaller than pointer */
#define HAVE_INT_LT_POINTER 1

/* Define to 1 if you have the `aio-lite' library (-laio-lite). */
/* #undef HAVE_LIBAIO_LITE */

/* Define to 1 if you have the `daos' library (-ldaos). */
/* #undef HAVE_LIBDAOS */

/* Define to 1 if you have the `daos_common' library (-ldaos_common). */
/* #undef HAVE_LIBDAOS_COMMON */

/* Define to 1 if you have the `dfs' library (-ldfs). */
/* #undef HAVE_LIBDFS */

/* Define to 1 if you have the `duns' library (-lduns). */
/* #undef HAVE_LIBDUNS */

/* Define to 1 if you have the `gurt' library (-lgurt). */
/* #undef HAVE_LIBGURT */

/* Define to 1 if you have the `uuid' library (-luuid). */
/* #undef HAVE_LIBUUID */

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the <linux/lustre/lustre_user.h> header file. */
/* #undef HAVE_LINUX_LUSTRE_LUSTRE_USER_H */

/* Define if long long is 64 bits */
#define HAVE_LONG_LONG_64 1

/* Define to 1 if you have the `lseek64' function. */
#define HAVE_LSEEK64 1

/* Define to 1 if you have the `lstat' function. */
#define HAVE_LSTAT 1

/* Define if LUSTRE_LOCKAHEAD is enabled. */
/* #undef HAVE_LUSTRE_LOCKAHEAD */

/* Define to 1 if you have the <lustre/lustre_user.h> header file. */
/* #undef HAVE_LUSTRE_LUSTRE_USER_H */

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the `memalign' function. */
#define HAVE_MEMALIGN 1

/* Define to 1 if you have the <minix/config.h> header file. */
/* #undef HAVE_MINIX_CONFIG_H */

/* Set if MPI-IO prototypes use const qualifier */
#define HAVE_MPIIO_CONST 1

/* Define to 1 if you have the <mpix.h> header file. */
/* #undef HAVE_MPIX_H */

/* Define if MPI combiners available */
/* #undef HAVE_MPI_COMBINERS */

/* Define if Darray is available */
#define HAVE_MPI_DARRAY_SUBARRAY 1

/* Define if generalized requests avaliable */
/* #undef HAVE_MPI_GREQUEST */

/* Define if MPI Info is available */
#define HAVE_MPI_INFO 1

/* Define if supports long long int */
#define HAVE_MPI_LONG_LONG_INT 1

/* Will always be 1 - OMPI has MPI_OFFSET */
#define HAVE_MPI_OFFSET 1

/* Define to 1 if you have the `MPI_Status_set_elements_x' function. */
#define HAVE_MPI_STATUS_SET_ELEMENTS_X 1

/* Define to 1 if you have the `MPI_Type_size_x' function. */
/* #undef HAVE_MPI_TYPE_SIZE_X */

/* Define if multiple weak symbols may be defined */
/* #undef HAVE_MULTIPLE_PRAGMA_WEAK */

/* Define to 1 if the system has the type `pan_fs_client_raidn_encoding_t'. */
/* #undef HAVE_PAN_FS_CLIENT_RAIDN_ENCODING_T */

/* Cray style weak pragma */
/* #undef HAVE_PRAGMA_CRI_DUP */

/* HP style weak pragma */
/* #undef HAVE_PRAGMA_HP_SEC_DEF */

/* Supports weak pragma */
/* #undef HAVE_PRAGMA_WEAK */

/* Define if PVFS_sys_create does not have layout parameter */
/* #undef HAVE_PVFS2_CREATE_WITHOUT_LAYOUT */

/* Define to 1 if you have the <pvfs2.h> header file. */
/* #undef HAVE_PVFS2_H */

/* Define if PVFS2_SUPER_MAGIC defined. */
/* #undef HAVE_PVFS2_SUPER_MAGIC */

/* Define to 1 if you have the `readlink' function. */
#define HAVE_READLINK 1

/* Define to 1 if you have the <signal.h> header file. */
/* #undef HAVE_SIGNAL_H */

/* Define to 1 if you have the `stat' function. */
#define HAVE_STAT 1

/* Define to 1 if you have the `statfs' function. */
#define HAVE_STATFS 1

/* Define if status set bytes available */
#define HAVE_STATUS_SET_BYTES 1

/* Define to 1 if you have the `statvfs' function. */
#define HAVE_STATVFS 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strdup' function. */
#define HAVE_STRDUP 1

/* Define to 1 if you have the `strerror' function. */
#define HAVE_STRERROR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if `aio_fildes' is a member of `struct aiocb'. */
/* #undef HAVE_STRUCT_AIOCB_AIO_FILDES */

/* Define to 1 if `aio_handle' is a member of `struct aiocb'. */
/* #undef HAVE_STRUCT_AIOCB_AIO_HANDLE */

/* Define to 1 if `aio_reqprio' is a member of `struct aiocb'. */
/* #undef HAVE_STRUCT_AIOCB_AIO_REQPRIO */

/* Define to 1 if `aio_sigevent' is a member of `struct aiocb'. */
/* #undef HAVE_STRUCT_AIOCB_AIO_SIGEVENT */

/* Define to 1 if `aio_whence' is a member of `struct aiocb'. */
/* #undef HAVE_STRUCT_AIOCB_AIO_WHENCE */

/* Define if struct statfs can be compiled */
#define HAVE_STRUCT_STATFS 1

/* Define to 1 if you have the <sys/aio.h> header file. */
/* #undef HAVE_SYS_AIO_H */

/* Define to 1 if you have the <sys/attr.h> header file. */
/* #undef HAVE_SYS_ATTR_H */

/* Define to 1 if you have the <sys/mount.h> header file. */
#define HAVE_SYS_MOUNT_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/statvfs.h> header file. */
#define HAVE_SYS_STATVFS_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/vfs.h> header file. */
#define HAVE_SYS_VFS_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `usleep' function. */
#define HAVE_USLEEP 1

/* Whether C compiler supports symbol visibility or not */
#define HAVE_VISIBILITY 1

/* Define to 1 if you have the <wchar.h> header file. */
#define HAVE_WCHAR_H 1

/* Attribute style weak pragma */
/* #undef HAVE_WEAK_ATTRIBUTE */

/* Define if weak symbols available */
/* #undef HAVE_WEAK_SYMBOLS */

/* Define if int smaller than pointer */
#define INT_LT_POINTER 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define if using MPICH */
/* #undef MPICH */

/* Define if using HP MPI */
/* #undef MPIHP */

/* hack to make ROMIO build without profiling */
#define MPIO_BUILD_PROFILING 1

/* Define if SGI MPI */
/* #undef MPISGI */

/* Define to "MPI_Aint" if MPI does not provide MPI_Count */
/* #undef MPI_Count */

/* Define if MPI_Offset is int */
/* #undef MPI_OFFSET_IS_INT */

/* Define if fsync needs a declaration */
/* #undef NEEDS_FSYNC_DECL */

/* Define if ftruncate needs a declaration */
/* #undef NEEDS_FTRUNCATE_DECL */

/* Define if l_start and l_len data should be cast as int */
/* #undef NEEDS_INT_CAST_WITH_FLOCK */

/* Define if lseek64 needs a declaration */
/* #undef NEEDS_LSEEK64_DECL */

/* Define if lstat needs a declaration */
/* #undef NEEDS_LSTAT_DECL */

/* Define if mpi_test needed */
/* #undef NEEDS_MPI_TEST */

/* Define if readlink needs a declaration */
/* #undef NEEDS_READLINK_DECL */

/* Define if strdup needs a declaration */
/* #undef NEEDS_STRDUP_DECL */

/* Define if usleep needs a declaration */
/* #undef NEEDS_USLEEP_DECL */

/* Define if no MPI type is contig */
/* #undef NO_MPI_SGI_type_is_contig */

/* Name of package */
#define PACKAGE "romio"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "discuss@mpich.org"

/* Define to the full name of this package. */
#define PACKAGE_NAME "ROMIO"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "ROMIO Open MPI"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "romio"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.mpich.org/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "Open MPI"

/* PE platform */
/* #undef PEPLATFORM */

/* Define for printing error messages */
/* #undef PRINT_ERR_MSG */

/* Define for ROMIO with DAOS */
/* #undef ROMIO_DAOS */

/* Define for ROMIO with GPFS */
/* #undef ROMIO_GPFS */

/* Define if AIO calls need file descriptor */
/* #undef ROMIO_HAVE_AIO_CALLS_NEED_FILEDES */

/* Define if aio_suspend needs two arguments */
/* #undef ROMIO_HAVE_AIO_SUSPEND_TWO_ARGS */

/* Define if statfs has f_fstypename */
/* #undef ROMIO_HAVE_STRUCT_STATFS_WITH_F_FSTYPENAME */

/* Define if statfs has f_type */
#define ROMIO_HAVE_STRUCT_STATFS_WITH_F_TYPE 1

/* defined if struct statvfs has a f_basetype member */
/* #undef ROMIO_HAVE_STRUCT_STATVFS_WITH_F_BASETYPE */

/* Define if struct stat has a st_fstype member */
/* #undef ROMIO_HAVE_STRUCT_STAT_WITH_ST_FSTYPE */

/* Define if AIO calls seem to work */
/* #undef ROMIO_HAVE_WORKING_AIO */

/* Define for ROMIO with IME */
/* #undef ROMIO_IME */

/* Define if compiling within MPICH */
/* #undef ROMIO_INSIDE_MPICH */

/* Define for ROMIO with LUSTRE */
/* #undef ROMIO_LUSTRE */

/* Define for ROMIO with NFS */
#define ROMIO_NFS 1

/* Define for ROMIO with PANFS */
/* #undef ROMIO_PANFS */

/* Define for ROMIO with PVFS2 */
/* #undef ROMIO_PVFS2 */

/* Define for ROMIO with QUOBYTEFS */
/* #undef ROMIO_QUOBYTEFS */

/* Define if run on Linux */
#define ROMIO_RUN_ON_LINUX 1

/* Define for ROMIO with TESTFS */
#define ROMIO_TESTFS 1

/* Define for ROMIO with UFS */
#define ROMIO_UFS 1

/* Define for ROMIO with XFS */
/* #undef ROMIO_XFS */

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Define if performing coverage tests */
/* #undef USE_COVERAGE */

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable general extensions on macOS.  */
#ifndef _DARWIN_C_SOURCE
# define _DARWIN_C_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable X/Open compliant socket functions that do not require linking
   with -lxnet on HP-UX 11.11.  */
#ifndef _HPUX_ALT_XOPEN_SOCKET_API
# define _HPUX_ALT_XOPEN_SOCKET_API 1
#endif
/* Identify the host operating system as Minix.
   This macro does not affect the system headers' behavior.
   A future release of Autoconf may stop defining this macro.  */
#ifndef _MINIX
/* # undef _MINIX */
#endif
/* Enable general extensions on NetBSD.
   Enable NetBSD compatibility extensions on Minix.  */
#ifndef _NETBSD_SOURCE
# define _NETBSD_SOURCE 1
#endif
/* Enable OpenBSD compatibility extensions on NetBSD.
   Oddly enough, this does nothing on OpenBSD.  */
#ifndef _OPENBSD_SOURCE
# define _OPENBSD_SOURCE 1
#endif
/* Define to 1 if needed for POSIX-compatible behavior.  */
#ifndef _POSIX_SOURCE
/* # undef _POSIX_SOURCE */
#endif
/* Define to 2 if needed for POSIX-compatible behavior.  */
#ifndef _POSIX_1_SOURCE
/* # undef _POSIX_1_SOURCE */
#endif
/* Enable POSIX-compatible threading on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-5:2014.  */
#ifndef __STDC_WANT_IEC_60559_ATTRIBS_EXT__
# define __STDC_WANT_IEC_60559_ATTRIBS_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-1:2014.  */
#ifndef __STDC_WANT_IEC_60559_BFP_EXT__
# define __STDC_WANT_IEC_60559_BFP_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-2:2015.  */
#ifndef __STDC_WANT_IEC_60559_DFP_EXT__
# define __STDC_WANT_IEC_60559_DFP_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-4:2015.  */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TS 18661-3:2015.  */
#ifndef __STDC_WANT_IEC_60559_TYPES_EXT__
# define __STDC_WANT_IEC_60559_TYPES_EXT__ 1
#endif
/* Enable extensions specified by ISO/IEC TR 24731-2:2010.  */
#ifndef __STDC_WANT_LIB_EXT2__
# define __STDC_WANT_LIB_EXT2__ 1
#endif
/* Enable extensions specified by ISO/IEC 24747:2009.  */
#ifndef __STDC_WANT_MATH_SPEC_FUNCS__
# define __STDC_WANT_MATH_SPEC_FUNCS__ 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable X/Open extensions.  Define to 500 only if necessary
   to make mbstate_t available.  */
#ifndef _XOPEN_SOURCE
/* # undef _XOPEN_SOURCE */
#endif


/* Define if weak symbols should be used */
/* #undef USE_WEAK_SYMBOLS */

/* Version number of package */
#define VERSION "Open MPI"

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Provide blksize_t if not available */
/* #undef blksize_t */

/* Define to `long int' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported only directly.  */
#define restrict __restrict__
/* Work around a bug in older versions of Sun C++, which did not
   #define __restrict__ or support _Restrict or __restrict__
   even though the corresponding Sun C compiler ended up with
   "#define restrict _Restrict" or "#define restrict __restrict__"
   in the previous line.  This workaround can be removed once
   we assume Oracle Developer Studio 12.5 (2016) or later.  */
#if defined __SUNPRO_CC && !defined __RESTRICT && !defined __restrict__
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

/* Define to "unsigned char" if sys/types.h does not define. */
/* #undef u_char */

/* Define to "unsigned int" if sys/types.h does not define. */
/* #undef u_int */

/* Define to "unsigned long" if sys/types.h does not define. */
/* #undef u_long */

/* Define to "unsigned short" if sys/types.h does not define. */
/* #undef u_short */


/* quash PACKAGE and PACKAGE_* vars, see MPICH top-level configure.ac for
 * more info */
#include "nopackage.h"

#endif /* !defined(ROMIOCONF_H_INCLUDED) */

