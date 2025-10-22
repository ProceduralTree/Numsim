#ifndef INCLUDE_MPLCONFIG_H
#define INCLUDE_MPLCONFIG_H 1

/* include/mplconfig.h. Generated automatically at end of configure. */
/* include/config.h.  Generated from config.h.in by configure.  */
/* include/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define cache-line size. */
#ifndef MPL_CACHELINE_SIZE
#define MPL_CACHELINE_SIZE 64
#endif

/* Defined the keyword for thread-local storage. */
#ifndef MPL_COMPILER_TLS
#define MPL_COMPILER_TLS _Thread_local
#endif

/* Define to 1 if MPL enables MPL_aligned_alloc. */
#ifndef MPL_DEFINE_ALIGNED_ALLOC
#define MPL_DEFINE_ALIGNED_ALLOC 1
#endif

/* Define if force compiler to always inline functions with
   MPL_STATIC_INLINE_PREFIX|SUFFIX */
/* #undef MPL_ENABLE_ALWAYS_INLINE */

/* Define to 1 if you have the `aligned_alloc' function. */
#ifndef MPL_HAVE_ALIGNED_ALLOC
#define MPL_HAVE_ALIGNED_ALLOC 1
#endif

/* Define to 1 if you have the <arpa/inet.h> header file. */
#ifndef MPL_HAVE_ARPA_INET_H
#define MPL_HAVE_ARPA_INET_H 1
#endif

/* Defined if backtrace() could be fully identified. */
#ifndef MPL_HAVE_BACKTRACE
#define MPL_HAVE_BACKTRACE 1
#endif

/* define if valgrind is old and/or broken compared to what we are expecting
   */
#ifndef MPL_HAVE_BROKEN_VALGRIND
#define MPL_HAVE_BROKEN_VALGRIND 1
#endif

/* Define to 1 if the compiler supports __builtin_expect. */
#ifndef MPL_HAVE_BUILTIN_EXPECT
#define MPL_HAVE_BUILTIN_EXPECT 1
#endif

/* Define to 1 if we have support for C11 atomic intrinsics */
#ifndef MPL_HAVE_C11_ATOMICS
#define MPL_HAVE_C11_ATOMICS 1
#endif

/* Define if C11 _Static_assert is supported. */
#ifndef MPL_HAVE_C11__STATIC_ASSERT
#define MPL_HAVE_C11__STATIC_ASSERT 1
#endif

/* Define to 1 if you have the `clock_getres' function. */
#ifndef MPL_HAVE_CLOCK_GETRES
#define MPL_HAVE_CLOCK_GETRES 1
#endif

/* Define to 1 if you have the `clock_gettime' function. */
#ifndef MPL_HAVE_CLOCK_GETTIME
#define MPL_HAVE_CLOCK_GETTIME 1
#endif

/* Define to 1 if you have the <ctype.h> header file. */
#ifndef MPL_HAVE_CTYPE_H
#define MPL_HAVE_CTYPE_H 1
#endif

/* Define if CUDA is available */
/* #undef MPL_HAVE_CUDA */

/* Define to 1 if you have the <dlfcn.h> header file. */
#ifndef MPL_HAVE_DLFCN_H
#define MPL_HAVE_DLFCN_H 1
#endif

/* Define to 1 if you have the <drd.h> header file. */
/* #undef MPL_HAVE_DRD_H */

/* Define to 1 if you have the <errno.h> header file. */
#ifndef MPL_HAVE_ERRNO_H
#define MPL_HAVE_ERRNO_H 1
#endif

/* Define to 1 if you have the <execinfo.h> header file. */
#ifndef MPL_HAVE_EXECINFO_H
#define MPL_HAVE_EXECINFO_H 1
#endif

/* Define to 1 if you have the `fdopen' function. */
#ifndef MPL_HAVE_FDOPEN
#define MPL_HAVE_FDOPEN 1
#endif

/* Define to 1 if the system has the `fallthrough' function attribute */
#ifndef MPL_HAVE_FUNC_ATTRIBUTE_FALLTHROUGH
#define MPL_HAVE_FUNC_ATTRIBUTE_FALLTHROUGH 1
#endif

/* Define if GNU __attribute__ is supported */
#ifndef MPL_HAVE_GCC_ATTRIBUTE
#define MPL_HAVE_GCC_ATTRIBUTE 1
#endif

/* Define to 1 if we have support for gcc __atomic intrinsics */
#ifndef MPL_HAVE_GCC_INTRINSIC_ATOMIC
#define MPL_HAVE_GCC_INTRINSIC_ATOMIC 1
#endif

/* Define to 1 if we have support for gcc __sync intrinsics */
#ifndef MPL_HAVE_GCC_INTRINSIC_SYNC
#define MPL_HAVE_GCC_INTRINSIC_SYNC 1
#endif

/* Define to 1 if you have the `gethrtime' function. */
/* #undef MPL_HAVE_GETHRTIME */

/* Define to 1 if you have the `getifaddrs' function. */
#ifndef MPL_HAVE_GETIFADDRS
#define MPL_HAVE_GETIFADDRS 1
#endif

/* Define to 1 if you have the `getpid' function. */
#ifndef MPL_HAVE_GETPID
#define MPL_HAVE_GETPID 1
#endif

/* Define to 1 if you have the `gettimeofday' function. */
#ifndef MPL_HAVE_GETTIMEOFDAY
#define MPL_HAVE_GETTIMEOFDAY 1
#endif

/* Define to 1 if you have the <helgrind.h> header file. */
/* #undef MPL_HAVE_HELGRIND_H */

/* Define to 1 if you have the <ifaddrs.h> header file. */
#ifndef MPL_HAVE_IFADDRS_H
#define MPL_HAVE_IFADDRS_H 1
#endif

/* Define to 1 if you have the `inet_ntop' function. */
#ifndef MPL_HAVE_INET_NTOP
#define MPL_HAVE_INET_NTOP 1
#endif

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef MPL_HAVE_INTTYPES_H
#define MPL_HAVE_INTTYPES_H 1
#endif

/* Define to 1 if you have the `abt' library (-labt). */
/* #undef MPL_HAVE_LIBABT */

/* Define to 1 if you have the `cuda' library (-lcuda). */
/* #undef MPL_HAVE_LIBCUDA */

/* Define to 1 if you have the `cudart' library (-lcudart). */
/* #undef MPL_HAVE_LIBCUDART */

/* Define to 1 if you have the `uti' library (-luti). */
/* #undef MPL_HAVE_LIBUTI */

/* Define to 1 if you have the `ze_loader' library (-lze_loader). */
/* #undef MPL_HAVE_LIBZE_LOADER */

/* Define to 1 if you have the `mach_absolute_time' function. */
/* #undef MPL_HAVE_MACH_ABSOLUTE_TIME */

/* Define if C99-style variable argument list macro functionality */
#ifndef MPL_HAVE_MACRO_VA_ARGS
#define MPL_HAVE_MACRO_VA_ARGS 1
#endif

/* Define to 1 if you have the <memcheck.h> header file. */
/* #undef MPL_HAVE_MEMCHECK_H */

/* Define to 1 if you have the <minix/config.h> header file. */
/* #undef MPL_HAVE_MINIX_CONFIG_H */

/* Define to 1 if you have the `mkstemp' function. */
#ifndef MPL_HAVE_MKSTEMP
#define MPL_HAVE_MKSTEMP 1
#endif

/* Define to 1 if you have the `mmap' function. */
#ifndef MPL_HAVE_MMAP
#define MPL_HAVE_MMAP 1
#endif

/* Define to 1 if you have the `munmap' function. */
#ifndef MPL_HAVE_MUNMAP
#define MPL_HAVE_MUNMAP 1
#endif

/* Define to 1 if we have support for Windows NT atomic intrinsics */
/* #undef MPL_HAVE_NT_INTRINSICS */

/* Define to 1 if you have the `posix_memalign' function. */
#ifndef MPL_HAVE_POSIX_MEMALIGN
#define MPL_HAVE_POSIX_MEMALIGN 1
#endif

/* Define to 1 if you have the <pthread.h> header file. */
#ifndef MPL_HAVE_PTHREAD_H
#define MPL_HAVE_PTHREAD_H 1
#endif

/* Define to 1 if you have the `pthread_mutexattr_setpshared' function. */
#ifndef MPL_HAVE_PTHREAD_MUTEXATTR_SETPSHARED
#define MPL_HAVE_PTHREAD_MUTEXATTR_SETPSHARED 1
#endif

/* Define to 1 if you have the `pthread_yield' function. */
/* #undef MPL_HAVE_PTHREAD_YIELD */

/* Define to 1 if you have the `putenv' function. */
#ifndef MPL_HAVE_PUTENV
#define MPL_HAVE_PUTENV 1
#endif

/* Define to 1 if you have the <sched.h> header file. */
#ifndef MPL_HAVE_SCHED_H
#define MPL_HAVE_SCHED_H 1
#endif

/* Define to 1 if you have the `sched_yield' function. */
#ifndef MPL_HAVE_SCHED_YIELD
#define MPL_HAVE_SCHED_YIELD 1
#endif

/* Define to 1 if you have the `select' function. */
#ifndef MPL_HAVE_SELECT
#define MPL_HAVE_SELECT 1
#endif

/* Define to 1 if you have the `shmat' function. */
/* #undef MPL_HAVE_SHMAT */

/* Define to 1 if you have the `shmctl' function. */
/* #undef MPL_HAVE_SHMCTL */

/* Define to 1 if you have the `shmdt' function. */
/* #undef MPL_HAVE_SHMDT */

/* Define to 1 if you have the `shmget' function. */
/* #undef MPL_HAVE_SHMGET */

/* Define to 1 if you have the `sleep' function. */
#ifndef MPL_HAVE_SLEEP
#define MPL_HAVE_SLEEP 1
#endif

/* Define to 1 if you have the `snprintf' function. */
#ifndef MPL_HAVE_SNPRINTF
#define MPL_HAVE_SNPRINTF 1
#endif

/* Define to 1 if you have the <stdarg.h> header file. */
#ifndef MPL_HAVE_STDARG_H
#define MPL_HAVE_STDARG_H 1
#endif

/* Define to 1 if stdbool.h conforms to C99. */
/* #undef MPL_HAVE_STDBOOL_H */

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef MPL_HAVE_STDINT_H
#define MPL_HAVE_STDINT_H 1
#endif

/* Define to 1 if you have the <stdio.h> header file. */
#ifndef MPL_HAVE_STDIO_H
#define MPL_HAVE_STDIO_H 1
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef MPL_HAVE_STDLIB_H
#define MPL_HAVE_STDLIB_H 1
#endif

/* Define to 1 if you have the `strdup' function. */
#ifndef MPL_HAVE_STRDUP
#define MPL_HAVE_STRDUP 1
#endif

/* Define to 1 if you have the `strerror' function. */
#ifndef MPL_HAVE_STRERROR
#define MPL_HAVE_STRERROR 1
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef MPL_HAVE_STRINGS_H
#define MPL_HAVE_STRINGS_H 1
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef MPL_HAVE_STRING_H
#define MPL_HAVE_STRING_H 1
#endif

/* Define to 1 if you have the <sys/mman.h> header file. */
#ifndef MPL_HAVE_SYS_MMAN_H
#define MPL_HAVE_SYS_MMAN_H 1
#endif

/* Define to 1 if you have the <sys/param.h> header file. */
#ifndef MPL_HAVE_SYS_PARAM_H
#define MPL_HAVE_SYS_PARAM_H 1
#endif

/* Define to 1 if you have the <sys/select.h> header file. */
#ifndef MPL_HAVE_SYS_SELECT_H
#define MPL_HAVE_SYS_SELECT_H 1
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef MPL_HAVE_SYS_STAT_H
#define MPL_HAVE_SYS_STAT_H 1
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef MPL_HAVE_SYS_TYPES_H
#define MPL_HAVE_SYS_TYPES_H 1
#endif

/* Define to 1 if you have the <sys/uio.h> header file. */
#ifndef MPL_HAVE_SYS_UIO_H
#define MPL_HAVE_SYS_UIO_H 1
#endif

/* Define to 1 if you have the <thread.h> header file. */
/* #undef MPL_HAVE_THREAD_H */

/* Define to 1 if you have the `thr_yield' function. */
/* #undef MPL_HAVE_THR_YIELD */

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef MPL_HAVE_UNISTD_H
#define MPL_HAVE_UNISTD_H 1
#endif

/* Define to 1 if you have the `usleep' function. */
#ifndef MPL_HAVE_USLEEP
#define MPL_HAVE_USLEEP 1
#endif

/* Define to 1 if you have the <valgrind/drd.h> header file. */
/* #undef MPL_HAVE_VALGRIND_DRD_H */

/* Define to 1 if you have the <valgrind.h> header file. */
/* #undef MPL_HAVE_VALGRIND_H */

/* Define to 1 if you have the <valgrind/helgrind.h> header file. */
/* #undef MPL_HAVE_VALGRIND_HELGRIND_H */

/* Define to 1 if you have the <valgrind/memcheck.h> header file. */
/* #undef MPL_HAVE_VALGRIND_MEMCHECK_H */

/* Define to 1 if you have the <valgrind/valgrind.h> header file. */
/* #undef MPL_HAVE_VALGRIND_VALGRIND_H */

/* Define to 1 if the system has the `aligned' variable attribute */
#ifndef MPL_HAVE_VAR_ATTRIBUTE_ALIGNED
#define MPL_HAVE_VAR_ATTRIBUTE_ALIGNED 1
#endif

/* Define to 1 if the system has the `used' variable attribute */
#ifndef MPL_HAVE_VAR_ATTRIBUTE_USED
#define MPL_HAVE_VAR_ATTRIBUTE_USED 1
#endif

/* Define to 1 if you have the <wchar.h> header file. */
#ifndef MPL_HAVE_WCHAR_H
#define MPL_HAVE_WCHAR_H 1
#endif

/* Define to 1 if you have the <windows.h> header file. */
/* #undef MPL_HAVE_WINDOWS_H */

/* Define to 1 if you have the `yield' function. */
/* #undef MPL_HAVE_YIELD */

/* Define if ZE is available */
/* #undef MPL_HAVE_ZE */

/* Define to 1 if the system has the type `_Bool'. */
#ifndef MPL_HAVE__BOOL
#define MPL_HAVE__BOOL 1
#endif

/* defined if the C compiler supports __typeof(variable) */
#ifndef MPL_HAVE___TYPEOF
#define MPL_HAVE___TYPEOF 1
#endif

/* Define which x86 cycle counter to use */
/* #undef MPL_LINUX86_CYCLE_CPUID_RDTSC32 */

/* Define which x86 cycle counter to use */
/* #undef MPL_LINUX86_CYCLE_CPUID_RDTSC64 */

/* Define which x86 cycle counter to use */
/* #undef MPL_LINUX86_CYCLE_RDTSC */

/* Define which x86 cycle counter to use */
/* #undef MPL_LINUX86_CYCLE_RDTSCP */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#ifndef MPL_LT_OBJDIR
#define MPL_LT_OBJDIR ".libs/"
#endif

/* Define if aligned_alloc needs a declaration */
/* #undef MPL_NEEDS_ALIGNED_ALLOC_DECL */

/* Define if fdopen needs a declaration */
/* #undef MPL_NEEDS_FDOPEN_DECL */

/* Define if mkstemp needs a declaration */
/* #undef MPL_NEEDS_MKSTEMP_DECL */

/* Define if pthread_mutexattr_settype needs a declaration */
/* #undef MPL_NEEDS_PTHREAD_MUTEXATTR_SETTYPE_DECL */

/* Define if putenv needs a declaration */
/* #undef MPL_NEEDS_PUTENV_DECL */

/* Define if snprintf needs a declaration */
/* #undef MPL_NEEDS_SNPRINTF_DECL */

/* Define if strdup needs a declaration */
/* #undef MPL_NEEDS_STRDUP_DECL */

/* Define if strerror needs a declaration */
/* #undef MPL_NEEDS_STRERROR_DECL */

/* Define if sys/time.h is required to get timer definitions */
/* #undef MPL_NEEDS_SYS_TIME_H */

/* Define if usleep needs a declaration */
/* #undef MPL_NEEDS_USLEEP_DECL */

/* Name of package */
#ifndef MPL_PACKAGE
#define MPL_PACKAGE "mpl"
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef MPL_PACKAGE_BUGREPORT
#define MPL_PACKAGE_BUGREPORT ""
#endif

/* Define to the full name of this package. */
#ifndef MPL_PACKAGE_NAME
#define MPL_PACKAGE_NAME "MPL"
#endif

/* Define to the full name and version of this package. */
#ifndef MPL_PACKAGE_STRING
#define MPL_PACKAGE_STRING "MPL 0.1"
#endif

/* Define to the one symbol short name of this package. */
#ifndef MPL_PACKAGE_TARNAME
#define MPL_PACKAGE_TARNAME "mpl"
#endif

/* Define to the home page for this package. */
#ifndef MPL_PACKAGE_URL
#define MPL_PACKAGE_URL ""
#endif

/* Define to the version of this package. */
#ifndef MPL_PACKAGE_VERSION
#define MPL_PACKAGE_VERSION "0.1"
#endif

/* set to the choice of the posix mutex */
#ifndef MPL_POSIX_MUTEX_NAME
#define MPL_POSIX_MUTEX_NAME MPL_POSIX_MUTEX_NATIVE
#endif

/* set to the name of the interprocess mutex package */
#ifndef MPL_PROC_MUTEX_PACKAGE_NAME
#define MPL_PROC_MUTEX_PACKAGE_NAME MPL_PROC_MUTEX_PACKAGE_POSIX
#endif

/* Define to an expression that will result in an error checking mutex type.
   */
/* #undef MPL_PTHREAD_MUTEX_ERRORCHECK_VALUE */

/* The size of `void *', as computed by sizeof. */
#ifndef MPL_SIZEOF_VOID_P
#define MPL_SIZEOF_VOID_P 8
#endif

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#ifndef MPL_STDC_HEADERS
#define MPL_STDC_HEADERS 1
#endif

/* set to the name of the thread package */
#ifndef MPL_THREAD_PACKAGE_NAME
#define MPL_THREAD_PACKAGE_NAME MPL_THREAD_PACKAGE_POSIX
#endif

/* Define if performing coverage tests */
/* #undef MPL_USE_COVERAGE */

/* Define to enable logging macros */
/* #undef MPL_USE_DBG_LOGGING */

/* Define to 1 if mutex-based synchronization is used */
/* #undef MPL_USE_LOCK_BASED_PRIMITIVES */

/* Define to enable memory tracing */
/* #undef MPL_USE_MEMORY_TRACING */

/* Define if we have sysv shared memory */
#ifndef MPL_USE_MMAP_SHM
#define MPL_USE_MMAP_SHM 1
#endif

/* Define to use nothing to yield processor */
#ifndef MPL_USE_NOTHING_FOR_YIELD
#define MPL_USE_NOTHING_FOR_YIELD 1
#endif

/* Define to 1 if no atomic primitives are used */
/* #undef MPL_USE_NO_ATOMIC_PRIMITIVES */

/* Define if use Windows shared memory */
/* #undef MPL_USE_NT_SHM */

/* Define to use sched_yield to yield processor */
/* #undef MPL_USE_SCHED_YIELD_FOR_YIELD */

/* Define to use select to yield processor */
/* #undef MPL_USE_SELECT_FOR_YIELD */

/* Define to use sleep to yield processor */
/* #undef MPL_USE_SLEEP_FOR_YIELD */

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


/* Define if we have sysv shared memory */
/* #undef MPL_USE_SYSV_SHM */

/* Define to use usleep to yield processor */
/* #undef MPL_USE_USLEEP_FOR_YIELD */

/* Define to use yield to yield processor */
/* #undef MPL_USE_YIELD_FOR_YIELD */

/* Version number of package */
#ifndef MPL_VERSION
#define MPL_VERSION "0.1"
#endif

/* Defined to return type of backtrace(). */
#ifndef _mpl_backtrace_size_t
#define _mpl_backtrace_size_t int
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported only directly.  */
#ifndef restrict
#define restrict __restrict__
#endif
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
/* once: INCLUDE_MPLCONFIG_H */
#endif
