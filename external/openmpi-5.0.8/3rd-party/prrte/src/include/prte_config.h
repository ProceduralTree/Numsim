/* src/include/prte_config.h.  Generated from prte_config.h.in by configure.  */
/* src/include/prte_config.h.in.  Generated from configure.ac by autoheader.  */

/* -*- c -*-
 *
 * Copyright (c) 2004-2005 The Trustees of Indiana University.
 *                         All rights reserved.
 * Copyright (c) 2004-2005 The Trustees of the University of Tennessee.
 *                         All rights reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2014      Intel, Inc. All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 *
 * Function: - OS, CPU and compiler dependent configuration
 */

#ifndef PRTE_CONFIG_H
#define PRTE_CONFIG_H

#include "prte_config_top.h"


/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* The normal alignment of `bool', in bytes. */
#define ALIGNOF_BOOL 1

/* The normal alignment of `double', in bytes. */
#define ALIGNOF_DOUBLE 8

/* The normal alignment of `int', in bytes. */
#define ALIGNOF_INT 4

/* The normal alignment of `int32_t', in bytes. */
#define ALIGNOF_INT32_T 4

/* The normal alignment of `int64_t', in bytes. */
#define ALIGNOF_INT64_T 8

/* The normal alignment of `long', in bytes. */
#define ALIGNOF_LONG 8

/* The normal alignment of `long long', in bytes. */
#define ALIGNOF_LONG_LONG 8

/* The normal alignment of `size_t', in bytes. */
#define ALIGNOF_SIZE_T 8

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define to 1 if you have the `asprintf' function. */
#define HAVE_ASPRINTF 1

/* Define to 1 if you have the declaration of `AF_INET6', and to 0 if you
   don't. */
#define HAVE_DECL_AF_INET6 1

/* Define to 1 if you have the declaration of `AF_UNSPEC', and to 0 if you
   don't. */
#define HAVE_DECL_AF_UNSPEC 1

/* Define to 1 if you have the declaration of `ethtool_cmd_speed', and to 0 if
   you don't. */
#define HAVE_DECL_ETHTOOL_CMD_SPEED 1

/* Define to 1 if you have the declaration of `PF_INET6', and to 0 if you
   don't. */
#define HAVE_DECL_PF_INET6 1

/* Define to 1 if you have the declaration of `PF_UNSPEC', and to 0 if you
   don't. */
#define HAVE_DECL_PF_UNSPEC 1

/* Define to 1 if you have the declaration of `RLIMIT_AS', and to 0 if you
   don't. */
#define HAVE_DECL_RLIMIT_AS 1

/* Define to 1 if you have the declaration of `RLIMIT_CORE', and to 0 if you
   don't. */
#define HAVE_DECL_RLIMIT_CORE 1

/* Define to 1 if you have the declaration of `RLIMIT_FSIZE', and to 0 if you
   don't. */
#define HAVE_DECL_RLIMIT_FSIZE 1

/* Define to 1 if you have the declaration of `RLIMIT_MEMLOCK', and to 0 if
   you don't. */
#define HAVE_DECL_RLIMIT_MEMLOCK 1

/* Define to 1 if you have the declaration of `RLIMIT_NOFILE', and to 0 if you
   don't. */
#define HAVE_DECL_RLIMIT_NOFILE 1

/* Define to 1 if you have the declaration of `RLIMIT_NPROC', and to 0 if you
   don't. */
#define HAVE_DECL_RLIMIT_NPROC 1

/* Define to 1 if you have the declaration of `RLIMIT_STACK', and to 0 if you
   don't. */
#define HAVE_DECL_RLIMIT_STACK 1

/* Define to 1 if you have the declaration of `SIOCETHTOOL', and to 0 if you
   don't. */
#define HAVE_DECL_SIOCETHTOOL 1

/* Define to 1 if you have the declaration of `__func__', and to 0 if you
   don't. */
#define HAVE_DECL___FUNC__ 1

/* Define to 1 if you have the <dirent.h> header file. */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <endian.h> header file. */
#define HAVE_ENDIAN_H 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `fileno_unlocked' function. */
#define HAVE_FILENO_UNLOCKED 1

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Define to 1 if you have the <hostLib.h> header file. */
/* #undef HAVE_HOSTLIB_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <ioLib.h> header file. */
/* #undef HAVE_IOLIB_H */

/* Define to 1 if you have the <libutil.h> header file. */
/* #undef HAVE_LIBUTIL_H */

/* Define to 1 if you have the <linux/ethtool.h> header file. */
#define HAVE_LINUX_ETHTOOL_H 1

/* Define to 1 if you have the <linux/sockios.h> header file. */
#define HAVE_LINUX_SOCKIOS_H 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <minix/config.h> header file. */
/* #undef HAVE_MINIX_CONFIG_H */

/* Define to 1 if you have the <netdb.h> header file. */
#define HAVE_NETDB_H 1

/* Define to 1 if you have the <netinet/in.h> header file. */
#define HAVE_NETINET_IN_H 1

/* Define to 1 if you have the <netinet/tcp.h> header file. */
#define HAVE_NETINET_TCP_H 1

/* Define to 1 if you have the <net/if.h> header file. */
#define HAVE_NET_IF_H 1

/* Define to 1 if you have the <net/uio.h> header file. */
/* #undef HAVE_NET_UIO_H */

/* Define to 1 if you have the <paths.h> header file. */
#define HAVE_PATHS_H 1

/* Define to 1 if you have the <poll.h> header file. */
#define HAVE_POLL_H 1

/* Define to 1 if you have the `printstack' function. */
/* #undef HAVE_PRINTSTACK */

/* Define to 1 if you have the `pthread_condattr_setpshared' function. */
#define HAVE_PTHREAD_CONDATTR_SETPSHARED 1

/* Define to 1 if you have the `pthread_mutexattr_setpshared' function. */
#define HAVE_PTHREAD_MUTEXATTR_SETPSHARED 1

/* Define to 1 if you have the `pthread_setaffinity_np' function. */
#define HAVE_PTHREAD_SETAFFINITY_NP 1

/* Define to 1 if you have the <pty.h> header file. */
#define HAVE_PTY_H 1

/* Define to 1 if you have the <pwd.h> header file. */
#define HAVE_PWD_H 1

/* Define to 1 if you have the <sched.h> header file. */
#define HAVE_SCHED_H 1

/* Define to 1 if you have the `setpgid' function. */
#define HAVE_SETPGID 1

/* Define to 1 if you have the `setsid' function. */
#define HAVE_SETSID 1

/* Define to 1 if `si_band' is a member of `siginfo_t'. */
#define HAVE_SIGINFO_T_SI_BAND 1

/* Define to 1 if `si_fd' is a member of `siginfo_t'. */
#define HAVE_SIGINFO_T_SI_FD 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if the system has the type `socklen_t'. */
#define HAVE_SOCKLEN_T 1

/* Define to 1 if you have the <sockLib.h> header file. */
/* #undef HAVE_SOCKLIB_H */

/* Define to 1 if you have the <stdatomic.h> header file. */
#define HAVE_STDATOMIC_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strsignal' function. */
#define HAVE_STRSIGNAL 1

/* Define to 1 if the system has the type `struct ethtool_cmd'. */
#define HAVE_STRUCT_ETHTOOL_CMD 1

/* Define to 1 if `speed_hi' is a member of `struct ethtool_cmd'. */
#define HAVE_STRUCT_ETHTOOL_CMD_SPEED_HI 1

/* Define to 1 if the system has the type `struct ifreq'. */
#define HAVE_STRUCT_IFREQ 1

/* Define to 1 if the system has the type `struct sockaddr_in'. */
#define HAVE_STRUCT_SOCKADDR_IN 1

/* Define to 1 if the system has the type `struct sockaddr_in6'. */
#define HAVE_STRUCT_SOCKADDR_IN6 1

/* Define to 1 if the system has the type `struct sockaddr_storage'. */
#define HAVE_STRUCT_SOCKADDR_STORAGE 1

/* Define to 1 if you have the `syslog' function. */
#define HAVE_SYSLOG 1

/* Define to 1 if you have the <syslog.h> header file. */
#define HAVE_SYSLOG_H 1

/* Define to 1 if you have the <sys/fcntl.h> header file. */
#define HAVE_SYS_FCNTL_H 1

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#define HAVE_SYS_IOCTL_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Whether or not we have the ptrace header */
#define HAVE_SYS_PTRACE_H 1

/* Define to 1 if you have the <sys/queue.h> header file. */
#define HAVE_SYS_QUEUE_H 1

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/uio.h> header file. */
#define HAVE_SYS_UIO_H 1

/* Define to 1 if you have the <sys/un.h> header file. */
#define HAVE_SYS_UN_H 1

/* Define to 1 if you have the <sys/utsname.h> header file. */
#define HAVE_SYS_UTSNAME_H 1

/* Define to 1 if you have the <sys/wait.h> header file. */
#define HAVE_SYS_WAIT_H 1

/* Define to 1 if you have the <termios.h> header file. */
#define HAVE_TERMIOS_H 1

/* Define to 1 if you have the <ucontext.h> header file. */
#define HAVE_UCONTEXT_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* whether unix byteswap routines -- htonl, htons, nothl, ntohs -- are
   available */
#define HAVE_UNIX_BYTESWAP 1

/* Define to 1 if you have the <util.h> header file. */
/* #undef HAVE_UTIL_H */

/* Define to 1 if you have the `vasprintf' function. */
#define HAVE_VASPRINTF 1

/* Define to 1 if you have the `vsnprintf' function. */
#define HAVE_VSNPRINTF 1

/* Define to 1 if you have the <wchar.h> header file. */
#define HAVE_WCHAR_H 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Whether or not we have apple */
#define OAC_HAVE_APPLE 0

/* Whether or not we have solaris */
#define OAC_HAVE_SOLARIS 0

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://github.com/openpmix/prrte/"

/* Define to the full name of this package. */
#define PACKAGE_NAME "prte"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "prte gitclone"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "prrte"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "gitclone"

/* PRTE architecture string */
#define PRTE_ARCH "x86_64-pc-linux-gnu"

/* Use C11 style atomics */
#define PRTE_ATOMIC_C11 1

/* Use GCC builtin style atomics */
#define PRTE_ATOMIC_GCC_BUILTIN 0

/* whether building on x86_64 platform */
#define PRTE_ATOMIC_X86_64 1

/* The compiler $lower which OMPI was built with */
#define PRTE_BUILD_PLATFORM_COMPILER_FAMILYID 0

/* The compiler $lower which OMPI was built with */
#define PRTE_BUILD_PLATFORM_COMPILER_VERSION 0

/* PRRTE underlying C compiler */
#define PRTE_CC "gcc"

/* Use static const char[] strings for C files */
#define PRTE_CC_USE_CONST_CHAR_IDENT 0

/* Use #ident strings for C files */
#define PRTE_CC_USE_IDENT 1

/* Use #pragma comment for C files */
#define PRTE_CC_USE_PRAGMA_COMMENT 

/* Use #pragma ident strings for C files */
#define PRTE_CC_USE_PRAGMA_IDENT 0

/* Capture the configure cmd line */
#define PRTE_CONFIGURE_CLI " \'--disable-option-checking\' \'--prefix=NONE\' \'--enable-prte-ft\' \'--with-proxy-version-string=5.0.8\' \'--with-proxy-package-name=Open MPI\' \'--with-proxy-bugreport=https://www.open-mpi.org/community/help/\' \'--disable-devel-check\' \'--enable-prte-prefix-by-default\' \'--with-libevent\' \'--disable-libevent-lib-checks\' \'--with-libevent-extra-libs=/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/libevent-2.1.12-stable-ompi/libevent_core.la /home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/libevent-2.1.12-stable-ompi/libevent_pthreads.la\' \'--disable-hwloc-lib-checks\' \'--with-hwloc-extra-libs=/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/hwloc-2.7.1/hwloc/libhwloc.la\' \'--disable-pmix-lib-checks\' \'--with-pmix-extra-libs=/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/openpmix/src/libpmix.la\' \'CPPFLAGS= -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/libevent-2.1.12-stable-ompi -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/libevent-2.1.12-stable-ompi/include -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/hwloc-2.7.1/include -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/hwloc-2.7.1/include -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/openpmix/include -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/openpmix/include -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/openpmix/ -I/home/proceduraltree/Projects/Numsim/external/openmpi-5.0.8/3rd-party/openpmix/\' \'PKG_CONFIG_PATH=/usr/lib64/openmpi/lib/pkgconfig\' \'--cache-file=/dev/null\' \'--srcdir=.\'"

/* Date when PMIx was built */
#define PRTE_CONFIGURE_DATE "Wed Oct 22 08:50:07 UTC 2025"

/* Hostname where PMIx was built */
#define PRTE_CONFIGURE_HOST "fedora"

/* User who built PMIx */
#define PRTE_CONFIGURE_USER "proceduraltree"

/* Whether C compiler supports atomic convenience variables in stdatomic.h */
#define PRTE_C_HAVE_ATOMIC_CONV_VAR 1

/* Whether C compiler supports __builtin_clz */
#define PRTE_C_HAVE_BUILTIN_CLZ 0

/* Whether C compiler supports __builtin_expect */
#define PRTE_C_HAVE_BUILTIN_EXPECT 0

/* Whether C compiler supports __builtin_prefetch */
#define PRTE_C_HAVE_BUILTIN_PREFETCH 0

/* Whether C compiler supports symbol visibility or not */
#define PRTE_C_HAVE_VISIBILITY 1

/* Whether C compiler supports __Atomic keyword */
#define PRTE_C_HAVE__ATOMIC 1

/* Whether C compiler supports __Generic keyword */
#define PRTE_C_HAVE__GENERIC 1

/* Whether C compiler supports _Static_assert keyword */
#define PRTE_C_HAVE__STATIC_ASSERT 1

/* Whether C compiler supports __Thread_local */
#define PRTE_C_HAVE__THREAD_LOCAL 1

/* Whether C compiler supports __thread */
#define PRTE_C_HAVE___THREAD 1

/* Command to detach from process being traced */
#define PRTE_DETACH PT_DETACH

/* Whether we have lt_dladvise or not */
#define PRTE_DL_LIBLTDL_HAVE_LT_DLADVISE 0

/* Whether we want developer-level debugging code or not */
#define PRTE_ENABLE_DEBUG 0

/* Whether we want to enable dlopen support */
#define PRTE_ENABLE_DLOPEN_SUPPORT 1

/* Disable getpwuid support (default: enabled) */
#define PRTE_ENABLE_GETPWUID 1

/* Enable IPv6 support, but only if the underlying system supports it */
#define PRTE_ENABLE_IPV6 0

/* Whether we should enable thread support within the PRTE code base */
#define PRTE_ENABLE_MULTI_THREADS 1

/* Whether user wants PTY support or not */
#define PRTE_ENABLE_PTY_SUPPORT 1

/* If built from a git repo */
/* #undef PRTE_GIT_REPO_BUILD */

/* Greek - alpha, beta, etc - release number of PMIx Reference Run-Time
   Environment */
#define PRTE_GREEK_VERSION "a1"

/* Whether your compiler has __attribute__ or not */
#define PRTE_HAVE_ATTRIBUTE 1

/* Whether your compiler has __attribute__ aligned or not */
#define PRTE_HAVE_ATTRIBUTE_ALIGNED 1

/* Whether your compiler has __attribute__ always_inline or not */
#define PRTE_HAVE_ATTRIBUTE_ALWAYS_INLINE 1

/* Whether your compiler has __attribute__ cold or not */
#define PRTE_HAVE_ATTRIBUTE_COLD 1

/* Whether your compiler has __attribute__ const or not */
#define PRTE_HAVE_ATTRIBUTE_CONST 1

/* Whether your compiler has __attribute__ deprecated or not */
#define PRTE_HAVE_ATTRIBUTE_DEPRECATED 1

/* Whether your compiler has __attribute__ deprecated with optional argument
   */
#define PRTE_HAVE_ATTRIBUTE_DEPRECATED_ARGUMENT 1

/* Whether your compiler has __attribute__ destructor or not */
#define PRTE_HAVE_ATTRIBUTE_DESTRUCTOR 1

/* Whether your compiler has __attribute__ extension or not */
#define PRTE_HAVE_ATTRIBUTE_EXTENSION 1

/* Whether your compiler has __attribute__ format or not */
#define PRTE_HAVE_ATTRIBUTE_FORMAT 1

/* Whether your compiler has __attribute__ format and it works on function
   pointers */
#define PRTE_HAVE_ATTRIBUTE_FORMAT_FUNCPTR 1

/* Whether your compiler has __attribute__ hot or not */
#define PRTE_HAVE_ATTRIBUTE_HOT 1

/* Whether your compiler has __attribute__ malloc or not */
#define PRTE_HAVE_ATTRIBUTE_MALLOC 1

/* Whether your compiler has __attribute__ may_alias or not */
#define PRTE_HAVE_ATTRIBUTE_MAY_ALIAS 1

/* Whether your compiler has __attribute__ noinline or not */
#define PRTE_HAVE_ATTRIBUTE_NOINLINE 1

/* Whether your compiler has __attribute__ nonnull or not */
#define PRTE_HAVE_ATTRIBUTE_NONNULL 1

/* Whether your compiler has __attribute__ noreturn or not */
#define PRTE_HAVE_ATTRIBUTE_NORETURN 1

/* Whether your compiler has __attribute__ noreturn and it works on function
   pointers */
#define PRTE_HAVE_ATTRIBUTE_NORETURN_FUNCPTR 1

/* Whether your compiler has __attribute__ no_instrument_function or not */
#define PRTE_HAVE_ATTRIBUTE_NO_INSTRUMENT_FUNCTION 1

/* Whether your compiler has __attribute__ optnone or not */
#define PRTE_HAVE_ATTRIBUTE_OPTNONE 0

/* Whether your compiler has __attribute__ packed or not */
#define PRTE_HAVE_ATTRIBUTE_PACKED 1

/* Whether your compiler has __attribute__ pure or not */
#define PRTE_HAVE_ATTRIBUTE_PURE 1

/* Whether your compiler has __attribute__ sentinel or not */
#define PRTE_HAVE_ATTRIBUTE_SENTINEL 1

/* Whether your compiler has __attribute__ unused or not */
#define PRTE_HAVE_ATTRIBUTE_UNUSED 1

/* Whether your compiler has __attribute__ visibility or not */
#define PRTE_HAVE_ATTRIBUTE_VISIBILITY 1

/* Whether your compiler has __attribute__ warn unused result or not */
#define PRTE_HAVE_ATTRIBUTE_WARN_UNUSED_RESULT 1

/* Whether your compiler has __attribute__ weak alias or not */
#define PRTE_HAVE_ATTRIBUTE_WEAK_ALIAS 1

/* whether backtrace_execinfo is found and available */
#define PRTE_HAVE_BACKTRACE_EXECINFO 1

/* whether ceil is found and available */
#define PRTE_HAVE_CEIL 1

/* Whether we have Clang __c11 atomic functions */
/* #undef PRTE_HAVE_CLANG_BUILTIN_ATOMIC_C11_FUNC */

/* Whether the PRTE DL framework is functional or not */
#define PRTE_HAVE_DL_SUPPORT 1

/* whether gethostbyname is found and available */
#define PRTE_HAVE_GETHOSTBYNAME 1

/* Whether we are building against libev */
#define PRTE_HAVE_LIBEV 0

/* Whether we are building against libevent */
#define PRTE_HAVE_LIBEVENT 1

/* Does ptrace have the Linux signature */
#define PRTE_HAVE_LINUX_PTRACE 0

/* If PTHREADS implementation supports PTHREAD_MUTEX_ERRORCHECK */
#define PRTE_HAVE_PTHREAD_MUTEX_ERRORCHECK 1

/* If PTHREADS implementation supports PTHREAD_MUTEX_ERRORCHECK_NP */
#define PRTE_HAVE_PTHREAD_MUTEX_ERRORCHECK_NP 1

/* whether shm_open_rt is found and available */
/* #undef PRTE_HAVE_SHM_OPEN_RT */

/* whether socket is found and available */
#define PRTE_HAVE_SOCKET 1

/* Whether or not we have stop-on-exec support */
#define PRTE_HAVE_STOP_ON_EXEC 1

/* whether yp_all_nsl is found and available */
#define PRTE_HAVE_YP_ALL_NSL 0

/* ident string for Open MPI */
#define PRTE_IDENT_STRING "3.0.11a1"

/* Major release number of PMIx Reference Run-Time Environment */
#define PRTE_MAJOR_VERSION 3

/* Whether or not we are using memory sanitizers */
#define PRTE_MEMORY_SANITIZERS 0

/* Minor release number of PMIx Reference Run-Time Environment */
#define PRTE_MINOR_VERSION 0

/* package/branding string for Open MPI */
#define PRTE_PACKAGE_STRING "Open MPI proceduraltree@fedora Distribution"

/* Whether or not we are using picky compiler settings */
#define PRTE_PICKY_COMPILERS 0

/* Minimum supported PMIx version */
#define PRTE_PMIX_MINIMUM_VERSION 0x00040204

/* Bugreport string to be returned by prte when in proxy mode */
#define PRTE_PROXY_BUGREPORT "https://www.open-mpi.org/community/help/"

/* Package name to be returned by prte when in proxy mode */
#define PRTE_PROXY_PACKAGE_NAME "Open MPI"

/* Version string to be returned by prte when in proxy mode */
#define PRTE_PROXY_VERSION_STRING "5.0.8"

/* Release date of PMIx Reference Run-Time Environment */
#define PRTE_RELEASE_DATE "Nov 24, 2018"

/* Release release number of PMIx Reference Run-Time Environment */
#define PRTE_RELEASE_VERSION 11

/* The repository version PMIx Reference Run-Time Environment */
#define PRTE_REPO_REV "2025-10-22"

/* Default value for mca_base_component_show_load_errors MCA variable */
#define PRTE_SHOW_LOAD_ERRORS_DEFAULT 0

/* Tarball filename version string of PMIx Reference Run-Time Environment */
#define PRTE_TARBALL_VERSION "gitclone"

/* Enable testbuild PLM launchers (default: disabled) */
#define PRTE_TESTBUILD_LAUNCHERS 0

/* Command for declaring that process expects to be traced by parent */
#define PRTE_TRACEME PT_TRACE_ME

/* Complete release number of PMIx Reference Run-Time Environment */
#define PRTE_VERSION "11"

/* Enable per-user config files */
#define PRTE_WANT_HOME_CONFIG_FILES 1

/* if want pretty-print stack trace feature */
#define PRTE_WANT_PRETTY_PRINT_STACKTRACE 1

/* Whether we want prte to effect "--prefix $prefix" by default */
#define PRTE_WANT_PRTE_PREFIX_BY_DEFAULT 1

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

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


/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#define YYTEXT_POINTER 1

/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif


#include "prte_config_bottom.h"
#endif /* PRTE_CONFIG_H */

