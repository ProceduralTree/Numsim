/* opal/include/opal_config.h.  Generated from opal_config.h.in by configure.  */
/* opal/include/opal_config.h.in.  Generated from configure.ac by autoheader.  */

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

#ifndef OPAL_CONFIG_H
#define OPAL_CONFIG_H

#include "opal_config_top.h"



/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* The normal alignment of `bool', in bytes. */
#define ALIGNOF_BOOL 1

/* The normal alignment of `char', in bytes. */
#define ALIGNOF_CHAR 1

/* The normal alignment of `double', in bytes. */
#define ALIGNOF_DOUBLE 8

/* The normal alignment of `double _Complex', in bytes. */
#define ALIGNOF_DOUBLE__COMPLEX 8

/* The normal alignment of `float', in bytes. */
#define ALIGNOF_FLOAT 4

/* The normal alignment of `float _Complex', in bytes. */
#define ALIGNOF_FLOAT__COMPLEX 4

/* The normal alignment of `int', in bytes. */
#define ALIGNOF_INT 4

/* The normal alignment of `int128_t', in bytes. */
/* #undef ALIGNOF_INT128_T */

/* The normal alignment of `int16_t', in bytes. */
#define ALIGNOF_INT16_T 2

/* The normal alignment of `int32_t', in bytes. */
#define ALIGNOF_INT32_T 4

/* The normal alignment of `int64_t', in bytes. */
#define ALIGNOF_INT64_T 8

/* The normal alignment of `int8_t', in bytes. */
#define ALIGNOF_INT8_T 1

/* The normal alignment of `long', in bytes. */
#define ALIGNOF_LONG 8

/* The normal alignment of `long double', in bytes. */
#define ALIGNOF_LONG_DOUBLE 16

/* The normal alignment of `long double _Complex', in bytes. */
#define ALIGNOF_LONG_DOUBLE__COMPLEX 16

/* The normal alignment of `long long', in bytes. */
#define ALIGNOF_LONG_LONG 8

/* The normal alignment of `opal_short_float_t', in bytes. */
#define ALIGNOF_OPAL_SHORT_FLOAT_T 2

/* The normal alignment of `short', in bytes. */
#define ALIGNOF_SHORT 2

/* The normal alignment of `short float', in bytes. */
/* #undef ALIGNOF_SHORT_FLOAT */

/* The normal alignment of `short float _Complex', in bytes. */
/* #undef ALIGNOF_SHORT_FLOAT__COMPLEX */

/* The normal alignment of `size_t', in bytes. */
#define ALIGNOF_SIZE_T 8

/* The normal alignment of `void *', in bytes. */
#define ALIGNOF_VOID_P 8

/* The normal alignment of `wchar_t', in bytes. */
#define ALIGNOF_WCHAR_T 4

/* The normal alignment of `__float128', in bytes. */
#define ALIGNOF___FLOAT128 16

/* Define to 1 if you have the `ABT_unit_get_thread' function. */
/* #undef HAVE_ABT_UNIT_GET_THREAD */

/* Define to 1 if you have the <aio.h> header file. */
#define HAVE_AIO_H 1

/* Define to 1 if you have the <alloca.h> header file. */
#define HAVE_ALLOCA_H 1

/* Define to 1 if you have the <arpa/inet.h> header file. */
#define HAVE_ARPA_INET_H 1

/* Define to 1 if you have the `asprintf' function. */
#define HAVE_ASPRINTF 1

/* Define to 1 if you have the <complex.h> header file. */
#define HAVE_COMPLEX_H 1

/* Define to 1 if you have the <crt_externs.h> header file. */
/* #undef HAVE_CRT_EXTERNS_H */

/* Define to 1 if you have the `dbm_open' function. */
/* #undef HAVE_DBM_OPEN */

/* Define to 1 if you have the `dbopen' function. */
/* #undef HAVE_DBOPEN */

/* Define to 1 if you have the <db.h> header file. */
/* #undef HAVE_DB_H */

/* Define to 1 if you have the declaration of `AF_INET6', and to 0 if you
   don't. */
#define HAVE_DECL_AF_INET6 1

/* Define to 1 if you have the declaration of `AF_UNSPEC', and to 0 if you
   don't. */
#define HAVE_DECL_AF_UNSPEC 1

/* Define to 1 if you have the declaration of `ethtool_cmd_speed', and to 0 if
   you don't. */
#define HAVE_DECL_ETHTOOL_CMD_SPEED 1

/* Define to 1 if you have the declaration of `FI_OPT_FI_HMEM_P2P', and to 0
   if you don't. */
/* #undef HAVE_DECL_FI_OPT_FI_HMEM_P2P */

/* Define to 1 if you have the declaration of `HWLOC_OBJ_OSDEV_COPROC', and to
   0 if you don't. */
/* #undef HAVE_DECL_HWLOC_OBJ_OSDEV_COPROC */

/* Define to 1 if you have the declaration of `open_memstream', and to 0 if
   you don't. */
#define HAVE_DECL_OPEN_MEMSTREAM 1

/* Define to 1 if you have the declaration of `PF_INET6', and to 0 if you
   don't. */
#define HAVE_DECL_PF_INET6 1

/* Define to 1 if you have the declaration of `PF_UNSPEC', and to 0 if you
   don't. */
#define HAVE_DECL_PF_UNSPEC 1

/* Define to 1 if you have the declaration of `PMIX_PACKAGE_RANK', and to 0 if
   you don't. */
/* #undef HAVE_DECL_PMIX_PACKAGE_RANK */

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

/* Define to 1 if you have the declaration of `ucm_test_events', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCM_TEST_EVENTS */

/* Define to 1 if you have the declaration of `ucm_test_external_events', and
   to 0 if you don't. */
/* #undef HAVE_DECL_UCM_TEST_EXTERNAL_EVENTS */

/* Define to 1 if you have the declaration of `UCP_ATOMIC_FETCH_OP_FAND', and
   to 0 if you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_FETCH_OP_FAND */

/* Define to 1 if you have the declaration of `UCP_ATOMIC_FETCH_OP_FOR', and
   to 0 if you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_FETCH_OP_FOR */

/* Define to 1 if you have the declaration of `UCP_ATOMIC_FETCH_OP_FXOR', and
   to 0 if you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_FETCH_OP_FXOR */

/* Define to 1 if you have the declaration of `ucp_atomic_op_nbx', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_OP_NBX */

/* Define to 1 if you have the declaration of `UCP_ATOMIC_POST_OP_AND', and to
   0 if you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_POST_OP_AND */

/* Define to 1 if you have the declaration of `UCP_ATOMIC_POST_OP_OR', and to
   0 if you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_POST_OP_OR */

/* Define to 1 if you have the declaration of `UCP_ATOMIC_POST_OP_XOR', and to
   0 if you don't. */
/* #undef HAVE_DECL_UCP_ATOMIC_POST_OP_XOR */

/* Define to 1 if you have the declaration of `UCP_ATTR_FIELD_MEMORY_TYPES',
   and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_ATTR_FIELD_MEMORY_TYPES */

/* Define to 1 if you have the declaration of `UCP_EP_ATTR_FIELD_TRANSPORTS',
   and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_EP_ATTR_FIELD_TRANSPORTS */

/* Define to 1 if you have the declaration of `ucp_ep_flush_nb', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_EP_FLUSH_NB */

/* Define to 1 if you have the declaration of `ucp_ep_flush_nbx', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_EP_FLUSH_NBX */

/* Define to 1 if you have the declaration of `ucp_get_nb', and to 0 if you
   don't. */
/* #undef HAVE_DECL_UCP_GET_NB */

/* Define to 1 if you have the declaration of `ucp_get_nbx', and to 0 if you
   don't. */
/* #undef HAVE_DECL_UCP_GET_NBX */

/* Define to 1 if you have the declaration of `UCP_MEM_MAP_SYMMETRIC_RKEY',
   and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_MEM_MAP_SYMMETRIC_RKEY */

/* Define to 1 if you have the declaration of `UCP_OP_ATTR_FLAG_MULTI_SEND',
   and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_OP_ATTR_FLAG_MULTI_SEND */

/* Define to 1 if you have the declaration of
   `UCP_PARAM_FIELD_ESTIMATED_NUM_PPN', and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_PARAM_FIELD_ESTIMATED_NUM_PPN */

/* Define to 1 if you have the declaration of `ucp_put_nb', and to 0 if you
   don't. */
/* #undef HAVE_DECL_UCP_PUT_NB */

/* Define to 1 if you have the declaration of `ucp_put_nbx', and to 0 if you
   don't. */
/* #undef HAVE_DECL_UCP_PUT_NBX */

/* Define to 1 if you have the declaration of `ucp_request_check_status', and
   to 0 if you don't. */
/* #undef HAVE_DECL_UCP_REQUEST_CHECK_STATUS */

/* Define to 1 if you have the declaration of `ucp_rkey_compare', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_RKEY_COMPARE */

/* Define to 1 if you have the declaration of `ucp_tag_recv_nbx', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_TAG_RECV_NBX */

/* Define to 1 if you have the declaration of `ucp_tag_send_nbr', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_TAG_SEND_NBR */

/* Define to 1 if you have the declaration of `ucp_tag_send_nbx', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCP_TAG_SEND_NBX */

/* Define to 1 if you have the declaration of `ucp_tag_send_sync_nbx', and to
   0 if you don't. */
/* #undef HAVE_DECL_UCP_TAG_SEND_SYNC_NBX */

/* Define to 1 if you have the declaration of
   `UCP_WORKER_ATTR_FIELD_ADDRESS_FLAGS', and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_WORKER_ATTR_FIELD_ADDRESS_FLAGS */

/* Define to 1 if you have the declaration of
   `UCP_WORKER_FLAG_IGNORE_REQUEST_LEAK', and to 0 if you don't. */
/* #undef HAVE_DECL_UCP_WORKER_FLAG_IGNORE_REQUEST_LEAK */

/* Define to 1 if you have the declaration of `ucp_worker_flush_nb', and to 0
   if you don't. */
/* #undef HAVE_DECL_UCP_WORKER_FLUSH_NB */

/* Define to 1 if you have the declaration of `UCS_MEMORY_TYPE_RDMA', and to 0
   if you don't. */
/* #undef HAVE_DECL_UCS_MEMORY_TYPE_RDMA */

/* Define to 1 if you have the declaration of `UCT_CB_FLAG_SYNC', and to 0 if
   you don't. */
/* #undef HAVE_DECL_UCT_CB_FLAG_SYNC */

/* Define to 1 if you have the declaration of `UCT_PROGRESS_THREAD_SAFE', and
   to 0 if you don't. */
/* #undef HAVE_DECL_UCT_PROGRESS_THREAD_SAFE */

/* Define to 1 if you have the declaration of `__func__', and to 0 if you
   don't. */
#define HAVE_DECL___FUNC__ 1

/* Define to 1 if you have the declaration of `__syscall', and to 0 if you
   don't. */
#define HAVE_DECL___SYSCALL 0

/* Define to 1 if you have the <dirent.h> header file. */
#define HAVE_DIRENT_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `dlsym' function. */
/* #undef HAVE_DLSYM */

/* Define to 1 if the system has the type `Dl_info'. */
/* #undef HAVE_DL_INFO */

/* Define to 1 if the system has the type `double _Complex'. */
#define HAVE_DOUBLE__COMPLEX 1

/* Define to 1 if you have the <endian.h> header file. */
#define HAVE_ENDIAN_H 1

/* Define to 1 if you have the <err.h> header file. */
#define HAVE_ERR_H 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the `execve' function. */
#define HAVE_EXECVE 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if the system has the type `float _Complex'. */
#define HAVE_FLOAT__COMPLEX 1

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Define to 1 if you have the `getpwuid' function. */
#define HAVE_GETPWUID 1

/* Define to 1 if you have the <glob.h> header file. */
/* #undef HAVE_GLOB_H */

/* Define to 1 if you have the `GNI_GetJobResInfo' function. */
/* #undef HAVE_GNI_GETJOBRESINFO */

/* Define to 1 if you have the <grp.h> header file. */
#define HAVE_GRP_H 1

/* Define to 1 if you have the `hcoll_context_free' function. */
/* #undef HAVE_HCOLL_CONTEXT_FREE */

/* Define to 1 if you have the <hostLib.h> header file. */
/* #undef HAVE_HOSTLIB_H */

/* Define to 1 if you have the `hwloc_topology_dup' function. */
/* #undef HAVE_HWLOC_TOPOLOGY_DUP */

/* Define to 1 if you have the <ieee754.h> header file. */
#define HAVE_IEEE754_H 1

/* Define to 1 if you have the <ifaddrs.h> header file. */
#define HAVE_IFADDRS_H 1

/* Define to 1 if the system has the type `int128_t'. */
/* #undef HAVE_INT128_T */

/* Define to 1 if the system has the type `int16_t'. */
#define HAVE_INT16_T 1

/* Define to 1 if the system has the type `int32_t'. */
#define HAVE_INT32_T 1

/* Define to 1 if the system has the type `int64_t'. */
#define HAVE_INT64_T 1

/* Define to 1 if the system has the type `int8_t'. */
#define HAVE_INT8_T 1

/* Define to 1 if the system has the type `intptr_t'. */
#define HAVE_INTPTR_T 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <ioLib.h> header file. */
/* #undef HAVE_IOLIB_H */

/* Define to 1 if you have the `isatty' function. */
#define HAVE_ISATTY 1

/* Define to 1 if you have the <libgen.h> header file. */
#define HAVE_LIBGEN_H 1

/* Define to 1 if you have the <libutil.h> header file. */
/* #undef HAVE_LIBUTIL_H */

/* Define to 1 if you have the <linux/ethtool.h> header file. */
#define HAVE_LINUX_ETHTOOL_H 1

/* Define to 1 if you have the <linux/kcmp.h> header file. */
#define HAVE_LINUX_KCMP_H 1

/* Define to 1 if you have the <linux/mman.h> header file. */
#define HAVE_LINUX_MMAN_H 1

/* Define to 1 if you have the <linux/sockios.h> header file. */
#define HAVE_LINUX_SOCKIOS_H 1

/* Define to 1 if the system has the type `long double'. */
#define HAVE_LONG_DOUBLE 1

/* Define to 1 if the system has the type `long double _Complex'. */
#define HAVE_LONG_DOUBLE__COMPLEX 1

/* Define to 1 if the system has the type `long long'. */
#define HAVE_LONG_LONG 1

/* Define to 1 if you have the <mach/mach_time.h> header file. */
/* #undef HAVE_MACH_MACH_TIME_H */

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <minix/config.h> header file. */
/* #undef HAVE_MINIX_CONFIG_H */

/* Define to 1 if you have the `mkfifo' function. */
#define HAVE_MKFIFO 1

/* Define to 1 if you have the `mmap' function. */
#define HAVE_MMAP 1

/* Define to 1 if you have the <mntent.h> header file. */
#define HAVE_MNTENT_H 1

/* Define to 1 if the system has the type `mode_t'. */
#define HAVE_MODE_T 1

/* Define to 1 if you have the <ndbm.h> header file. */
/* #undef HAVE_NDBM_H */

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

/* Define to 1 if you have the `on_exit' function. */
#define HAVE_ON_EXIT 1

/* Define to 1 if the system has the type `opal_short_float_complex_t'. */
#define HAVE_OPAL_SHORT_FLOAT_COMPLEX_T 1

/* Define to 1 if the system has the type `opal_short_float_t'. */
#define HAVE_OPAL_SHORT_FLOAT_T 1

/* Define to 1 if you have the `openpty' function. */
#define HAVE_OPENPTY 1

/* Define to 1 if you have the <paths.h> header file. */
#define HAVE_PATHS_H 1

/* Define to 1 if you have the `pipe' function. */
#define HAVE_PIPE 1

/* Define to 1 if you have the <poll.h> header file. */
#define HAVE_POLL_H 1

/* Define to 1 if you have the `posix_memalign' function. */
#define HAVE_POSIX_MEMALIGN 1

/* Define to 1 if you have the `preadv' function. */
#define HAVE_PREADV 1

/* Define to 1 if you have the `printstack' function. */
/* #undef HAVE_PRINTSTACK */

/* have PSM2_LIB_REFCOUNT_CAP in psm2.h */
/* #undef HAVE_PSM2_LIB_REFCOUNT_CAP */

/* Define to 1 if you have the `pthread_condattr_setpshared' function. */
#define HAVE_PTHREAD_CONDATTR_SETPSHARED 1

/* Define to 1 if you have the <pthread.h> header file. */
#define HAVE_PTHREAD_H 1

/* Define to 1 if you have the `pthread_mutexattr_setpshared' function. */
#define HAVE_PTHREAD_MUTEXATTR_SETPSHARED 1

/* Define to 1 if the system has the type `ptrdiff_t'. */
#define HAVE_PTRDIFF_T 1

/* Define to 1 if you have the `ptsname' function. */
#define HAVE_PTSNAME 1

/* Define to 1 if you have the <pty.h> header file. */
#define HAVE_PTY_H 1

/* Define to 1 if you have the <pwd.h> header file. */
#define HAVE_PWD_H 1

/* Define to 1 if you have the `pwritev' function. */
#define HAVE_PWRITEV 1

/* Define to 1 if you have the <rdma/fi_ext.h> header file. */
/* #undef HAVE_RDMA_FI_EXT_H */

/* Define to 1 if you have the `regcmp' function. */
/* #undef HAVE_REGCMP */

/* Define to 1 if you have the `regexec' function. */
#define HAVE_REGEXEC 1

/* Define to 1 if you have the <regex.h> header file. */
#define HAVE_REGEX_H 1

/* Define to 1 if you have the `regfree' function. */
#define HAVE_REGFREE 1

/* Define to 1 if you have the <sched.h> header file. */
#define HAVE_SCHED_H 1

/* Define to 1 if you have the `sem_init' function. */
#define HAVE_SEM_INIT 1

/* Define to 1 if you have the `sem_open' function. */
#define HAVE_SEM_OPEN 1

/* Define to 1 if you have the `setenv' function. */
#define HAVE_SETENV 1

/* Define to 1 if you have the `setpgid' function. */
#define HAVE_SETPGID 1

/* Define to 1 if you have the `setsid' function. */
#define HAVE_SETSID 1

/* Define to 1 if you have the <shlwapi.h> header file. */
/* #undef HAVE_SHLWAPI_H */

/* Define to 1 if the system has the type `short float'. */
/* #undef HAVE_SHORT_FLOAT */

/* Define to 1 if the system has the type `short float _Complex'. */
/* #undef HAVE_SHORT_FLOAT__COMPLEX */

/* Define to 1 if `si_band' is a member of `siginfo_t'. */
#define HAVE_SIGINFO_T_SI_BAND 1

/* Define to 1 if `si_fd' is a member of `siginfo_t'. */
#define HAVE_SIGINFO_T_SI_FD 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `socketpair' function. */
#define HAVE_SOCKETPAIR 1

/* Define to 1 if the system has the type `socklen_t'. */
#define HAVE_SOCKLEN_T 1

/* Define to 1 if you have the <sockLib.h> header file. */
/* #undef HAVE_SOCKLIB_H */

/* Define to 1 if the system has the type `ssize_t'. */
#define HAVE_SSIZE_T 1

/* Define to 1 if you have the `statfs' function. */
#define HAVE_STATFS 1

/* Define to 1 if you have the `statvfs' function. */
#define HAVE_STATVFS 1

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

/* Define to 1 if you have the <stropts.h> header file. */
/* #undef HAVE_STROPTS_H */

/* Define to 1 if you have the `strsignal' function. */
#define HAVE_STRSIGNAL 1

/* Define to 1 if `d_type' is a member of `struct dirent'. */
#define HAVE_STRUCT_DIRENT_D_TYPE 1

/* Define to 1 if the system has the type `struct ethtool_cmd'. */
#define HAVE_STRUCT_ETHTOOL_CMD 1

/* Define to 1 if `speed_hi' is a member of `struct ethtool_cmd'. */
#define HAVE_STRUCT_ETHTOOL_CMD_SPEED_HI 1

/* Define to 1 if the system has the type `struct fi_ops_mem_monitor'. */
/* #undef HAVE_STRUCT_FI_OPS_MEM_MONITOR */

/* Define to 1 if the system has the type `struct ifreq'. */
#define HAVE_STRUCT_IFREQ 1

/* Define to 1 if `ifr_hwaddr' is a member of `struct ifreq'. */
#define HAVE_STRUCT_IFREQ_IFR_HWADDR 1

/* Define to 1 if `ifr_mtu' is a member of `struct ifreq'. */
#define HAVE_STRUCT_IFREQ_IFR_MTU 1

/* Define to 1 if the system has the type `struct sockaddr_in'. */
#define HAVE_STRUCT_SOCKADDR_IN 1

/* Define to 1 if the system has the type `struct sockaddr_in6'. */
#define HAVE_STRUCT_SOCKADDR_IN6 1

/* Define to 1 if `sa_len' is a member of `struct sockaddr'. */
/* #undef HAVE_STRUCT_SOCKADDR_SA_LEN */

/* Define to 1 if the system has the type `struct sockaddr_storage'. */
#define HAVE_STRUCT_SOCKADDR_STORAGE 1

/* Define to 1 if `f_fstypename' is a member of `struct statfs'. */
/* #undef HAVE_STRUCT_STATFS_F_FSTYPENAME */

/* Define to 1 if `f_type' is a member of `struct statfs'. */
#define HAVE_STRUCT_STATFS_F_TYPE 1

/* Define to 1 if `f_basetype' is a member of `struct statvfs'. */
/* #undef HAVE_STRUCT_STATVFS_F_BASETYPE */

/* Define to 1 if `f_fstypename' is a member of `struct statvfs'. */
/* #undef HAVE_STRUCT_STATVFS_F_FSTYPENAME */

/* Define to 1 if `tv_nsec' is a member of `struct timespec'. */
#define HAVE_STRUCT_TIMESPEC_TV_NSEC 1

/* Define to 1 if you have the `syscall' function. */
/* #undef HAVE_SYSCALL */

/* Define to 1 if you have the `sysconf' function. */
#define HAVE_SYSCONF 1

/* Define to 1 if you have the <syslimits.h> header file. */
/* #undef HAVE_SYSLIMITS_H */

/* Define to 1 if you have the `syslog' function. */
#define HAVE_SYSLOG 1

/* Define to 1 if you have the <syslog.h> header file. */
#define HAVE_SYSLOG_H 1

/* Define to 1 if you have the <sys/fcntl.h> header file. */
#define HAVE_SYS_FCNTL_H 1

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#define HAVE_SYS_IOCTL_H 1

/* Define to 1 if you have the <sys/ipc.h> header file. */
#define HAVE_SYS_IPC_H 1

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H 1

/* Define to 1 if you have the <sys/mount.h> header file. */
#define HAVE_SYS_MOUNT_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/prctl.h> header file. */
#define HAVE_SYS_PRCTL_H 1

/* Define to 1 if you have the <sys/queue.h> header file. */
#define HAVE_SYS_QUEUE_H 1

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/select.h> header file. */
#define HAVE_SYS_SELECT_H 1

/* Define to 1 if you have the <sys/shm.h> header file. */
#define HAVE_SYS_SHM_H 1

/* Define to 1 if you have the <sys/socket.h> header file. */
#define HAVE_SYS_SOCKET_H 1

/* Define to 1 if you have the <sys/sockio.h> header file. */
/* #undef HAVE_SYS_SOCKIO_H */

/* Define to 1 if you have the <sys/statfs.h> header file. */
#define HAVE_SYS_STATFS_H 1

/* Define to 1 if you have the <sys/statvfs.h> header file. */
#define HAVE_SYS_STATVFS_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/synch.h> header file. */
/* #undef HAVE_SYS_SYNCH_H */

/* Define to 1 if you have the <sys/syscall.h> header file. */
#define HAVE_SYS_SYSCALL_H 1

/* Define to 1 if you have the <sys/sysctl.h> header file. */
/* #undef HAVE_SYS_SYSCTL_H */

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/tree.h> header file. */
/* #undef HAVE_SYS_TREE_H */

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <sys/uio.h> header file. */
#define HAVE_SYS_UIO_H 1

/* Define to 1 if you have the <sys/un.h> header file. */
#define HAVE_SYS_UN_H 1

/* Define to 1 if you have the <sys/utsname.h> header file. */
#define HAVE_SYS_UTSNAME_H 1

/* Define to 1 if you have the <sys/vfs.h> header file. */
#define HAVE_SYS_VFS_H 1

/* Define to 1 if you have the <sys/wait.h> header file. */
#define HAVE_SYS_WAIT_H 1

/* Define to 1 if you have the <TargetConditionals.h> header file. */
/* #undef HAVE_TARGETCONDITIONALS_H */

/* Define to 1 if you have the `tcgetpgrp' function. */
#define HAVE_TCGETPGRP 1

/* Define to 1 if you have the <termios.h> header file. */
#define HAVE_TERMIOS_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the `ucc_comm_free' function. */
/* #undef HAVE_UCC_COMM_FREE */

/* Define to 1 if you have the <ucontext.h> header file. */
#define HAVE_UCONTEXT_H 1

/* have memory types attribute */
/* #undef HAVE_UCP_ATTR_MEMORY_TYPES */

/* Define to 1 if the system has the type `ucp_request_param_t'. */
/* #undef HAVE_UCP_REQUEST_PARAM_T */

/* have ucp_tag_send_nbr() */
/* #undef HAVE_UCP_TAG_SEND_NBR */

/* have worker address attribute */
/* #undef HAVE_UCP_WORKER_ADDRESS_FLAGS */

/* Define to 1 if the system has the type `uint128_t'. */
/* #undef HAVE_UINT128_T */

/* Define to 1 if the system has the type `uint16_t'. */
#define HAVE_UINT16_T 1

/* Define to 1 if the system has the type `uint32_t'. */
#define HAVE_UINT32_T 1

/* Define to 1 if the system has the type `uint64_t'. */
#define HAVE_UINT64_T 1

/* Define to 1 if the system has the type `uint8_t'. */
#define HAVE_UINT8_T 1

/* Define to 1 if the system has the type `uintptr_t'. */
#define HAVE_UINTPTR_T 1

/* Define to 1 if you have the <ulimit.h> header file. */
#define HAVE_ULIMIT_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* whether unix byteswap routines -- htonl, htons, nothl, ntohs -- are
   available */
#define HAVE_UNIX_BYTESWAP 1

/* Define to 1 if you have the `usleep' function. */
#define HAVE_USLEEP 1

/* Define to 1 if you have the <util.h> header file. */
/* #undef HAVE_UTIL_H */

/* Define to 1 if you have the <utmp.h> header file. */
#define HAVE_UTMP_H 1

/* Define to 1 if you have the <valgrind/valgrind.h> header file. */
/* #undef HAVE_VALGRIND_VALGRIND_H */

/* Define to 1 if you have the `vasprintf' function. */
#define HAVE_VASPRINTF 1

/* Define to 1 if you have the `vsnprintf' function. */
#define HAVE_VSNPRINTF 1

/* Define to 1 if you have the `vsyslog' function. */
#define HAVE_VSYSLOG 1

/* Define to 1 if you have the `waitpid' function. */
#define HAVE_WAITPID 1

/* Define to 1 if you have the <wchar.h> header file. */
#define HAVE_WCHAR_H 1

/* is xpmem.h available */
/* #undef HAVE_XPMEM_H */

/* Define to 1 if the system has the type `_Float16'. */
#define HAVE__FLOAT16 1

/* Define to 1 if you have the `_NSGetEnviron' function. */
/* #undef HAVE__NSGETENVIRON */

/* Define to 1 if you have the `__clear_cache' function. */
#define HAVE___CLEAR_CACHE 1

/* Define to 1 if you have the `__curbrk' function. */
#define HAVE___CURBRK 1

/* Define to 1 if the system has the type `__float128'. */
#define HAVE___FLOAT128 1

/* Define to 1 if the system has the type `__int128'. */
#define HAVE___INT128 1

/* Define to 1 if you have the `__malloc_initialize_hook' function. */
/* #undef HAVE___MALLOC_INITIALIZE_HOOK */

/* Define to 1 if you have the `__munmap' function. */
/* #undef HAVE___MUNMAP */

/* Define to 1 if you have the `__syscall' function. */
/* #undef HAVE___SYSCALL */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Custom matching engine to use in pml/ob1 */
#define MCA_PML_OB1_CUSTOM_MATCHING MCA_PML_OB1_CUSTOM_MATCHING_NONE

/* Complete set of command line arguments given to ROMIOs configure script */
#define MCA_io_romio341_COMPLETE_CONFIGURE_FLAGS " FROM_OMPI=yes CC='gcc' CFLAGS='-O3 -DNDEBUG  -finline-functions -mcx16 -D__EXTENSIONS__' CPPFLAGS='' FFLAGS='' LDFLAGS='' --enable-shared --disable-static   --disable-aio --disable-weak-symbols --enable-strict --disable-f77 --disable-f90 ac_cv_lib_cuda_cuMemGetAddressRange=no ac_cv_lib_cudart_cudaStreamSynchronize=no"

/* Set of user-defined configure flags given to ROMIOs configure script via
   --with-io-romio-flags */
#define MCA_io_romio341_USER_CONFIGURE_FLAGS ""

/* Header to include for memcpy implementation */
#define MCA_memcpy_IMPLEMENTATION_HEADER "opal/mca/memcpy/base/memcpy_base_default.h"

/* Header to include for parts of the memory implementation */
#define MCA_memory_IMPLEMENTATION_HEADER "opal/mca/memory/base/empty.h"

/* Defined to 1 if ompi:mtl should use direct calls instead of components */
#define MCA_ompi_mtl_DIRECT_CALL 0

/* name of component to use for direct calls, if MCA_ompi_mtl_DIRECT_CALL is 1
   */
#define MCA_ompi_mtl_DIRECT_CALL_COMPONENT 

/* Header ompi:mtl includes to be direct called */
#define MCA_ompi_mtl_DIRECT_CALL_HEADER ""

/* Defined to 1 if ompi:part should use direct calls instead of components */
#define MCA_ompi_part_DIRECT_CALL 0

/* name of component to use for direct calls, if MCA_ompi_part_DIRECT_CALL is
   1 */
#define MCA_ompi_part_DIRECT_CALL_COMPONENT 

/* Header ompi:part includes to be direct called */
#define MCA_ompi_part_DIRECT_CALL_HEADER ""

/* Defined to 1 if ompi:pml should use direct calls instead of components */
#define MCA_ompi_pml_DIRECT_CALL 0

/* name of component to use for direct calls, if MCA_ompi_pml_DIRECT_CALL is 1
   */
#define MCA_ompi_pml_DIRECT_CALL_COMPONENT 

/* Header ompi:pml includes to be direct called */
#define MCA_ompi_pml_DIRECT_CALL_HEADER ""

/* Defined to 1 if opal:smsc should use direct calls instead of components */
#define MCA_opal_smsc_DIRECT_CALL 0

/* name of component to use for direct calls, if MCA_opal_smsc_DIRECT_CALL is
   1 */
#define MCA_opal_smsc_DIRECT_CALL_COMPONENT 

/* Header opal:smsc includes to be direct called */
#define MCA_opal_smsc_DIRECT_CALL_HEADER ""

/* Defined to 1 if oshmem:memheap should use direct calls instead of
   components */
#define MCA_oshmem_memheap_DIRECT_CALL 0

/* name of component to use for direct calls, if
   MCA_oshmem_memheap_DIRECT_CALL is 1 */
#define MCA_oshmem_memheap_DIRECT_CALL_COMPONENT 

/* Header oshmem:memheap includes to be direct called */
#define MCA_oshmem_memheap_DIRECT_CALL_HEADER ""

/* Defined to 1 if oshmem:spml should use direct calls instead of components
   */
#define MCA_oshmem_spml_DIRECT_CALL 0

/* name of component to use for direct calls, if MCA_oshmem_spml_DIRECT_CALL
   is 1 */
#define MCA_oshmem_spml_DIRECT_CALL_COMPONENT 

/* Header oshmem:spml includes to be direct called */
#define MCA_oshmem_spml_DIRECT_CALL_HEADER ""

/* Header to include for threads implementation */
#define MCA_threads_base_include_HEADER "opal/mca/threads/pthreads/threads_pthreads_threads.h"

/* Header to include for mutex implementation */
#define MCA_threads_mutex_base_include_HEADER "opal/mca/threads/pthreads/threads_pthreads_mutex.h"

/* Header to include for tsd implementation */
#define MCA_threads_tsd_base_include_HEADER "opal/mca/threads/pthreads/threads_pthreads_tsd.h"

/* Header to include for timer implementation */
#define MCA_timer_IMPLEMENTATION_HEADER "opal/mca/timer/linux/timer_linux.h"

/* Macro that is set to 1 when CUDA-aware support is configured in and 0 when
   it is not */
#define MPIX_CUDA_AWARE_SUPPORT 0

/* Macro that is set to 1 when ROCM-aware support is configured in and 0 when
   it is not */
#define MPIX_ROCM_AWARE_SUPPORT 0

/* Maximum value for an MPI_Count */
#define MPI_COUNT_MAX 0x7fffffffffffffffll

/* Whether we want to check MPI parameters always, never, or decide at
   run-time */
#define MPI_PARAM_CHECK ompi_mpi_param_check

/* MPI Standard Minor version number */
#define MPI_SUBVERSION 1

/* MPI Standard Major version number */
#define MPI_VERSION 3

/* Alignment of Fortran CHARACTER */
#define OMPI_ALIGNMENT_FORTRAN_CHARACTER 1

/* Alignment of Fortran COMPLEX */
#define OMPI_ALIGNMENT_FORTRAN_COMPLEX 4

/* Alignment of Fortran COMPLEX*16 */
#define OMPI_ALIGNMENT_FORTRAN_COMPLEX16 8

/* Alignment of Fortran COMPLEX*32 */
#define OMPI_ALIGNMENT_FORTRAN_COMPLEX32 16

/* Alignment of Fortran COMPLEX*4 */
#define OMPI_ALIGNMENT_FORTRAN_COMPLEX4 0

/* Alignment of Fortran COMPLEX*8 */
#define OMPI_ALIGNMENT_FORTRAN_COMPLEX8 4

/* Alignment of Fortran DOUBLE COMPLEX */
#define OMPI_ALIGNMENT_FORTRAN_DOUBLE_COMPLEX 8

/* Alignment of Fortran DOUBLE PRECISION */
#define OMPI_ALIGNMENT_FORTRAN_DOUBLE_PRECISION 8

/* Alignment of Fortran INTEGER */
#define OMPI_ALIGNMENT_FORTRAN_INTEGER 4

/* Alignment of Fortran INTEGER*1 */
#define OMPI_ALIGNMENT_FORTRAN_INTEGER1 1

/* Alignment of Fortran INTEGER*16 */
#define OMPI_ALIGNMENT_FORTRAN_INTEGER16 0

/* Alignment of Fortran INTEGER*2 */
#define OMPI_ALIGNMENT_FORTRAN_INTEGER2 2

/* Alignment of Fortran INTEGER*4 */
#define OMPI_ALIGNMENT_FORTRAN_INTEGER4 4

/* Alignment of Fortran INTEGER*8 */
#define OMPI_ALIGNMENT_FORTRAN_INTEGER8 8

/* Alignment of Fortran LOGICAL */
#define OMPI_ALIGNMENT_FORTRAN_LOGICAL 4

/* Alignment of Fortran LOGICAL*1 */
#define OMPI_ALIGNMENT_FORTRAN_LOGICAL1 1

/* Alignment of Fortran LOGICAL*2 */
#define OMPI_ALIGNMENT_FORTRAN_LOGICAL2 2

/* Alignment of Fortran LOGICAL*4 */
#define OMPI_ALIGNMENT_FORTRAN_LOGICAL4 4

/* Alignment of Fortran LOGICAL*8 */
#define OMPI_ALIGNMENT_FORTRAN_LOGICAL8 8

/* Alignment of Fortran REAL */
#define OMPI_ALIGNMENT_FORTRAN_REAL 4

/* Alignment of Fortran REAL*16 */
#define OMPI_ALIGNMENT_FORTRAN_REAL16 16

/* Alignment of Fortran REAL*2 */
#define OMPI_ALIGNMENT_FORTRAN_REAL2 0

/* Alignment of Fortran REAL*4 */
#define OMPI_ALIGNMENT_FORTRAN_REAL4 4

/* Alignment of Fortran REAL*8 */
#define OMPI_ALIGNMENT_FORTRAN_REAL8 8

/* The level of fortran bindings to be built */
#define OMPI_BUILD_FORTRAN_BINDINGS 3

/* OMPI underlying C++ compiler */
#define OMPI_CXX "g++"

/* whether we support Grequest extensions */
#define OMPI_ENABLE_GREQUEST_EXTENSIONS 0

/* whether we want MPI-1.x support */
#define OMPI_ENABLE_MPI1_COMPAT 0

/* Underlying Fortran compiler */
#define OMPI_FC "gfortran"

/* Absolute path to the underlying Fortran compiler found by configure */
#define OMPI_FC_ABSOLUTE "/usr/bin/gfortran"

/* Whether the mpif.h interface supports the MPI_SIZEOF interface or not */
#define OMPI_FORTRAN_BUILD_SIZEOF 1

/* Whether fortran symbols are all caps or not */
#define OMPI_FORTRAN_CAPS 0

/* Whether fortran symbols have a trailing double underscore or not */
#define OMPI_FORTRAN_DOUBLE_UNDERSCORE 0

/* How many bytes the mpi_f08 TYPE(MPI_<foo>) handles will be aligned to */
#define OMPI_FORTRAN_F08_HANDLE_ALIGNMENT 4

/* How many bytes the mpi_f08 TYPE(MPI_<foo>) handles will be */
#define OMPI_FORTRAN_F08_HANDLE_SIZE 4

/* Max handle value for fortran MPI handles, effectively min(INT_MAX, max
   fortran INTEGER value) */
#define OMPI_FORTRAN_HANDLE_MAX 2147483647

/* For mpi-f08-interfaces-callbacks.f90 and ompi_info: whether the compiler
   supports the "abstract" keyword or not */
#define OMPI_FORTRAN_HAVE_ABSTRACT 1

/* For ompi/mpi/fortran/use-mpi-f08/blah.F90 and blah.h and ompi_info: whether
   the compiler supports the "asynchronous" keyword or not */
#define OMPI_FORTRAN_HAVE_ASYNCHRONOUS 1

/* Whether the compiler supports Fortran ATTRIBUTES DEPRECATED or not */
#define OMPI_FORTRAN_HAVE_ATTR_DEPRECATED 1

/* For ompi_info: Whether the compiler supports all forms of BIND(C) that we
   need */
#define OMPI_FORTRAN_HAVE_BIND_C 1

/* For ompi_info: Whether the compiler supports SUBROUTINE ... BIND(C) or not
   */
#define OMPI_FORTRAN_HAVE_BIND_C_SUB 1

/* For ompi_info: Whether the compiler supports TYPE, BIND(C) or not */
#define OMPI_FORTRAN_HAVE_BIND_C_TYPE 1

/* For ompi_info: Whether the compiler supports TYPE, BIND(C, NAME="name") or
   not */
#define OMPI_FORTRAN_HAVE_BIND_C_TYPE_NAME 1

/* For ompi/mpi/fortran/use-mpi-f08/blah.F90 and blah.h and ompi_info: whether
   the compiler supports c_funloc or not */
#define OMPI_FORTRAN_HAVE_C_FUNLOC 1

/* For ompi_info: Whether the Fortran compiler supports the Fortran 2008
   "assumed rank" syntax or not */
#define OMPI_FORTRAN_HAVE_F08_ASSUMED_RANK 1

/* Whether the Fortran compiler supports ignore TKR functionality or not */
#define OMPI_FORTRAN_HAVE_IGNORE_TKR 1

/* Whether the compiler supports INTERFACE or not */
#define OMPI_FORTRAN_HAVE_INTERFACE 1

/* For ompi_info: Whether the compiler supports ISO_C_BINDING or not */
#define OMPI_FORTRAN_HAVE_ISO_C_BINDING 1

/* Whether the compiler supports ISO_FORTRAN_ENV or not */
#define OMPI_FORTRAN_HAVE_ISO_FORTRAN_ENV 1

/* For ompi_info: whether the Fortran compiler supports optional arguments or
   not */
#define OMPI_FORTRAN_HAVE_OPTIONAL_ARGS 1

/* For mpif-status.h, mpi-f08-types.f90 and ompi_info: whether the compiler
   supports the "private" keyword or not (used in TYPE(MPI_Status)) */
#define OMPI_FORTRAN_HAVE_PRIVATE 1

/* For ompi/mpi/fortran/use-mpi-f08/blah.F90 and blah.h and ompi_info: whether
   the compiler supports the "procedure" keyword or not */
#define OMPI_FORTRAN_HAVE_PROCEDURE 1

/* Whether the compiler supports STORAGE_SIZE on relevant types */
#define OMPI_FORTRAN_HAVE_STORAGE_SIZE 1

/* For ompi/mpi/fortran/use-mpi-f08/blah.F90 and blah.h and ompi_info: whether
   the compiler supports "USE ... ONLY" notation properly or not */
#define OMPI_FORTRAN_HAVE_USE_ONLY 1

/* Pre declaration for FORTRAN ignore parameter TKR behavior */
#define OMPI_FORTRAN_IGNORE_TKR_PREDECL "!GCC$ ATTRIBUTES NO_ARG_CHECK ::"

/* Type declaration for FORTRAN ignore parameter TKR behavior */
#define OMPI_FORTRAN_IGNORE_TKR_TYPE 

/* Max dimension rank of Fortran arrays */
#define OMPI_FORTRAN_MAX_ARRAY_RANK 15

/* Whether we are building support for the mpif.h bindings or not */
#define OMPI_FORTRAN_MPIFH_BINDINGS 1

/* Whether the mpi_f08 implementation is using wrapper routines ("bad" Fortran
   compiler) or weak symbols ("good" Fortran compiler) for the F08 interface
   definition implementations */
#define OMPI_FORTRAN_NEED_WRAPPER_ROUTINES 1

/* Whether fortran symbols have no trailing underscore or not */
#define OMPI_FORTRAN_PLAIN 0

/* Whether fortran symbols have a trailing underscore or not */
#define OMPI_FORTRAN_SINGLE_UNDERSCORE 1

/* The number or Fortran INTEGER in MPI Status */
#define OMPI_FORTRAN_STATUS_SIZE 6

/* Whether we are building support for the "use mpif08" bindings or not */
#define OMPI_FORTRAN_USEMPIF08_BINDINGS 3

/* Whether we are building support for the "use mpi" bindings or not */
#define OMPI_FORTRAN_USEMPI_BINDINGS 2

/* Fortran value for LOGICAL .TRUE. value */
#define OMPI_FORTRAN_VALUE_TRUE 1

/* Greek - alpha, beta, etc - release number of Open MPI */
#define OMPI_GREEK_VERSION ""

/* Whether we want sparse process groups */
#define OMPI_GROUP_SPARSE 0

/* Whether we have Fortran CHARACTER or not */
#define OMPI_HAVE_FORTRAN_CHARACTER 1

/* Whether we have Fortran COMPLEX or not */
#define OMPI_HAVE_FORTRAN_COMPLEX 1

/* Whether we have Fortran COMPLEX*16 or not */
#define OMPI_HAVE_FORTRAN_COMPLEX16 1

/* Whether we have Fortran COMPLEX*32 or not */
#define OMPI_HAVE_FORTRAN_COMPLEX32 1

/* Whether we have Fortran COMPLEX*4 or not */
#define OMPI_HAVE_FORTRAN_COMPLEX4 0

/* Whether we have Fortran COMPLEX*8 or not */
#define OMPI_HAVE_FORTRAN_COMPLEX8 1

/* Whether we have Fortran DOUBLE COMPLEX or not */
#define OMPI_HAVE_FORTRAN_DOUBLE_COMPLEX 1

/* Whether we have Fortran DOUBLE PRECISION or not */
#define OMPI_HAVE_FORTRAN_DOUBLE_PRECISION 1

/* Whether we have Fortran INTEGER or not */
#define OMPI_HAVE_FORTRAN_INTEGER 1

/* Whether we have Fortran INTEGER*1 or not */
#define OMPI_HAVE_FORTRAN_INTEGER1 1

/* Whether we have Fortran INTEGER*16 or not */
#define OMPI_HAVE_FORTRAN_INTEGER16 0

/* Whether we have Fortran INTEGER*2 or not */
#define OMPI_HAVE_FORTRAN_INTEGER2 1

/* Whether we have Fortran INTEGER*4 or not */
#define OMPI_HAVE_FORTRAN_INTEGER4 1

/* Whether we have Fortran INTEGER*8 or not */
#define OMPI_HAVE_FORTRAN_INTEGER8 1

/* Whether we have Fortran LOGICAL or not */
#define OMPI_HAVE_FORTRAN_LOGICAL 1

/* Whether we have Fortran LOGICAL*1 or not */
#define OMPI_HAVE_FORTRAN_LOGICAL1 1

/* Whether we have Fortran LOGICAL*2 or not */
#define OMPI_HAVE_FORTRAN_LOGICAL2 1

/* Whether we have Fortran LOGICAL*4 or not */
#define OMPI_HAVE_FORTRAN_LOGICAL4 1

/* Whether we have Fortran LOGICAL*8 or not */
#define OMPI_HAVE_FORTRAN_LOGICAL8 1

/* Whether we have Fortran REAL or not */
#define OMPI_HAVE_FORTRAN_REAL 1

/* Whether we have Fortran REAL*16 or not */
#define OMPI_HAVE_FORTRAN_REAL16 1

/* Whether we have Fortran REAL*2 or not */
#define OMPI_HAVE_FORTRAN_REAL2 0

/* Whether we have Fortran REAL*4 or not */
#define OMPI_HAVE_FORTRAN_REAL4 1

/* Whether we have Fortran REAL*8 or not */
#define OMPI_HAVE_FORTRAN_REAL8 1

/* Whether or not PRRTE is available */
#define OMPI_HAVE_PRRTE 1

/* Fortrn KIND number for CHARACTER */
#define OMPI_KIND_FORTRAN_CHARACTER C_SIGNED_CHAR

/* Fortrn KIND number for COMPLEX */
#define OMPI_KIND_FORTRAN_COMPLEX C_FLOAT_COMPLEX

/* Fortrn KIND number for COMPLEX*16 */
#define OMPI_KIND_FORTRAN_COMPLEX16 C_DOUBLE_COMPLEX

/* Fortrn KIND number for COMPLEX*32 */
#define OMPI_KIND_FORTRAN_COMPLEX32 C_LONG_DOUBLE_COMPLEX

/* Fortrn KIND number for COMPLEX*4 */
#define OMPI_KIND_FORTRAN_COMPLEX4 0

/* Fortrn KIND number for COMPLEX*8 */
#define OMPI_KIND_FORTRAN_COMPLEX8 C_FLOAT_COMPLEX

/* Fortrn KIND number for DOUBLE COMPLEX */
#define OMPI_KIND_FORTRAN_DOUBLE_COMPLEX C_DOUBLE_COMPLEX

/* Fortrn KIND number for DOUBLE PRECISION */
#define OMPI_KIND_FORTRAN_DOUBLE_PRECISION C_DOUBLE

/* Fortrn KIND number for INTEGER */
#define OMPI_KIND_FORTRAN_INTEGER C_INT

/* Fortrn KIND number for INTEGER*1 */
#define OMPI_KIND_FORTRAN_INTEGER1 C_SIGNED_CHAR

/* Fortrn KIND number for INTEGER*16 */
#define OMPI_KIND_FORTRAN_INTEGER16 0

/* Fortrn KIND number for INTEGER*2 */
#define OMPI_KIND_FORTRAN_INTEGER2 C_SHORT

/* Fortrn KIND number for INTEGER*4 */
#define OMPI_KIND_FORTRAN_INTEGER4 C_INT

/* Fortrn KIND number for INTEGER*8 */
#define OMPI_KIND_FORTRAN_INTEGER8 C_LONG_LONG

/* Fortrn KIND number for LOGICAL */
#define OMPI_KIND_FORTRAN_LOGICAL C_INT

/* Fortrn KIND number for LOGICAL*1 */
#define OMPI_KIND_FORTRAN_LOGICAL1 C_SIGNED_CHAR

/* Fortrn KIND number for LOGICAL*2 */
#define OMPI_KIND_FORTRAN_LOGICAL2 C_SHORT

/* Fortrn KIND number for LOGICAL*4 */
#define OMPI_KIND_FORTRAN_LOGICAL4 C_INT

/* Fortrn KIND number for LOGICAL*8 */
#define OMPI_KIND_FORTRAN_LOGICAL8 C_LONG_LONG

/* Fortrn KIND number for REAL */
#define OMPI_KIND_FORTRAN_REAL C_FLOAT

/* Fortrn KIND number for REAL*16 */
#define OMPI_KIND_FORTRAN_REAL16 C_LONG_DOUBLE

/* Fortrn KIND number for REAL*2 */
#define OMPI_KIND_FORTRAN_REAL2 0

/* Fortrn KIND number for REAL*4 */
#define OMPI_KIND_FORTRAN_REAL4 C_FLOAT

/* Fortrn KIND number for REAL*8 */
#define OMPI_KIND_FORTRAN_REAL8 C_DOUBLE

/* Major release number of Open MPI */
#define OMPI_MAJOR_VERSION 5

/* AVX supported in the current build */
#define OMPI_MCA_OP_HAVE_AVX 1

/* AVX2 supported in the current build */
#define OMPI_MCA_OP_HAVE_AVX2 1

/* AVX512 supported in the current build */
#define OMPI_MCA_OP_HAVE_AVX512 1

/* NEON supported in the current build */
/* #undef OMPI_MCA_OP_HAVE_NEON */

/* NEON FP supported in the current build */
/* #undef OMPI_MCA_OP_HAVE_NEON_FP */

/* SSE3 supported in the current build */
#define OMPI_MCA_OP_HAVE_SSE3 1

/* SSE4.1 supported in the current build */
#define OMPI_MCA_OP_HAVE_SSE41 1

/* SVE supported in the current build */
/* #undef OMPI_MCA_OP_HAVE_SVE */

/* SVE not supported */
#define OMPI_MCA_OP_SVE_EXTRA_FLAGS 0

/* Minor release number of Open MPI */
#define OMPI_MINOR_VERSION 0

/* MPI Extensions included in libmpi */
#define OMPI_MPIEXT_COMPONENTS "affinity, cuda, ftmpi, rocm, shortfloat"

/* Type of MPI_Aint */
#define OMPI_MPI_AINT_TYPE ptrdiff_t

/* Size of the MPI_Count datatype */
#define OMPI_MPI_COUNT_SIZE 8

/* Type of the MPI_Count datatype */
#define OMPI_MPI_COUNT_TYPE long long

/* Size of the MPI_Offset */
#define OMPI_MPI_OFFSET_SIZE 8

/* Type of MPI_Offset */
#define OMPI_MPI_OFFSET_TYPE long long

/* Enable flow control for Portals4 MTL */
#define OMPI_MTL_PORTALS4_FLOW_CONTROL 1

/* MPI datatype corresponding to MPI_Offset */
#define OMPI_OFFSET_DATATYPE MPI_LONG_LONG

/* Whether we want to check MPI parameters never or possible (an integer
   constant) */
#define OMPI_PARAM_CHECK 1

/* Index into endpoint array for BML */
#define OMPI_PROC_ENDPOINT_TAG_BML 0

/* Maximum number of endpoint entries to be attached to an ompi_proc_t */
#define OMPI_PROC_ENDPOINT_TAG_MAX 1

/* Index into endpoint array for MTL */
/* #undef OMPI_PROC_ENDPOINT_TAG_MTL */

/* Index into endpoint array for PML */
/* #undef OMPI_PROC_ENDPOINT_TAG_PML */

/* Index into endpoint array for PORTALS4 */
/* #undef OMPI_PROC_ENDPOINT_TAG_PORTALS4 */

/* Index into endpoint array for UCX */
/* #undef OMPI_PROC_ENDPOINT_TAG_UCX */

/* Path to prterun */
/* #undef OMPI_PRTERUN_PATH */

/* Whether Fortran REAL*16 matches the bit format of the equivalent C type */
#define OMPI_REAL16_MATCHES_C 1

/* Release date of Open MPI */
#define OMPI_RELEASE_DATE "May 30, 2025"

/* Release release number of Open MPI */
#define OMPI_RELEASE_VERSION 8

/* The repository version Open MPI */
#define OMPI_REPO_REV "v5.0.8"

/* Size of Fortran CHARACTER */
#define OMPI_SIZEOF_FORTRAN_CHARACTER 1

/* Size of Fortran COMPLEX */
#define OMPI_SIZEOF_FORTRAN_COMPLEX 8

/* Size of Fortran COMPLEX*16 */
#define OMPI_SIZEOF_FORTRAN_COMPLEX16 16

/* Size of Fortran COMPLEX*32 */
#define OMPI_SIZEOF_FORTRAN_COMPLEX32 32

/* Size of Fortran COMPLEX*4 */
#define OMPI_SIZEOF_FORTRAN_COMPLEX4 0

/* Size of Fortran COMPLEX*8 */
#define OMPI_SIZEOF_FORTRAN_COMPLEX8 8

/* Size of Fortran DOUBLE COMPLEX */
#define OMPI_SIZEOF_FORTRAN_DOUBLE_COMPLEX 16

/* Size of Fortran DOUBLE PRECISION */
#define OMPI_SIZEOF_FORTRAN_DOUBLE_PRECISION 8

/* Size of Fortran INTEGER */
#define OMPI_SIZEOF_FORTRAN_INTEGER 4

/* Size of Fortran INTEGER*1 */
#define OMPI_SIZEOF_FORTRAN_INTEGER1 1

/* Size of Fortran INTEGER*16 */
#define OMPI_SIZEOF_FORTRAN_INTEGER16 16

/* Size of Fortran INTEGER*2 */
#define OMPI_SIZEOF_FORTRAN_INTEGER2 2

/* Size of Fortran INTEGER*4 */
#define OMPI_SIZEOF_FORTRAN_INTEGER4 4

/* Size of Fortran INTEGER*8 */
#define OMPI_SIZEOF_FORTRAN_INTEGER8 8

/* Size of Fortran LOGICAL */
#define OMPI_SIZEOF_FORTRAN_LOGICAL 4

/* Size of Fortran LOGICAL*1 */
#define OMPI_SIZEOF_FORTRAN_LOGICAL1 1

/* Size of Fortran LOGICAL*2 */
#define OMPI_SIZEOF_FORTRAN_LOGICAL2 2

/* Size of Fortran LOGICAL*4 */
#define OMPI_SIZEOF_FORTRAN_LOGICAL4 4

/* Size of Fortran LOGICAL*8 */
#define OMPI_SIZEOF_FORTRAN_LOGICAL8 8

/* Size of Fortran REAL */
#define OMPI_SIZEOF_FORTRAN_REAL 4

/* Size of Fortran REAL*16 */
#define OMPI_SIZEOF_FORTRAN_REAL16 16

/* Size of Fortran REAL*2 */
#define OMPI_SIZEOF_FORTRAN_REAL2 0

/* Size of Fortran REAL*4 */
#define OMPI_SIZEOF_FORTRAN_REAL4 4

/* Size of Fortran REAL*8 */
#define OMPI_SIZEOF_FORTRAN_REAL8 8

/* Tarball filename version string of Open MPI */
#define OMPI_TARBALL_VERSION "5.0.8"

/* Whether or not we are using the internal PRRTE */
#define OMPI_USING_INTERNAL_PRRTE 1

/* Complete release number of Open MPI */
#define OMPI_VERSION "5.0.8"

/* do we want java mpi bindings */
#define OMPI_WANT_JAVA_BINDINGS 0

/* Enable warnings when using deprecated MPI functions */
#define OMPI_WANT_MPI_INTERFACE_WARNING 1

/* if the peruse interface should be enabled */
#define OMPI_WANT_PERUSE 0

/* CFLAGS to pass through the wrapper compilers */
#define OMPI_WRAPPER_CFLAGS " "

/* CXXFLAGS to pass through the wrapper compilers */
#define OMPI_WRAPPER_CXXFLAGS " "

/* FCFLAGS to pass through the wrapper compilers */
#define OMPI_WRAPPER_FCFLAGS ""

/* LDFLAGS to pass through the wrapper compilers */
#define OMPI_WRAPPER_LDFLAGS "-L${libdir}  -Wl,-rpath -Wl,${libdir} -Wl,--enable-new-dtags"

/* LIBS to pass through the wrapper compilers */
#define OMPI_WRAPPER_LIBS "-lmpi"

/* Alignment of bool */
#define OPAL_ALIGNMENT_BOOL 1

/* Alignment of char */
#define OPAL_ALIGNMENT_CHAR 1

/* Alignment of double */
#define OPAL_ALIGNMENT_DOUBLE 8

/* Alignment of double _Complex */
#define OPAL_ALIGNMENT_DOUBLE_COMPLEX 8

/* Alignment of float */
#define OPAL_ALIGNMENT_FLOAT 4

/* Alignment of float _Complex */
#define OPAL_ALIGNMENT_FLOAT_COMPLEX 4

/* Alignment of int */
#define OPAL_ALIGNMENT_INT 4

/* Alignment of int128_t */
/* #undef OPAL_ALIGNMENT_INT128 */

/* Alignment of int16_t */
#define OPAL_ALIGNMENT_INT16 2

/* Alignment of int32_t */
#define OPAL_ALIGNMENT_INT32 4

/* Alignment of int64_t */
#define OPAL_ALIGNMENT_INT64 8

/* Alignment of int8_t */
#define OPAL_ALIGNMENT_INT8 1

/* Alignment of long */
#define OPAL_ALIGNMENT_LONG 8

/* Alignment of long double */
#define OPAL_ALIGNMENT_LONG_DOUBLE 16

/* Alignment of long double _Complex */
#define OPAL_ALIGNMENT_LONG_DOUBLE_COMPLEX 16

/* Alignment of long long */
#define OPAL_ALIGNMENT_LONG_LONG 8

/* Alignment of opal_short_float_t */
#define OPAL_ALIGNMENT_OPAL_SHORT_FLOAT_T 2

/* Alignment of short */
#define OPAL_ALIGNMENT_SHORT 2

/* Alignment of short float */
/* #undef OPAL_ALIGNMENT_SHORT_FLOAT */

/* Alignment of short float _Complex */
/* #undef OPAL_ALIGNMENT_SHORT_FLOAT_COMPLEX */

/* Alignment of size_t */
#define OPAL_ALIGNMENT_SIZE_T 8

/* Alignment of void * */
#define OPAL_ALIGNMENT_VOID_P 8

/* Alignment of wchar_t */
#define OPAL_ALIGNMENT_WCHAR 4

/* Alignment of __float128 */
#define OPAL_ALIGNMENT___FLOAT128 16

/* set to 1 if word-size integers must be aligned to word-size padding to
   prevent bus errors */
#define OPAL_ALIGN_WORD_SIZE_INTEGERS 0

/* OMPI architecture string */
#define OPAL_ARCH "x86_64-pc-linux-gnu"

/* Enable flow control for Portals4 BTL */
#define OPAL_BTL_PORTALS4_FLOW_CONTROL 1

/* define to 1 if usnic BTL unit tests are enabled, 0 otherwise */
#define OPAL_BTL_USNIC_UNIT_TESTS 0

/* The compiler $lower which OMPI was built with */
#define OPAL_BUILD_PLATFORM_COMPILER_FAMILYID 1

/* The compiler $lower which OMPI was built with */
#define OPAL_BUILD_PLATFORM_COMPILER_VERSION 983553

/* OMPI underlying C compiler */
#define OPAL_CC "gcc"

/* If kcmp is available */
#define OPAL_CMA_KCMP_AVAIL 1

/* Need CMA syscalls defined */
#define OPAL_CMA_NEED_SYSCALL_DEFS 0

/* Whether we have CUDA GDR support available */
#define OPAL_CUDA_GDR_SUPPORT 0

/* Whether we have CUDA cuPointerGetAttributes function available */
#define OPAL_CUDA_GET_ATTRIBUTES 

/* Whether we want cuda device pointer support */
#define OPAL_CUDA_SUPPORT 0

/* Whether we have CUDA CU_POINTER_ATTRIBUTE_SYNC_MEMOPS support available */
#define OPAL_CUDA_SYNC_MEMOPS 

/* Whether we have CU_MEM_LOCATION_TYPE_HOST_NUMA support available */
#define OPAL_CUDA_VMM_SUPPORT 0

/* Whether C compiler supports GCC style inline assembly */
#define OPAL_C_GCC_INLINE_ASSEMBLY 1

/* Whether C compiler support atomic convenience variables in stdatomic.h */
#define OPAL_C_HAVE_ATOMIC_CONV_VAR 1

/* Whether C compiler supports atomic operations on _Atomic variables without
   warnings */
#define OPAL_C_HAVE_ATOMIC_SUPPORT_FOR__ATOMIC 1

/* Whether C compiler supports __builtin_clz */
#define OPAL_C_HAVE_BUILTIN_CLZ 1

/* Whether C compiler supports __builtin_expect */
#define OPAL_C_HAVE_BUILTIN_EXPECT 1

/* Whether C compiler supports __builtin_prefetch */
#define OPAL_C_HAVE_BUILTIN_PREFETCH 1

/* Whether C compiler supports symbol visibility or not */
#define OPAL_C_HAVE_VISIBILITY 1

/* Whether C compiler supports __Atomic keyword */
#define OPAL_C_HAVE__ATOMIC 1

/* Whether C compiler supports __Generic keyword */
#define OPAL_C_HAVE__GENERIC 1

/* Whether C compiler support _Static_assert keyword */
#define OPAL_C_HAVE__STATIC_ASSERT 1

/* Whether C compiler supports __Thread_local */
#define OPAL_C_HAVE__THREAD_LOCAL 1

/* Whether C compiler supports __thread */
#define OPAL_C_HAVE___THREAD 1

/* Whether we have lt_dladvise or not */
#define OPAL_DL_LIBLTDL_HAVE_LT_DLADVISE 0

/* Whether we want developer-level debugging code or not */
#define OPAL_ENABLE_DEBUG 0

/* Enable fault tolerance general components and logic */
#define OPAL_ENABLE_FT 1

/* Enable fault tolerance MPI ULFM components and logic */
#define OPAL_ENABLE_FT_MPI 1

/* Disable getpwuid support (default: enabled) */
#define OPAL_ENABLE_GETPWUID 1

/* Enable features required for heterogeneous support */
#define OPAL_ENABLE_HETEROGENEOUS_SUPPORT 0

/* Enable IPv6 support, but only if the underlying system supports it */
#define OPAL_ENABLE_IPV6 0

/* Whether we want the memory profiling or not */
#define OPAL_ENABLE_MEM_DEBUG 0

/* Whether we want the memory profiling or not */
#define OPAL_ENABLE_MEM_PROFILE 0

/* Whether we want BTL progress threads enabled */
#define OPAL_ENABLE_PROGRESS_THREADS 0

/* Whether user wants PTY support or not */
#define OPAL_ENABLE_PTY_SUPPORT 1

/* Whether we want developer-level timing framework or not */
#define OPAL_ENABLE_TIMING 0

/* Greek - alpha, beta, etc - release number of Open Portable Access Layer */
#define OPAL_GREEK_VERSION ""

/* Whether your compiler has __attribute__ or not */
#define OPAL_HAVE_ATTRIBUTE 1

/* Whether your compiler has __attribute__ aligned or not */
#define OPAL_HAVE_ATTRIBUTE_ALIGNED 1

/* Whether your compiler has __attribute__ always_inline or not */
#define OPAL_HAVE_ATTRIBUTE_ALWAYS_INLINE 1

/* Whether your compiler has __attribute__ cold or not */
#define OPAL_HAVE_ATTRIBUTE_COLD 1

/* Whether your compiler has __attribute__ const or not */
#define OPAL_HAVE_ATTRIBUTE_CONST 1

/* Whether your compiler has __attribute__ constructor or not */
#define OPAL_HAVE_ATTRIBUTE_CONSTRUCTOR 1

/* Whether your compiler has __attribute__ deprecated or not */
#define OPAL_HAVE_ATTRIBUTE_DEPRECATED 1

/* Whether your compiler has __attribute__ deprecated with optional argument
   */
#define OPAL_HAVE_ATTRIBUTE_DEPRECATED_ARGUMENT 1

/* Whether your compiler has __attribute__ destructor or not */
#define OPAL_HAVE_ATTRIBUTE_DESTRUCTOR 1

/* Whether your compiler has __attribute__ error or not */
#define OPAL_HAVE_ATTRIBUTE_ERROR 1

/* Whether your compiler has __attribute__ extension or not */
#define OPAL_HAVE_ATTRIBUTE_EXTENSION 1

/* Whether your compiler has __attribute__ format or not */
#define OPAL_HAVE_ATTRIBUTE_FORMAT 1

/* Whether your compiler has __attribute__ format and it works on function
   pointers */
#define OPAL_HAVE_ATTRIBUTE_FORMAT_FUNCPTR 1

/* Whether your compiler has __attribute__ hot or not */
#define OPAL_HAVE_ATTRIBUTE_HOT 1

/* Whether your compiler has __attribute__ malloc or not */
#define OPAL_HAVE_ATTRIBUTE_MALLOC 1

/* Whether your compiler has __attribute__ may_alias or not */
#define OPAL_HAVE_ATTRIBUTE_MAY_ALIAS 1

/* Whether your compiler has __attribute__ noinline or not */
#define OPAL_HAVE_ATTRIBUTE_NOINLINE 1

/* Whether your compiler has __attribute__ nonnull or not */
#define OPAL_HAVE_ATTRIBUTE_NONNULL 1

/* Whether your compiler has __attribute__ noreturn or not */
#define OPAL_HAVE_ATTRIBUTE_NORETURN 1

/* Whether your compiler has __attribute__ noreturn and it works on function
   pointers */
#define OPAL_HAVE_ATTRIBUTE_NORETURN_FUNCPTR 1

/* Whether your compiler has __attribute__ no_instrument_function or not */
#define OPAL_HAVE_ATTRIBUTE_NO_INSTRUMENT_FUNCTION 1

/* Whether your compiler has __attribute__ optnone or not */
#define OPAL_HAVE_ATTRIBUTE_OPTNONE 0

/* Whether your compiler has __attribute__ packed or not */
#define OPAL_HAVE_ATTRIBUTE_PACKED 1

/* Whether your compiler has __attribute__ pure or not */
#define OPAL_HAVE_ATTRIBUTE_PURE 1

/* Whether your compiler has __attribute__ sentinel or not */
#define OPAL_HAVE_ATTRIBUTE_SENTINEL 1

/* Whether your compiler has __attribute__ unused or not */
#define OPAL_HAVE_ATTRIBUTE_UNUSED 1

/* Whether your compiler has __attribute__ visibility or not */
#define OPAL_HAVE_ATTRIBUTE_VISIBILITY 1

/* Whether your compiler has __attribute__ warn unused result or not */
#define OPAL_HAVE_ATTRIBUTE_WARN_UNUSED_RESULT 1

/* Whether your compiler has __attribute__ weak alias or not */
#define OPAL_HAVE_ATTRIBUTE_WEAK_ALIAS 1

/* whether backtrace_execinfo is found and available */
#define OPAL_HAVE_BACKTRACE_EXECINFO 1

/* whether qsort is broken or not */
#define OPAL_HAVE_BROKEN_QSORT 0

/* Whether C11 atomic compare swap is both supported and lock-free on 128-bit
   values */
/* #undef OPAL_HAVE_C11_CSWAP_INT128 */

/* whether ceil is found and available */
#define OPAL_HAVE_CEIL 1

/* whether clock_gettime is found and available */
#define OPAL_HAVE_CLOCK_GETTIME 1

/* Whether the processor supports the cmpxchg16b instruction */
/* #undef OPAL_HAVE_CMPXCHG16B */

/* whether dirname is found and available */
#define OPAL_HAVE_DIRNAME 1

/* Whether the OPAL DL framework is functional or not */
#define OPAL_HAVE_DL_SUPPORT 1

/* whether fbtl_posix is found and available */
#define OPAL_HAVE_FBTL_POSIX 1

/* Whether the __atomic builtin atomic compare swap is both supported and
   lock-free on 128-bit values */
#define OPAL_HAVE_GCC_BUILTIN_CSWAP_INT128 0

/* whether gethostbyname is found and available */
#define OPAL_HAVE_GETHOSTBYNAME 1

/* Do not use outside of mpi.h. Define to 1 if the system has the type `long
   long'. */
#define OPAL_HAVE_LONG_LONG 1

/* whether openpty is found and available */
#define OPAL_HAVE_OPENPTY 1

/* If PTHREADS implementation supports PTHREAD_MUTEX_ERRORCHECK */
#define OPAL_HAVE_PTHREAD_MUTEX_ERRORCHECK 1

/* If PTHREADS implementation supports PTHREAD_MUTEX_ERRORCHECK_NP */
#define OPAL_HAVE_PTHREAD_MUTEX_ERRORCHECK_NP 1

/* Whether we have SA_RESTART in <signal.h> or not */
#define OPAL_HAVE_SA_RESTART 1

/* whether sched_yield is found and available */
#define OPAL_HAVE_SCHED_YIELD 1

/* whether sharedfp_sm is found and available */
#define OPAL_HAVE_SHAREDFP_SM 1

/* whether shmem_posix is found and available */
#define OPAL_HAVE_SHMEM_POSIX 1

/* whether socket is found and available */
#define OPAL_HAVE_SOCKET 1

/* Whether or not we have solaris */
#define OPAL_HAVE_SOLARIS 0

/* Whether the __sync builtin atomic compare and swap supports 128-bit values
   */
#define OPAL_HAVE_SYNC_BUILTIN_CSWAP_INT128 1

/* Do not use outside of mpi.h. Define to 1 if you have the <sys/synch.h>
   header file. */
/* #undef OPAL_HAVE_SYS_SYNCH_H */

/* Do not use outside of mpi.h. Define to 1 if you have the <sys/time.h>
   header file. */
#define OPAL_HAVE_SYS_TIME_H 1

/* Whether we have __va_copy or not */
#define OPAL_HAVE_UNDERSCORE_VA_COPY 1

/* Whether we have va_copy or not */
#define OPAL_HAVE_VA_COPY 1

/* Whether we have weak symbols or not */
#define OPAL_HAVE_WEAK_SYMBOLS 1

/* Define to 1 ifyou have the declaration of _SC_NPROCESSORS_ONLN, and to 0
   otherwise */
#define OPAL_HAVE__SC_NPROCESSORS_ONLN 1

/* ident string for Open MPI */
#define OPAL_IDENT_STRING "5.0.8"

/* Major release number of Open Portable Access Layer */
#define OPAL_MAJOR_VERSION 5

/* Maximum length of datarep strings (default is 128) */
#define OPAL_MAX_DATAREP_STRING 128

/* Maximum length of error strings (default is 256) */
#define OPAL_MAX_ERROR_STRING 256

/* Maximum length of info keys (default is 36) */
#define OPAL_MAX_INFO_KEY 36

/* Maximum length of info vals (default is 256) */
#define OPAL_MAX_INFO_VAL 256

/* Maximum length of object names (default is 64) */
#define OPAL_MAX_OBJECT_NAME 64

/* Maximum length of port names (default is 1024) */
#define OPAL_MAX_PORT_NAME 1024

/* Maximum length of processor names (default is 256) */
#define OPAL_MAX_PROCESSOR_NAME 256

/* Maximum length of pset name lens (default is 512) */
#define OPAL_MAX_PSET_NAME_LEN 512

/* Maximum length of stringtag lens (default is 1024) */
#define OPAL_MAX_STRINGTAG_LEN 1024

/* MCA cmd line identifier */
#define OPAL_MCA_CMD_LINE_ID "mca"

/* MCA prefix string for envars */
#define OPAL_MCA_PREFIX "OMPI_MCA_"

/* Whether any opal memory mca components were found */
#define OPAL_MEMORY_HAVE_COMPONENT 0

/* Minor release number of Open Portable Access Layer */
#define OPAL_MINOR_VERSION 0

/* check if FI_HMEM_ROCR avaiable in fi_hmem_iface */
/* #undef OPAL_OFI_HAVE_FI_HMEM_ROCR */

/* check if iface avaiable in fi_mr_attr */
/* #undef OPAL_OFI_HAVE_FI_MR_IFACE */

/* check if pci data is available in ofi */
/* #undef OPAL_OFI_PCI_DATA_AVAILABLE */

/* package/branding string for Open MPI */
#define OPAL_PACKAGE_STRING "Open MPI proceduraltree@fedora Distribution"

/* Log base 2 of the maximum size in bytes of a memory descriptor. Set to 0 if
   MD can bind all of memory. */
#define OPAL_PORTALS4_MAX_MD_SIZE 0

/* Log base 2 of the maximum size in bytes of the user virtual address space.
   Set to 0 if MD can bind all of memory. */
#define OPAL_PORTALS4_MAX_VA_SIZE 0

/* Release date of Open Portable Access Layer */
#define OPAL_RELEASE_DATE "May 30, 2025"

/* Release release number of Open Portable Access Layer */
#define OPAL_RELEASE_VERSION 8

/* The repository version Open Portable Access Layer */
#define OPAL_REPO_REV "v5.0.8"

/* Disable ROCm support */
#define OPAL_ROCM_SUPPORT 0

/* Whether we have shared memory support for mmap or not */
#define OPAL_SHMEM_MMAP 1

/* Whether we have shared memory support for POSIX or not */
#define OPAL_SHMEM_POSIX 1

/* Whether we have shared memory support for SYSV or not */
#define OPAL_SHMEM_SYSV 1

/* Whether $opal_short_float_type needs __extension__ keyword */
#define OPAL_SHORT_FLOAT_NEEDS_EXTENSION 0

/* User-selected alternate C type of short float (used to redefine
   opal_short_float_t in opal_bottom.h) */
#define OPAL_SHORT_FLOAT_TYPE _Float16

/* Default value for mca_base_component_show_load_errors MCA variable */
#define OPAL_SHOW_LOAD_ERRORS_DEFAULT "all"

/* Do not use outside of mpi.h. Define to 1 if you have the ANSI C header
   files. */
#define OPAL_STDC_HEADERS 1

/* Tarball filename version string of Open Portable Access Layer */
#define OPAL_TARBALL_VERSION "5.0.8"

/* Whether to use assembly-coded atomics for atomics implementation */
#define OPAL_USE_ASM_ATOMICS 0

/* Whether to use C11 atomics for atomics implementation */
#define OPAL_USE_C11_ATOMICS 0

/* Whether to use GCC-style built-in atomics for atomics implementation */
#define OPAL_USE_GCC_BUILTIN_ATOMICS 1

/* Whether or not we are using the internal PMIx */
#define OPAL_USING_INTERNAL_PMIX 1

/* Complete release number of Open Portable Access Layer */
#define OPAL_VERSION "5.0.8"

/* Enable per-user config files */
#define OPAL_WANT_HOME_CONFIG_FILES 1

/* if the memory and buffer checking should be enabled */
#define OPAL_WANT_MEMCHECKER 0

/* if want pretty-print stack trace feature */
#define OPAL_WANT_PRETTY_PRINT_STACKTRACE 1

/* Greek - alpha, beta, etc - release number of Open SHMEM */
#define OSHMEM_GREEK_VERSION ""

/* Major release number of Open SHMEM */
#define OSHMEM_MAJOR_VERSION 5

/* Minor release number of Open SHMEM */
#define OSHMEM_MINOR_VERSION 0

/* Whether we want to check OSHMEM parameters always or never */
#define OSHMEM_PARAM_CHECK 1

/* Release date of Open SHMEM */
#define OSHMEM_RELEASE_DATE "May 30, 2025"

/* Release release number of Open SHMEM */
#define OSHMEM_RELEASE_VERSION 8

/* The repository version Open SHMEM */
#define OSHMEM_REPO_REV "v5.0.8"

/* Whether user wants OSHMEM in compatibility mode or not */
#define OSHMEM_SPEC_COMPAT 1

/* Whether we have shared memory support for mmap or not */
#define OSHMEM_SSHMEM_MMAP 1

/* Whether we have shared memory support for SYSV or not */
#define OSHMEM_SSHMEM_SYSV 1

/* Tarball filename version string of Open SHMEM */
#define OSHMEM_TARBALL_VERSION "5.0.8"

/* Complete release number of Open SHMEM */
#define OSHMEM_VERSION "5.0.8"

/* CFLAGS to pass through the wrapper compilers */
#define OSHMEM_WRAPPER_CFLAGS " "

/* CXXFLAGS to pass through the wrapper compilers */
#define OSHMEM_WRAPPER_CXXFLAGS " "

/* FCFLAGS to pass through the wrapper compilers */
#define OSHMEM_WRAPPER_FCFLAGS ""

/* LDFLAGS to pass through the wrapper compilers */
#define OSHMEM_WRAPPER_LDFLAGS "-L${libdir}  -Wl,-rpath -Wl,${libdir} -Wl,--enable-new-dtags"

/* LIBS to pass through the wrapper compilers */
#define OSHMEM_WRAPPER_LIBS "-loshmem -lmpi"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://www.open-mpi.org/community/help/"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Open MPI"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Open MPI 5.0.8"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "openmpi"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "5.0.8"

/* The size of `atomic_int', as computed by sizeof. */
#define SIZEOF_ATOMIC_INT 4

/* The size of `atomic_llong', as computed by sizeof. */
#define SIZEOF_ATOMIC_LLONG 8

/* The size of `atomic_long', as computed by sizeof. */
#define SIZEOF_ATOMIC_LONG 8

/* The size of `atomic_short', as computed by sizeof. */
#define SIZEOF_ATOMIC_SHORT 2

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `double _Complex', as computed by sizeof. */
#define SIZEOF_DOUBLE__COMPLEX 16

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `float _Complex', as computed by sizeof. */
#define SIZEOF_FLOAT__COMPLEX 8

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long double', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE 16

/* The size of `long double _Complex', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE__COMPLEX 32

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `opal_short_float_complex_t', as computed by sizeof. */
#define SIZEOF_OPAL_SHORT_FLOAT_COMPLEX_T 4

/* The size of `opal_short_float_t', as computed by sizeof. */
#define SIZEOF_OPAL_SHORT_FLOAT_T 2

/* The size of `pid_t', as computed by sizeof. */
#define SIZEOF_PID_T 4

/* The size of `ptrdiff_t', as computed by sizeof. */
#define SIZEOF_PTRDIFF_T 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `short float', as computed by sizeof. */
/* #undef SIZEOF_SHORT_FLOAT */

/* The size of `short float _Complex', as computed by sizeof. */
/* #undef SIZEOF_SHORT_FLOAT__COMPLEX */

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* The size of `ssize_t', as computed by sizeof. */
#define SIZEOF_SSIZE_T 8

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* The size of `wchar_t', as computed by sizeof. */
#define SIZEOF_WCHAR_T 4

/* The size of `_Bool', as computed by sizeof. */
#define SIZEOF__BOOL 1

/* The size of `__float128', as computed by sizeof. */
#define SIZEOF___FLOAT128 16

/* If the software-based performance counters capability should be enabled. */
#define SPC_ENABLE 0

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Check if float128 and float32(64,128)_complex dt are available in ucc. */
/* #undef UCC_HAVE_COMPLEX_AND_FLOAT128_DT */

/* Whether to use the legacy Solaris munmap prototype or not */
/* #undef USE_SOLARIS_LEGACY_MUNMAP_PROTOTYPE */

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

/* Whether the wrapper compilers add rpath flags by default */
#define WRAPPER_RPATH_SUPPORT "runpath"

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#define YYTEXT_POINTER 1

/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif

/* A bogus type that allows us to have sentinel type values that are still
   valid */
#define ompi_fortran_bogus_type_t int

/* C type corresponding to Fortran CHARACTER */
#define ompi_fortran_character_t char

/* C type corresponding to Fortran COMPLEX*16 */
/* #undef ompi_fortran_complex16_t */

/* C type corresponding to Fortran COMPLEX*32 */
/* #undef ompi_fortran_complex32_t */

/* C type corresponding to Fortran COMPLEX*4 */
/* #undef ompi_fortran_complex4_t */

/* C type corresponding to Fortran COMPLEX*8 */
/* #undef ompi_fortran_complex8_t */

/* C type corresponding to Fortran COMPLEX */
/* #undef ompi_fortran_complex_t */

/* C type corresponding to Fortran DOUBLE COMPLEX */
/* #undef ompi_fortran_double_complex_t */

/* C type corresponding to Fortran DOUBLE PRECISION */
#define ompi_fortran_double_precision_t double

/* C type corresponding to Fortran INTEGER*16 */
#define ompi_fortran_integer16_t 

/* C type corresponding to Fortran INTEGER*1 */
#define ompi_fortran_integer1_t char

/* C type corresponding to Fortran INTEGER*2 */
#define ompi_fortran_integer2_t short

/* C type corresponding to Fortran INTEGER*4 */
#define ompi_fortran_integer4_t int

/* C type corresponding to Fortran INTEGER*8 */
#define ompi_fortran_integer8_t long long

/* C type corresponding to Fortran INTEGER */
#define ompi_fortran_integer_t int

/* C type corresponding to Fortran LOGICAL*1 */
#define ompi_fortran_logical1_t char

/* C type corresponding to Fortran LOGICAL*2 */
#define ompi_fortran_logical2_t short

/* C type corresponding to Fortran LOGICAL*4 */
#define ompi_fortran_logical4_t int

/* C type corresponding to Fortran LOGICAL*8 */
#define ompi_fortran_logical8_t long long

/* C type corresponding to Fortran LOGICAL */
#define ompi_fortran_logical_t int

/* C type corresponding to Fortran REAL*16 */
#define ompi_fortran_real16_t long double

/* C type corresponding to Fortran REAL*2 */
#define ompi_fortran_real2_t ompi_fortran_bogus_type_t

/* C type corresponding to Fortran REAL*4 */
#define ompi_fortran_real4_t float

/* C type corresponding to Fortran REAL*8 */
#define ompi_fortran_real8_t double

/* C type corresponding to Fortran REAL */
#define ompi_fortran_real_t float

/* User-selected alternate C type of short float _Complex */
#define opal_short_float_complex_t _Float16 [2]

/* User-selected alternate C type of short float */
#define opal_short_float_t _Float16


#include "opal_config_bottom.h"
#endif /* OPAL_CONFIG_H */

