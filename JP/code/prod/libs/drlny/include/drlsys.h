/************************************************************************
 * Module:	DRL - DR C Utilities
 * Submodule:	
 * File:	drlsys.h
 * Function:	Standard include for platform dependent systems calls
 *		MUST BE INCLUDED BEFORE ANY SYSTEM HEADER FILE
 * Author:	Christian Daher
 * Revision:	$Header$
 ************************************************************************/
#ifndef	_drlsys_H
#define	_drlsys_H

/* include for system calls */


#if defined(UNIX)
# ifndef _INCLUDE_POSIX_SOURCE
#  define _INCLUDE_POSIX_SOURCE
# endif
# ifndef _INCLUDE_XOPEN_SOURCE
#  define _INCLUDE_XOPEN_SOURCE
# endif
# ifndef _INCLUDE_HPUX_SOURCE
#  define _INCLUDE_HPUX_SOURCE
# endif


# define _XOPEN_SOURCE
    /* SunOS does not like it */
/*#   define _POSIX_SOURCE*/

# ifndef _ALL_SOURCE
#  define _ALL_SOURCE
# endif

# undef KERNEL
# include <errno.h>
# include <unistd.h>
# include <stdio.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <sys/times.h>
# include <sys/param.h>
# include <sys/utsname.h>
# include <sys/ioctl.h>
# include <sys/ipc.h>
# include <sys/shm.h>
# include <signal.h>
# include <fcntl.h>

# include <math.h>

#if !defined(SOLARIS) 
# ifndef CLOCKS_PER_SEC
#  define CLOCKS_PER_SEC  1000000L
# endif
#endif

# define Open	open
# define Close	close
# define Read	read
# define Write	write
# define Umask	umask
# define Creat	creat


/* tcp/ip */
# if defined(DRL_NET)
#  include <sys/socket.h>
#  include <netinet/in.h>
#  include <arpa/inet.h>
# endif /*DRL_NET*/



#elif defined(_WINDLL) || defined(WIN32) || defined(_WIN32)

# include <fcntl.h>
# include <sys\stat.h>
# include <sys\types.h>
# include <io.h>

# define Open	_open
# define Close	_close
# define Read	_read
# define Write	_write
# define Umask	_umask
# define Creat	_creat

# define O_APPEND	_O_APPEND
# define O_CREAT	_O_CREAT
# define O_EXCL		_O_EXCL
# define O_RDONLY	_O_RDONLY
# define O_WRONLY	_O_WRONLY
# define O_RDWR		_O_RDWR
# define O_TRUNC	_O_TRUNC


/* tcp/ip */
# if defined(DRL_NET)
#  include <winsock.h>
#  include <memory.h>
#  include <process.h>
# endif /*DRL_NET*/


#endif	/*UNIX,_WINDLL,etc.*/


/*
 * Other include
 */

# include <time.h>

#endif	/* _drlsys_H */
