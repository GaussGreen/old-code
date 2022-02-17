#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*
#define COMPUTER_MAC2MS
*/
#if defined(COMPUTER_HP93C) || defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C)
#ifndef _INCLUDE_POSIX_SOURCE
#define _INCLUDE_POSIX_SOURCE
#endif
#endif

#include <errno.h>
#ifdef IMSL_MACHINE_SUN
#include <sys/file.h>
#include <sys/mman.h>
#endif

#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
#include <types.h>
#include <stat.h>
#include <file.h>
#elif defined COMPUTER_MAC2MS
#include <unix.h>
#include <Files.h>
#include <Types.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#ifdef IMSL_MACHINE_NT
#include <io.h>
#include <windows.h>
#else
#include <sys/types.h>
#include <unistd.h>
#endif

typedef unsigned short int   UINT16;

#define     MAX_FILES       3

#define     ID              29432
#define     ENV_NAME        "IMSLERRPATH"
#define     MAX_FILENAME    1024
#define     MAX_NAME        32

    /*  SEPARATOR symbol between directories in path */
    /*  END_OF_ENV an environment variable is ended by something in this set */
    /*  PATH is a list of directories searched for error message file */

#ifdef IMSL_MACHINE_DOS
							    /* DOS */
#define     SEPARATOR       ';'
#define     END_OF_ENV      "\\;"
#define     PATH            ".\\;\\imslclib;\\lib\\;c:\\imslclib;c:\\lib\\"
#endif

#ifdef IMSL_MACHINE_NT
							    /* NT */
#define     SEPARATOR       ';'
#define     END_OF_ENV      ";"
#define     PATH            ".;$LIB;/lib;c:/lib"
#endif

#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
							     /* VAX/VMS */
#define     SEPARATOR       ','
#define     END_OF_ENV
#define     PATH            "[],SYS$LOGIN:,SYS$SYSDEVICE:[IMSLCBASE]"
#endif
#ifdef COMPUTER_MAC2MS
							    /* Mac */
#define     SEPARATOR       ','
#define     END_OF_ENV      ";"
#define     PATH            "LXXX MB:Cmath:work:"
#endif

#ifndef     END_OF_ENV
							     /* UNIX */
#define     SEPARATOR       ':'
#define     END_OF_ENV      "/:"
#define     PATH            "./:$HOME/:/usr/lib/:/usr/local/lib/"
#endif

		/* names of the files containing error messages */
static Mchar        *lv_filename[MAX_FILES] = {
    "imslerr.bin", "egerr.bin", "XXXXXX.bin"
};

#ifndef O_BINARY
#define O_BINARY            0
#endif

static UINT16  *PROTO(l_malloc_16,(Mint n));
static void     PROTO(l_open_error_file,(Mlong code));
static Mint     PROTO(l_find,(Mlong code));
static Mint     PROTO(l_parse_path,(Mint (*proc)(Mchar*), Mchar* path));
static Mint     PROTO(l_look_for_file,(Mchar *filename));
static int      PROTO(l_compare,(UINT16 *p, UINT16*q));
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE) || defined(COMPUTER_MAC2MS)
static void     PROTO(l_swab,(Mchar *from, Mchar *to, Mint nbytes));
#endif

#if defined(COMPUTER_MAC2MS)
#define         READ(BUF,NINT,SIZE)     \
		    if (__read (lv_cur->unit,(void*)(BUF),sizeof(SIZE)*(NINT)) \
			!= sizeof(SIZE)*(NINT)) goto READ_ERROR
#elif defined(IMSL_MACHINE_NT)
#define         READ(BUF,NINT,SIZE)     \
		    if (_read (lv_cur->unit,(void*)(BUF),sizeof(SIZE)*(NINT)) \
			!= sizeof(SIZE)*(NINT)) goto READ_ERROR
#else
#define         READ(BUF,NINT,SIZE)     \
		    if (read (lv_cur->unit,(void*)(BUF),sizeof(SIZE)*(NINT)) \
			!= sizeof(SIZE)*(NINT)) goto READ_ERROR
#endif                  
			
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE) || defined(COMPUTER_MAC2MS)
#define         SWAP16(FROM,TO,NINTS) \
		  l_swab((Mchar*)(FROM),(Mchar*)(TO),sizeof(UINT16)*(NINTS))
#elif defined(IMSL_MACHINE_NT)
#define         SWAP16(FROM,TO,NINTS) \
		    _swab((Mchar*)(FROM),(Mchar*)(TO),sizeof(UINT16)*(NINTS))
#else
#define         SWAP16(FROM,TO,NINTS) \
		    swab((Mchar*)(FROM),(Mchar*)(TO),sizeof(UINT16)*(NINTS))
#endif

typedef struct {
    UINT16   id;
    UINT16   version;
    UINT16   n_messages;
    UINT16   text_size;
    UINT16   max_message;
    UINT16   max_name;
    UINT16   shift;
    UINT16   unused;
} Header;

typedef enum {
    open_never_tried,
    open_failed,
    open_ok
} Open_status;

typedef struct {
    Open_status     open_status;
    Mint            unit;
/*  OLD COMMENT: Need to fix this for Mac2ms */
/*
** I put this back. Apparently, it was removed because there 
** was a problem (originally) with this definition when
** MAC2MS was defined. However, it is needed for UNIX compilation.
*/

#if !defined(MAC2MS)
    struct stat     stat_buf;
#endif

    Header          header;
    UINT16          *map;
    UINT16          *start;
    Mlong           offset;
#ifdef IMSL_MACHINE_SUN
    Mchar           *data;
#endif
    Mint            swap;
} State;

static Mint         lv_file_num_max = -1;
static State        **lv_state = NULL;
static State        *lv_cur;
static Mchar        lv_special_message[64];
static Mint         lv_file_number;
static Mchar        lv_name[MAX_NAME];
static Mchar        *lv_read_error;


#ifdef ANSI
static UINT16   *l_malloc_16(Mint n)
#else
static UINT16   *l_malloc_16(n)
    Mint        n;
#endif
{
    UINT16      *space;
    space = (UINT16*)imsl_malloc(sizeof(UINT16)*n);
    if (space == NULL) lv_read_error = "Out-of-space";
    return space;
}

static Mint l_check_header()
{
    if (lv_cur->header.id == ID) {
	lv_cur->swap = 0;
	return 0;
    }
    if (lv_cur->header.id != ID) {
	UINT16  swap_id;
	SWAP16(&lv_cur->header.id, &swap_id, 1);
	if (swap_id == ID) {
	    Header      work;
	    SWAP16(&lv_cur->header, &work, sizeof(Header)/sizeof(UINT16));
	    lv_cur->header = work;
	    lv_cur->swap = 1;
	    return 0;
	}
    }
    lv_read_error = "File is not in imslerr.bin format";
    return 1;
}

#ifdef ANSI
static Mint l_parse_path(Mint (*proc)(Mchar*), Mchar *path)
#else
static Mint l_parse_path(proc, path)
    Mint        (*proc)();
    Mchar       *path;
#endif
{
    Mchar       filename[MAX_FILENAME];
    Mchar       *next;
    Mchar       *end;
    Mint        len;
#ifdef IMSL_MACHINE_NT
    Mchar       search_path[MAX_FILENAME];
    Mchar       *filepart;      /* a dummy variable */

    /*
     *  Use native Win32 services to search the directories for us.
     *  This also has the advantage of being able to have multiple
     *  directories in an environment variable and does not require
     *  any trailing directory delimiters.
     */
    search_path[0] = '\0';
    for (end = next = path;  end != NULL; ) {
	end = strchr(next, SEPARATOR);
	len = strlen(next);
	if (end != NULL)
	    len = end - next;
	else
	    len = strlen(next);
	if (*next == '$') {
	    Mchar   *end_env_name;
	    /*Mchar   *env;*/
	    Mchar   l_env[MAX_PATH];

	    if ((end_env_name = strchr(next, SEPARATOR)) != NULL)
		*end_env_name = '\0';
	    /* if ((env = getenv(next+1)) != NULL) */
	    if (GetEnvironmentVariable(next+1, l_env, MAX_PATH))
		strcat(search_path, l_env);
	    else {
		next += len + 1;
		continue;
	    }
	} else
	    strncat(search_path, next, len);
	if (end != NULL) {
	    strcat(search_path, END_OF_ENV);
	    next += len + 1;
	}
    }
    if (SearchPath(search_path, imsl_error_param.name, NULL, MAX_FILENAME,
		   filename, &filepart))
	if (proc(filename) == 0) return 0;
#else
    for (end = next = path;  end != NULL; ) {
	end = strchr(next, SEPARATOR);
	len = strlen(next);
	if (end != NULL)
	    len = end - next;
	else
	    len = strlen(next);
#if !defined(COMPUTER_VAX) && !defined(COMPUTER_ALFAC_IEEE)
	if (*next == '$') {
	    Mchar   *end_env_name;
	    Mchar   *env;
	    Mint    left;
	    Mint    len_env_name;

	    end_env_name = strpbrk(next, END_OF_ENV);
	    len_env_name = end_env_name - next - 1;
	    left         = end - next - len_env_name - 1;
	    strncpy(filename, next+1, len_env_name);
	    *(filename+len_env_name) = '\0';
	    env = getenv(filename);
	    len = strlen(env);
	    strcpy(filename, env);
	    strncpy(filename+len, end_env_name, left);
	    len  += left;
	    next += len_env_name + left + 2;
	} else
#endif
	{
	    strncpy(filename, next, len);
	    next += len + 1;
	}
	strcpy(filename+len, imsl_error_param.name);
	if (proc(filename) == 0) return 0;
    }
#endif
    lv_read_error = "File not found";
    return 1;
}

#ifdef ANSI
static int l_compare (UINT16 *p, UINT16*q)
#else
static int l_compare (p, q)
    UINT16   *p;
    UINT16   *q;
#endif
{
    if (*p <  *q) return -1;
    if (*p >  *q) return  1;
    return  0;
}

#ifdef ANSI
static Mint l_look_for_file(Mchar *filename)
#else
static Mint l_look_for_file(filename)
    Mchar       *filename;
#endif
{
/*
#if defined(COMPUTER_MAC2MS)
    OSErr status;
    FSSpec spec;
    FInfo fndrInfo;
    spec.name=filename;
    status = FSpGetFInfo(*spec,fndrInfo);
#endif
*/
    Mint        status=0;
/* Need to fix this for Mac2ms
    Mint        status;
    status = stat(filename, &lv_cur->stat_buf);
*/

#if defined(COMPUTER_MAC2MS)
    status=0;    
/*     if (status == 0) lv_cur->unit = __open(filename, 0); */
    if (status == 0) lv_cur->unit = __open(filename, O_RDONLY|O_BINARY);
#elif defined(IMSL_MACHINE_NT)
    if (status == 0) lv_cur->unit = _open(filename, O_RDONLY|O_BINARY, 0);
#else    /* UNIX flavors (I put this back in::Jim M.) */
    status = stat(filename, &lv_cur->stat_buf);
    if (status == 0) lv_cur->unit = open(filename, O_RDONLY|O_BINARY, 0);
#endif
    return status;
}

#ifdef ANSI
static void l_open_error_file(Mlong code)
#else
static void l_open_error_file(code)
    Mlong       code;
#endif
{
    Mint        n_messages;
    Mlong       modulo;
#if defined(COMPUTER_MAC2MS)
    FILE *f;
#endif
#ifdef IMSL_MACHINE_NT
    Mint        print_type;
    static Mchar        l_env[MAX_PATH];
#endif

    if (imsl_error_param.path == NULL)
#ifndef IMSL_MACHINE_NT
#if defined(COMPUTER_MAC2MS)
#else
	imsl_error_param.path = getenv(ENV_NAME);
#endif
#else
	/* Microslop doesn't like to use UNIX style functions.
	   Since they don't allow getenv()/putenv() to work
	   correctly with the corresponding Win32 functions,
	   we can't use getenv() here.  Maybe they'll fix this
	   deficiency in the future. */
	if (GetEnvironmentVariable(ENV_NAME, l_env, MAX_PATH))
	    imsl_error_param.path = l_env;
#endif
#if defined(COMPUTER_MAC2MS)
    if ((f=fopen("imslerr.bin","r"))==NULL) {
	if ((f=fopen("imslerr.bin alias","r"))==NULL) { 
/*
	    printf("Cannot find imslerr.bin or imslerr.bin alias in the working folder.\n");
	    printf("Please make a copy or make an alias of imslerr.bin in your working folder.\n");
*/
	    goto READ_ERROR;
	    }
	 fclose(f);
	 l_look_for_file("imslerr.bin alias");
	 }
    else {
       fclose(f);
       l_look_for_file("imslerr.bin");
    }      
/*    if (l_look_for_file("imslerr.bin"))
	    goto READ_ERROR;     */
#else
    if (imsl_error_param.path == NULL){
        /** must take copy as l_parse_path modifies supplied string - optimised
            code puts strings in read-only memory */
        imsl_error_param.path = (char*)malloc(strlen(PATH)+1);
        if (imsl_error_param.path == NULL){
            goto READ_ERROR;
        }
        strcpy(imsl_error_param.path, PATH);
    }
    if (l_parse_path(l_look_for_file, imsl_error_param.path))
	goto READ_ERROR;
#endif
#ifdef IMSL_MACHINE_SUN
    lv_cur->data = mmap(0, lv_cur->stat_buf.st_size, PROT_READ, MAP_PRIVATE,
		  lv_cur->unit, 0);
    close(lv_cur->unit);
    if ((int)lv_cur->data == -1) {
	lv_read_error = "Error in mmap";
	goto READ_ERROR;
    }
    memcpy((Mchar*)&lv_cur->header, lv_cur->data, sizeof(Header));
    if (l_check_header()) goto READ_ERROR;
    n_messages    = lv_cur->header.n_messages;
    lv_cur->map   = (UINT16*)(lv_cur->data+sizeof(Header));
    lv_cur->start = (UINT16*)(lv_cur->data+sizeof(Header)
					  +sizeof(UINT16)*n_messages);
    if (lv_cur->swap) {
	UINT16      *s_map;
	UINT16      *s_start;
	s_map   = l_malloc_16(n_messages);
	s_start = l_malloc_16(n_messages+1);
	if (s_map==NULL || s_start==NULL) goto READ_ERROR;
	SWAP16 (lv_cur->map,   s_map,   n_messages);
	SWAP16 (lv_cur->start, s_start, n_messages+1);
	lv_cur->map   = s_map;
	lv_cur->start = s_start;
    }
#else
    READ (&lv_cur->header, 1, Header);
    if (l_check_header()) goto READ_ERROR;
    n_messages    = lv_cur->header.n_messages;
    lv_cur->map   = l_malloc_16(n_messages);
    lv_cur->start = l_malloc_16(n_messages+1);
    if (lv_cur->map==NULL || lv_cur->start==NULL) goto READ_ERROR;
    if (lv_cur->swap) {
	UINT16  *work = l_malloc_16(n_messages+1);
	if (work==NULL) goto READ_ERROR;
	READ   (work, n_messages, UINT16);
	SWAP16 (work, lv_cur->map, n_messages);
	READ   (work, n_messages+1, UINT16);
	SWAP16 (work, lv_cur->start, n_messages+1);
	imsl_free(work);
    } else {
	READ (lv_cur->map,   n_messages,   UINT16);
	READ (lv_cur->start, n_messages+1, UINT16);
    }
#endif
    lv_cur->offset = sizeof(Header) + 2*sizeof(UINT16)*n_messages + 1;
    modulo = 1L << lv_cur->header.shift;
    if (lv_cur->offset & (modulo-1) != 0)
	lv_cur->offset += modulo - (lv_cur->offset & (modulo-1));
    lv_cur->open_status = open_ok;
    return;

    READ_ERROR:
#ifdef IMSL_MACHINE_NT
	imsl_error_options(IMSL_GET_PRINT_TYPE, &print_type, 0);
	if (print_type == IMSL_WINDOWS_PRINT_PROC) {
	    Mchar       process_name[MAX_PATH];
	    Mchar       *ptr;
	    Mchar       *msg = (Mchar*)malloc((strlen(imsl_error_param.name)+
					       strlen(lv_read_error) + 21) *
					      sizeof(Mchar));

	    sprintf(msg, "Error in reading %s\n%s.\n", imsl_error_param.name,
		    lv_read_error);

	    process_name[0] = '\0';
	    GetModuleFileName(NULL, process_name, MAX_PATH);
	    /*  Remove any path from the process name.  */
	    if ((ptr = strrchr(process_name, '\\')) == NULL)
		ptr = process_name;
	    else
		ptr++;

        /* 
	    MessageBox(NULL, msg, ptr, MB_ICONINFORMATION|MB_OK);
        */

	    free(msg);
	} else
#endif
#if defined(COMPUTER_MAC2MS)
/*
    fprintf(stderr,"Cannot find imslerr.bin or imslerr.bin alias in the working folder.\n");
    fprintf(stderr,"Please make a copy or make an alias of imslerr.bin in your working folder.\n");
*/
#else
    /*
    ** I commented this call to fprintf because, for Excel Addins, making a 
    ** to fprintf inside the spreadsheet could crash Excel. On Unix, if 
    ** fprintf is called inside of Kapital, it may cause problems. 
    **
    ** Jim Mitsiopoulos
    */ 
/*
	fprintf(stderr, "Error in reading %s\n%s.\n", imsl_error_param.name,
	    lv_read_error);
*/
#endif
	lv_cur->open_status = open_failed;
}

#ifdef ANSI
static Mint l_find(Mlong code)
#else
static Mint l_find(code)
    Mlong           code;
#endif
{
    UINT16          key;
    UINT16          *pindex;

    lv_file_number = code / 100000L;
    imsl_error_param.name = lv_filename[lv_file_number];

    if (lv_file_number > lv_file_num_max) {
	Mint    k;
	lv_state = (lv_state==NULL) ? (State**)malloc(sizeof(State))
	    : (State**) realloc((Mvoid*)lv_state,(lv_file_num_max+1)*sizeof(State));
	if (lv_state == NULL) return (-2);
	for (k=lv_file_num_max+1; k<=lv_file_number; k++)
	    lv_state[k] = NULL;
	lv_file_num_max = lv_file_number + 1;
    }
    if (lv_state[lv_file_number] == NULL) {
	lv_state[lv_file_number] = (State*)imsl_malloc(sizeof(State));
	if (lv_state[lv_file_number] == NULL) return (-2);
	lv_state[lv_file_number]->open_status = open_never_tried;
    }
    lv_cur = lv_state[lv_file_number];

    if (lv_cur->open_status == open_never_tried)
	l_open_error_file(code);
    if (lv_cur->open_status != open_ok)
	return (-1);

    key    = (UINT16)(code%100000L);
#if defined(COMPUTER_HP93C) || defined(COMPUTER_SUN4) || defined(COMPUTER_SUN5) || defined(COMPUTER_RTXLXS) || defined(COMPUTER_SGRUXS) || defined(COMPUTER_HP97C) || defined(COMPUTER_ALFAC_IEEE) || defined(IMSL_MACHINE_NT) || defined(COMPUTER_LINUX)
    pindex = (UINT16*)bsearch ((char*)&key, (char*)lv_cur->map,
	lv_cur->header.n_messages, sizeof(UINT16), (int (*)())l_compare);
#else
    pindex = (UINT16*)bsearch ((char*)&key, (char*)lv_cur->map,
	lv_cur->header.n_messages, sizeof(UINT16), l_compare);
#endif
    if (pindex == NULL)
	return (-1);
    return (pindex - lv_cur->map);
}

#ifdef ANSI
Mchar *imsl_find_message(Mlong code)
#else
Mchar *imsl_find_message(code)
    Mlong           code;
#endif
{
    static Mint     serial    = 0;
#ifndef IMSL_MACHINE_SUN
    static Mchar    *message  = NULL;
    Mlong           len;
#endif
    static Mlong    last_code = 0;
#ifdef COMPUTER_DECOSF
    unsigned int    off;
#else
    unsigned long   off;
#endif

    if (code != last_code)  serial = l_find(code);
    if (serial == -1) {
	sprintf(lv_special_message, "Error code %ld.", code);
	return lv_special_message;
    } else if (serial == -2) {
	lv_read_error = "Out-of-memory";
	goto READ_ERROR;
    }
    off = lv_cur->start[serial];
    off = off << lv_cur->header.shift;
#ifdef IMSL_MACHINE_SUN
    return (Mchar*)(lv_cur->data+lv_cur->offset+off);
#else
    if (message == NULL) {
	Mint    size = lv_cur->header.max_name + lv_cur->header.max_message + 2;
	message = (Mchar*)imsl_malloc(size);
	if (message == NULL) {
	    lv_read_error = "Out-of-memory";
	    goto READ_ERROR;
	}
    }
    if (code == last_code)  return message;

    off += lv_cur->offset;
#if defined(COMPUTER_MAC2MS)
    if (__lseek(lv_cur->unit, off, 0) == -1) {
#elif defined(IMSL_MACHINE_NT)
    if (_lseek(lv_cur->unit, off, 0) == -1) {
#else 
    if (lseek(lv_cur->unit, off, 0) == -1) {
#endif
	lv_read_error = "Cannot seek to correct location";
	goto READ_ERROR;
    }
    len = lv_cur->start[serial+1] - lv_cur->start[serial] + 1L;
    len = len << lv_cur->header.shift;
    READ(message, len, Mchar);
    last_code = code;
    return message;
#endif
    READ_ERROR:
	sprintf(lv_special_message,
	    "Error in reading %s for error message %d.\n%s.\n",
	    imsl_error_param.name, code, lv_read_error);
	lv_cur->open_status = open_failed;
	return lv_special_message;
}

#ifdef ANSI
Mchar *imsl_find_name(Mlong code)
#else
Mchar *imsl_find_name(code)
    Mlong           code;
#endif
{
    Mchar           *message;

    message = imsl_find_message(code);
    if (lv_cur->open_status == open_ok)
	strcpy(lv_name, message + strlen(message) + 1);
    else
	sprintf(lv_name, "%ld", code);
    return lv_name;
}

#ifdef ANSI
void imsl_ermes(Mint type, Mlong code)
#else
void imsl_ermes(type, code)
    Mint            type;
    Mlong           code;
#endif
{
    imsl_e1mes(type, code, imsl_find_message(code));
}

#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE) || defined(COMPUTER_MAC2MS)
#ifdef ANSI
static void l_swab(Mchar *from, Mchar *to, Mint nbytes)
#else
static void l_swab(from, to, nbytes)
	      Mchar     *from;
	      Mchar     *to;
	      Mint      nbytes;
#endif
{
   UINT16 *src = NULL, *dst = NULL;
   Mint i, nwords;
   nwords = (nbytes / 2);
   src = (UINT16 *)from;
   dst = (UINT16 *)to;
   for(i = 0; i < nwords; i++)  {
       *dst = ( (*src << 8) | (*src >> 8) );
       *src++;
       *dst++;
       }
}
#endif
