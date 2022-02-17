/*
#define COMPUTER_MAC2MS
#define Mac2ms
*/ 
#ifdef _IBMR2
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE
#endif
#endif

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
#include <file.h>
#elif defined(COMPUTER_MAC2MS)
#include <unix.h>
#else
#if defined(COMPUTER_HP98C)
#ifndef _INCLUDE_HPUX_SOURCE
#define _INCLUDE_HPUX_SOURCE
#endif
#ifndef _INCLUDE_POSIX_SOURCE
#define _INCLUDE_POSIX_SOURCE
#endif
#endif

#if defined(COMPUTER_HP97C)
# ifndef _INCLUDE_POSIX_SOURCE
#  define _INCLUDE_POSIX_SOURCE
# endif
#endif

#ifdef WIN32
#include <io.h>
#include <stdlib.h>
#endif

#include <fcntl.h>
#endif

typedef unsigned short int   UINT16;
#ifdef COMPUTER_DECOSF
#define Mint            int
#else
typedef long int        Mint;
#endif

#ifdef COMPUTER_DECOSF
int         atol();
#else
long        atol();
#endif
void        *malloc();

#define     ID              29432
#define     VERSION             1

#define     MAX_LINE_LENGTH  1024
#define     MAX_MESSAGES     2000
#define     MAX_NAME_LEN       31

#define     MAX(A,B)        ((A)>(B)) ? (A) : (B)
#define     NEW_STRING(S)   strcpy(malloc(strlen(S)+1),S)

#ifndef O_BINARY
#define O_BINARY            0
#endif

struct {
    UINT16   id;
    UINT16   version;
    UINT16   n_messages;
    UINT16   text_size;
    UINT16   max_message;
    UINT16   max_name;
    UINT16   shift;
    UINT16   unused;
} header;

FILE        *data_file;
Mint        binary_file;
FILE        *include_file;
/* FILE     *log_file; */

Mint        shift = 0;
Mint        mask;
Mint        modulo;
Mint        skip;
Mint        n_messages = 0;
Mint        version = VERSION;
char        filename[128];
char        *lc_basename = "imsl";
char        *uc_basename = "IMSL";
char        *error_text = NULL;
Mint        last_error_number = 0;
UINT16      *map = NULL;
UINT16      *start;
char        name[2*MAX_NAME_LEN];
/*
char        *text;
unsigned char       *text;
*/
char        *text;
unsigned char       *text1;
char        *line = NULL;
char        *input_line = NULL;
Mint        total_size;
Mint        text_size;
Mint        max_message = 0;
Mint        max_name = 0;
Mint        len_message;
Mint        len_name;
char        *cp;
char        *cq;
Mint        k;
Mint        max_messages    = MAX_MESSAGES;
Mint        max_text        = 65536;
Mint        max_line_length = MAX_LINE_LENGTH;

static void l_alloc_memory()
{
    map   = (UINT16*)malloc(max_messages*sizeof(UINT16));
    start = (UINT16*)malloc(max_messages*sizeof(UINT16));
    error_text = text =  (char *) calloc(max_text, sizeof(char));
   /*  this is a CW4 compiler bug
    line  = (char *) realloc(line, max_line_length*sizeof(char));
    input_line = (char *) realloc(input_line, max_line_length*sizeof(char));
    
			printf("...in l_alloc_memory()... , line=%s\n",line); */
    if (map==NULL || start==NULL || text==NULL || line==NULL ||
	input_line==NULL) {
	fprintf(stderr, "Not enough memory.\n");
	exit(1);
    }
}

				/* trim trailing blanks from a line */
static void strip_trailing_blanks(line)
    char        *line;
{
    int         len = strlen(line);

    if (line[len-1] == '\n') len--;
    while (line[len-1] == ' ') len--;
    line[len] = '\n';
    line[len+1] = 0;
}


/*  Copy string  dest  to  src  and expand C escape sequences.  */
static void copy_and_expand(dest, src)
    char        *dest;
    char        *src;
{
    char        nextch;
    char        ch;
    int         k;

    strip_trailing_blanks(src);

    for (; *src != '\0'; src++)
	if (*src == '\\') {
	    switch (nextch = *++src) {
		case 'e':   ch = '\033';    break;
		case 'a':   ch = '\a';      break;
		case 'b':   ch = '\b';      break;
		case 'f':   ch = '\f';      break;
		case 'n':   ch = '\n';      break;
		case 'r':   ch = '\r';      break;
		case 't':   ch = '\t';      break;
		case 'v':   ch = '\v';      break;
		default:    ch = nextch;    break;
		case '0':       /* octal */
		    ch = 0;
		    for (k=0; k<3; k++) {
			nextch = *++src;
			if (nextch<'0' || nextch>'7') {--src; break;}
			ch = 8*ch + (nextch-'0');
		    }
		    break;
		case 'x':       /* hex */
		    ch = 0;
		    for (k=0; k<2; k++) {
			nextch = *++src;
			if (!isxdigit(nextch)) {--src; break;}
			if (isdigit(nextch))
			    ch = 16*ch + (nextch-'0');
			else if (isupper(nextch))
			    ch = 16*ch + (nextch-'A'+10);
			else
			    ch = 16*ch + (nextch-'a'+10);
		    }
		    break;
	    }
	    *dest++ = ch;
	}
	else
	    *dest++ = *src;
    *dest = 0;
}

static char *copy_to_uc(lc_basename)
    char        *lc_basename;
{
    char        *uc_basename;
    char        *ch;

    uc_basename = NEW_STRING(lc_basename);
    for (ch=uc_basename;  *ch != 0;  ch++)
	if ('a' <= *ch  &&  *ch <= 'z')
	    *ch += 'A' - 'a';
    return uc_basename;
}
main(argc, argv)
    int         argc;
    char        *argv[];
{
#if defined(COMPUTER_MAC2MS)   
    fpos_t pos;
#endif

/*
MaxApplZone();
MoreMasters();
MoreMasters();
MoreMasters();
MoreMasters();
MoreMasters();
MoreMasters();
MoreMasters();
MoreMasters();
MoreMasters();
*/
    if (argc == 2) {
	lc_basename = argv[1];
	uc_basename = copy_to_uc(lc_basename);
    }  
		    /* open  *err.dat  */
    sprintf(filename, "%serr.dat", lc_basename);
/*    printf("filename= %s.\n",filename); */
    data_file = fopen(filename, "r");
    if (data_file == NULL) {
	fprintf(stderr, "Cannot open input file %s.\n", filename);
	exit(1);
    }
		    /* open  *err.h  */
    sprintf(filename, "%serr.h", lc_basename);
    include_file = fopen (filename, "w");
    if (include_file == NULL) {
	fprintf (stderr, "Cannot open  %s  for writing.\n", filename);
	exit(1);
    }
/*    sprintf(filename, "%serr.log", lc_basename);
    log_file = fopen (filename, "w");
    if (log_file == NULL) {
	fprintf (stderr, "Cannot open  %s  for writing.\n", filename);
	exit(1);
    } */
    fprintf (include_file, "typedef enum {\n");

    modulo = 1L << shift;
    mask   = modulo - 1;
    line       = (char *) malloc(max_line_length*sizeof(char));
    input_line = (char *) malloc(max_line_length*sizeof(char));
#if defined(COMPUTER_MAC2MS)   
    fgetpos(data_file, &pos); 
#endif
/*            printf("1.............\n");
    gets(filename); */
/*    while (fgets(input_line, max_line_length, data_file) != NULL) { */
    while (fgets(input_line, max_line_length, data_file) != NULL) {    
    /* ????? */
/*   printf(" input line: %s\n",input_line); */
/*          fprintf (log_file, "%s\n", input_line); */
	copy_and_expand(line, input_line);
/*        printf("... line=%s\n",line); */
	if (line[0] == '#'){
	#if defined(COMPUTER_MAC2MS)   
    fgetpos(data_file, &pos); 
    #endif

	 continue;      /* comments are flagged by # */
	}       
	if (line[0] == '%') {
	    char        argument[128];
/*                            printf("0.............\n");
    gets(filename); */
	    if (error_text != NULL) {
		fprintf(stderr,
		    "Command must be at the top of the file.\n%s\n", line);
		exit(1);
	    }
	    /*    printf("...before  call sscanf(line+1..., line=%s\n",line); */
	    sscanf(line+1, "%s %s", name, argument);
	    /*    printf("...after call sscanf(line+1..., line=%s\n",line); */
	    if (strcmp(name,"shift") == 0) {
		shift  = atol(argument);
		modulo = 1L << shift;
		mask   = modulo - 1;
		max_text = 1L << (16+shift);
		continue;
	    } else if (strcmp(name,"max_text") == 0) {
		continue;
	    } else if (strcmp(name,"max_messages") == 0) {
		max_messages = atol(argument);
		continue;
	    } else if (strcmp(name,"max_line_length") == 0) {
		max_line_length = atol(argument);
		continue;
	    } else if (strcmp(name,"version") == 0) {
		version = atol(argument);
		continue;
	    }
	    fprintf(stderr, "Unknown command %s\n", line);
	    exit(1);
	}
		    /* allocate memory */
		    /*  printf("...before call l_alloc_memory()... , line=%s\n",line); */
	if (error_text == NULL) l_alloc_memory();
		    /*  printf("...after call l_alloc_memory()... , line=%s\n",line); */
	if (error_text - text > max_text - 100) {
	    fprintf(stderr, "Too much error message text.  ");
	    fprintf(stderr, "Use directive  \"%%shift %d\" in %serr.dat\n",
		    shift+1, lc_basename);
	    exit(1);
	}
	/* printf("...before call isspace(line[0]) , line=%s\n",line); */
	if (isspace(line[0])) { /* continuation line */
	 /*   printf("continuation line=%s\n",line); */
	    for (cp=line; *cp==' ' || *cp=='\t'; cp++);
	    error_text -= len_name + 2;
	    strcpy(error_text, cp+1);
	    error_text  += strlen(cp+1);
	    len_message += strlen(cp+1);
	    max_message  = MAX(max_message, len_message);
	    *(error_text-1) = (char) (NULL);
	    strcpy(error_text, name);
	    error_text += len_name + 1;
	  /*   printf("....error_text=%s\n",error_text); */
	} else {
	    Mint        off;

	    if (n_messages == max_messages) {
		fprintf(stderr, "Too many error messages.  ");
		fprintf(stderr, "Use directive  \"%%max_messages %d\" in %serr.dat\n",
		    2*max_messages, lc_basename);
		exit(1);
	    }
/*  printf("...before call atoi(line), line=%s\n",line); */
	    map[n_messages] = atoi(line);       /* get error msg number */
/*    printf("line=%s\n",line); */
/*      printf(".....map[n_messages]=map[%d]=%d\n",n_messages,map[n_messages]); */
	    cp = strpbrk(line, " \t");          /* skip over number */
	    while (*cp==' ' || *cp=='\t') cp++; /* skip to number */
	    cq = strpbrk(cp, " \t");            /* skip to end of name */
	    *cq = (char)NULL;
	    strcpy(name, cp);                   /* copy name */
	    cp = cq + 1;
	    while (*cp==' ' || *cp=='\t') cp++; /* skip to message */
	    off  = error_text - text;
	    skip = off & mask;
	    if (skip != 0) {
		off        += modulo - skip;
		error_text += modulo - skip;
	    }
	    off = off >> shift;
	    start[n_messages] = (UINT16) off;
	    strcpy(error_text, cp);
	    len_message = strlen(cp);
	    max_message = MAX(max_message, len_message);
	    error_text += len_message;
	    *(error_text-1) = (char)NULL;
	    len_name = strlen(name);
	    if (len_name > MAX_NAME_LEN) {
		fprintf(stderr, "Name %s is more than %d characters long.\n",
			name, MAX_NAME_LEN);
		fprintf(stderr, "      %*s=column %d\n", MAX_NAME_LEN, "|",
			MAX_NAME_LEN+1);
	    }
	    strcpy(error_text, name);       /* save name */
	    max_name = MAX(max_name,len_name);
	    fprintf (include_file, "\t%-32s =%7d,\n", name, map[n_messages]);
	    error_text += len_name + 1;
	    *(error_text-1) = (char)NULL;
	    if (last_error_number >= map[n_messages]) {
/*          printf("last_error_number=%d\n",last_error_number);
	    printf("map[n_messages]=map[%d]=%d\n",n_messages,map[n_messages]); */
		fprintf(stderr, "%s%s %d is\n%s\n",
		    "Error messages must be in strictly increasing order, ",
		    "but following message",
		    last_error_number, line);
		exit(1);
	    }
	    last_error_number = map[n_messages];
	    n_messages++;
	    #if defined(COMPUTER_MAC2MS)   
    fgetpos(data_file, &pos); 
	#endif


	  }
/*            printf("1.............\n");
    gets(filename); */
    }
/*        printf("2.............\n");
    gets(filename); */

    start[n_messages] = (UINT16)(error_text-text);

    header.id          = ID;
    header.version     = version;
    header.n_messages  = n_messages;
    text_size          = error_text - text;
    header.text_size   = text_size >> header.shift;
    header.max_name    = max_name;
    header.max_message = max_message;
    header.shift       = shift;
    total_size         = sizeof(header)
		       + sizeof(UINT16)*(2*header.n_messages+1)
		       + text_size;
    printf("%s%6d\n%s%6d\n%s%6d\n%s%6d\n%s%6d\n%s%6d\n",
	"Number of error messages:            ", header.n_messages,
	"Number of char's in error messages:  ", text_size,
	"Max name length                      ", header.max_name,
	"Max message length                   ", header.max_message,
	"Size of errors.bin:                  ", total_size,
	"Shift:                               ", header.shift);

    sprintf(filename, "%serr.bin", lc_basename);
#if defined(SOLARIS) || defined(SUN) || defined(HP_UX) || defined(AIX) || defined(LINUX)
    binary_file = open(filename, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY,
	(int)0666);
#elif defined(WIN32)
    binary_file = _open(filename, O_WRONLY | O_CREAT | O_TRUNC | O_BINARY,
	(int)0666);
#else
	 binary_file = __open(filename,
	(int)0666);
#endif
	
#if defined(SOLARIS) || defined(SUN) || defined(HP_UX) || defined(AIX) || defined(LINUX)
    write (binary_file, (char*)&header, sizeof(header));
    write (binary_file, (char*)map,     sizeof(UINT16)*header.n_messages);
    write (binary_file, (char*)start,   sizeof(UINT16)*(header.n_messages+1));
    write (binary_file, text,           text_size);
#elif defined(WIN32)
    _write (binary_file, (char*)&header, sizeof(header));
    _write (binary_file, (char*)map,     sizeof(UINT16)*header.n_messages);
    _write (binary_file, (char*)start,   sizeof(UINT16)*(header.n_messages+1));
    _write (binary_file, text,           text_size);
#else
    __write (binary_file, (unsigned char*)&header, sizeof(header));
    __write (binary_file, (unsigned char*)map,     sizeof(UINT16)*header.n_messages);
    __write (binary_file, (unsigned char*)start,   sizeof(UINT16)*(header.n_messages+1));
    text1=(unsigned char *)(text);
    __write (binary_file, text1,                text_size);
#endif


#if defined(SOLARIS) || defined(SUN) || defined(HP_UX) || defined(AIX) || defined(LINUX)
    close (binary_file);
#elif defined(WIN32) 
    _close (binary_file); 
#else
    __close (binary_file);
#endif
    
/*    printf("3.............\n");
    gets(filename); */
			/* finish  *err.h  */
    fprintf (include_file, "\t%s_%-32s\n} %c%s_code;\n", uc_basename,
	     "LAST_ERROR_MSG", uc_basename[0], lc_basename+1);
    fclose(include_file);
/*    fclose(log_file); */
    exit(0);
}


