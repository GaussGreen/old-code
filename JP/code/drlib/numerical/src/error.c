#define ERROR_C
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef VLD
#include "imsl_vld.c"
#endif


static char lv_copyright[] =
    "Copyright 1990 by IMSL, Inc.  All Rights Reserved.";

#define MAX_STACK	200

static VA_LIST_HACK	    PROTO(l_error_options,(Mint code1, va_list argptr));
static Mint         PROTO(l_check_type,(Mint type, Mint first));
static void	    l_traceback();
static void         PROTO(l_error_print,(Imsl_error, Mlong, Mchar*, Mchar*));
#ifdef IMSL_MACHINE_NT
static void         PROTO(l_windows_error_print,(Imsl_error, Mlong, Mchar*, Mchar*));
static void         PROTO(l_console_error_print,(Imsl_error, Mlong, Mchar*, Mchar*));
#endif
static void         PROTO(l_parse_message,(Mchar *input));
static void	    l_alloc_message();
static void	    l_error_stop();
static void         l_error_init();

typedef struct {
    unsigned char   type;		/* error type number */
    unsigned char   user;		/* 1 if in user code */
    Mlong           code;		/* error code number */
    Mchar	    *routine;		/* routine name */
} stack_t;

			    /* Type 8 is for errors internal to the error
			       handler system. */
    /*   0,     1,     2,     3,     4,    5,      6,     7,    8  */
static Mchar lv_print_status[9] =
    {    0,     0,     0,     1,     1,    1,      1,     1,    1};
static Mchar lv_stop_status[9]  =
    {    0,     0,     0,     0,     1,    1,      0,     0,    1};
static Mchar lv_trace_status[9] =
    {    0,     0,     0,     0,     0,    0,      0,     0,    1};
static Mchar lv_immead_status[9]=
    {    0,     0,     0,     0,     0,    1,      1,     1,    1};

static Mint         lv_error_init = 0;
static Mlong        lv_error_n1rcd[2] = {0, 0};
static Mchar	    *lv_message = NULL;
static Mchar        *lv_end_message;
static Mint         lv_max_message_len = 256;
static stack_t	    lv_stack[MAX_STACK] = {{0, 1, 0, "USER"}};
static stack_t	    *lv_top = lv_stack;
static FILE	    *lv_error_file = NULL;
static Imsl_error_print_proc lv_error_print = l_error_print;
static Mint	    lv_error_print_type = IMSL_DEFAULT_PRINT_PROC;
static void	    PROTO((*lv_at_checksum),(Mint, Mint, Mchar*)) = NULL;
static Mint	    lv_full_traceback = 0;
static Mchar	    *lv_type_name[] = {
			"INTERNAL", "NOTE    ", "ALERT   ",
			"WARNING ", "FATAL   ", "TERMINAL",
			"WARNING_IMMEDIATE",  "FATAL_IMMEDIATE",
			"INTERNAL"
		    };

#ifdef ANSI
static void l_error_print(Imsl_error type, Mlong code, Mchar *routine,
                          Mchar *messagep)
#else
static void l_error_print(type, code, routine, messagep)
    Imsl_error	type;
    Mlong	code;
    Mchar	*routine;
    Mchar	*messagep;
#endif
{
    Mint	max_len = 78;
    Mint	line_len;
    Mint	stars_len;
    Mint	last_line;
    Mchar	header[256];
    Mchar	*stars = "***         ";
    Mchar	savech;
    Mchar	*split;
    Mchar	*end;
    Mchar	*newlinep;
    Mchar       *called_routine;
    stack_t     *p_stack;

    if (messagep[0] == '\0') return;
    	/* find the routine called just after USER mode */
    for (p_stack=lv_top;  !p_stack->user;  p_stack--);
    called_routine = (p_stack+1)->routine;

#if defined(COMPUTER_PMXUX) || defined(COMPUTER_HP93C) || defined(COMPUTER_DECFRS)
    if ((int)type < 5)
	sprintf(header, "\n*** %s Error %s from %s.  ", lv_type_name[(int)type],
	    imsl_find_name(code), called_routine);
    else
	sprintf(header, "\n*** %s Error from %s.  ", lv_type_name[(int)type],
	    called_routine);
#else
    if (type < 5)
	sprintf(header, "\n*** %s Error %s from %s.  ", lv_type_name[type],
	    imsl_find_name(code), called_routine);
    else
	sprintf(header, "\n*** %s Error from %s.  ", lv_type_name[type],
	    called_routine);
#endif
    stars_len = strlen(stars);
    line_len  = strlen(header);

    imsl_umach(3, &lv_error_file);
    fprintf(lv_error_file, header);
    if (line_len >= max_len) {
	fprintf(lv_error_file, "\n%s ", stars);
	line_len = stars_len;
    }
    while (1) {
	last_line = (Mint)strlen(messagep)+line_len < max_len;
	   /* Put a null at the end of the longest part left that will fit */
	if (!last_line) {
	    end    = messagep + (max_len - line_len);
	    savech = *end;
	    *end   = '\0';
	} else
	    savech = '\0';
	if (!last_line) {
		/* Look back for a blank */
	    split = strrchr(messagep, ' ');
	        /* Restore the character overwritten by the null */
	    if (savech != '\0') *end = savech;
		/* Print the part that fits */
	    if (split != NULL)
                *split = '\0';
            else {    /* NULL if not enough room for one word */
                fprintf(lv_error_file, "\n%s ", stars);
                line_len = stars_len;
                continue;
            }
	}
	    /* Check for newline in the message */
	if ((newlinep=strchr(messagep, '\n')) != NULL) {
	    *split = ' ';	/* Print up to newline instead of blank */
	    split  = newlinep;	/* pointed to by end */
	    *split = '\0';
	    last_line = 0;
	}
	fprintf(lv_error_file, "%s\n", messagep);
	if (last_line) break;
	if (split != NULL) *split = ' ';
	    /* Update message to point past where we have printed */
	messagep = split;
	    /* Print start of the next line */
	fprintf(lv_error_file, "%s", stars);
	line_len = stars_len;
    }
    fputc('\n', lv_error_file);
}

#ifdef IMSL_MACHINE_NT
#include <windows.h>
#ifdef ANSI
static void l_windows_error_print(Imsl_error type, Mlong code, Mchar *routine,
                          Mchar *messagep)
#else
static void l_windows_error_print(type, code, routine, messagep)
    Imsl_error	type;
    Mlong	code;
    Mchar	*routine;
    Mchar	*messagep;
#endif
{
    Mchar	type_name[32];
    Mchar	*error_text;
    Mint	error_length;
    Mchar       *called_routine;
    Mchar	process_name[MAX_PATH];
    Mchar	*ptr;
    stack_t     *p_stack;

    if (messagep[0] == '\0') return;
    	/* find the routine called just after USER mode */
    for (p_stack=lv_top;  !p_stack->user;  p_stack--);
    called_routine = (p_stack+1)->routine;

    strcpy(type_name, lv_type_name[type]);
    /*  Remove trailing blanks from the error type string.  */
    if ((ptr = strchr(type_name, ' ')) != NULL)
	*ptr = '\0';

    if (type < 5) {
	error_length = strlen(type_name) +
		       strlen(imsl_find_name(code)) +
		       strlen(called_routine) +
		       strlen(messagep) + 16;
	error_text = (Mchar*)malloc(error_length*sizeof(Mchar));
	sprintf(error_text, "%s Error %s from %s.  %s", type_name,
	    imsl_find_name(code), called_routine, messagep);
    } else {
	error_length = strlen(type_name) +
		       strlen(called_routine) +
		       strlen(messagep) + 16;
	error_text = (Mchar*)malloc(error_length*sizeof(Mchar));
	sprintf(error_text, "%s Error from %s.  %s", type_name,
	    called_routine, messagep);
    }

    process_name[0] = '\0';
    GetModuleFileName(NULL, process_name, MAX_PATH);
    /*  Remove any path from the process name.  */
    if ((ptr = strrchr(process_name, '\\')) == NULL)
	ptr = process_name;
    else
	ptr++;

    MessageBox(NULL, error_text, ptr, lv_stop_status[type] ?
	       MB_ICONSTOP|MB_OK|MB_SETFOREGROUND :
	       MB_ICONINFORMATION|MB_OK|MB_SETFOREGROUND);

    free(error_text);
}

#ifdef ANSI
static void l_console_error_print(Imsl_error type, Mlong code, Mchar *routine,
                          Mchar *messagep)
#else
static void l_console_error_print(type, code, routine, messagep)
    Imsl_error	type;
    Mlong	code;
    Mchar	*routine;
    Mchar	*messagep;
#endif
{
    Mint	max_len = 78;
    Mint	line_len;
    Mint	stars_len;
    Mint	last_line;
    Mchar	header[256];
    Mchar	*stars = "***         ";
    Mchar	savech;
    Mchar	*split;
    Mchar	*end;
    Mchar	*newlinep;
    Mchar       *called_routine;
    stack_t     *p_stack;
    HANDLE	console_handle = GetStdHandle(STD_OUTPUT_HANDLE);

    if (!console_handle) {
	/* If no console is available, create a default one. */
	imsl_create_console();
	console_handle = GetStdHandle(STD_OUTPUT_HANDLE);
    }

    if (messagep[0] == '\0') return;
    	/* find the routine called just after USER mode */
    for (p_stack=lv_top;  !p_stack->user;  p_stack--);
    called_routine = (p_stack+1)->routine;

    if (type < 5)
	sprintf(header, "\r\n*** %s Error %s from %s.  ", lv_type_name[type],
	    imsl_find_name(code), called_routine);
    else
	sprintf(header, "\r\n*** %s Error from %s.  ", lv_type_name[type],
	    called_routine);

    stars_len = strlen(stars);
    line_len  = strlen(header);

    WriteConsole(console_handle, header, strlen(header), NULL, NULL);
    if (line_len >= max_len) {
	sprintf(header, "\r\n%s ", stars);
	WriteConsole(console_handle, header, strlen(header), NULL, NULL);
	line_len = stars_len;
    }
    while (1) {
	last_line = (Mint)strlen(messagep)+line_len < max_len;
	   /* Put a null at the end of the longest part left that will fit */
	if (!last_line) {
	    end    = messagep + (max_len - line_len);
	    savech = *end;
	    *end   = '\0';
	} else
	    savech = '\0';
	if (!last_line) {
		/* Look back for a blank */
	    split = strrchr(messagep, ' ');
	        /* Restore the character overwritten by the null */
	    if (savech != '\0') *end = savech;
		/* Print the part that fits */
	    if (split != NULL)
                *split = '\0';
            else {    /* NULL if not enough room for one word */
                sprintf(header, "\r\n%s ", stars);
                WriteConsole(console_handle, header, strlen(header), NULL, NULL);
                line_len = stars_len;
                continue;
            }
	}
	    /* Check for newline in the message */
	if ((newlinep=strchr(messagep, '\n')) != NULL) {
	    *split = ' ';	/* Print up to newline instead of blank */
	    split  = newlinep;	/* pointed to by end */
	    *split = '\0';
	    last_line = 0;
	}
	sprintf(header, "%s\r\n", messagep);
	WriteConsole(console_handle, header, strlen(header), NULL, NULL);
	if (last_line) break;
	if (split != NULL) *split = ' ';
	    /* Update message to point past where we have printed */
	messagep = split;
	    /* Print start of the next line */
	sprintf(header, "%s", stars);
	WriteConsole(console_handle, header, strlen(header), NULL, NULL);
	line_len = stars_len;
    }
    sprintf(header, "\r\n");
    WriteConsole(console_handle, header, strlen(header), NULL, NULL);
}

#ifdef ANSI
void imsl_create_console(void)
#else
void imsl_create_console()
#endif
{
    FreeConsole();
    AllocConsole();
    /* Add any console attributes here */
}
#endif	/* IMSL_MACHINE_NT */

void imsl_e1sti(n, value)
    Mint	n;
    Mint	value;
{
    if (0<n && n<10) imsl_error_param.i[n] = value;
}

void imsl_e1str(n, value)
    Mint	n;
    Mdouble	value;
{
    if (0<n && n<10) imsl_error_param.d[n] = value;
}

void imsl_e1std(n, value)
    Mint	n;
    Mdouble	value;
{
    if (0<n && n<10) imsl_error_param.d[n] = value;
}

void imsl_e1stc(n, value)
    Mint	n;
    Mf_complex	value;
{
    Md_complex	z;

    z.re = (Mdouble)value.re;
    z.im = (Mdouble)value.im;
    if (0<n && n<10) imsl_error_param.z[n] = z;
}

void imsl_e1stz(n, value)
    Mint	n;
    Md_complex	value;
{
    if (0<n && n<10) imsl_error_param.z[n] = value;
}

void imsl_e1stl(n, value)
    Mint	n;
    Mchar	*value;
{
    n = abs(n);		/* n<0 means retain blanks. */
			/* HOW SHOULD n<0 BE HANDLED ???????? */
    if (0<n && n<10) imsl_error_param.l[n] = value;
}

#ifdef ANSI
static void l_parse_message(Mchar *input)
#else
static void l_parse_message(input)
    Mchar	*input;
#endif
{
    Mint	ind;
    Mchar	*mp;
    Mchar	*ip;
    Mchar	*np;
    Mchar	number[512];
    Mchar	format[128];
    Mint	len;

    for (mp=lv_message, ip=input;  *ip!='\0'; ) {
        if (mp+sizeof(number) >= lv_end_message) {
            len = mp - lv_message;
            l_alloc_message();
	    mp = lv_message + len;
        }

	if (*ip=='%' && *(ip+1)=='(' && isdigit(*(ip+3))) {
	    if (*(ip+4) == '%') { /* syntax %(R2%20.10f) */
		len = strchr(ip+4,')') - (ip+4);
		strncpy(format, ip+4, len);
		format[len] = 0;
	    } else		/* syntax %(R2) */
		len = 0;
	    ip += 2;
	    ind = *(ip+1) - '0';

	    switch (*ip) {
		case 'i':  case 'I':
   		    sprintf(number, (len==0) ? "%d" : format,
			    imsl_error_param.i[ind]);
		    break;
		case 'r':  case 'R':
		case 'd':  case 'D':
		case 'f':  case 'F':
		    sprintf(number, (len==0) ? "%e" : format,
			    imsl_error_param.d[ind]);
		    break;
		case 'c':  case 'C':
		case 'z':  case 'Z':
		    sprintf(number, (len==0) ? "(%e,%e)" : format,
			    imsl_error_param.z[ind].re,
                            imsl_error_param.z[ind].im);
		    break;
		case 'l':  case 'L':
		case 's':  case 'S':
		    sprintf(number, (len==0) ? "%s" : format,
			    imsl_error_param.l[ind]);
		    break;
		default:	/* Error in format */
		    sprintf(number, "%%(%c%c)", *ip, *(ip+1));
		    break;
	    }
	    for (np=number; *np!='\0'; ) *mp++ = *np++;
	    ip += len + 3;
	} else if (*ip=='%' && *(ip+1)=='/') {	/*  %/ means newline	*/
	    *mp++ = '\n';
	    ip += 2;
	} else
	    *mp++ = *ip++;
    }
    *mp = '\0';
}


#ifdef ANSI
void imsl_e1mes(Mint arg_type, Mlong arg_code, Mchar *arg_message)
#else
void imsl_e1mes(arg_type, arg_code, arg_message)
    Mint	arg_type;
    Mlong	arg_code;
    Mchar	*arg_message;
#endif
{
    if (arg_type<-1 || arg_type>8) {
	imsl_e1psh("imsl_e1mes");
	sprintf(lv_message,
	    "Error type must be -1,...,8, but type = %d.  Message = \"%s\"",
	    arg_type, arg_message);
	(*lv_error_print)(8, IMSL_BAD_ERROR_TYPE, lv_top->routine, lv_message);
        if (lv_trace_status[8]) l_traceback();
	l_error_stop();
    }
		/* type -1 and code -1 changes the message while retaining
		   the old type and code */
    if (arg_type == -1  &&  arg_code == -1) {
	arg_type = lv_top->type;
	arg_code = lv_top->code;
    }
		/* type -1 and code > 0 changes the code */
    if (arg_type == -1  &&  arg_code > 0) {
	lv_top->code = arg_code;
	goto RETURN;
    }
		/* type 0 and code 0 clears the error message */
    if (arg_type == 0  &&  arg_code == 0) {
	lv_message[0] = '\0';
	lv_top->type = lv_top->code = 0;
	goto RETURN;
    }

    lv_top->type = arg_type;
    lv_top->code = arg_code;

		/* For f77 compatiblity, treat the string " " as if it
		   were NULL */
    if (arg_message!=NULL && arg_message[0]==' ' && arg_message[1]=='\0')
	arg_message = NULL;
    if (arg_message != NULL)
	l_parse_message(arg_message);
    else
	lv_message[0] = '\0';
    if (lv_immead_status[lv_top->type]) {
	if (lv_at_checksum != NULL)
	     (*lv_at_checksum)(lv_top->type, lv_top->code, lv_message);
	if (lv_print_status[lv_top->type]) {
	    (*lv_error_print)((Mint)lv_top->type, (Mint)lv_top->code, lv_top->routine,
			      lv_message);
	    if (lv_trace_status[lv_top->type]) l_traceback();
	}
	lv_message[0] = '\0';
    }
		/* stop at once on internal errors */
    if (lv_top->type == 8) l_error_stop();

    RETURN:
	imsl_error_n1rty[0] = lv_top->type;
        lv_error_n1rcd[0] = lv_top->code;
}

static void l_alloc_message()
{
    Mint    dot_len = 3;
    if (lv_message == NULL) {
        lv_message = (Mchar*) malloc(lv_max_message_len*sizeof(Mchar));
            /* force these 3 chars before error message. */
            /* needed for compatability with f77 error checksum routine */
        if (lv_message != NULL) strcpy(lv_message, ".  ");
    }
    else {
        lv_max_message_len *= 2;
        lv_message = (Mchar*) realloc(lv_message-dot_len,
                                    lv_max_message_len*sizeof(Mchar));
    }
    if (lv_message == NULL) {
        lv_message = "Out of memory in error handler.";
        (*lv_error_print)(8, IMSL_OUT_OF_MEMORY, lv_top->routine, lv_message);
        if (lv_trace_status[8]) l_traceback();
        l_error_stop();
    }
    lv_message += dot_len;
    lv_end_message = lv_message + lv_max_message_len - dot_len;
}


static void l_error_init()
{
#ifdef HIGHLAND_LM
    char	feature[10];
#ifdef COMPUTER_DECOSF
    int 	t  = time(NULL);
#else
    long	t  = time(NULL);
#endif
    int		k  = t%100;

    feature[1] = 'M';
    feature[3] = 'T';
    feature[5] = 0;
    feature[2] = 'A';
    feature[4] = 'H';
    feature[0] = 'C';
		  /* This 'c' flags that imsl_highland_init */
		  /* should return an encoded value based on last 2 nulls, */
		  /* set based on the time of day */
    feature[9] = 0;
    feature[6] = 'c';
    feature[7] = (k/10)%10 + '0';
    feature[8] = k%10 + '0';
    if (imsl_highland_init(feature, CMATH_VERSION) != (13*k*k+8*k-6)) {
	lv_top = NULL;
	imsl_machine.f[2] = imsl_machine.f[3] = -1;
	imsl_machine.d[2] = imsl_machine.d[3] = -1;
    }
    feature[6] = 'x';
    imsl_highland_init(feature, CMATH_VERSION);
#endif
    imsl_signal();
#ifdef VLD
    imsl_cbase_vlds();
#endif
    l_alloc_message();
    imsl_error_param.path = NULL;
    imsl_error_param.name = NULL;
#ifdef VLD
    imsl_cbase_vlck();
#endif
    lv_error_init = 1;
}


#ifdef ANSI
void imsl_e1psh(Mchar *routine)
#else
void imsl_e1psh(routine)
    Mchar   *routine;
#endif
{
#ifdef WAVE
    if (lv_error_init) {
#ifdef ANSI
	extern void wave_product_name(char *, char *);
#else
	extern void wave_product_name();
#endif
	char prod_name[10];

	prod_name[0] = 'a';
	prod_name[9] = 0;
	prod_name[2] = 'v';
	prod_name[7] = 'g';
	prod_name[5] = 't';
	prod_name[3] = 'a';
	prod_name[6] = 'a';
	prod_name[4] = 'n';
	prod_name[1] = 'd';
	prod_name[8] = 'e';
	wave_product_check(prod_name, routine);
    }
#endif
    if (!lv_error_init) l_error_init();
    if (lv_top == lv_stack) {
       lv_message[0]   = '\0';
    }
    lv_top++;
#ifdef HIGHLAND_LM
    imsl_highland_check();
#endif
    if (lv_top == lv_stack+MAX_STACK) {
	imsl_ermes(8, IMSL_NESTED_TOO_DEEP);
	l_error_stop();
    }
    lv_top->type    = 0;
    lv_top->code    = 0;
    lv_top->user    = 0;
    lv_top->routine = routine;
    imsl_error_n1rty[0] = 0;
    imsl_error_n1rty[1] = 0;
    lv_error_n1rcd[0] = 0;
    lv_error_n1rcd[1] = 0;
}

#ifdef ANSI
void imsl_e1pop(Mchar *routine)
#else
void imsl_e1pop(routine)
    Mchar   *routine;
#endif
{
    if (imsl_error_param.signal_set) {  /* if signal then pop stack till we get a match */
        imsl_error_param.signal_set = 0;
	while (strcmp(routine,lv_top->routine) != 0) lv_top--;
    }
		/* Check that name on stack matches routine */
    if (routine!=lv_top->routine && strcmp(routine,lv_top->routine)!=0) {
	imsl_e1stl(1, routine);
	imsl_e1stl(2, lv_top->routine);
	imsl_ermes(8, IMSL_ERROR_STACK_MISMATCH);
    }
		/* Print & stop when popping back to USER level */
    if ((lv_top-1)->user) {
	if (lv_at_checksum != NULL)
	     (*lv_at_checksum)(lv_top->type, lv_top->code, lv_message);
	if (lv_message[0] != 0 && lv_print_status[lv_top->type]) {
	    (*lv_error_print)((Mint)lv_top->type, (Mint)lv_top->code,
		lv_top->routine, lv_message);
	    if (lv_trace_status[lv_top->type]) l_traceback();
	}
       if (lv_stop_status[lv_top->type])  l_error_stop();
    }
		/* Propogate back type and code */
    if ((lv_top-1)->type <= lv_top->type) {
	(lv_top-1)->type = lv_top->type;
	(lv_top-1)->code = lv_top->code;
	       /* Pop */
        imsl_error_n1rty[0] = lv_top->type;
        lv_error_n1rcd[0] = lv_top->code;
    }
    imsl_error_n1rty[1] = lv_top->type;
    lv_error_n1rcd[1]   = lv_top->code;
    lv_top--;
}

	/* User mode if flag is ON,  IMSL mode if flag is OFF */
#ifdef ANSI
void imsl_e1usr(Mchar *flag)
#else
void imsl_e1usr(flag)
    Mchar   *flag;
#endif
{
    lv_top->user = (flag[1]=='N') || (flag[1]=='n');
}

static void l_traceback()
{
    stack_t *p;
    stack_t *actual_error = NULL;
    Mchar   *m_format = "Here is a traceback of the calls in reverse order.";
    Mchar   *h_format = "  Error Type        Error Code              Routine";
    Mchar   *s_format = "  ----------        ----------              -------";
    Mchar   *l_format = " %5s%-9s   %-25s %s\n";

    imsl_umach(3,&lv_error_file);
    fprintf(lv_error_file, "%s\n%s\n%s\n", m_format, h_format, s_format);
    for (p=lv_top; p>=lv_stack; p--) {
	if (p->type != 0) actual_error = p;
	if (p==lv_stack || (p-1)->user || lv_full_traceback) {
	    if (actual_error == NULL || p->user) {
	        fprintf(lv_error_file, l_format, "", "", "", p->routine);
	    } else {
	        fprintf(lv_error_file, l_format, "IMSL_",
			lv_type_name[actual_error->type],
		        imsl_find_name(actual_error->code),
			p->routine);
	    }
	    actual_error = NULL;
	}
    }
}

static void l_error_stop()
{
  /*
  **  Removed call to exit so that CMATH will NOT exit 
  ** on internal errors.

    exit((Mint)lv_top->type);

  */
}


#ifdef ANSI
void imsl_e1pos(Mint type, Mint *print, Mint *stop)
#else
void imsl_e1pos(type, print, stop)
    Mint    type;
    Mint    *print;
    Mint    *stop;
#endif
{
    Mint    k;

    if (type<-8 || type>8) {
	imsl_e1psh("imsl_e1pos");
	imsl_e1sti(1, -8);
	imsl_e1sti(2, 8);
	imsl_e1sti(3, type);
	imsl_e1stl(1, "type");
	imsl_ermes(IMSL_TERMINAL, IMSL_INTEGER_OUT_OF_RANGE);
	imsl_e1pop("imsl_e1pos");
	return;
    }
    if (type == 0) {
	for (k=1; k<=7; k++) {
	    if (*print>=0) lv_print_status[k] = *print;
	    if (*stop>=0)  lv_stop_status[k] = *stop;
	}
    } else if (type>0) {
	if (*print>0) lv_print_status[type] = *print;
	if (*stop>0)  lv_stop_status[type] = *stop;
    } else {
	*print = lv_print_status[-type];
	*stop  = lv_stop_status[-type];
    }
}


#ifdef ANSI
Mint imsl_n1rty(Mint level)
#else
Mint imsl_n1rty(level)
    Mint    level;
#endif
{
    if (level<0 || level>1) {
	imsl_e1sti (1, 0);
	imsl_e1sti (2, 1);
	imsl_e1sti (3, level);
	imsl_e1stl (1, "level");
	imsl_ermes(8, IMSL_ILLEGAL_INTEGER_2);
    }
    return imsl_error_n1rty[level];
}


#ifdef ANSI
Mlong imsl_n1rcd(Mint level)
#else
Mlong imsl_n1rcd(level)
    Mint    level;
#endif
{
    if (level<0 || level>1) {
	imsl_e1sti (1, 0);
	imsl_e1sti (2, 1);
	imsl_e1sti (3, level);
	imsl_e1stl (1, "level");
	imsl_ermes (8, IMSL_ILLEGAL_INTEGER_2);
    }
    return lv_error_n1rcd[level];
}

#ifdef ANSI
void imsl_set_at_checksum(void(*checksum)(Mint, Mint, Mchar*))
#else
void imsl_set_at_checksum(checksum)
    void    (*checksum)();
#endif
{
    lv_at_checksum = checksum;
}


Imsl_error imsl_error_type()
{
    return lv_top->type;
}


Mlong imsl_error_code()
{
    return lv_error_n1rcd[1];
}


#ifdef ANSI
void imsl_error_options(Mint code1, ...)
#else
void imsl_error_options(code1, va_alist)
    Mint    code1;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, code1);
    imsl_e1psh("imsl_error_options");
    IMSL_CALL(l_error_options(code1, argptr));
    va_end(argptr);
    imsl_e1pop("imsl_error_options");
}


#ifdef ANSI
static VA_LIST_HACK l_error_options(Mint code1, va_list argptr)
#else
static VA_LIST_HACK l_error_options(code1, argptr)
    Mint        code1;
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number  = 0;
    Mint	    type;
    Mint            setting;
    Mint            *ans;
    FILE            *file;
    FILE            **pfile;
    Mchar           *path;
    Mchar           *name;
    Mint	    *return_print_type;
    Mint	    print_type;

    code = code1;
    while (code != 0) {
	if (arg_number > 0) code = va_arg(argptr, Mint);
	arg_number++;
	switch (code) {
	    case IMSL_SET_PRINT:
		arg_number += 2;
		type = va_arg(argptr, Mint);
		setting = va_arg(argptr, Mint);
                if (l_check_type(type,1)) lv_print_status[type] = setting;
		break;
	    case IMSL_SET_STOP:
		arg_number += 2;
		type = va_arg(argptr, Mint);
		setting = va_arg(argptr, Mint);
                if (l_check_type(type,1)) lv_stop_status[type] = setting;
		break;
	    case IMSL_SET_TRACEBACK:
		arg_number += 2;
		type = va_arg(argptr, Mint);
		setting = va_arg(argptr, Mint);
                if (l_check_type(type,1)) lv_trace_status[type] = setting;
		break;
	    case IMSL_GET_PRINT:
		arg_number += 2;
		type = va_arg(argptr, Mint);
		ans  = va_arg(argptr, Mint*);
                if (l_check_type(type,1)) *ans = lv_print_status[type];
		break;
	    case IMSL_GET_STOP:
		arg_number += 2;
		type = va_arg(argptr, Mint);
		ans  = va_arg(argptr, Mint*);
                if (l_check_type(type,1)) *ans = lv_stop_status[type];
		break;
	    case IMSL_GET_TRACEBACK:
		arg_number += 2;
		type = va_arg(argptr, Mint);
		ans  = va_arg(argptr, Mint*);
                if (l_check_type(type,1)) *ans = lv_trace_status[type];
		break;
            case IMSL_SET_ERROR_FILE:
                arg_number++;
                file = va_arg(argptr, FILE*);
                imsl_umach(-3, &file);
                break;
            case IMSL_GET_ERROR_FILE:
                arg_number++;
                pfile = va_arg(argptr, FILE**);
                imsl_umach(3, pfile);
                break;
            case IMSL_ERROR_PRINT_PROC:
                arg_number++;
                lv_error_print = va_arg(argptr, Imsl_error_print_proc);
		if (lv_error_print == NULL) lv_error_print = l_error_print;
                break;
	    case IMSL_GET_PRINT_TYPE:
                arg_number++;
                return_print_type = va_arg(argptr, Mint*);
                *return_print_type = lv_error_print_type;
		break;
#ifdef IMSL_MACHINE_NT
	    case IMSL_SET_PRINT_TYPE:
                arg_number++;
                print_type = va_arg(argptr, Mint);
                switch (print_type) {
                    case IMSL_WINDOWS_PRINT_PROC:

#if defined(NT_PRINT_CONSOLE)

			lv_error_print = l_console_error_print;
			lv_error_print_type = IMSL_CONSOLE_PRINT_PROC;
			break;

#else

			lv_error_print = l_windows_error_print;
			lv_error_print_type = IMSL_WINDOWS_PRINT_PROC;
			break;

#endif
                    case IMSL_CONSOLE_PRINT_PROC:
			lv_error_print = l_console_error_print;
			lv_error_print_type = IMSL_CONSOLE_PRINT_PROC;
			break;
                    default:
			lv_error_print = l_error_print;
                        lv_error_print_type = IMSL_DEFAULT_PRINT_PROC;
                        break;
                }
		break;
#endif	/* IMSL_MACHINE_NT */
            case IMSL_ERROR_MSG_PATH:
                arg_number++;
                path = va_arg(argptr, Mchar*);
                                /* copy path into new string */
                imsl_error_param.path = strcpy(malloc((int) strlen(path)+1),path);
                break;
            case IMSL_ERROR_MSG_NAME:
                arg_number++;
                name = va_arg(argptr, Mchar*);
                                /* copy path into new string */
                imsl_error_param.name = strcpy(malloc((int) strlen(name)+1),name);
                break;
            case IMSL_FULL_TRACEBACK:
		arg_number++;
                lv_full_traceback = va_arg(argptr, Mint);
                break;
	    case 0:
		break;
	    default:
		imsl_e1sti (1, code);
		imsl_e1sti (2, arg_number);
		imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
                goto RETURN;
		break;
	}
    }
    RETURN:
    return (argptr);
}

#ifdef ANSI
static Mint l_check_type(Mint type, Mint first)
#else
static Mint l_check_type(type, first)
    Mint    type;
    Mint    first;
#endif
{
    if (type<first || type>5) {
        imsl_e1sti(1, first);
        imsl_e1sti(2, 5);
        imsl_e1sti(3, type);
        imsl_e1stl(1, "type");
        imsl_ermes(IMSL_TERMINAL, IMSL_INTEGER_OUT_OF_RANGE);
        return 0;
    }
    return 1;
}
