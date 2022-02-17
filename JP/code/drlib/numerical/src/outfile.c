#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK PROTO(l_output_file,(Mint,va_list));

#ifdef ANSI
void imsl_output_file(Mint code1, ...)
#else
void imsl_output_file(code1, va_alist)
    Mint    code1;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, code1);

    imsl_e1psh("imsl_output_file");
    IMSL_CALL(l_output_file(code1, argptr));
    va_end(argptr);
    imsl_e1pop("imsl_output_file");
}


#ifdef ANSI
static VA_LIST_HACK l_output_file(Mint code1, va_list argptr)
#else
static VA_LIST_HACK l_output_file(code1, argptr)
    Mint        code1;
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number = 0;
    FILE            *file;
    FILE            **pfile;

    code = code1;
    while (code != 0) {
	if (arg_number > 0) code = va_arg(argptr, Mint);
	arg_number++;
	switch (code) {
             case IMSL_SET_OUTPUT_FILE:
                arg_number++;
                file = va_arg(argptr, FILE*);
                imsl_umach(-2, &file);
                break;
             case IMSL_SET_ERROR_FILE:
                arg_number++;
                file = va_arg(argptr, FILE*);
                imsl_umach(-3, &file);
                break;
             case IMSL_GET_OUTPUT_FILE:
                arg_number++;
                pfile = va_arg(argptr, FILE**);
                imsl_umach(2, pfile);
                break;
             case IMSL_GET_ERROR_FILE:
                arg_number++;
                pfile = va_arg(argptr, FILE**);
                imsl_umach(3, pfile);
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
