/* <nag.h>
 *
 * Copyright 1990 Numerical Algorithms Group
 *
 * General include file for all NAG C software.
 *
 * Malcolm Cohen, NAG Ltd., Oxford, U.K., 1990.
 *
 * Mark 1, 1990.
 *
 * Mark 2, revised 
 *    The defined types and error names have been placed into
 *    separate files, nag_types.h and  nag_errlist.h, respectively.
 *    SPD, October 91.
 *    Reworked handling of implementation specifics.  AJS feb92.
 *
 * Purpose : To define a macro that is implementation specific so that
 *           variant code may be easily inserted by implementors and carried
 *           foreward into future releases.
 *           A number of other macros are also defined to indicate the
 *           compliance (or otherwise) to particular common characteristics.
 *           These macros are always set to their default state in the
 *           introduction.  There is a later section where the implementation-
 *           specific macro controls changes to the defaults.
 * Mark 5 revised. IER-2134 (Feb 1998).
 */

#ifndef NAG_H
#define NAG_H

/* The type of system should be defined here.
   It is defined by the "system" and "compiler" components of the NAG product
   code for the library (c.f. a00aact.c).  Some examples are :
       APOA_NAG_IMP for Apollo SR 10.3
       IBPA_NAG_IMP for the microsoft PC compiler
       IBPC_NAG_IMP for the Zortech PC compiler
       IBPE_NAG_IMP for the Borland PC compiler
       SU3A_NAG_IMP for Sun 3

   These macros are used in preprocessor directives in other parts of the code
   to introduce the implementation specific variations.  It is only anticipated
   that variations will be required in nag.h , nag_stddef.h , nag_stdlib.h ,
   nag_string.h and x05aact.c.
   */
/*#define APOA_NAG_IMP  Apollo SR 10.3 */
/*#define D31A_NAG_IMP  DEC RISC Ultrix */
/*#define DVVA_NAG_IMP  DEC VAX/VMS */
/*#define DVVG_NAG_IMP  DEC VAX/VMS G_float option */
/*#define IBPA_NAG_IMP  MS/PC-DOS PCs , Microsoft compiler */
/*#define IBPA_DLL_NAG_IMP  MS-Windows DLL, Microsoft compiler */
/*#define IBPC_NAG_IMP  MS/PC-DOS PCs , Zortech compiler */
/*#define IBPE_NAG_IMP  MS/PC-DOS & MS/WINDOWS PCs , Borland compiler */
/*#define SU3A_NAG_IMP  Sun 3 bundled C compiler */
/*#define SU4A_NAG_IMP  Sun 4 (Sparc) */
/*#define SG4A_NAG_IMP  Silicon Graphics IRIX */
#define IBPA_STL32_NAG_IMP /* MS-Windows 32 bit Static Library, Microsoft compiler */
/* #define IBPA_DLL32_NAG_IMP  MS-Windows 32 bit DLL, Microsoft compiler */
/*#define IBPA_STL32_NAG_IMP MS-Windows 32 bit Static Library, Microsoft compiler */
/*#define DAUA_NAG_IMP  /* DEC Alpha, Digital UNIX */

/* The following 14 sections define the default settings for macros concerned
   with common characteristics that are shared by many compilers. These
   settings are listed below and a macro is defined for each 
   setting. If for a particular compiler, a setting is not appropriate,
   then it should be corrected later, in the implementation specific section
   of this file by undefining the corresponding macro. For example, we assume
   that most C compilers support function prototypes, hence NAG_PROTO is
   defined. If the compiler happens to be the bundled SUN C compiler which
   does not support function prototypes, then this macro is undefined, see below.
   */

/* NAG_PROTO 
   This is used to indicate that the compiler supports function prototypes
   as  specified in the ANSI standard. For the older K&R compilers which do
   not support function prototypes undefine NAG_PROTO */
#define NAG_PROTO

/* NAG_FABS
   This macro is defined on the assumption that the function fabs is available
   in the math library. If this is not so and the compiler does not offer
   anything better then use the macro ABS which is defined later in this file
   by undefining NAG_FABS.
   */
#define NAG_FABS

/* NAG_EXIT
   The default action for p01acc is to terminate the program by the use of
   exit.  On some systems it may be more appropriate to terminate with abort,
   perhaps causing a core dump.  However, any parent process (or script)
   may wish to continue after testing the exit status of the program.
   If it is approriate to use abort rather than exit then undefine NAG_EXIT.  */
#define NAG_EXIT

/* NAG_YES_STDLIB
   Most compilers provide the include file stdlib.h containining declarations
   for the functions abort , atoi , free , exit , atol , calloc , malloc ,
   realloc , alloca. If this is not the case then undefine stdlib.h
   */
#define NAG_YES_STDLIB

  /* NAG_YES_STRING 
   Most compilers provide the include file string.h, however some older BSD
   derived, non-ANSI compilers provide strings.h instead. If strings.h is 
   provided then undefine NAG_YES_STRING
   */
#define NAG_YES_STRING

  /*
    NAG_YES_STDDEF
    Most compilers provide the include file stddef.h containining definitions
    for ptrdiff_t and size_t. If these definitions are not provided then undefine
    NAG_YES_STDDEF.
    */
#define NAG_YES_STDDEF

/* NAG_VOID_STAR
   Newer ANSI based compilers provide the generic pointer void *. If the
   target system  does not allow the  type void *, then undefine NAG_VOID_STAR.
   */
#define NAG_VOID_STAR

/* NAG_IND_STR
  Newer ANSI based compilers provide the string functions strchr and strrchr.
  Some older compilers provide index instead of strchr and rindex instead of
  strrchr. If this is the case with the target system then undefine NAG_IND_STR.
  */
#define NAG_IND_STR

/* NAG_NO_INC_MEMORY_DOT_H
   The functions memcmp , memcpy are normally declared string(s).h. If on
   target system they are declared in memory.h then undefine
   NAG_NO_INC_MEMORY_DOT_H */
#define NAG_NO_INC_MEMORY_DOT_H

/* NAG_NO_BC_MEM_FUN
   The functions memcpy and memcmp are meant to manipulate objects as
   character arrays. Some BSD derived systems provide bcopy and bcmp
   as either additional to or instead of memcpy and memcmp. If it
   is felt that it is desirable to use bcopy and bcmp then undefine
   NAG_NO_BC_MEM_FUN
   */
#define NAG_NO_BC_MEM_FUN

/* NAG_ARRAY_SIZE
   In the stringent test program for c06fuc, a number of large arrays
   are being used. This might prove too large, especially for the PC
   type machines. In this case undefine NAG_ARRAY_SIZE.
   */
#define NAG_ARRAY_SIZE

/* NAG_TIME_DOT_H
   x05aac function is used to set up a non-repeatable seed for the
   random number generator function g05ccc. This is achieved by setting
   the seed to a value derived from the system time by calling an appropriate
   function.
   Ansi derived compilers call the function time with the argument of type
   time_t and include files time.h and sys/types.h.
   If this is not the case undefine NAG_TIME_DOT_H. System specific
   code may then have to be inserted in x05aact.c
*/
#define NAG_TIME_DOT_H

/* NAG_YES_MALLOC
   The C dynamic memory allocator, malloc returns size_t whose maximum
   value may not be large enough for practical applications, e.g the
   Microsoft C compiler whose size_t has the maximum value of 65536. In
   this case the compiler may provide an alternative function, such as
   halloc. See nag_stdlib.h for the definition of halloc.
*/
#define NAG_YES_MALLOC

/* NAG_LINT
   The Apollo version of lint is not ANSI compatible. For lint to work
   successfully, it seems necessary to assume a 'pure' bsd system.
*/

/* NAG_CHAR_BOOL
Most compilers seem to be happy to accept 
        typedef char Boolean;
If yours insists on defining 
        typedef signed char Boolean;
then
undef NAG_CHAR_BOOL.
If it is preferred to define Boolean as an int rather than
char, then undef NAG_CHAR_BOOL and define NAG_INT_BOOL.
*/
#define NAG_CHAR_BOOL

/* NAG_THREAD_SAFE
   If the compiler supports POSIX threads, then this should be defined to
   produce thread safe versions of the routines g05cac, g05cay and g05caz.
*/
/* #define NAG_THREAD_SAFE */

/* These are defined to nothing, they are used only for Windows code, 16 and 32 bit */
#define NAGW_S
#define NAG_FPE
#define NAG_FP
#define NAG_FAR
#define NAG_HUGE 
#define NAG_DLL_IMP
#define NAG_DLL_EXP
#define NAG_DLL_EXPIMP
#define NAG_CALL

/* The following sections override the default settings of the above 14
   macros, where it is necessary for a particular implementation.
*/

#ifdef IBPA_NAG_IMP
#undef NAG_EXIT
#undef NAG_ARRAY_SIZE
#undef NAG_YES_MALLOC
#undef NAG_CHAR_BOOL
#define NAG_INT_BOOL

#else
#ifdef IBPA_DLL_NAG_IMP
#define IBPA_NAG_IMP
#undef NAG_EXIT
#undef NAG_ARRAY_SIZE
#undef NAG_YES_MALLOC
#undef NAG_CHAR_BOOL
#define NAG_INT_BOOL
#include <windows.h>
#define NAG_WIN
#define NAG_DLL
#undef NAG_FPE
#define NAG_FPE FAR PASCAL _export
#undef NAG_FP
#define NAG_FP FAR PASCAL
#undef NAG_FAR
#define NAG_FAR FAR
#undef NAG_HUGE 
#define NAG_HUGE huge

#else
#ifdef IBPA_DLL32_NAG_IMP
#undef NAG_EXIT
#undef NAG_ARRAY_SIZE
#undef NAG_CHAR_BOOL
#define NAG_INT_BOOL
#define NAG_WIN
#define NAG_DLL
#undef NAG_DLL_EXPIMP
#undef NAG_CALL
#define NAG_CALL __stdcall
#define NAG_MICROSOFT_THREAD_SAFE 
#ifdef  NAG_DLL_EXPORT
#define NAG_DLL_EXPIMP __declspec(dllexport)
#else
#define NAG_DLL_EXPIMP __declspec(dllimport)
#endif

#else
#ifdef IBPA_STL32_NAG_IMP
#undef NAG_EXIT
#undef NAG_ARRAY_SIZE
#undef NAG_CHAR_BOOL
#define NAG_INT_BOOL
#define NAG_WIN
#define NAG_DLL
#undef NAG_CALL
#define NAG_CALL __stdcall


#else
#ifdef IBPE_NAG_IMP
#undef NAG_ARRAY_SIZE
#undef NAG_YES_MALLOC
extern unsigned  _floatconvert;
#pragma extref _floatconvert

#else
#ifdef SU3A_NAG_IMP
#undef NAG_PROTO
#undef NAG_EXIT
#undef NAG_NO_BC_MEM_FUN
#undef NAG_IND_STR

#else
#ifdef SU4A_NAG_IMP
#undef NAG_PROTO
#undef NAG_EXIT
#undef NAG_NO_BC_MEM_FUN
#undef NAG_IND_STR

#else
#ifdef DVVA_NAG_IMP
#undef NAG_TIME_DOT_H

#else
#ifdef NAG_LINT
#undef NAG_PROTO
#undef NAG_VOID_STAR
/*
#undef NAG_YES_STRING
#undef NAG_YES_STDLIB
#undef __STDC__
*/

#else
#ifdef SG4A_NAG_IMP
#undef NAG_CHAR_BOOL

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif


/* ****** End of anticipated implementation specifics ***** */


#ifndef NULLFN
#define NULLFN (void(*)())0
#endif
#ifndef NULLDFN
#define NULLDFN (double(*)())0
#endif

/*
#ifndef NULLFN2
#define NULLFN2 (void(NAG_FPE *)())0
#endif
#ifndef NULLDFN2
#define NULLDFN2 (double(NAG_FPE*)())0
#endif
*/

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define NAGERR_DEFAULT (NagError NAG_HUGE *)0
#define NAGUSER_DEFAULT (Nag_User NAG_HUGE *)0
#define NAGCOMM_NULL (Nag_Comm NAG_HUGE *)0
#define NAGMESG_DEFAULT (Nag_Mesg NAG_HUGE *)0
#define E01_DEFAULT (Nag_E01_Opt NAG_HUGE *)0
#define E04_DEFAULT (Nag_E04_Opt NAG_HUGE *)0
#define G13_DEFAULT (Nag_G13_Opt NAG_HUGE *)0
#define H02_DEFAULT (Nag_H02_Opt NAG_HUGE *)0
#define TRANSF_ORDER_DEFAULT (Nag_TransfOrder NAG_HUGE *)0

#define MPS_NAME_LEN 8

#define NAG_ERROR_BUF_LEN 512
#define NAG_FILE_LEN 512

/* INIT_FAIL initialises a fail structure with print set to FALSE */
#define INIT_FAIL(E) (E).code = NE_NOERROR, (E).print = 0, \
(E).handler = 0, (E).errnum = 0

/* SET_FAIL initialises a fail structure with print set to TRUE */
#define SET_FAIL(E) (E).code = NE_NOERROR, (E).print = TRUE, \
(E).handler = 0, (E).errnum = 0

#define INIT_MESG(M) (M).code = NM_NO_MESG, (M).print = FALSE, \
(M).res_file = (char NAG_HUGE *)0, (M).name = (char NAG_HUGE *)0, (M).opt_name = (char NAG_HUGE *)0, \
(M).print_mesg = 0, (M).init_mesg1 = IDUMMY, (M).init_mesg2 = INIT2DUMMY

#define INIT_STREAM(S) (S).set = FALSE, (S).st_out = FALSE, \
(S).st_err = FALSE, (S).st_in = FALSE, (S).chapter = (char NAG_HUGE *)0, \
(S).error = NE_NOERROR

/* d02rac use, define two or four functions as null pointers in the d02rac
 * call. Saves us having to design a limited option setting facility for
 * this routine
 */
#define NULL_4_FUN NULLFN, NULLFN, NULLFN, NULLFN
#define NULL_2_FUN NULLFN, NULLFN

/* Magic numbers to denote initialisation state */
#define RDUMMY -11111.0
#define IDUMMY -11111
#define INIT2DUMMY -23456


#define Vprintf (void)printf
#define Vfprintf (void)fprintf
#define Vsprintf (void)sprintf
#define Vscanf (void)scanf
#define Vfscanf (void)fscanf
#define Vstrcpy (void)strcpy


#define ABS(x) (x>=0 ? x : -(x)) /* Absolute value */

#ifdef NAG_FABS
#define FABS(x) fabs(x)
#else
#define FABS(x) ABS(x)
#endif

#define SIGN(x,y) (y>0 ? ABS(x) : -ABS(x)) /* Sign transfer */

#define MAX(x,y) (x>=y ? x : y) /* Maximum of two arguments */

#define MIN(x,y) (x<y ? x : y) /* Minimum of two arguments */

/* Round real number to nearest integer, returned as a double: valid for all X */
#define DROUND(X) floor((X)+0.5)

/* Round real number to nearest integer: valid for all X below MAXINT */
#define ROUND(X) (Integer)(floor((X)+0.5))

/* Defines for Complex Values */

/* Square of the modulus of a Complex number */
#define SQZABS(Z) ((Z).re*(Z).re + (Z).im*(Z).im)

/* Complex conjugate of itself */
#define CONJ(A) ((A).im = -(A).im)

/* A gets the Complex conjugate of B */
#define VCONJ(A,B) ((A).re = (B).re, (A).im = -((B).im))

/* Complex multiplication. X = Y*Z */
#define ZMULT(X,Y,Z) { double zz_t; zz_t = (Y).re*(Z).re - (Y).im*(Z).im; \
                                  (X).im = (Y).re*(Z).im + (Y).im*(Z).re, \
                                  (X).re = zz_t; }                  

#include <nag_types.h>
#include <nag_errlist.h>
#include <nag_names.h>

#endif /* not NAG_H */
