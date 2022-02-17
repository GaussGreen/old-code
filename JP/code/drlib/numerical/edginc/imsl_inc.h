#define ANSI

#ifndef IMSL_INC_H

#define IMSL_INC_H



#include "imsl.h"



/*

 *      ANSI is defined using -DANSI on the compile line if function prototypes

 *      are used.

 *

 *      COMPUTER_?????? is defined using the -D option on the compile line.

 *      An attempt to define the variable is also made in this file.

 *      Nonstandard COMPUTER_?????? names are

 *              GSUN3    GNU compiler on a Sun 3

 *              GSUN4    GNU compiler on a Sun 4

 *              TC       Turbo C on a PC
 
 *              IRIX_N32 SGI IRIX n32 ABI

 *      Here ?????? is the standard IMSL name for the environment,

 *      e.g., SUN, PMXUX, etc.

 *      Based on this variable, preprocessor

 *      variables of the form IMSL_MACHINE_?????? are defined in this

 *      header file.

 *       Variabled                              Defined if

 *      IMSL_MACHINE_SUN3               Sun 3 hardware, any compiler

 *      IMSL_MACHINE_SUN4               Sun 4 (SPARC) hardware, any compiler

 *      IMSL_MACHINE_SUN                Any Sun hardware

 *      IMSL_MACHINE_MC68020            MC68020 hardware, any vendor, OS, compiler

 *      IMSL_MACHINE_MIPS               MIPS hardware, any vendor, OS, compiler

 *      IMSL_MACHINE_VAX                VAX hardware, any OS or compiler

 *      IMSL_MACHINE_DECOSF             ALPHA hardware, any OSF or compiler

 *      IMSL_MACHINE_80X86              Intel 8086, 80286, 80386, etc., hardware

 *      IMSL_MACHINE_DOS                DOS, any compiler

 *      IMSL_MACHINE_NT                 NT, any compiler

 *      IMSL_MACHINE_GNU                GNU compiler, any hardware or OS

 *      IMSL_MACHINE_IEEE               IEEE binary arithmetic is used.

 *      IMSL_MACHINE_LITTLE_ENDIAN      Little endian integers are used.

 *

 *      COMPUTER_?????? =>              IMSL_MACHINE_?????

 *      SUN3            =>      SUN3, SUN, MC68020

 *      SUN4            =>      SUN4, SUN

 *      GSUN3           =>      SUN3, SUN, MC68020, GNU

 *      GSUN4           =>      SUN4, SUN, GNU

 *      VAX             =>      VAX

 *      MIPS            =>      MIPS
 
 *      IRIX_N32        =>      MIPS

 *      PMXUX           =>      MIPS

 *      RTXLXS          =>      IEEE

 *      TC              =>      80X86, DOS

 *

 *      Based on the above IMSL_MACHINE_????? definitions, the following

 *      IMSL_MACHINE_????? definitions are made.  These depend on the hardware,

 *      regardless of the compiler or the OS.

 *

 *      IEEE            <=      !VAX

 *      LITTLE_ENDIAN   <=      80X86, COMPUTER_MIPS, VAX

 *

 *      Based on the above IMSL_MACHINE_????? definitions, the following

 *      IMSL_OS_????? definitions are made.  These depend on the OS,

 *      regardless of the compiler or the hardware.

 *

 *      SHORT_FILENAMES <=      DOS

 */







    /* define COMPUTER_???? based on other defined variables */

#if (defined(sun) && defined(mc68000)) && !defined(__GNUC__) && !defined(COMPUTER_SUN)

#define COMPUTER_SUN

#endif

/*
** Had to put and extra check for SOLARIS because "sparc" is defined 
** for both SunOS 4... and Sun 5... (Solaris).
*/

#if defined(sparc) && !defined(__GNUC__) && !defined(COMPUTER_SUN4) && !defined(SOLARIS)

#define COMPUTER_SUN4

#endif

/*
** Definition for IRIX (added for port)
*/
#if defined(IRIX)
#define COMPUTER_IRIX_N32
#endif

/*
** Definition for LINUX (added for port)
*/
#if defined(LINUX)
#define COMPUTER_LINUX
#endif


/*
** Define a machine for compiling under SunOS 5.5.1 (Solaris).
** The user need only set define SOLARIS.
*/

#if defined(SOLARIS)
#   if !defined(COMPUTER_SUN5)
#       define COMPUTER_SUN5
#   endif
#
#   if defined(__GNUC__)
#      if !defined(COMPUTER_SLRSGCC)
#       define COMPUTER_SLRSGCC
#      endif
#   else
#      if !defined(COMPUTER_SLRSACC) 
#       define COMPUTER_SLRSACC
#      endif
#   endif
#endif

#if (defined(sun) && defined(mc68000)) && defined(__GNUC__) && !defined(COMPUTER_GSUN)

#define COMPUTER_GSUN

#endif



#if defined(sparc) && defined(__GNUC__) && !defined(COMPUTER_GSUN4)  && !defined(COMPUTER_SLRSGCC)

#define COMPUTER_GSUN4

#endif



#if defined(vax) && !defined(COMPUTER_VAX)

#define COMPUTER_VAX

#endif



#if defined(vax) && !defined(COMPUTER_VAXG) && !defined(COMPUTER_VAX)

#define COMPUTER_VAXG

#endif



    /* define IMSL_MACHINE_???? based on definitions of COMPUTER_???? */



#if defined(COMPUTER_APLC)

#define IMSL_MACHINE_APLC

#endif



#if defined(COMPUTER_SUN) || defined(COMPUTER_GSUN)

#define IMSL_MACHINE_SUN3

#endif



#if defined(COMPUTER_SUN4) || defined(COMPUTER_GSUN4)

#define IMSL_MACHINE_SUN4

#endif


#if defined(COMPUTER_SUN5)

#define IMSL_MACHINE_SUN5

#endif



#if defined(IMSL_MACHINE_SUN3) || defined(IMSL_MACHINE_SUN4) || defined(IMSL_MACHINE_SUN5)

#define IMSL_MACHINE_SUN

#endif



#if defined(COMPUTER_GSUN) || defined(COMPUTER_GSUN4) || defined(COMPUTER_SLRSGCC)

#define IMSL_MACHINE_GNU

#endif



#if defined(IMSL_MACHINE_SUN3)

#define IMSL_MACHINE_MC68020

#endif



#if defined(COMPUTER_HP98C)

#define IMSL_MACHINE_HP98C

#endif



#if defined(COMPUTER_HP97C)

#define IMSL_MACHINE_HP97C

#endif



#if defined(COMPUTER_VAX)

#define IMSL_MACHINE_VAX

#endif



#if defined(COMPUTER_DECOSF)

#define IMSL_MACHINE_DECOSF

#endif



#if defined(COMPUTER_VAXG)

#define IMSL_MACHINE_VAX

#endif



#if defined(COMPUTER_PMXUX) || defined(COMPUTER_MIPNT) || defined(COMPUTER_IRIX_N32)

#define IMSL_MACHINE_MIPS

#endif



#if defined(COMPUTER_TC)

#define IMSL_MACHINE_DOS

#endif


#if defined(COMPUTER_LINUX)

#define IMSL_MACHINE_LINUX

#endif


#if defined(IMSL_MACHINE_DOS) || defined(COMPUTER_NTLNT) || defined(IMSL_MACHINE_LINUX)

#define IMSL_MACHINE_80X86

#endif



#if defined(COMPUTER_ALFAC_IEEE)

#if !defined(COMPUTER_ALFAC)

#define COMPUTER_ALFAC

#endif

#endif



#if defined(COMPUTER_NTLNT) || defined(COMPUTER_ALFNT) || defined(COMPUTER_MIPNT)

#define IMSL_MACHINE_NT

#endif


#if defined(COMPUTER_ALFAC) || defined(COMPUTER_ALFNT)

#define IMSL_MACHINE_ALPHA

#endif



#if !defined(IMSL_MACHINE_VAX)

#define IMSL_MACHINE_IEEE

#endif



#if defined(IMSL_MACHINE_VAX) || defined(IMSL_MACHINE_DECOSF) || defined(IMSL_MACHINE_MIPS) || defined(IMSL_MACHINE_80X86) || defined(IMSL_MACHINE_NT)

#define IMSL_MACHINE_LITTLE_ENDIAN

#endif



#if defined(IMSL_MACHINE_DOS)

#define IMSL_OS_SHORT_FILENAMES

#endif



#ifdef COMPUTER_RTXLXS

#define _POSIX_SOURCE

#endif


/*
** Always use the macro definition of nint for Solaris ACC 
** compilations. (Jim M.)
** Always use the macro definition of nint for ALL platforms
** (Charles Irwin). 
*/

#define  nint(x) ((x) < 0 ? (int)((x)-.5) : (int)((x)+.5))

#ifndef STDIO_H__

#define STDIO_H__

#include <stdio.h>

#endif



#ifndef STRING_H__

#define STRING_H__

#include <string.h>

#endif



#ifndef MATH_H__

#define MATH_H__

#include <math.h>

#endif



#ifndef SETJMP_H__

#define SETJMP_H__

#include <setjmp.h>

#endif



#ifndef CTYPE_H__

#define CTYPE_H__

#include <ctype.h>

#ifdef IMSL_MACHINE_NT

#undef isalpha

#undef isupper

#undef islower

#undef isdigit

#undef isxdigit

#undef isspace

#undef ispunct

#undef isalnum

#undef isprint

#undef isgraph

#undef iscntrl

#undef tolower

#undef toupper

#endif

#endif



#ifndef LIMITS_H__

#define LIMITS_H__

#include <limits.h>

#endif



#ifdef ANSI

#ifndef STDARG_H__

#define STDARG_H__

#include <stdarg.h>

#endif

#else

#ifndef VARARGS_H__

#define VARARGS_H__

#include <varargs.h>

#endif

#endif



/*

#if defined(COMPUTER_RTXLXS) || defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C)

extern     double   hypot(double,double);

#else

#if defined(COMPUTER_HP93C)

extern     double   hypot();

#endif

#endif

*/

#ifdef IMSL_MACHINE_NT

#define hypot   _hypot

#else

#ifdef ANSI

extern     double   hypot(double,double);

#endif

#endif



/******************************************

#if defined(ANSI) || defined(COMPUTER_PMXUX)

void        *malloc();

void        *realloc();

#else

char        *malloc();

char        *realloc();

#endif

int         free();

****************************************/

/**** ABOVE BLOCK REPLACED BY THE FOLLOWING ******/

#if defined(ANSI) || defined(IMSL_MACHINE_SUN)

#ifndef STDLIB_H__

#define STDLIB_H__

#include <stdlib.h>

#endif

#endif

#if defined(COMPUTER_PMXUX)

char        *malloc();

char        *realloc();

void         free();

#endif



#ifdef DOUBLE

#define  E1PSH(FNAME,DNAME)     imsl_e1psh(DNAME)

#define  E1POP(FNAME,DNAME)     imsl_e1pop(DNAME)

#else

#define  E1PSH(FNAME,DNAME)     imsl_e1psh(FNAME)

#define  E1POP(FNAME,DNAME)     imsl_e1pop(FNAME)

#endif



#ifdef ANSI

#define VA_START(VA_LIST,LASTFIX)  va_start(VA_LIST,LASTFIX)

#else /* not ANSI */

#define VA_START(VA_LIST,LASTFIX)  va_start(VA_LIST)

#endif /* ANSI */




#ifndef LINUX64
/* On AMD64 va_list is an array and you cannot copy/return it. So instead we return void. */
#define VA_LIST_HACK va_list
#define IMSL_CALL(L_PROC) \
    if (setjmp(imsl_error_param.jmpbuf_environ[imsl_got_environ++]) == 0) {\
	imsl_set_signal(1);             \
	argptr = L_PROC;                \
	imsl_got_environ--;             \
    } else                              \
	imsl_ermes(IMSL_TERMINAL, IMSL_MAJOR_VIOLATION);\
    imsl_set_signal(0);
#else /* LINUX64 */
#define VA_LIST_HACK void
#define IMSL_CALL(L_PROC) \
    if (setjmp(imsl_error_param.jmpbuf_environ[imsl_got_environ++]) == 0) {\
	imsl_set_signal(1);             \
	L_PROC;				\
	imsl_got_environ--;             \
    } else                              \
	imsl_ermes(IMSL_TERMINAL, IMSL_MAJOR_VIOLATION);\
    imsl_set_signal(0);
#endif


#ifndef BLAS_H__

#define BLAS_H__

#include "blas.h"

#endif



#if defined(COMPUTER_PMXUX)

typedef char        Mvoid;

#else

#if defined(COMPUTER_VAX) || defined(COMPUTER_VAXG) || defined(COMPUTER_ALFAC_IEEE) || defined(COMPUTER_DECOSF)

#define Mvoid       void

#else

typedef void            Mvoid;

#endif

#endif



#if (defined(COMPUTER_DECOSF) && defined(CMATH_MINT_TO_LONG))

#define Mchar           char

#define Mint            long

#define Muint           unsigned int

#define Mlong           long

#define Mdouble         double

#elif defined(COMPUTER_DECOSF)

#define Mchar           char

#define Mint            int

#define Muint           unsigned int

#define Mlong           int

#define Mdouble         double

#else

typedef char            Mchar;

typedef int             Mint;

typedef unsigned int    Muint;

typedef long            Mlong;

typedef double          Mdouble;

#endif

typedef d_complex       Md_complex;

#ifndef DOUBLE

#if defined(COMPUTER_DECOSF)

#define Mfloat          float

#else

typedef float           Mfloat;

#endif

typedef f_complex       Mf_complex;

typedef Imsl_f_ppoly    Mf_ppoly;

typedef Imsl_f_spline   Mf_spline;

typedef ip_f_node       Mf_node;

typedef Imsl_f_radial_basis_fit Mf_radial_basis_fit;

typedef f_radial_fcn    Mf_radial_fcn;

typedef f_constraint_struct     Mf_constraint_struct;

typedef Imsl_f_sparse_elem      Mf_sparse_elem;

typedef Imsl_f_numeric_factor   Mf_numeric_factor;

typedef Imsl_c_sparse_elem      Mc_sparse_elem;

typedef Imsl_c_numeric_factor   Mc_numeric_factor;

typedef Imsl_f_sparse_list_element Mf_sparse_list_element;

typedef Imsl_c_sparse_list_element Mc_sparse_list_element;

typedef Imsl_f_header_element   Mf_header_element;

typedef Imsl_c_header_element   Mc_header_element;

typedef Imsl_f_memory_list      Mf_memory_list;

typedef Imsl_c_memory_list      Mc_memory_list;

typedef Imsl_f_sparse_lu_factor Mf_sparse_lu_factor;

typedef Imsl_c_sparse_lu_factor Mc_sparse_lu_factor;

#else /* DOUBLE */

#if defined(COMPUTER_DECOSF)

#define Mfloat          double

#else

typedef double          Mfloat;

#endif

typedef d_complex       Mf_complex;

typedef Imsl_d_ppoly    Mf_ppoly;

typedef Imsl_d_spline   Mf_spline;

typedef ip_d_node       Mf_node;

typedef Imsl_d_radial_basis_fit Mf_radial_basis_fit;

typedef d_radial_fcn    Mf_radial_fcn;

typedef d_constraint_struct     Mf_constraint_struct;

typedef Imsl_d_sparse_elem      Mf_sparse_elem;

typedef Imsl_d_numeric_factor   Mf_numeric_factor;

typedef Imsl_z_sparse_elem      Mc_sparse_elem;

typedef Imsl_z_numeric_factor   Mc_numeric_factor;

typedef Imsl_d_sparse_list_element Mf_sparse_list_element;

typedef Imsl_z_sparse_list_element Mc_sparse_list_element;

typedef Imsl_d_header_element   Mf_header_element;

typedef Imsl_z_header_element   Mc_header_element;

typedef Imsl_d_memory_list      Mf_memory_list;

typedef Imsl_z_memory_list      Mc_memory_list;

typedef Imsl_d_sparse_lu_factor Mf_sparse_lu_factor;

typedef Imsl_z_sparse_lu_factor Mc_sparse_lu_factor;

#endif /* DOUBLE */



#ifdef DOUBLE



#ifdef WAVE_RENAME

#include "renamefd.h"



#ifdef saxpy

#undef saxpy

#define saxpy(N,A,X,INCX,Y,INCY)        AXPY(double,N,A,X,INCX,Y,INCY)

#endif



#ifdef scopy

#undef scopy

#define scopy(N,X,INCX,Y,INCY)  COPY(double,N,X,INCX,Y,INCY)

#endif



#ifdef sswap

#undef sswap

#define sswap(N,X,INCX,Y,INCY)  SWAP(double,N,X,INCX,Y,INCY)

#endif



#ifdef sscal

#undef sscal

#define sscal(N,A,X,INCX)       SCAL(double,N,A,X,INCX)

#endif



#ifdef sset

#undef sset

#define sset(N,A,X,INCX)        SET(double,N,A,X,INCX)

#endif



#ifdef sadd

#undef sadd

#define sadd(N,A,X,INCX)        ADD(double,N,A,X,INCX)

#endif



#define Imsl_f_ppoly                    Imsl_d_ppoly

#define Imsl_f_spline                   Imsl_d_spline



#else



#include "imsl_fd.h"



#endif /* WAVE_RENAME */

#endif /*   DOUBLE  */



#ifndef ERROR_C

extern

#endif

struct {

    Mint        signal_set;

    Mint        i[10];

    Mdouble     d[10];

    Md_complex  z[10];

    Mchar      *l[10];

    Mchar      *path;

    Mchar      *name;

    jmp_buf    jmpbuf_environ[100];

} imsl_error_param;





#ifdef ERROR_C

Mint                imsl_got_environ;

jmp_buf             imsl_environ;

Mint                imsl_error_n1rty[2];

#ifdef IMSL_MACHINE_VAX

struct {

    float     f[8];

    double    d[8];

} imsl_machine = {

/* it's a VAX using the /NOG_FLOAT option */

#if defined(COMPUTER_VAX) && !defined(COMPUTER_VAXG)

    2.93941e-39, /* amach */ /************************************/

    1.70141e38,              /* THIS SHOULD BE REPLACED BY EXACT */

    5.96184e-8,              /* VALUES, AS BELOW.                */

    1.19237e-7,              /************************************/

    0.301029995663981195,

    1.70101e38,

    1.70102e38,

   -1.70102e38,

    2.93941e-39,    /* dmach */

    1.70141e38,

    1.38810e-17,

    2.77620e-17,

    0.3010299956639811952137388947245,

    1.70101e38,

    1.70102e38,

   -1.70102e38

#else  /* it's a VAX using the /G_FLOAT option */

       /* values obtained from Mike A.'s const.f program */

    0.29387359e-38,

    0.17014117e+39,

    0.59604645e-07,

    0.11920929e-06,

    0.30103001e+00,

    0.17014116e+39,

    0.17014117e+39,

   -0.17014117e+39,

    0.56e-308,            /* dmach  */

    0.898846567431e+308,

    0.111022302462516e-15,

    0.222044604925031e-15,

    0.301029995663981e+00,

    0.898845710224272e+308,

    0.898846567431e+308,

   -0.898846567431e+308

#endif

};

#else  /* not VAX */

struct {

#if (defined(COMPUTER_DECOSF) && defined(CMATH_MINT_TO_LONG))

    int         f[8];

    int         d[16];

#elif defined(COMPUTER_DECOSF)

    Mint         f[8];

    Mint         d[16];

#elif defined(LINUX64) && defined(IMSL_MACHINE_LITTLE_ENDIAN)
  int          f[8];
  long int     d[8];
#else

    long int     f[8];

    long int     d[16];

#endif

} imsl_machine = {

    /* amach */
  /*
    .f[]={
    GSL_FLOAT_MIN_NORMAL,
    GSL_FLOAT_MAX_NORMAL,
    2^-24,
    2^-23,
    IEEE Log10(2),
    IEEE Quiet NaN,
    IEEE Positive Infinity,
    IEEE Negatie Infinity}
  */
    
    
    
  0x00800000, 0x7F7FFFFF, 0x33800000, 0x34000000, 0x3e9a209b,

   0x7FFFFFFF, 0x7F800000, 0xFF800000,

    /* dmach */

#ifdef IMSL_MACHINE_LITTLE_ENDIAN
#ifdef LINUX64

    0x0010000000000000UL,  0x7FEFFFFFFFFFFFFFUL,  0x3ca0000000000000UL,
    0x3cb0000000000000UL, /* 0x3FD34413509F79FFUL */  0x3fd344133fd34413UL /* BUG?*/,  0x7FFFFFFFFFFFFFFFUL,
    0x7FF0000000000000UL,  0xFFF0000000000000UL
#else

    0x00000000, 0x00100000,  0xFFFFFFFF, 0x7FEFFFFF,  0x00000000, 0x3ca00000,

    0x00000000, 0x3cb00000,  0x3fd34413 /* bug?*/ , 0x3fd34413,  0xFFFFFFFF, 0x7FFFFFFF,

    0x00000000, 0x7FF00000,  0x00000000, 0xFFF00000
#endif
#else

    0x00100000, 0x00000000,  0x7FEFFFFF, 0xFFFFFFFF,  0x3ca00000, 0x00000000,

    0x3cb00000, 0x00000000,  0x3fd34413, 0x3fd34413 /* bug */,  0x7FFFFFFF, 0xFFFFFFFF,

    0x7FF00000, 0x00000000,  0xFFF00000, 0x00000000

#endif  /* IMSL_MACHINE_LITTLE_ENDIAN */

};

#endif  /* VAX */

#endif  /* ERROR_C */





	/* definition of DOCNT is from f_rt.h */

	/* support for DO and arithmetic IF translations */

/* formula:  DOCNT(ini,tst,inc) = max( (long)((tst-ini+inc)/inc), 0 ) */

#define DOCNT(i,t,n)    (_d_l=n, (_d_m=(t-(i)+_d_l)/_d_l) > 0 ? _d_m : 0L )



#define dble(_Z)    (_Z).re

#define dimag(_Z)   (_Z).im



#ifdef sun

#define imsl_ifnan     isnan

#else

#if defined(COMPUTER_VAX) || defined(COMPUTER_VAXG) || defined(COMPUTER_ALFAC_IEEE)



#ifdef DOUBLE

#define imsl_ifnan(P)   ((P)==(imsl_dmach(6)))

#else

#define imsl_ifnan(P)   ((P)==(imsl_amach(6)))

#endif



#else



#if defined(COMPUTER_APLC) || defined(IMSL_MACHINE_80X86) || defined(IMSL_MACHINE_NT) || defined(COMPUTER_PMXUX)

Mint    IMSL_PROTO(imsl_ifnan,(Mfloat));

#else

# ifdef imsl_ifnan

#  undef imsl_ifnan

# endif 

#define imsl_ifnan(P)   ((P)!=(P))

#endif

#endif

#endif



#define mod(a,b)            ((a) % (b))

/*

#define sign(P,Q)           (((Q)<F_ZERO) ? -(P) : (P))

*/

#define sign(P,Q)           (((Q)<F_ZERO) ? -fabs(P) : fabs(P))

#define imsl_toupper(c)     ((c)-'a'+'A')



#ifdef powl

#undef powl

#endif

#define powl(X,K)           imsl_fi_power(X,K)



	/* Define names for common float constants */

#ifdef IMSL_MACHINE_NT



#define F_ZERO  0.

#define F_ONE   1.

#define F_TWO   2.

#define F_THREE 3.

#define F_FOUR  4.

#define F_FIVE  5.

#define F_SIX   6.

#define F_SEVEN 7.

#define F_EIGHT 8.

#define F_NINE  9.

#define F_TEN   10.

#define F_HALF  0.5

#define C_ZERO  imsl_cf_convert(0.,0.)

#define C_ONE   imsl_cf_convert(1.,0.)



#else



#ifdef ERROR_C

float  imsl_F_NUMBER[12] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 0.5};

double imsl_D_NUMBER[12] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 0.5};

f_complex imsl_C_NUMBER[2] = {{0.0,0.0}, {1.0,0.0}};

d_complex imsl_Z_NUMBER[2] = {{0.0,0.0}, {1.0,0.0}};

#else

extern float     imsl_F_NUMBER[12];

extern double    imsl_D_NUMBER[12];

extern f_complex imsl_C_NUMBER[2];

extern d_complex imsl_Z_NUMBER[2];

#endif



#ifndef DOUBLE

#define F_ZERO  imsl_F_NUMBER[0]

#define F_ONE   imsl_F_NUMBER[1]

#define F_TWO   imsl_F_NUMBER[2]

#define F_THREE imsl_F_NUMBER[3]

#define F_FOUR  imsl_F_NUMBER[4]

#define F_FIVE  imsl_F_NUMBER[5]

#define F_SIX   imsl_F_NUMBER[6]

#define F_SEVEN imsl_F_NUMBER[7]

#define F_EIGHT imsl_F_NUMBER[8]

#define F_NINE  imsl_F_NUMBER[9]

#define F_TEN   imsl_F_NUMBER[10]

#define F_HALF  imsl_F_NUMBER[11]

#define C_ZERO  imsl_C_NUMBER[0]

#define C_ONE   imsl_C_NUMBER[1]

#else

#define F_ZERO  imsl_D_NUMBER[0]

#define F_ONE   imsl_D_NUMBER[1]

#define F_TWO   imsl_D_NUMBER[2]

#define F_THREE imsl_D_NUMBER[3]

#define F_FOUR  imsl_D_NUMBER[4]

#define F_FIVE  imsl_D_NUMBER[5]

#define F_SIX   imsl_D_NUMBER[6]

#define F_SEVEN imsl_D_NUMBER[7]

#define F_EIGHT imsl_D_NUMBER[8]

#define F_NINE  imsl_D_NUMBER[9]

#define F_TEN   imsl_D_NUMBER[10]

#define F_HALF  imsl_D_NUMBER[11]

#define C_ZERO  imsl_Z_NUMBER[0]

#define C_ONE   imsl_Z_NUMBER[1]

#endif



#endif



/*   #if !defined(CCOMPLEX_C)  */

/*   for now, we are not using the inline complex macros, JFB */

#if 0



#ifdef ERROR_C

Mint                    imsl_CS_k = 0;

f_complex               imsl_CS[10];

Mint                    imsl_ZS_k = 0;

d_complex               imsl_ZS[10];

#else /* ERROR_C */

extern Mint             imsl_CS_k;

extern f_complex        imsl_CS[10];

extern Mint             imsl_ZS_k;

extern d_complex        imsl_ZS[10];

#endif /* ERROR_C */



#define CS0     imsl_CS[imsl_CS_k]

#define CS1     imsl_CS[imsl_CS_k-1]

#define CS2     imsl_CS[imsl_CS_k-2]

#define CS0pp   imsl_CS[imsl_CS_k++]

#define ZS0     imsl_ZS[imsl_ZS_k]

#define ZS1     imsl_ZS[imsl_ZS_k-1]

#define ZS2     imsl_ZS[imsl_ZS_k-2]

#define ZS0pp   imsl_ZS[imsl_ZS_k++]



#define imsl_cz_convert(Q) (ZS0 = (Q), \

	CS0.re = (float)ZS0.re, CS0.im = (float)ZS0.im, CS0)

#define imsl_zc_convert(Q) (CS0 = (Q), \

	ZS0.re = (double)CS0.re, ZS0.im = (double)CS0.im, ZS0)





#ifndef DOUBLE

#define imsl_c_eq(P,Q)          (CS0pp = (P), imsl_CS[imsl_CS_k] = (Q), \

	imsl_CS_k--, (CS0.re==imsl_CS[imsl_CS_k+1].re) && (CS0.im==imsl_CS[imsl_CS_k+1].im))

#define imsl_c_neq(P,Q)         (CS0pp = (P), imsl_CS[imsl_CS_k] = (Q), \

	imsl_CS_k--, CS0.re != imsl_CS[imsl_CS_k+1].re || CS0.im != imsl_CS[imsl_CS_k+1].im)

#define imsl_cf_convert(P,Q)    (CS0.re = (P), CS0.im = (Q), CS0)

#define imsl_fc_convert(Q)      (CS0 = (Q), (float)CS0.re)

#define imsl_c_aimag(Q)         (CS0 = (Q), (float)CS0.im)

#define imsl_c_conjg(Q)         (CS0 = (Q), CS0.im = -CS0.im, CS0)

#define imsl_c_neg(Q)           (CS0 = (Q), CS0.re = -CS0.re, CS0.im = -CS0.im, CS0)

#define imsl_c_abs1(P)          (CS0 = (P), (float)(fabs(CS0.re)+fabs(CS0.im)))

#define imsl_c_abs(P)           (CS0 = (P), (float)hypot(CS0.re, CS0.im))

#define imsl_c_add(P,Q)         (CS0pp = (P), CS0 = (Q), \

	CS1.re += CS0.re, CS1.im += CS0.im, imsl_CS[--imsl_CS_k])

#define imsl_c_sub(P,Q)  ( CS0pp = (P), CS0 = (Q), \

	CS1.re -= CS0.re, CS1.im -= CS0.im, imsl_CS[--imsl_CS_k])

#define imsl_c_mul(P,Q)  ( \

	CS0pp = (P), CS0pp = (Q), \

	CS0.re = CS2.re*CS1.re - CS2.im*CS1.im, \

	CS0.im = CS2.re*CS1.im + CS2.im*CS1.re, \

	imsl_CS_k -= 2, imsl_CS[imsl_CS_k+2])

#else /* DOUBLE */

#define imsl_c_eq(P,Q)          (ZS0pp = (P), imsl_ZS[imsl_ZS_k] = (Q), \

	imsl_ZS_k--, (ZS0.re==imsl_ZS[imsl_ZS_k+1].re) && (ZS0.im==imsl_ZS[imsl_ZS_k+1].im))

#define imsl_c_neq(P,Q)         (ZS0pp = (P), imsl_ZS[imsl_ZS_k] = (Q), \

	imsl_ZS_k--, ZS0.re != imsl_ZS[imsl_ZS_k+1].re || ZS0.im != imsl_ZS[imsl_ZS_k+1].im)

#define imsl_cf_convert(P,Q)    (ZS0.re = (P), ZS0.im = (Q), ZS0)

#define imsl_fc_convert(Q)      (ZS0 = (Q), ZS0.re)

#define imsl_c_aimag(Q)         (ZS0 = (Q), ZS0.im)

#define imsl_c_conjg(Q)         (ZS0 = (Q), ZS0.im = -ZS0.im, ZS0)

#define imsl_c_neg(Q)           (ZS0 = (Q), ZS0.re = -ZS0.re, ZS0.im = -ZS0.im, ZS0)

#define imsl_c_abs(P)           (ZS0 = (P), hypot(ZS0.re, ZS0.im))

#define imsl_c_abs1(P)          (ZS0 = (P), fabs(ZS0.re)+fabs(ZS0.im))

#define imsl_c_add(P,Q)  (ZS0pp = (P), ZS0 = (Q), \

	ZS1.re += ZS0.re, ZS1.im += ZS0.im, imsl_ZS[--imsl_ZS_k])

#define imsl_c_sub(P,Q)  ( ZS0pp = (P), ZS0 = (Q), \

	ZS1.re -= ZS0.re, ZS1.im -= ZS0.im, imsl_ZS[--imsl_ZS_k])

#define imsl_c_mul(P,Q)  ( \

	ZS0pp = (P), ZS0pp = (Q), \

	ZS0.re = ZS2.re*ZS1.re - ZS2.im*ZS1.im, \

	ZS0.im = ZS2.re*ZS1.im + ZS2.im*ZS1.re, \

	imsl_ZS_k -= 2, imsl_ZS[imsl_ZS_k+2])

#endif  /* DOUBLE */

#endif  /* CCOMPLEX_C   (inline complex) */





#ifndef ERROR_C

extern Mint         imsl_got_environ;

extern jmp_buf      imsl_environ;

extern Mint         imsl_error_n1rty[2];

#define  imsl_e1sti(P,Q)    imsl_error_param.i[P] = Q

#define  imsl_e1str(P,Q)    imsl_error_param.d[P] = Q

#define  imsl_e1std(P,Q)    imsl_error_param.d[P] = Q

#ifdef DOUBLE

#define  imsl_e1stc(P,Q)    imsl_error_param.z[P] = Q;

#else

#define  imsl_e1stc(P,Q)    imsl_error_param.z[P] = imsl_zc_convert(Q);

#endif

#define  imsl_e1stz(P,Q)    imsl_error_param.z[P] = Q

#define  imsl_e1stl(P,Q)    imsl_error_param.l[P] = Q

# ifdef imsl_n1rty

#  undef imsl_n1rty

# endif

#define  imsl_n1rty(P)      imsl_error_n1rty[P]

#ifdef IMSL_MACHINE_NT

Mfloat  IMSL_PROTO(imsl_amach,(Mint));

#ifndef DOUBLE

Mdouble IMSL_PROTO(imsl_dmach,(Mint));

#else

#define imsl_amach imsl_dmach

#endif

#undef imsl_n1rty

Mint    IMSL_PROTO(imsl_n1rty,(Mint));

#else   /* IMSL_MACHINE_NT */

#undef imsl_amach

#undef imsl_dmach

#ifdef DOUBLE

#define imsl_amach(P)       (imsl_machine.d[(P)-1])

#else

#define imsl_amach(P)       (imsl_machine.f[(P)-1])

#endif

#define imsl_dmach(P)       (imsl_machine.d[(P)-1])

#endif  /* IMSL_MACHINE_NT */

extern struct {

    float     f[8];

    double    d[8];

} imsl_machine;

#endif /* ERROR_C */



#define PROTO(PP,QQ)        IMSL_PROTO(PP,QQ)



#ifndef iposdif

#define iposdif(A_, B_)     ((A_) > (B_) ? (A_) - (B_) : 0)

#endif





	/* Chapter 1 --- Linear System */

void        PROTO(imsl_l2trg,(Mint, Mfloat*, Mint, Mfloat*, Mint, Mint*, Mfloat*));

void        PROTO(imsl_lfsrg,(Mint, Mfloat*, Mint, Mint*, Mfloat*, Mint*, Mfloat*));

void        PROTO(imsl_l2rrr,(Mint*, Mint*, Mfloat*, Mint*, Mint*, Mint*, Mfloat*, Mint*,

			      Mfloat*, Mfloat*, Mfloat*));

void        PROTO(imsl_l2qrr,(Mint*, Mint*, Mfloat*, Mint *, Mfloat *,

			       Mfloat *, Mfloat *, Mfloat *, Mint *,

			       Mfloat *, Mfloat *, Mint *, Mfloat *));

void        PROTO(imsl_lqrsl,(Mint *, Mint *, Mfloat *, Mint *, Mfloat *,

			       Mfloat *, Mint *, Mfloat *, Mfloat *, Mfloat *,

			       Mfloat *, Mfloat *));

void        PROTO(imsl_cgeru,(Mint*, Mint*, Mf_complex*, Mf_complex*,

			       Mint*, Mf_complex*, Mint*, 

			       Mf_complex*, Mint*));



	/* Chapter 2 --- Eigenvalues */



	/* Chapter 3 --- Interpolation and Approximation */

void        PROTO(imsl_c2int,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mfloat*, Mint[]));

Mfloat      PROTO(imsl_csder,(Mint*, Mfloat*, Mint*, Mfloat[], Mfloat*));

Mfloat      PROTO(imsl_ppder,(Mint, Mfloat, Mint, Mint, Mfloat[], Mfloat*));

void        PROTO(imsl_p3der,(Mint, Mint, Mfloat*, Mfloat, Mint*));

Mfloat      PROTO(imsl_csitg,(Mfloat*, Mfloat*, Mint*, Mfloat[], Mfloat*));

Mfloat      PROTO(imsl_ppitg,(Mfloat, Mfloat, Mint, Mint, Mfloat[], Mfloat*));



void        PROTO(imsl_c1sor,(Mint, Mfloat[], Mfloat[], Mfloat[], Mfloat*, Mint, Mint[]));

void        PROTO(imsl_c2dec,(Mint*, Mfloat[], Mfloat[], Mint*, Mfloat*, Mint*,

			      Mfloat*, Mfloat[], Mfloat*, Mint[]));





Mfloat      PROTO(imsl_b2der,(Mint*, Mfloat*, Mint*, Mfloat[],

			      Mint*, Mfloat[],Mfloat[], Mfloat[], Mfloat[]));

Mfloat      PROTO(imsl_b3der,(Mint*, Mfloat*, Mint*, Mfloat[],

			      Mint*, Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_b4der,(Mfloat[], Mint*, Mfloat*, Mint*, Mint*));

void        PROTO(imsl_b2int,(Mint*, Mfloat[], Mfloat[], Mint*, Mfloat[], Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[], Mint[]));

void        PROTO(imsl_b3int,(Mint*, Mfloat[], Mint*));

void        PROTO(imsl_b4int,(Mfloat[], Mint*, Mfloat*, Mint*, Mfloat[],

			      Mfloat[], Mfloat[]));

void        PROTO(imsl_b5int,(Mint*, Mfloat[], Mfloat[], Mint*, Mfloat[], Mfloat[],

			      Mfloat*, Mint*, Mfloat[], Mint[]));

Mfloat      PROTO(imsl_b2itg,(Mfloat*, Mfloat*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat[]));

Mfloat      PROTO(imsl_b3itg,(Mfloat*, Mfloat*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat[]));

Mfloat      PROTO(imsl_b4itg,(Mfloat*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat[], Mfloat[], Mint*));

void        PROTO(imsl_b5itg,(Mfloat[], Mint*, Mfloat*, Mint*, Mint*));

void        PROTO(imsl_b2nak,(Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat[],Mint[]));

void        PROTO(imsl_c1not,(Mchar*, Mchar*, Mint*, Mfloat[],

			      Mint*, Mfloat[]));

void        PROTO(imsl_b2opk,(Mint*, Mfloat[], Mint*, Mfloat[],

			      Mint*, Mfloat[], Mint[]));

void        PROTO(imsl_b3opk,(Mint*, Mfloat[], Mint*, Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat[],

			      Mint*, Mfloat[], Mint[]));

void        PROTO(imsl_b4opk,(Mint*, Mint*, Mfloat*, Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[],

			      Mfloat[], Mfloat[], Mfloat*, Mint*));

void        PROTO(imsl_crbrb,(Mint*, Mfloat*, Mint*, Mint*, Mint*,

			      Mfloat*, Mint*, Mint*, Mint*));

void        PROTO(imsl_l2lrb,(Mint*, Mfloat*, Mint*, Mint*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat[], Mint[], Mfloat[]));

void        PROTO(imsl_l2crb,(Mint*, Mfloat*, Mint*, Mint*, Mint*, Mfloat*,

			      Mint*, Mint[], Mfloat*, Mfloat[]));

void        PROTO(imsl_l2trb,(Mint*, Mfloat*, Mint*, Mint*, Mint*,

			      Mfloat*, Mint*, Mint[], Mfloat[]));

void        PROTO(imsl_lfsrb,(Mint*, Mfloat*, Mint*, Mint*,

			      Mint*, Mint[], Mfloat[], Mint*, Mfloat[]));

void        PROTO(imsl_nr1rb,(Mint*, Mfloat*, Mint*, Mint*, Mint*, Mfloat*));

void        PROTO(imsl_stbsv,(Mchar*, Mchar*, Mchar*,

			      Mint*, Mint*, Mfloat*, Mint*, Mfloat*, Mint*));

Mfloat      PROTO(imsl_b22dr,(Mint*, Mint*, Mfloat*, Mfloat*,

			      Mint*, Mint*, Mfloat[],

			      Mfloat[], Mint*, Mint*, Mfloat*, Mfloat[]));

void        PROTO(imsl_b22gd, (Mint*, Mint*, Mint*, Mfloat[],

			      Mint*, Mfloat[], Mint*, Mint*,

			      Mfloat[], Mfloat[],  Mint*,

			      Mint*, Mfloat*, Mfloat*, Mint*,

			      Mint[], Mint[], Mfloat*, Mfloat*,

			      Mfloat*, Mfloat*, Mfloat*, Mfloat*));

void        PROTO(imsl_b32dr,(Mchar*, Mint*, Mint*));

Mfloat      PROTO(imsl_b22ig,(Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*, Mint*,

			      Mfloat[], Mfloat[], Mint*,

			      Mint*, Mfloat*, Mfloat[]));

Mfloat      PROTO(imsl_b32ig,(Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*, Mint*,

			      Mfloat[], Mfloat[], Mint*,

			      Mint*, Mfloat*, Mfloat[], Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_b22in,(Mint*, Mfloat[], Mint*, Mfloat[], Mfloat*, Mint*,

			      Mint*, Mint*, Mfloat[], Mfloat[], Mfloat[],

			      Mfloat[], Mint[]));

void        PROTO(imsl_b32in,(Mchar*, Mint*, Mfloat[], Mint*));

void        PROTO(imsl_b42in,(Mchar*, Mint*, Mint*, Mfloat[], Mfloat*,

			      Mint*, Mint*, Mfloat[], Mfloat*, Mfloat*,

			      Mint*, Mfloat[], Mfloat[], Mfloat[],

			      Mint[]));

void        PROTO(imsl_b2ls2,(Mint*, Mfloat[], Mint*, Mfloat[],

			      Mfloat[], Mint*, Mint*, Mint*,

			      Mfloat[], Mfloat[], Mint*, Mint*,

			      Mfloat[], Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_b3ls2,(Mint*, Mfloat[], Mint*, Mfloat[],

			      Mint*, Mfloat[], Mchar*, Mint,

			      Mint*));

void        PROTO(imsl_b4ls2,(Mint*, Mfloat[], Mint*, Mfloat[],

			      Mfloat*, Mint*, Mint*, Mint*,

			      Mfloat[], Mfloat[], Mint*, Mint*,

			      Mfloat[], Mfloat[], Mfloat*, Mfloat*,

			      Mfloat[], Mfloat*, Mfloat*, Mfloat[]));

void        PROTO(imsl_b5ls2,(Mint*, Mfloat[], Mfloat[], Mfloat[],

			      Mint*, Mfloat[], Mint*, Mfloat[],

			      Mfloat*, Mfloat[], Mint*));

void        PROTO(imsl_b6ls2,(Mint*, Mint*, Mint*, Mfloat*,

			      Mfloat*, Mfloat*));

void        PROTO(imsl_b7ls2,(Mint*));

void        PROTO(imsl_b2lsq,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mint*,

			      Mfloat[], Mint*, Mfloat[], Mfloat[], Mfloat[],

			      Mfloat[], Mfloat[], Mint[]));

void        PROTO(imsl_b3lsq,(Mint*, Mint*, Mfloat[], Mfloat[], Mint *));

void        PROTO(imsl_b4lsq,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mint*,

			      Mfloat[], Mint*, Mfloat[], Mfloat*, Mfloat[]));

void        PROTO(imsl_b5lsq,(Mfloat*, Mint*, Mint*));

void        PROTO(imsl_b6lsq,(Mfloat*, Mint*, Mint*, Mfloat[]));

void        PROTO(imsl_b2vls,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mint*, Mint*,

			      Mfloat[], Mfloat[], Mfloat[], Mfloat*, Mint[], Mfloat[]));

void        PROTO(imsl_b3vls,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mint*, Mint*,

			      Mfloat[], Mfloat[], Mfloat[], Mfloat*, Mfloat*, Mfloat[],

			      Mfloat*, Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_b4vls,(Mint*, Mint*, Mfloat[], Mfloat[], Mfloat[], Mint*,

			      Mfloat[], Mint*, Mfloat[], Mfloat[], Mfloat[], Mfloat[],

			      Mfloat[], Mfloat[]));

void        PROTO(imsl_b5vls,(Mfloat[], Mfloat*, Mint*, Mint*, Mfloat[], Mint*, Mfloat*));

Mfloat      PROTO(imsl_b6vls,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mint*, Mfloat[],

			      Mint*, Mfloat[], Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_b7vls,(Mint*, Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*, Mfloat*));

void        PROTO(imsl_b8vls,(Mfloat*, Mint*, Mint*, Mfloat[]));

void        PROTO(imsl_b2cpp,(Mint*, Mfloat[], Mint*, Mfloat[], Mint*, Mfloat[],

			      Mfloat[], Mfloat[]));

void        PROTO(imsl_b3cpp,(Mint*, Mfloat[], Mint*, Mfloat[], Mint*, Mfloat[],

			      Mfloat*, Mfloat[], Mfloat[], Mfloat[], Mfloat*));

void        PROTO(imsl_f2lsq,(Mfloat (*f) (Mint, Mfloat), Mint*, Mint*, Mint*,

			      Mfloat[], Mfloat[], Mint*,

			      Mfloat[], Mfloat[], Mfloat*, Mfloat[]));

void        PROTO(imsl_c2scv,(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mint *iequal, Mfloat break_[], Mfloat *cscoef,

			      Mfloat weight[], Mfloat *wk, Mfloat ywk[], Mint iwk[]));

void        PROTO(imsl_c3scv,(Mfloat xdata[], Mfloat *avh, Mint iwk[], Mfloat dfcopy[], Mfloat *avdy,

			      Mint *ndata, Mfloat ydata[],Mfloat *cscoef, Mfloat *r, Mfloat *t));

void        PROTO(imsl_c4scv,(Mfloat xdata[], Mfloat *avh, Mfloat dfcopy[], Mint *ndata,

			      Mfloat *rho, Mfloat *p, Mfloat *q, Mfloat *fun,

			      Mfloat *var, Mfloat stat[], Mfloat ydata[], Mfloat *cscoef,

			      Mfloat *r, Mfloat *t, Mfloat u[], Mfloat v[]));

void        PROTO(imsl_c5scv,(Mfloat xdata[], Mfloat *avh, Mfloat dfcopy[], Mint *ndata, Mfloat *p,

			      Mfloat *q, Mfloat ydata[], Mfloat *cscoef, Mfloat u[], Mfloat v[]));

void        PROTO(imsl_c2smh,(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat weight[], Mfloat *smpar, Mfloat break_[],

			      Mfloat *cscoef, Mfloat wk[], Mint iwk[]));

void        PROTO(imsl_c3smh,(Mint *ndata, Mfloat fdata[], Mfloat weight[], Mfloat *smpar, Mfloat break_[], Mfloat *cscoef,

			      Mfloat r[], Mfloat r1[], Mfloat r2[], Mfloat t[], Mfloat u[], Mfloat v[], Mint iwk[]));

Mfloat      PROTO(imsl_c4smh,(Mint k, Mfloat x));

void        PROTO(imsl_b21gd,(Mint*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mint*, Mfloat[], Mfloat[],

			      Mfloat*, Mfloat[], Mint[],

			      Mfloat[], Mfloat[], Mfloat[]));

#if 0

void        PROTO(imsl_c21gd,(Mint*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat*, Mfloat[],

			      Mint[], Mfloat[], Mfloat[]));

#endif

void        PROTO(imsl_c21gd,(Mint*, Mint*, Mfloat[], Mint*,

			      Mfloat[], Mfloat**, Mfloat[],

			      Mint[], Mfloat[], Mfloat[]));



	/* Chapter 4 --- Intergration and Differentiation */

void        PROTO(imsl_q3awo,(Mfloat (*) (Mfloat), Mfloat*, Mfloat*,

				Mfloat*, Mint*, Mfloat*,

				Mfloat*, Mint*, Mint*,

				Mint*, Mfloat*, Mfloat*,

				Mint*, Mint*, Mfloat[],

				Mfloat[], Mfloat[], Mfloat[],

				Mint[], Mint[], Mint*, Mint*,

				Mfloat*));

void        PROTO(imsl_q4awo,(Mint*, Mfloat[], Mfloat*, Mfloat*,

				Mfloat[], Mint*));

void        PROTO(imsl_q7awo,(Mfloat[], Mfloat[], Mfloat[],

				Mfloat[]));

void        PROTO(imsl_q8awo,(Mfloat (*)(Mfloat), Mfloat (*)(Mfloat*,

				Mfloat*, Mfloat*, Mfloat*, Mfloat*,

				Mint*), Mfloat*, Mfloat*, Mfloat*,

				Mfloat*, Mint*, Mfloat*, Mfloat*,

				Mfloat*, Mfloat*, Mfloat*, Mfloat*));



void        PROTO(imsl_q9ag,(Mfloat (*)(Mfloat), Mfloat*, Mfloat*,

				Mfloat*, Mfloat*, Mfloat*, Mfloat*));

void        PROTO(imsl_q10g,(Mint*, Mint*, Mint*, Mfloat*, Mfloat[],

				Mint[], Mint*));

void        PROTO(imsl_q4ng,(Mfloat*, Mfloat*, Mfloat*));



void        PROTO(imsl_q3agi,(Mfloat (*) (Mfloat), Mfloat*, Mint*,

				Mfloat*, Mfloat*, Mint*,

				Mfloat*, Mfloat*, Mint*,

				Mint*, Mfloat[], Mfloat[],

				Mfloat[], Mfloat[], Mint[],

				Mint*));



void        PROTO(imsl_g2rcf,(Mint*, Mfloat[], Mfloat[], Mint*, Mfloat[],

			      Mfloat[], Mfloat[], Mfloat[]));

Mfloat      PROTO(imsl_g3rcf,(Mfloat*, Mint*, Mfloat[], Mfloat[]));

void        PROTO(imsl_g4rcf,(Mint*, Mfloat[], Mfloat[], Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_g2rul,(Mint*, Mint*, Mfloat*, Mfloat*,

			      Mint*, Mfloat[], Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_reccf,(Mint*, Mint*, Mfloat*, Mfloat*, Mfloat[], Mfloat[]));



	/* Chapter 5 --- Differential Equations */



	/* Chapter 6 --- Transforms */



void        PROTO(imsl_f2trf,(Mint*, Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_f2trb,(Mint*, Mfloat*, Mfloat*, Mfloat*));



	/* Chapter 7 --- Nonlinear Equations */



	/* Chapter 8 --- Optimization */

void        PROTO(imsl_fdgrd,(Mfloat (*fcn) (Mint, Mfloat[], Mfloat*), Mint,

			      Mfloat[], Mfloat[], Mfloat*, Mfloat*, Mfloat[]));

void        PROTO(imsl_cdgrd,(Mfloat (*fcn) (Mint, Mfloat[], Mfloat*), Mint,

			      Mfloat[], Mfloat[], Mfloat*, Mfloat[]));

void        PROTO(imsl_d2prs,(Mint, Mint, Mfloat*, Mint, Mfloat*,

		      Mfloat*, Mfloat*, Mint*, Mfloat*, Mfloat*,

		      Mint, Mint*, Mint*, Mfloat*, Mfloat*, Mfloat*,

		      Mfloat*, Mint*, Mint));

void        PROTO(imsl_q2rog, (Mint, Mint, Mint, Mfloat*, Mint, Mfloat*,

			Mfloat*, Mfloat*, Mint, Mfloat*, Mfloat*,

			Mint*, Mint*, Mfloat*, Mfloat*));



	/* Chapter 9 --- Basic Matrix/Vector Operations */

	    /* Matrix copy */

void        PROTO(imsl_crgrg,(Mint, Mfloat*, Mint, Mfloat*, Mint));

void        PROTO(imsl_csfrg,(Mint*, Mfloat*, Mint*));

void        PROTO(imsl_i_m1ran,(Mint, Mint, Mint*, Mint*));

void        PROTO(imsl_f_m1ran,(Mint, Mint, Mfloat*, Mfloat*));

void        PROTO(imsl_c_m1ran,(Mint, Mint, Mf_complex*, Mf_complex*));

void        PROTO(imsl_trnrr,(Mint, Mint, Mfloat*, Mint, Mint, Mint, Mfloat*, Mint));

void        PROTO(imsl_trncr,(Mint, Mint, Mf_complex*, Mint, Mint, Mint, Mf_complex*, Mint));

void        PROTO(imsl_svrgp,(Mint, Mfloat[], Mfloat[], Mint[]));

void        PROTO(imsl_svrgn,(Mint, Mfloat*, Mfloat*));

void        PROTO(imsl_svrbn,(Mint*, Mfloat*, Mfloat*));

void        PROTO(imsl_prime,(Mint, Mint*, Mint*, Mint*, Mint*));

		/* BLAS 1 */

Mint        PROTO(imsl_isamax,(Mint, Mfloat*, Mint));

void        PROTO(imsl_scopy,(Mint, Mfloat[], Mint, Mfloat[], Mint));

Mfloat      PROTO(imsl_sdot,(Mint, Mfloat*, Mint, Mfloat*, Mint));

Mfloat      PROTO(imsl_sxyz,(Mint, Mfloat[], Mint, Mfloat[], Mint, Mfloat[], Mint));

Mfloat      PROTO(imsl_sasum,(Mint, Mfloat*, Mint));

Mfloat      PROTO(imsl_snrm2,(Mint n, Mfloat *sx, Mint incx));

Mfloat      PROTO(imsl_scnrm2,(Mint*, Mf_complex*, Mint*));

void        PROTO(imsl_svcal,(Mint, Mfloat, Mfloat[], Mint, Mfloat[], Mint));

void        PROTO(imsl_srotm,(Mint,Mfloat[],Mint,Mfloat[],Mint,Mfloat[]));

void        PROTO(imsl_srotmg,(Mfloat*,Mfloat*,Mfloat*,Mfloat*,Mfloat[]));

void        PROTO(imsl_srot,(Mint, Mfloat[], Mint, Mfloat[], Mint, Mfloat,

			     Mfloat));

void        PROTO(imsl_srotg,(Mfloat*, Mfloat*, Mfloat*, Mfloat*));

Mf_complex  PROTO(imsl_cdotc,(Mint*, Mf_complex*, Mint*, Mf_complex*,

			Mint*));

Mf_complex  PROTO(imsl_cdotu,(Mint*, Mf_complex*, Mint*, Mf_complex*,

			Mint*));

Mfloat      PROTO(imsl_scasum,(Mint*, Mf_complex*, Mint*));

void        PROTO(imsl_cadd,(Mint*, Mf_complex*, Mf_complex*, Mint*));

void        PROTO(imsl_ccgcg,(Mint*, Mf_complex*, Mint*, Mf_complex*,

			Mint*));

Mint        PROTO(imsl_icamax,(Mint*, Mf_complex*, Mint*));

void        PROTO(imsl_cset,(Mint*, Mf_complex*, Mf_complex*, Mint*));

void        PROTO(imsl_cgemv,(Mchar*, unsigned, Mint*, Mint*,

			Mf_complex*, Mf_complex*, Mint*, Mf_complex*,

			Mint*, Mf_complex*, Mf_complex*, Mint*));

void        PROTO(imsl_caxpy,(Mint*, Mf_complex*, Mf_complex*, Mint*,

			Mf_complex*, Mint*));

void        PROTO(imsl_ccbcb, (Mint*, Mf_complex*, Mint*, Mint*,

		Mint*, Mf_complex*, Mint*, Mint*, Mint*));

void        PROTO(imsl_ctbsv, (char*, Mint, char*, Mint,

		char*, Mint, Mint*, Mint*, Mf_complex*,

		Mint*, Mf_complex[], Mint*));

void        PROTO (imsl_ctrsv, (Mchar*, unsigned,

				 Mchar*, unsigned,

				 Mchar*, unsigned,

				 Mint*, Mf_complex*, Mint*,

				 Mf_complex*, Mint*));

		/* BLAS 2 */

void        PROTO(imsl_strsv,(Mchar*, Mchar*, Mchar*, Mint, Mfloat*, Mint, Mfloat*, Mint));

void        PROTO(imsl_sgemv,(Mchar*, Mint, Mint*, Mint*, Mfloat*, Mfloat*,

			Mint*, Mfloat*, Mint*, Mfloat*, Mfloat*, Mint*));

void        PROTO(imsl_sger,(Mint, Mint, Mfloat, Mfloat*, Mint, Mfloat*,

			Mint, Mfloat*, Mint));

void        PROTO(imsl_cgerc, (Mint*, Mint*, Mf_complex*, Mf_complex[],

			       Mint*, Mf_complex[], Mint*, Mf_complex[], Mint*));

void        PROTO(imsl_ssymv,(Mchar*, Mint, Mint*, Mfloat*, Mfloat*,

		Mint*, Mfloat*, Mint*, Mfloat*, Mfloat*, Mint*));

void        PROTO(imsl_ssyr2,(Mchar*, Mint, Mint*, Mfloat*, Mfloat*,

		Mint*, Mfloat*, Mint*, Mfloat*, Mint*));

void        PROTO(imsl_ccopy,(Mint*, Mf_complex*, Mint*, Mf_complex*, Mint*));

void        PROTO(imsl_cscal,(Mint*, Mf_complex*, Mf_complex*, Mint*));

void        PROTO(imsl_csscal,(Mint*, Mfloat*, Mf_complex*, Mint*));

void        PROTO(imsl_cswap,(Mint*, Mf_complex*, Mint*, Mf_complex*, Mint*));



	/* Chapter 10 --- Statisitics */

void        PROTO(imsl_r2ivn,(Mint, Mint, Mint, Mfloat*, Mint, Mint, Mint, Mint[],

			Mint, Mint[], Mint, Mint, Mint, Mfloat, Mfloat*, Mint,

			Mfloat*, Mint, Mfloat[], Mint*, Mfloat*, Mfloat*,

			Mint, Mint*, Mfloat[], Mfloat[], Mfloat[]));

void        PROTO(imsl_c1wfr,(Mint, Mint, Mfloat*, Mint, Mint, Mint, Mint, Mint, Mfloat, Mint*, Mfloat*, Mfloat*, Mint*));

void        PROTO(imsl_rcovb,(Mint, Mfloat*, Mint, Mfloat, Mfloat*, Mint));



		/* Random Number Routines */

void        PROTO(imsl_rnun,(Mint, Mfloat[]));

void        PROTO(imsl_rnses,(Mfloat[]));

void        PROTO(imsl_r1int,(Mint));

void        PROTO(imsl_r1clk,(Mint*));



	/* Chapter ?? --- Special Functions */

Mfloat      PROTO(imsl_csevl,(Mfloat, Mfloat*, Mint));

Mdouble     PROTO(imsl_dcsevl,(Mdouble, Mdouble*, Mint));

Mint        PROTO(imsl_inits,(Mfloat*, Mint, Mfloat));

Mfloat      PROTO(imsl_r9lgmc,(Mfloat));

Mfloat      PROTO(imsl_bsi0e,(Mfloat));

Mfloat      PROTO(imsl_bsi1e,(Mfloat));

Mfloat      PROTO(imsl_bsk0e,(Mfloat));

Mfloat      PROTO(imsl_bsk1e,(Mfloat));

Mfloat      PROTO(imsl_betin,(Mfloat,Mfloat,Mfloat));

void        PROTO(imsl_c3is,(Mf_complex*, Mfloat*, Mint*, Mf_complex[],

		Mf_complex[], Mf_complex[], Mf_complex[],

		Mint*));



	/* Chapter 11 --- Utilities */

void        PROTO(imsl_umach,(Mint, FILE**));

Mint        PROTO(imsl_imach,(Mint));

	    /* Error checking routines */

void        PROTO(imsl_c12ile,(Mint, Mchar*, Mint, Mchar*, Mint*));

void        PROTO(imsl_c12ile,(Mint, Mchar*, Mint, Mchar*, Mint*));

void        PROTO(imsl_c1dim,(Mint, Mint, Mchar*, Mint, Mchar*, Mint*));

void        PROTO(imsl_c1iarg,(Mint, Mchar*, Mint, Mint, Mint*));

void        PROTO(imsl_c1ind,(Mint, Mint, Mchar*, Mint, Mchar*, Mint*));

void        PROTO(imsl_c1ge0,(Mfloat, Mchar*, Mint*));

void        PROTO(imsl_c1r,(Mint,Mfloat[],Mint,Mint*));

Mfloat      PROTO(imsl_a1ot,(Mint, Mfloat[], Mint, Mfloat[], Mint));

void        PROTO(imsl_c1div,(Mfloat, Mfloat, Mfloat*));

void        PROTO(imsl_g1aov,(Mfloat, Mfloat, Mfloat, Mfloat, Mfloat, Mfloat[]));

Mint        PROTO(imsl_i_vmax,(Mint, ...));



		/* Writing Routines */

void        PROTO(imsl_i_wrirl,(Mchar*,Mint,Mint,Mint[],Mint,Mint,Mchar*,Mchar**,

			Mchar**,Mint,Mint));

void        PROTO(imsl_f_wrrrl,(Mchar*,Mint,Mint,Mfloat[],Mint,Mint,Mchar*,Mchar**,

			Mchar**,Mint,Mint));

void        PROTO(imsl_c_wrcrl,(Mchar*,Mint,Mint,Mf_complex[],Mint,Mint,Mchar*,Mchar**,

			Mchar**,Mint,Mint));



void        PROTO(imsl_wrtrl,(Mchar*, Mint));

void        PROTO(imsl_null_pointer,(Mchar*,Mint,void*));

void        PROTO(imsl_write_controller,(Mint*, Mint, Mint, Mint, Mchar**,Mchar**,

			Mint, Mint, Mchar*, Mchar**, Mint*, Mint, Mint*,Mint*,

			Mint*, Mint*, Mint*, Mint, Mint, Mint, Mint,Mint,Mint*,

			Mint, Mint, Mint, Mchar*, Mint*, Mint));

Mint        PROTO(imsl_write_initialize,(Mint*, Mint*, Mint*, Mint*, Mint*,

			Mint*, Mchar**, Mchar**,Mint*,Mint*,Mint,Mint*,Mchar*,

			Mint*, Mint*, Mint, Mint,Mint*,Mint*,Mint*,Mint,Mint*));

void        PROTO(imsl_write_title,(Mchar*, Mint, Mint*, Mchar*, Mint));

void        PROTO(imsl_write_labels,(Mint, Mchar**,Mchar*,Mchar*,Mint*,Mint,Mint,

			Mint, Mint, Mint, Mint,Mint, Mint, Mchar*, Mint, Mint));

Mchar      *PROTO(imsl_write_conversion,(Mchar*, Mint*, Mchar*, Mchar*, Mint*));

void        PROTO(imsl_c1nter,(Mint, Mint*, Mchar*));

void        PROTO(imsl_w5rrl_f,(Mint, Mint, Mchar**, Mchar*, Mint*,Mint*,Mchar**,

			Mint*, Mint));

void        PROTO(imsl_w6rrl,(Mchar*, Mint, Mchar*, Mchar*, Mint*));

Mchar      *PROTO(imsl_w7rrl,(Mint, Mchar*));

void        PROTO(imsl_w8rrl,(Mint, Mint, Mchar**, Mint, Mint, Mint, Mint, Mint,

			Mint*, Mchar*, Mint*));

void        PROTO(imsl_w12rl,(Mint, Mint, Mint, Mint, Mint*, Mint*));

void        PROTO(imsl_w1opt,(Mint, Mint*));

void        PROTO(imsl_write_line,(Mint, Mchar*));

Mchar      *PROTO(imsl_w1iss,(Mfloat*, Mchar*, Mint));

Mchar      *PROTO(imsl_fmtx,(Mfloat*, Mint, Mint));

void        PROTO(imsl_write_format,(Mint, Mint, Mint, Mint, Mchar*, Mint*,Mchar*,

			Mchar*, Mint*, Mint*));



	/* Error Handler and Workspace */

#ifdef USE_IMSL_MALLOC

void       *PROTO(imsl_calloc,(Mint,Mint));

void       *PROTO(imsl_malloc,(Mint));

void       *PROTO(imsl_realloc,(void*,Mint));

void        PROTO(imsl_free,(void*));

#elif defined(USE_ALIB_MALLOC)

#if !defined(_MEMORY_H)
#    include <memory.h>   /* memset() */
#endif


void *PROTO(imsl_alib_malloc, (long)); 
void  PROTO(imsl_alib_free, (void *)); 

#if defined(imsl_malloc)
#   undef   imsl_malloc
#endif
#define  imsl_malloc   imsl_alib_malloc

#if defined(imsl_free)
#   undef   imsl_free
#endif
#define  imsl_free     imsl_alib_free

#if defined(imsl_calloc)
#  undef imsl_calloc
#endif
#define imsl_calloc    calloc

#if defined(imsl_realloc)
#  undef imsl_realloc
#endif
#define imsl_realloc   realloc


/* USE_ALIB_MALLOC */                                 
#else

# ifdef imsl_calloc

#  undef imsl_calloc

# endif

#define imsl_calloc     calloc

# ifdef imsl_malloc

#  undef imsl_malloc

# endif

#define imsl_malloc     malloc

# ifdef imsl_realloc

#  undef imsl_realloc

# endif

#define imsl_realloc    realloc

# ifdef imsl_free

#  undef imsl_free

# endif

#define imsl_free       free

#endif

void        PROTO(imsl_signal,());

void        PROTO(imsl_set_signal,(Mint));

void        PROTO(imsl_e1psh,(Mchar*));

void        PROTO(imsl_e1pop,(Mchar*));

void        PROTO(imsl_e1mes,(Mint,Mlong,Mchar*));

void        PROTO(imsl_ermes,(Mint type, Mlong code));

void        PROTO(imsl_ercode,(Mlong code));

void        PROTO(imsl_e1usr,(Mchar*));

void        PROTO(imsl_e1pos,(Mint,Mint*,Mint*));

Mlong       PROTO(imsl_n1rcd,(Mint));

Mint        PROTO(imsl_n1rnof,(Mint));

Mchar      *PROTO(imsl_find_name,(Mlong));

Mchar      *PROTO(imsl_find_message,(Mlong code));

#ifdef IMSL_MACHINE_NT

void        PROTO(imsl_create_console,(void));

#endif



Mint        PROTO(imsl_l1ame,(Mchar*, Mint, Mchar*, Mint));



#endif  /* IMSL_INC_H */

