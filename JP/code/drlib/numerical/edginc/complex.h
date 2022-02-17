/* complex.h - definitions to support FORTRAN translations to C 

			made by FOR_C (TM) from COBALT BLUE 

         copyright Lightfoot & Associates, Inc., 1988

                ALL RIGHTS RESERVED

 */



#ifndef COMPLEX_H

#define COMPLEX_H	/* only define these objects once per file! */



/* Define the macro switch:

		INLN_COMPLEX

	if you wish to use macros expanding to inline code (instead of

	functions) for the complex conversion functions.

	Beware of duplicate names!

 */



/* NOTE: f_rt.h should be included in front of this file so that any

	configuration macro switches are properly defined! */



#ifdef imsl_c_abs

#undef imsl_c_abs

#endif



typedef struct {

	float re;

	float im;

	}	f_complex;

#define COMPLEX Mf_complex



typedef struct {

	double re;

	double im;

	}	d_complex;

#define COMPLEX16 d_complex



/* access to the real & imaginary components of the complex nos. */

#define ZRE(z)	(z).re

#define ZIM(z)	(z).im



	/* real & imag. components for I/O purposes */

#define ZADR(z)		&(z).re, &(z).img		/* used for input */

#define ZPAIR(z)	(z).re, (z).im			/* used for output */



	/* complex function prototypes */



f_complex PROTO(imsl_c_neg,(f_complex ));

d_complex PROTO(imsl_z_neg,(d_complex ));

f_complex PROTO(imsl_c_add,(f_complex, f_complex ));

d_complex PROTO(imsl_z_add,(d_complex, d_complex ));

f_complex PROTO(imsl_c_sub,(f_complex, f_complex ));

d_complex PROTO(imsl_z_sub,(d_complex, d_complex ));

f_complex PROTO(imsl_c_mul,(f_complex, f_complex ));

d_complex PROTO(imsl_z_mul,(d_complex, d_complex ));

f_complex PROTO(imsl_c_div,(f_complex, f_complex ));

d_complex PROTO(imsl_z_div,(d_complex, d_complex ));

int       PROTO(imsl_c_eq,(f_complex, f_complex ));

int       PROTO(imsl_z_eq,(d_complex, d_complex ));

f_complex PROTO(imsl_cz_convert,(d_complex ));

d_complex PROTO(imsl_zc_convert,(f_complex ));



		/* f_complex Library Functions */



float     PROTO(imsl_c_cftof,( f_complex ));

double	  PROTO(imsl_z_ctof,( d_complex ));

float     PROTO(imsl_c_aimag,( f_complex ));

double    PROTO(imsl_z_aimag,( d_complex ));

float     PROTO(imsl_fc_convert,(f_complex));

double    PROTO(imsl_dz_convert,(d_complex));

f_complex PROTO(imsl_cf_convert,(float, float));

d_complex PROTO(imsl_zd_convert,(double,double));

f_complex PROTO(imsl_c_conjg,(f_complex ));

d_complex PROTO(imsl_z_conjg,(d_complex ));

float     PROTO(imsl_c_arg,(f_complex ));

double    PROTO(imsl_z_arg,(d_complex));

f_complex PROTO(imsl_c_sqrt,(f_complex ));

d_complex PROTO(imsl_z_sqrt,(d_complex ));

f_complex PROTO(imsl_c_log,(f_complex ));

d_complex PROTO(imsl_z_log,(d_complex ));

f_complex PROTO(imsl_c_exp,(f_complex ));

d_complex PROTO(imsl_z_exp,(d_complex ));

f_complex PROTO(imsl_c_sin,(f_complex ));

d_complex PROTO(imsl_z_sin,(d_complex ));

f_complex PROTO(imsl_c_cos,(f_complex ));

d_complex PROTO(imsl_z_cos,(d_complex ));

float     PROTO(imsl_c_abs,(f_complex ));

double    PROTO(imsl_z_abs,(d_complex ));



			/* * * psuedo functions * * */



#define itocf(i,j)	imsl_zd_convert((double)i, (double)j)

#define ltocf(i,j)	imsl_zd_convert((double)i, (double)j)

#define cftoi( z )	(int)imsl_c_cftof(z)

#ifdef COMPUTER_DECOSF

#define cftol( z )	(int)imsl_c_cftof(z)

#else

#define cftol( z )	(long)imsl_c_cftof(z)

#endif



#define itoc(i,j)	imsl_zd_convert((double)i, (double)j)

#define ltoc(i,j)	imsl_zd_convert((double)i, (double)j)

#define ctoi( z )	(int)imsl_z_ctof(z)

#ifdef COMPUTER_DECOSF

#define ctol( z )	(int)imsl_z_ctof(z)

#else

#define ctol( z )	(long)imsl_z_ctof(z)

#endif



#ifdef INLN_COMPLEX	/* then use the following macro versions */

#ifndef 	EXTERNAL

#define 	EXTERNAL extern

#endif

	/* NOTE: EXTERNAL is defined as nothing in f_cmplx.c to define the

			the f_complex externals in that file.  This K&R approach avoids 

			problems with certain older Librarians in use. */



EXTERNAL f_complex _z_;	/* temp for complex macros, (to avoid side effects) */

#define	xxx_imsl_c_aimag(z)	(_z_=z, _z_.im)	/* return float */

#define	xxx_imsl_c_cftof(z)	(_z_=z, _z_.re)	/* return float */

#define	xxx_imsl_c_ftocf(r,i)   ((_z_.re=r, _z_.im=i), _z_)	/* return complex */

#define	xxx_imsl_c_conjg(z)	((_z_=z, _z_.im=-_z_.im), _z_) /* return complex */

#define	xxx_imsl_c_neg(z)	((_z_=z, (_z_.re=-_z_.re,_z_.im=-_z_.im)), _z_)



EXTERNAL d_complex _dz_; /* temp for complex macros, (to avoid side effects) */

#define xxx_imsl_z_aimag(z)	(_dz_=z, _dz_.im)

#define xxx_imsl_z_ctof(z)	(_dz_=z, _dz_.re)

#define xxx_imsl_z_ftoc(r,i)    ((_dz_.re=r, _dz_.im=i), _dz_)

#define xxx_imsl_z_conjg(z)	((_dz_=z, _dz_.im=-_dz_.im), _dz_)

#define xxx_imsl_z_neg(z)	((_dz_=z, (_dz_.re=-_dz_.re,_dz_.im=-_dz_.im)), _dz_)

#endif	/* of INLN_COMPLEX */



#endif	/* on COMPLEX_H */

