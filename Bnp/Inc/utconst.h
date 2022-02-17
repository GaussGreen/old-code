/* ===============================================================
   FILENAME : utconst.h

   PURPOSE:   useful constants
   =============================================================== */

#ifndef UTCONST_H
#define UTCONST_H

/*
#include "math.h"
*/
#include "stdio.h"
#include "stdlib.h"
#include "string.h"    
#include "float.h"    

/*------------------------platform dependant stuff-----------------------*/
#if defined(ALPHA) || defined(unix) || defined(_windows) || defined (_WIN32)

#define EXTERN  extern
#define EXTERND

#else

#define EXTERN  noshare extern
#define EXTERND noshare

#endif

/*------------------------system wide constants--------------------------*/

#define SRT_PI 				3.141592653589793
#define SQRT_TWO            1.414213562373095
#define SQRT_TWO_INV  		0.707106781186548     /**  1/sqrt(2)  **/
#define INV_SQRT_TWO_PI  	0.3989422804014327    /**  1/sqrt(2*Pi)  **/
#define INV_SQRT_PI  	    0.5641895835477563    /**  1/sqrt(Pi)  **/

/* Date constants */
#define	MIN_DATE			30317	/* 01-Jan-1983 */
#define	MAX_DATE			54789	/* 01-Jan-2050 */
#define MAX_NFP				100		/* Max of full periods */

/*----------------------------definitions----------------------------*/

#define EPS    				DBL_EPSILON
#define SRT_EPSILON    		1.0e-10
#define DTOL(x) 			(long)(floor(x+EPS)+ (x < 0 ? -EPS : EPS))

#define	YEARS_IN_DAY		(1.0/365.0)
#define	DAYS_IN_YEAR		365.0

#define SRTBUFSZ			32

#define MEMORY_ERR			-444444.44
#define UNKNOWN_GREEK		-333333.33

/*------------------------system wide types------------------------------*/

#endif
