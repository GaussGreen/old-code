/*********************************************************************************
 * T_CURVE.H 
 * crx london t_curves
 *
 ********************************************************************************/

#ifndef __T_CURVE_H__
#define __T_CURVE_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "l_date.h"

#define  MAX_ZRATE          1000          /* Max nb of elemts in input date arr */

/* Type of interpolation for T_CURVE and Z_CURVE -------------------------------*/
#define     SRM3_LINEAR_INTERP      0L
#define     SRM3_FLATFWD_INTERP     1L
#define     SRM3_LINEAR_DISC_INTERP 2L
    
/*********************************************************************************
 * Curves
 ********************************************************************************/
typedef struct
{
    long interpType;
    void *data;     /* data will be a BASIC_CURVE structure
                       see paryield.c for details                               */
    long today;     
} Z_CURVE;

typedef struct
{
    /* Base date information ---------------------------------------------------*/
    long    Today;                        /* Today's date                       */
    int     SpotDays;                     /* Spot days                          */
    long    ValueDate;                    /* Value date                         */
  
    /* zero curve --------------------------------------------------------------*/
    char    CurveDCC;
    char    CurveFreq;
    int     NbZero;                       /* Number of zeros                    */
    long    ZeroDate[MAX_ZRATE];          /* Zero maturity dates                */
    double  Zero[MAX_ZRATE];              /* Zero rates (Annual ACT/365 basis)  */
  
    /* Interpolation Type ------------------------------------------------------*/
    long InterpType;                      /* defaulted to SRM3_LINEAR_INTERP 
                                                 by all the curve loaders       */
  
} T_CURVE;

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
