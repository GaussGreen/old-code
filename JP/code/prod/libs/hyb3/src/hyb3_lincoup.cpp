/****************************************************************************/
/*      Bootstrapping of yield curve according to the linear coupon credo   */
/****************************************************************************/
/*      LINCOUP.c                                                           */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cupslib.h"



/*****  Yldconv  ************************************************************/
/*
*       Convert Money Market rates into Annual ACT/365 zero rates.
*/
int     Yldconv (       double  *Zero,                                          /* (O) Zeros at money market points 	        */
                        long    *ZeroDate,                                      /* (O) zero dates 	                        */
                        int     NbMMPoint,                                      /* (I) Number of money market points 	        */
                        int     *MMPoint,                                       /* (I) Money Market points in months 	        */
                        double  *MMYield,                                       /* (I) Money Market rates 	                */
                        long    ValueDate,                                      /* (I) Value date for the benchmark instruments */
                        int     MMB)                                            /* (I) Money Market basis (360 or 365) 	        */
{
        long
                Period;                                                         /* Number of days in the money market period */
        int
                i,
                status = FAILURE;                                               /* Error status = FAILURE initially */


        Zero[0]   = pow (1. + MMYield[0] / 100. / (double) MMB, 365.) - 1.;     /* O/N Money Market rate */
        ZeroDate[0] = Nxtday (ValueDate, (long) 1);

        for (i = 1; i < NbMMPoint; i++)                                         /* Other money market points */
        {
                ZeroDate[i] = Nxtmth (ValueDate, MMPoint[i], (long) 1);
                
                Period = Daysact (ValueDate, ZeroDate[i]);
                
                Zero[i] = pow (1. + MMYield[i] / 100. * Period / (double) MMB, 365. / (double) Period) - 1.;
                                        
        }  /* for i */

        status = SUCCESS;

        return (status);

}  /* Yldconv */



/*****  Coupint  ************************************************************/
/*
*       Linear interpolation of the swap benchmark yields by 6-month steps.
*/
int     Coupint (       double  *SwapCurve,                                     /* (O) Array of interpolated and benchmark swap yields */
                        int     NbSwapPoint,                                    /* (I) Number of swap points 	                       */
                        int     *SwapPoint,                                     /* (I) Swap points in years 	                       */
                        double  *SwapYield)                                     /* (I) Input swap yield                                */
{
        double
                yield;                                                          /* Current interpolated swap yield 	                 */
        int	                                                                /*	                                                 */
                term,                                                           /* Current interpolated swap yield maturity in 1/2 years */
                i,                                                              /*	                                                 */
                status = FAILURE;                                               /* Error status = FAILURE initially 	                 */


        for (i = 0; i < NbSwapPoint; i++)                                       /* We fill the swap yield array with the benchmark yields */
                SwapCurve [2*SwapPoint[i]] = SwapYield[i];


        term = 5;                                                               /* First point to interpolate: 2.5 year. It is assumed */
                                                                                /* the first swap is a 2 year swap.                    */       
        for (i = 0; i < NbSwapPoint - 1; i++)                                       
        {
                do
                {
                        linterp (       (double) term,
                                        &yield,
                                        2. * SwapPoint[i],
                                        2. * SwapPoint[i+1],
                                        SwapYield[i],
                                        SwapYield[i+1]);

                        SwapCurve[term] = yield;

                        term++;

                } while (term < 2 * SwapPoint[i+1]);

                term++;                                                         /* Skip the benchmark point */

        }  /* for i */


        for ( ;term <= 2 * MAXMAT; term++)	                                /* If the last swap has a maturity shorter than MAXMAT years */
                SwapCurve[term] = SwapYield[NbSwapPoint - 1];	                /* we use a flat curve after this last swap maturity.        */
                

        status = SUCCESS;

        return (status);

}  /* Coupint */



/*****  Spotrate  ***********************************************************/
/*
*       Bootstrapping of the swap yield curve. Fill the array Zero and
*       Period with the zeros and their maturity at 6 month interval. The
*       zero rates are output with an Annual ACT/365 basis.
*/
int     Spotrate (      double  *Zero,                                          /* (O) Zero rates at 6m interval                           */
                        long    *ZeroDate,                                      /* (O) Zero dates 	                                   */
                        double  *SwapCurve,                                     /* (I) Array of interpolated and benchmark swap yields 	   */
                        int     NbMMPoint,                                      /* (I) Number of money market points 	                   */
                        char    AoS,                                            /* (I) 'A'nnual or 'S'emi-Annual frequency 	           */
                        long    ValueDate,                                      /* (I) Value date for the benchmark instruments            */
                        char    *YBase)                                         /* (I) Year basis for the benchmark swaps ("ACT" or "365") */
{
        double
                *Zero6m,                                                        /* Zero at 6m interval 	                                                */
                *ZeroPeriod,                                                    /* Corresponding maturity in ACT/365 from value date 	                */
                *CouponPeriod,                                                  /* Length of swap coupon period according to the swap convention        */
                LastZero,                                                       /* Zero maturing on last coupon date of current benchmard swap 	        */
                Annuity;                                                        /* Price of the annuity, i.e. sum of the current benchmark swap coupons */
        long	                                                                /* 	                                                                */
                *Zero6mDate,                                                    /* Maturity of the zeros                                                */
                Date,                                                           /* Maturity of the current zero                                         */
                StartDate,                                                      /* Start of coupon period 	                                        */
                EndDate;                                                        /* End of coupon period                                                 */
        int                                                                     /* 	                                                                */
                F,                                                              /* Frequency of swap payment as an integer 	                        */
                Step,                                                           /* Index step 	                                                        */
                i, j, k,                                                        /* 	                                                                */
                status = FAILURE;                                               /* Error status = FAILURE initially 	                                */
        char	                                                                /* 	                                                                */
                ErrorMsg[MAXBUFF];                                              /* Error message to be sent to nrerror */
                

        Zero6m       = dvector (0, 2*MAXMAT);                                           
        Zero6mDate   = lvector (0, 2*MAXMAT);                                           
        ZeroPeriod   = dvector (0, 2*MAXMAT);
        CouponPeriod = dvector (0, 2*MAXMAT);                                           
                
        
        F    = Conv_AoS (AoS);
        Step = (F == 1) ? 2 : 1;
                
        for (i = 1; i <= 2; i++)                                                /* Find the 6m and 1yr zero rates */
        {        
                Date = Nxtmth (ValueDate, (long) i * 6, (long) 1);
                
                k = NbMMPoint - 1;                                              /* Find the 6m and 1yr money market input */
                
                while ((ZeroDate[k] > Date) && (k > 0))
                        k--;
                
                if (ZeroDate[k] != Date)                                        /* In case not 6m or 1yr money market was input */
                {
                        nrerror ("Could not find 6m / 1yr money market rate! (Spotrate)");
                        goto FREE_MEM_AND_RETURN;
        
                }  /* if */                                                
                                                
                Zero6m[i]     = Zero[k];
                Zero6mDate[i] = ZeroDate[k];
                ZeroPeriod[i] = Daysact (ValueDate, Date) / 365.;

        }  /* for i */                                          

                        
        for (i = 1; i <= 2 * MAXMAT; i++)                                       /* Calculate swap coupon period according to the input convention */
        {
                StartDate = Nxtmth (ValueDate, (long) i * 6 - 12 / F, (long) 1);
                EndDate   = Nxtmth (ValueDate, (long) i * 6, (long) 1);
                
                if (YBase[2] == '5')                                            /* ACT/365 convention */
                        CouponPeriod[i] = Daysact (StartDate, EndDate) / 365.;
                else if (YBase[2] == '0')                                       /* ACT/360 convention */                        
                        CouponPeriod[i] = Daysact (StartDate, EndDate) / 360.;
                else                                                            /* ACT/ACT convention */
                        CouponPeriod[i] = 1. / F;
        
        }  /* for i */        
        

        if (AoS == 'S')                                                         /* For semi-annual curve we have to interpolate the 18m swap rate */
        {
                for (i = 1, Annuity = 0.; i <= 2; i++)
                        Annuity += CouponPeriod[i] / pow (1. + Zero6m[i], ZeroPeriod[i]);
        
                LastZero = 1. / pow (1. + Zero6m[2], ZeroPeriod[2]);
                        
                SwapCurve[3] = (100. * (1. - LastZero) / Annuity + SwapCurve[4]) / 2.;  /* 18m rate is the linear interpolation of 1yr and 2yr rate */
                
        }  /* if */        


        for (i = Step + 2; i <= 2 * MAXMAT; i += Step)                          /* If F = 1 we bootstrap annual swaps, if F =2, semi-annual swaps */
        {
                for (j = Step, Annuity = 0.; j < i; j += Step)                                                             
                        Annuity += SwapCurve[i] * CouponPeriod[j] / pow (1. + Zero6m[j], ZeroPeriod[j]);
                                 
                Date = Nxtmth (ValueDate, (long) i * 6, (long) 1);              /* Maturity date of the current zero */
                                                
                ZeroPeriod[i] = Daysact (ValueDate, Date) / 365.;               
                Zero6mDate[i] = Date;

                Zero6m[i] = pow ((100. + SwapCurve[i] * CouponPeriod[i]) / (100. - Annuity), 1. / ZeroPeriod[i]) -1.;  	/* (1+zero rate)^(days/365) =       */
                                                                                                                        /* (100 + coupon) / (100 - annuity) */
                if (Zero6m[i] < ERROR)	                                        /* Problem in bootstrapping */
                {
                        sprintf (ErrorMsg, "Problem in bootstrapping swap #%d in Spotrate (.)", i);
                        nrerror (ErrorMsg);
                        goto FREE_MEM_AND_RETURN;

                }  /* if */

        }  /* for i */

                
        if (AoS == 'A')                                                         /* For annual curve we have to interpolate zero at 18m, 30m, ... */
        {
                for (i = 3; i < 2 * MAXMAT; i += 2)
                {
                        Date = Nxtmth (ValueDate, (long) i * 6, (long) 1);
                        
                        ZeroPeriod[i] = Daysact (ValueDate, Date) / 365.;       /* Maturity of the current zero */
                        Zero6mDate[i] = Date;

                        linterp (       ZeroPeriod[i],                          /* Interpolation of the zero maturing 6m before the current one */
                                        &(Zero6m[i]),
                                        ZeroPeriod[i-1],
                                        ZeroPeriod[i+1],
                                        Zero6m[i-1],
                                        Zero6m[i+1]);

                }  /* for i */
        }  /* if */


        for (i = 3; i <= 2 * MAXMAT; i++)                                       /* Copy the zero curve back in the original array */
        {                                                                       /* taking into account the NbMMPoint offset.      */
                Zero[i - 3 + NbMMPoint] = Zero6m[i];
                ZeroDate[i - 3 + NbMMPoint] = Zero6mDate[i];
                
        }  /* for i */


        status = SUCCESS;

        FREE_MEM_AND_RETURN:

        free_dvector (CouponPeriod, 0, 2*MAXMAT);
        free_dvector (ZeroPeriod,   0, 2*MAXMAT);
        free_lvector (Zero6mDate,   0, 2*MAXMAT);
        free_dvector (Zero6m,       0, 2*MAXMAT);
        
        return (status);

}  /* Spotrate */



/*****  CoB_Spotrate  *******************************************************/
/*
*       Bootstrapping of the swap yield curve with cost of basis. Assume
*	that the floating leg frequency is the same as the fixed leg which
*	is obviously not true.
*/
int     CoB_Spotrate (  double  *Zero,                                          /* (O) Zero rates at 6m interval                           */
                        long    *ZeroDate,                                      /* (O) Zero dates 	                                   */
                        int     NbDZero,                                        /* (I) Number of zeros in the discount zero curve          */
                        double  *DZero,                                         /* (I) Discount zero rates at 6m interval 	           */
                        long    *DZeroDate,                                     /* (I) Discount zero dates 	                           */
                        double  *SwapCurve,                                     /* (I) Array of interpolated and benchmark swap yields 	   */
                        int     NbMMPoint,                                      /* (I) Number of money market points 	                   */
                        char    AoS,                                            /* (I) 'A'nnual or 'S'emi-Annual frequency 	           */
                        long    ValueDate,                                      /* (I) Value date for the benchmark instruments            */
                        char    *YBase)                                         /* (I) Year basis for the benchmark swaps ("ACT" or "365") */
{
        double
                *Zero6m,                                                        /* Zero rate at 6m interval 	                                        */
                *DZero6m,                                                       /* Discount zero price at 6m interval 	                                */
                *Libor,	                                                        /* Libor at 6m interval	                                                */
                *ZeroPeriod,                                                    /* Corresponding maturity in ACT/365 from value date 	                */
                *CouponPeriod,                                                  /* Length of swap coupon period according to the swap convention        */
                ZeroRate,	                                                /*	                                                                */
                ZeroPrice,                                                      /* Zero maturing on last coupon date of current benchmard swap 	        */
                Annuity;                                                        /* Price of the annuity, i.e. sum of the current benchmark swap coupons */
        long	                                                                /* 	                                                                */
                *Zero6mDate,                                                    /* Maturity of the zeros                                                */
                Date,                                                           /* Maturity of the current zero                                         */
                StartDate,                                                      /* Start of coupon period 	                                        */
                EndDate;                                                        /* End of coupon period                                                 */
        int                                                                     /* 	                                                                */
                F,                                                              /* Frequency of swap payment as an integer 	                        */
                Step,                                                           /* Index step 	                                                        */
                i, j, k,                                                        /* 	                                                                */
                status = FAILURE;                                               /* Error status = FAILURE initially 	                                */
        char	                                                                /* 	                                                                */
                ErrorMsg[MAXBUFF];                                              /* Error message to be sent to nrerror */
                

        Zero6m       = dvector (0, 2*MAXMAT);                                           
        DZero6m      = dvector (0, 2*MAXMAT);                                           
        Libor        = dvector (0, 2*MAXMAT);                                           
        Zero6mDate   = lvector (0, 2*MAXMAT);                                           
        ZeroPeriod   = dvector (0, 2*MAXMAT);
        CouponPeriod = dvector (0, 2*MAXMAT);                                           
                

        Zero6m[0]     = 0.;
        Zero6mDate[0] = ValueDate;
        ZeroPeriod[0] = 0.;
        
                
        F    = Conv_AoS (AoS);
        Step = (F == 1) ? 2 : 1;
                
        for (i = 1; i <= 2; i++)                                                /* Find the 6m and 1yr zero rates */
        {        
                Date = Nxtmth (ValueDate, (long) i * 6, 1L);
                
                k = NbMMPoint - 1;                                              /* Find the 6m and 1yr money market input */
                
                while ((ZeroDate[k] > Date) && (k > 0))
                        k--;
                
                if (ZeroDate[k] != Date)                                        /* In case not 6m or 1yr money market was input */
                {
                        nrerror ("Could not find 6m / 1yr money market rate! (CoB_Spotrate)");
                        goto FREE_MEM_AND_RETURN;
        
                }  /* if */                                                
                                                
                Zero6m[i]     = Zero[k];
                Zero6mDate[i] = ZeroDate[k];
                ZeroPeriod[i] = Daysact (ValueDate, Date) / 365.;

        }  /* for i */                                          


        for (i = Step; i <= 2; i += Step)                          	        /* Calculate 6m and 1y Libor */
        {
                Libor[i]  = pow (1. + Zero6m[i], ZeroPeriod[i]);
                Libor[i] /= pow (1. + Zero6m[i-Step], ZeroPeriod[i-Step]);
                Libor[i] -= 1.;
                Libor[i] *= 100.;	                                        /* Libor payment in % */
                                                                       
        }  /* for i */                                                         

                        
        for (i = 1; i <= 2 * MAXMAT; i++)                                       /* Calculate coupon period according to payment convention */
        {
                StartDate = Nxtmth (ValueDate, (long) i * 6 - 12 / F, 1L);
                EndDate   = Nxtmth (ValueDate, (long) i * 6, 1L);
                
                if (YBase[2] == '5')                                            /* ACT/365 convention */
                        CouponPeriod[i] = Daysact (StartDate, EndDate) / 365.;
                else if (YBase[2] == '0')                                       /* ACT/360 convention */                        
                        CouponPeriod[i] = Daysact (StartDate, EndDate) / 360.;
                else                                                            /* ACT/ACT convention */
                        CouponPeriod[i] = 1. / F;

        
                ZeroPeriod[i] = Daysact (ValueDate, EndDate) / 365.;
                Zero6mDate[i] = EndDate;

                
                if (Get_Zero(&ZeroRate,	                                /* Interpolate discount curve every 6m */
                             &ZeroPrice,
                             NbDZero,
                             DZero,
                             DZeroDate,
                             ValueDate,
                             EndDate) == FAILURE)
                {
                        nrerror ("Problem interpolating discount curve! (CoB_Spotrate)");
                        goto FREE_MEM_AND_RETURN;
        
                }  /* if */                                                
                
                DZero6m[i] = ZeroPrice;
        
        }  /* for i */        
        

        if (AoS == 'S')                                                         /* For semi-annual curve we have to interpolate the 18m swap rate */
        {
                for (i = 1, Annuity = 0.; i <= 2; i++)
                        Annuity += CouponPeriod[i] / pow (1. + Zero6m[i], ZeroPeriod[i]);
        
                ZeroPrice = 1. / pow (1. + Zero6m[2], ZeroPeriod[2]);
                        
                SwapCurve[3] = (100. * (1. - ZeroPrice) / Annuity + SwapCurve[4]) / 2.;	/* 18m rate is the linear interpolation of 1yr and 2yr rate */
                
        }  /* if */        

        
        for (i = Step + 2; i <= 2 * MAXMAT; i += Step)                          /* If F = 1 we bootstrap annual swaps, if F =2, semi-annual swaps */
        {
                for (j = Step, Annuity = 0.; j < i; j += Step)
                {
                        Annuity += (SwapCurve[i] * CouponPeriod[j] - Libor[j]) * DZero6m[j];
                        
                }  /* for j */
                                 
                Libor[i] = SwapCurve[i] * CouponPeriod[i] + Annuity / DZero6m[i];
                                 
                if (Libor[i] < ERROR)	                                        /* Problem in bootstrapping */
                {
                        sprintf (ErrorMsg, "Problem in bootstrapping swap #%d in Spotrate (.)", i);
                        nrerror (ErrorMsg);
                        goto FREE_MEM_AND_RETURN;

                }  /* if */

                ZeroPrice = pow (1. + Zero6m[i-Step], ZeroPeriod[i-Step]) * (1. + Libor[i] / 100.);
                
                Zero6m[i] = pow (ZeroPrice, 1. / ZeroPeriod[i]) -1.;
                
        }  /* for i */

                
        if (AoS == 'A')                                                         /* For annual curve we have to interpolate zero at 18m, 30m, ... */
        {
                for (i = 3; i < 2 * MAXMAT; i += 2)
                {
                        linterp (       ZeroPeriod[i],                          /* Interpolation of the zero maturing 6m before the current one */
                                        &(Zero6m[i]),
                                        ZeroPeriod[i-1],
                                        ZeroPeriod[i+1],
                                        Zero6m[i-1],
                                        Zero6m[i+1]);

                }  /* for i */
        }  /* if */


        for (i = 3; i <= 2 * MAXMAT; i++)                                       /* Copy the zero curve back in the original array */
        {                                                                       /* taking into account the NbMMPoint offset.      */
                Zero[i - 3 + NbMMPoint] = Zero6m[i];
                ZeroDate[i - 3 + NbMMPoint] = Zero6mDate[i];
                
        }  /* for i */


        status = SUCCESS;

        FREE_MEM_AND_RETURN:

        free_dvector (CouponPeriod, 0, 2*MAXMAT);
        free_dvector (ZeroPeriod,   0, 2*MAXMAT);
        free_lvector (Zero6mDate,   0, 2*MAXMAT);
        free_dvector (Libor,        0, 2*MAXMAT);
        free_dvector (DZero6m,      0, 2*MAXMAT);
        free_dvector (Zero6m,       0, 2*MAXMAT);
        
        return (status);

}  /* CoB_Spotrate */



/*****  Zero_6m  ************************************************************/
/*
*       Calculate zeros and their corresponding maturity. This routine could
*       be bypassed by simply inputing these zero rates, their number and
*       their maturities. This is base on linear interpolation of the swap
*       yield curve.
*/
int     Zero_6m (       T_CURVE         *t_curve,                               /* (O) Structure of zero curve data 	*/
                        TERM_DATA       *term_data)                             /* (I) Structure of term structure data */
{
        double
                *SwapCurve;                                                     /* Array of interpolated and benchmark swap yields */
        int        	                                                        /*	                                           */
                status = FAILURE;                                               /* Error status = FAILURE initially 	           */


        SwapCurve = dvector (0, 2*MAXMAT);                                      /* We have to interpolate swap rates up to 50 years */

        t_curve->ValueDate = term_data->ValueDate;
        t_curve->NbZero = 2 * MAXMAT - 2 + term_data->NbMMPoint;                /* The zeros are: money markets and zeros at 6m interval from 2yr onwards */

        if (Yldconv (   t_curve->Zero,                                          /* Conversion of the money market rates. */
                        t_curve->ZeroDate,
                        term_data->NbMMPoint,
                        term_data->MMPoint,
                        term_data->MMYield,
                        t_curve->ValueDate,
                        term_data->MMB) == FAILURE)
        {
                goto FREE_MEM_AND_RETURN;
        
        }  /* if */                                                

        if (Coupint (   SwapCurve,                                              /* Interpolation of Swap yields by 6mo increments */
                        term_data->NbSwapPoint,                                 
                        term_data->SwapPoint,
                        term_data->SwapYield) == FAILURE)                                                    
        {
                goto FREE_MEM_AND_RETURN;
        
        }  /* if */                                                

        if (Spotrate (  t_curve->Zero,                                          /* Bootstapping of swap curve to get zeros at 6mo intervals */
                        t_curve->ZeroDate,                                      
                        SwapCurve,                                              
                        term_data->NbMMPoint,
                        term_data->AoS,
                        t_curve->ValueDate,
                        term_data->YearBase) == FAILURE)
        {
                goto FREE_MEM_AND_RETURN;
        
        }  /* if */                                                


        status = SUCCESS;

        FREE_MEM_AND_RETURN:

        free_dvector (SwapCurve, 0, 2*MAXMAT);

        return (status);

}  /* Zero_6m */



/*****  CoB_Zero_6m  ********************************************************/
/*
*       Calculate zeros and their corresponding maturity. This routine could
*       be bypassed by simply inputing these zero rates, their number and
*       their maturities. This is base on linear interpolation of the swap
*       yield curve.
*/
int     CoB_Zero_6m (   T_CURVE         *t_curve,                               /* (O) Structure of zero curve data          */
                        T_CURVE         *t_curve_d,                             /* (I) Structure of discount zero curve data */
                        TERM_DATA       *term_data)                             /* (I) Structure of term structure data      */
{
        double
                *SwapCurve;                                                     /* Array of interpolated and benchmark swap yields */
        int        
                status = FAILURE;                                               /* Error status = FAILURE initially */


        SwapCurve = dvector (0, 2*MAXMAT);                                      /* We have to interpolate swap rates up to 50 years */

        t_curve->ValueDate = term_data->ValueDate;
        t_curve->NbZero = 2 * MAXMAT - 2 + term_data->NbMMPoint;                /* The zeros are: money markets and zeros at 6m interval from 2yr onwards */

        if (Yldconv (   t_curve->Zero,                                          /* Conversion of the money market rates. */
                        t_curve->ZeroDate,
                        term_data->NbMMPoint,
                        term_data->MMPoint,
                        term_data->MMYield,
                        t_curve->ValueDate,
                        term_data->MMB) == FAILURE)
        {
                goto FREE_MEM_AND_RETURN;
        
        }  /* if */                                                

        if (Coupint (   SwapCurve,                                              /* Interpolation of Swap yields by 6mo increments */
                        term_data->NbSwapPoint,                                 
                        term_data->SwapPoint,
                        term_data->SwapYield) == FAILURE)                                                    
        {
                goto FREE_MEM_AND_RETURN;
        
        }  /* if */                                                

        if (CoB_Spotrate (  	t_curve->Zero,                                  /* Bootstapping of swap curve with cost of basis */
                                t_curve->ZeroDate,
                                t_curve_d->NbZero,
                                t_curve_d->Zero,
                                t_curve_d->ZeroDate,
                                SwapCurve,        
                                term_data->NbMMPoint,
                                term_data->AoS,
                                t_curve->ValueDate,
                                term_data->YearBase) == FAILURE)
        {
                goto FREE_MEM_AND_RETURN;
        
        }  /* if */                                                


        status = SUCCESS;

        FREE_MEM_AND_RETURN:

        free_dvector (SwapCurve, 0, 2*MAXMAT);

        return (status);

}  /* CoB_Zero_6m */
