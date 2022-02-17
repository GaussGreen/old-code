/****************************************************************************/
/*      Calculation routines customised for a turbo swap                    */
/****************************************************************************/
/*      TURBO.C                                                             */
/****************************************************************************/


/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"


/*****  Hyb3_Turbo_Pmt  *********************************************************/
/**                                                                          
         Price of one turbo payment.                                         
                                                                             
         Assumes dimensions as follows:                                      
             FX rates :              3-Dim                                   
             Domestic index:         2-Dim                                   
             Foreign index times FX: 3-Dim                                   
             Zero to payment:        2-Dim                                   
                                                                           */
int   Hyb3_Turbo_Pmt(TSLICE      TurboPtr,         /**< (O) The turbo payment     */
                /* Domestic portion */       
                TSLICE      DCouponRatePtr,   /**< (I) Domestic coupon rate  */
                double      DPrincForCoupon,  /**< (I) Domestic principal    */
                double      DDCFrac,          /**< (I) Domestic day count fr */
                char        DSimpOrComp,      /**< (I) Pmt simple or comp    */
                /* Foreign portion */
                TSLICE      FXFCouponRatePtr, /**< (I) For'gn cp X FX rate   */
                double      FPrincForCoupon,  /**< (I) Foreign principal     */
                double      FDCFrac,          /**< (I) Foreign day count fr  */
                char        FSimpOrComp,      /**< (I) Simple or compounding */
                /* Common for all  */
                double      FloorRate,        /**< (I) Floor level as rate   */
                double      CapRate,          /**< (I) Cap level as rate     */
                double      PrincForCapFloor, /**< (I) Princ for cap amount  */
                double      DCFracCapFloor,   /**< (I) Dcf for cap amount    */
                                              
                TSLICE      DZeroToPmtPtr,    /**< (I) Dom zero to turbo pmt */
                int         ResetForAll,      /**< (I) Reset flag for turbo  */                
                char        ResetEqualsPayDate, /**< (I) ='Y' if coincide    */

                int         t,                /**< (I) Current time point    */
                int         T,                /**< (I) Total nb of timepoints*/
                int         DCurve,           /**< (I) Curve to discount on  */
                HYB3_DEV_DATA   *dev_data,         /**< (I) Internal data         */
                HYB3_TREE_DATA  *tree_data)        /**< (I) Tree data             */
{


    /* Locals for addressing the slices */
    double    *Turbo;
    double    *FRate;  
    double    *Zero;   
    double    *DRate;  

    double
            DPmt,  /* Domestic payment */
            FPmt,  /* Foreign payment  */
            TurboPmt,
            FloorAmt,
            CapAmt,
            RateAux;

    int
            Top1,      Bottom1,        /* Limits of the tree (1rst dimension) */
           *Top2,     *Bottom2,
          **Top3,    **Bottom3,

            offset,
            i, j, k,                   /* Node indices                        */
            status = FAILURE;          /* Error status = FAILURE initially    */


    /* If not resetting/calculating the turbo payment
       simply DEV the existing turbo payment slice and return */  
    if (!ResetForAll)
    {
        if (Hyb3_Dev(TurboPtr,         
                t,
                T,
                DCurve,
                DISC_3D_CUPS,
                dev_data,
                tree_data) != SUCCESS)
        {
            goto RETURN;
        }
        return SUCCESS;
    }


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];	                


    /* If reset required, calculate turbo payment slice */
    if (ResetForAll)
    {

        /* caps and floors are applied to the payment values rather than the payment
           rates as floor/Cap rates are simple, turbo rates may be simple or compounding */
        FloorAmt = FloorRate * PrincForCapFloor * DCFracCapFloor;
        CapAmt   = CapRate   * PrincForCapFloor * DCFracCapFloor;


        for (i = Bottom1; i <= Top1; i ++)
        {
            
            /* Prepare for accessing 2-D slices */
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DRate = (double *)DCouponRatePtr + offset;
            Zero  = (double *)DZeroToPmtPtr + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {

                /* Prepare for accessing 3-D slices */
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                FRate = (double *)FXFCouponRatePtr + offset;
                Turbo = (double *)TurboPtr + offset;

                /* Domestic */
                if (DSimpOrComp ==  'C')
                {
                    RateAux = pow(1.0 + DRate[j] , DDCFrac) - 1.0;
                }
                else
                {
                    RateAux =  DRate[j] * DDCFrac;
                }
                DPmt = RateAux * DPrincForCoupon;


                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {

                    /* Foreign */
                    if (FSimpOrComp ==  'C')
                    {
                        RateAux = pow(1.0 + FRate[k] , FDCFrac) - 1.0;
                    }
                    else
                    {
                        RateAux =  FRate[k] * FDCFrac;
                    }
                    FPmt = RateAux * FPrincForCoupon;

                    TurboPmt = MAXMIN(DPmt + FPmt, CapAmt, FloorAmt);

                    /* No need to discount if coupon reset and pay dates coincide */
                    if (ResetEqualsPayDate == 'Y')
                    {
                        Turbo[k] = TurboPmt;
                    }
                    else
                    {
                        Turbo[k] = TurboPmt * Zero[j];
                    }

                }  /* For k */

            } /* For j */

        } /* For i */
    
        
    }  /* if it is reset for all payments */

  

    status = SUCCESS;

  RETURN:

    return (status);

}  /* Hyb3_Turbo_Pmt */
 

      
/*****  Hyb3_FiTurbo_Pmt  *********************************************************/
/**                                                                            
         Computes PV one FiTurbo cash-flow defined as:                        
                COLLAR(FX)*COLLAR(FLegRate) + COLLAR(DLegRate)                
         
         where FLegRate and DLegRate can be either domestic (2D Slice)
         or foreign (1D Slice) rates. 
         FX is assumed to be a 3D slice. 
                                                                              
  Note: In comparison,the standard turbo Pmt is such that                    
          -FLegRate is foreign only.                                         
          -DLegRate is domestic only                                          
          -There is a single COLLAR applied to the whole cash-flow         */
   
                                                                         
int   Hyb3_FiTurbo_Pmt(  TSLICE      FiTPtr,      /**< (O) FiT(urbo) payment      */
                /* Domestic Leg portion */       
                TSLICE      DLegRatePtr,     /**< (I) Domestic Leg cpn rate */
                int         DLegRateDim,     /**< (I) Dom leg rate dim      */  
                double      DPrincForCoupon, /**< (I) Domestic principal    */
                double      DDCFrac,         /**< (I) Domestic day count fr */
                double      DCapRate,        /**< (I) Domestic Cap rate     */
                double      DFloorRate,      /**< (I) Domestic Floor rate   */  
                /* Foreign Leg portion */
                TSLICE      FLegRatePtr,     /**< (I) Foreign Leg cpn rate  */
                int         FLegRateDim,     /**< (I) Foriegn Leg rate dim  */
                double      FPrincForCoupon, /**< (I) Foreign principal     */
                double      FDCFrac,         /**< (I) Foreign day count fr  */
                double      FCapRate,        /**< (I) Foreign Cap rate      */
                double      FFloorRate,      /**< (I) Foreign Floor rate    */
                
                TSLICE      FXRatePtr,       /**< (I) FX rate               */
                double      FXCap,           /**< (I) FX Cap                */
                double      FXFloor,         /**< (I) FX Floor              */
                                              
                TSLICE      DZeroToPmtPtr,   /**< (I) Dom zero to FiT pmt   */                
                int         ArrearsForAll,   /**< (I) Arrears flag for turbo*/
                int         t,               /**< (I) Current time point    */
                HYB3_TREE_DATA  *tree_data)       /**< (I) Tree data             */
{
    /* Locals for addressing the slices */

    double    *ZeroL=NULL;   
    double    *DRateL=NULL,*FRateL=NULL,*FXRateL=NULL;  
    double    *FiTurboL=NULL;  

    double
            DomAux,
            ForAux,
            FXAux;
            

    int
            Top1,      Bottom1,        /* Limits of the tree (1rst dimension) */
           *Top2,     *Bottom2,
          **Top3,    **Bottom3,

            offset,
            i, j, k,                   /* Node indices                        */
            status = FAILURE;          /* Error status = FAILURE initially    */

    
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];	                

    if (DLegRateDim !=1 && DLegRateDim != 2)
    {
        DR_Error("Slice dimensions must be 1 or 2 (FiT_Pmt)\n");
        goto RETURN;
    }

    if (FLegRateDim !=1 && FLegRateDim != 2)
    {
        DR_Error("Slice dimensions must be 1 or 2 (FiT_Pmt)\n");
        goto RETURN;
    }
   
       
    if (DLegRateDim == 1 && FLegRateDim == 1)
    {
        
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        DRateL = (double*) DLegRatePtr + offset;
        FRateL = (double*) FLegRatePtr + offset;
        for (i = Bottom1; i <= Top1; i++)
        {
            
            DomAux = MAXMIN(DRateL[i],DCapRate,DFloorRate);
            DomAux *= DDCFrac * DPrincForCoupon;

            ForAux = MAXMIN(FRateL[i],FCapRate,FFloorRate);
            ForAux *= FDCFrac * FPrincForCoupon;

            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);                
            ZeroL  = (double*) DZeroToPmtPtr + offset;
            for (j = Bottom2[i] ; j <= Top2[i] ; j++)
            {
                offset   = Hyb3_Node_Offset(3,i,j,t,tree_data);
                FXRateL  = (double*)FXRatePtr + offset;
                FiTurboL = (double*)FiTPtr + offset;
                for(k = Bottom3[i][j] ; k <= Top3[i][j]; k++)
                {

                    FXAux       = MAXMIN(FXRateL[k], FXCap,FXFloor);
                    FiTurboL[k] = DomAux + ForAux * FXAux;
                    if (!ArrearsForAll)
                    {
                        FiTurboL[k] *= ZeroL[j];
                     }
                        
                }/* for k */
            }/* for j */
        }/* for i */

    }/* DLeg=1dim and FLeg= 1 dim */
        
    else if (DLegRateDim == 1 && FLegRateDim == 2)
    {

        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        DRateL = (double*) DLegRatePtr + offset;
        for (i = Bottom1 ; i<= Top1 ; i++)
        {
            DomAux = MAXMIN(DRateL[i],DCapRate, DFloorRate);
            DomAux *= DDCFrac * DPrincForCoupon;
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            FRateL = (double*)FLegRatePtr + offset;
            ZeroL  = (double*)DZeroToPmtPtr + offset;
            for (j = Bottom2[i]; j <= Top2[i] ; j++)
            {

                ForAux   = MAXMIN(FRateL[j],FCapRate,FFloorRate);
                ForAux   *= FDCFrac * FPrincForCoupon;
                offset   = Hyb3_Node_Offset(3,i,j,t,tree_data);
                FXRateL  = (double*)FXRatePtr + offset;
                FiTurboL = (double*)FiTPtr + offset;
                for (k = Bottom3[i][j] ; k <= Top3[i][j] ; k++)
                {
                    FXAux       = MAXMIN(FXRateL[k], FXCap,FXFloor);
                    FiTurboL[k] = DomAux + ForAux * FXAux;
                    if (!ArrearsForAll)
                    {
                        FiTurboL[k] *= ZeroL[j];
                    }                        
                }/* for k */
            }/* for j */
        }/* for i */
    } /* Dleg = 1dim, FLeg = 2dim */

    else if (DLegRateDim == 2 && FLegRateDim == 1)
    {
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        FRateL = (double*) FLegRatePtr + offset;
        for (i = Bottom1 ; i <= Top1 ; i++)
        {
            ForAux = MAXMIN(FRateL[i],FCapRate, FFloorRate);
            ForAux *= FDCFrac * FPrincForCoupon;
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DRateL = (double*)DLegRatePtr + offset;
            ZeroL  = (double*)DZeroToPmtPtr + offset;
            for (j = Bottom2[i]; j <= Top2[i] ; j++)
            {
                DomAux   = MAXMIN(DRateL[j],DCapRate,DFloorRate);
                DomAux   *= DDCFrac * DPrincForCoupon;
                offset   = Hyb3_Node_Offset(3,i,j,t,tree_data);
                FXRateL  = (double*)FXRatePtr + offset;
                FiTurboL = (double*)FiTPtr + offset;
                for (k = Bottom3[i][j] ; k<= Top3[i][j] ; k++)
                {
                    FXAux       = MAXMIN(FXRateL[k], FXCap,FXFloor);
                    FiTurboL[k] = DomAux + ForAux * FXAux;
                    if (!ArrearsForAll)
                    {
                       FiTurboL[k] *= ZeroL[j];
                    }                        
                }/* for k */
            }/* for j */
        }/* for i */
    }/* DLeg = 2dim, FLeg = 1dim */

    else if (DLegRateDim == 2 && FLegRateDim == 2)
    {
        for (i=Bottom1 ; i<=Top1 ; i++)
        {
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DRateL = (double*) DLegRatePtr + offset;
            FRateL = (double*) FLegRatePtr + offset;
            ZeroL  = (double*)DZeroToPmtPtr + offset;
            for (j = Bottom2[i]; j <= Top2[i] ; j++)
            {
                DomAux = MAXMIN(DRateL[j],DCapRate,DFloorRate);
                DomAux *= DDCFrac * DPrincForCoupon;
                ForAux = MAXMIN(FRateL[j],FCapRate,FFloorRate);
                ForAux *= FDCFrac * FPrincForCoupon;

                offset   = Hyb3_Node_Offset(3,i,j,t,tree_data);
                FXRateL  = (double*)FXRatePtr + offset;
                FiTurboL = (double*)FiTPtr + offset;
                for (k = Bottom3[i][j] ; k <= Top3[i][j] ; k++)
                {
                    FXAux       = MAXMIN(FXRateL[k], FXCap,FXFloor);
                    FiTurboL[k] = DomAux + ForAux * FXAux;
                    if (!ArrearsForAll)
                    {
                        FiTurboL[k] *= ZeroL[j];
                    }                        
                }/* for k */
            }/* for j */
        }/* for i */
    }/* Dleg = 2; FLeg=2dim */



    status = SUCCESS;

RETURN:
    return (status);
}/* Hyb3_FiTurbo_Pmt */

/*****  Hyb3_Turbo_Price  *********************************************************/
/**                                                                             
         Price of a turbo swap including different amortisation schedules      
         for the domestic and foreign payments.                                
                                                                               
         Assumes dimensions as follows:                                        
             FX rates :              3-Dim                                     
             Domestic index:         2-Dim                                     
             Foreign index times FX: 3-Dim                                     
             Zero to payment:        2-Dim                                     
                                                                               
                                                                             */
int   Hyb3_Turbo_Price(TSLICE      TurboPtr,         /**< (O) The turbo             */
                  /* Domestic payments */       
                  TSLICE      DCouponRatePtr,   /**< (I) Domestic coupon rate  */
                  double      DPrincForCoupon,  /**< (I) Domestic principal    */
                  double      DDCFrac,          /**< (I) Domestic day count fr */
                  char        DSimpOrComp,      /**< (I) Pmt simple or comp    */
                  /* Foreign payments */
                  TSLICE      FXFCouponRatePtr, /**< (I) For'gn cp X FX rate   */
                  double      FPrincForCoupon,  /**< (I) Foreign principal     */
                  double      FDCFrac,          /**< (I) Foreign day count fr  */
                  char        FSimpOrComp,      /**< (I) Simple or compounding */
                  /* Common for all  */
                  double      FloorRate,        /**< (I) Floor level as rate   */
                  double      CapRate,          /**< (I) Cap level as rate     */
                  double      PrincForCapFloor, /**< (I) Princ for cap amount  */
                  double      DCFracCapFloor,   /**< (I) Dcf for cap amount    */
                                                
                  TSLICE      DZeroToPmtPtr,    /**< (I) Dom zero to turbo pmt */
                  int         ResetForAll,      /**< (I) Reset flag for turbo  */
                  char        ArrearsForAll,    /**< (I) Arrears flag for turbo*/

                  int         DPrincipalFlag,   /**< (I) Flag for dom principal*/
                  double      DPrincipal,       /**< (I) Domestic princ pmt    */
                  int         FPrincipalFlag,   /**< (I) Flag for foreign princ*/
                  double      FPrincipal,       /**< (I) Foreign princ pmt     */
                  TSLICE      FxSpotPtr,        /**< (I) FX rates at time point*/
                  int         t,                /**< (I) Current time point    */
                  int         T,                /**< (I) Total nb of timepoints*/
                  int         DCurve,           /**< (I) Curve to discount on  */
                  HYB3_DEV_DATA   *dev_data,         /**< (I) Internal data         */
                  HYB3_TREE_DATA  *tree_data)        /**< (I) Tree data             */
{
        /* Locals for addressing the slices */
        double    *Turbo;
        double    *FRate;  
        double    *FXRate; 


        double    *Zero;   
        double    *DRate;  



        double
                DPmt,  /* Domestic payment */
                FPmt,  /* Foreign payment  */
                TurboPmt,
                FloorAmt,
                CapAmt,
                RateAux;

        int
                Top1,      Bottom1,        /* Limits of the tree (1rst dimension) */
               *Top2,     *Bottom2,
              **Top3,    **Bottom3,

                offset,
                i, j, k,                   /* Node indices                        */
                status = FAILURE;          /* Error status = FAILURE initially    */

        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];	                


    /* Discounted expected value function */        
    if (Hyb3_Dev(TurboPtr,         
            t,
            T,
            DCurve,
            DISC_3D_CUPS,
            dev_data,
            tree_data) != SUCCESS)
    {
        goto RETURN;
    }
            
  

        /* Manage principal payments for foreign and domestic */
        if (FPrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    Turbo = (double *)TurboPtr + offset;
                    FXRate = (double *)FxSpotPtr + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Turbo[k] += FPrincipal * FXRate[k];

                    }  /* for k */
                }
            }	
        }  /* if */

        if (DPrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    Turbo = (double *)TurboPtr + offset;
                    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Turbo[k] += DPrincipal;

                    }  /* for k */	
                }
            }
        } 




        /* If reset, add turbo payments */
        if (ResetForAll)
        {



            FloorAmt = FloorRate * PrincForCapFloor * DCFracCapFloor;
            CapAmt   = CapRate   * PrincForCapFloor * DCFracCapFloor;


            for (i = Bottom1; i <= Top1; i ++)
            {
                
                /* Prepare for accessing 2-D slices */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                DRate = (double *)DCouponRatePtr + offset;
                Zero  = (double *)DZeroToPmtPtr + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {

                    /* Prepare for accessing 3-D slices */
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    FRate = (double *)FXFCouponRatePtr + offset;
                    Turbo = (double *)TurboPtr + offset;

                    /* Domestic */
                    if (DSimpOrComp ==  'C')
                    {
                        RateAux = pow(1.0 + DRate[j] , DDCFrac) - 1.0;
                    }
                    else
                    {
                        RateAux =  DRate[j] * DDCFrac;
                    }
                    DPmt = RateAux * DPrincForCoupon;


                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {

                        /* Foreign */
                        if (FSimpOrComp ==  'C')
                        {
                            RateAux = pow(1.0 + FRate[k] , FDCFrac) - 1.0;
                        }
                        else
                        {
                            RateAux =  FRate[k] * FDCFrac;
                        }
                        FPmt = RateAux * FPrincForCoupon;

                        TurboPmt = MAXMIN(DPmt + FPmt, CapAmt, FloorAmt);


                        if (ArrearsForAll == 'Y')
                        {
                            Turbo[k] += TurboPmt;
                        }
                        else
                        {
                            Turbo[k] += TurboPmt * Zero[j];
                        }

                    }  /* For k */

                } /* For j */

            } /* For i */
        
            
        }  /* if it is reset for all payments */

  

        status = SUCCESS;

      RETURN:

        return (status);



}  /* Hyb3_Turbo_Price */
 

/*****  Hyb3_FwdTurbo_t  **********************************************************/
/**
 *      Calculates the  forward swap price on the lattice minus an offset
 *      level.  In other  words, this  function calculates  the intrinsic
 *      value of an option on the turbo fwd swap struck at the given off-
 *      set level.
 *
 *      Assumptions are  made concerning the  dimensions of the variables
 *      being priced (see Hyb3_Turbo_Price above).
 *
 *      In addition,  the  foreign coupon rate must  be pre-multiplied by
 *      the FX rate.   Both the foreign and the domestic coupon rates are
 *      assumed to be DEV'd outside the code in case the arrears turbo is
 *      being priced.
 *
 *
 */


 int   Hyb3_FwdTurbo_t(TSLICE     FwdTurboPtr,         /**< (O) The forward         */
                  TSLICE     UlTurboPtr,          /**< (I) Underlying, no stubs*/

                  TSLICE     DCouponRatePtr,      /**< (I) Domestic rate       */
                  double     DPrincForCouponPaid, /**< (I) Domestic principal  */
                  double     DDCFracPaid,         /**< (I) Day count fraction  */
                  char       DSimpOrComp,         /**< (I) Simp or compounding */

                  TSLICE     FXFCouponRatePtr,    /**< (I) For rate X FX rate  */
                  double     FPrincForCouponPaid, /**< (I) Foreign principal   */
                  double     FDCFracPaid,         /**< (I) Foreign day count   */
                  char       FSimpOrComp,         /**< (I) Simp or compounding */

                  double     FloorRate,           /**< (I) Floor level as rate */
                  double     CapRate,             /**< (I) Cap level as rate   */
                  double     PrincForCapFloor,    /**< (I) Princ for cap amount*/
                  double     DCFracCapFloor,      /**< (I) Dcc for cap amount  */

                  TSLICE     DZeroToPmtPtr,       /**< (I) Domestic zero to pmt*/
                  int        ResetFlag,           /**< (I) Turbo reset flag    */
                  long       AccStartPaid,        /**< (I) Acc st of curr cp   */
                  char       DPmtDayCount,        /**< (I) Dom day count conv  */
                  char       FPmtDayCount,        /**< (I) For day count conv  */
                  char       CapFloorDayCount,    /**< (I) Cap/fl daycount conv*/
                  
                  int        ExerFlag,            /**< (I) True if exer date   */
                  double     Strike,              /**< (I) Strike amount       */
                  char       OptionStubConv,      /**< (I) Bond, swap or none  */
                  char       Arrears,             /**< (I) Arrears flag        */

                  long       CurrentDate,         /**< (I) Date at timepoint   */
                  int        t,                   /**< (I) Curr timepoint      */
                  int        T,                   /**< (I) Total nb of t'points*/
                  int       DCurve,               /**< (I) Discount curve      */
                  HYB3_DEV_DATA  *dev_data,            /**< (I) Internal data       */
                  HYB3_TREE_DATA *tree_data)           /**< (I) Tree data           */
{

    /* Locals for slice addressing */
    double     *FwdTurbo;
    double     *UlTurbo;   
    double     *DRate;     
    double     *FRate;

    double     *Zero;      

    double  
            RateAux,

            TurboPmtAcc,
            FPmtAcc,       
            DPmtAcc,

            TurboPmtStub,
            FPmtStub,       
            DPmtStub,

            TurboPmtFull,
            FPmtFull,       
            DPmtFull,

            FloorAmtAcc,
            CapAmtAcc,

            FloorAmtStub,
            CapAmtStub,

            FloorAmt,
            CapAmt,


            FDayCountAcc,
            DDayCountAcc,
            CapFloorDayCountAcc;



    int    
        

              Top1,   Bottom1,   /* Tree limits (1rst dim)             */
             *Top2,  *Bottom2,   /* Tree limits (2nd dim)              */
            **Top3, **Bottom3,   /* Tree limits (3rd dim)              */
            i, j, k,             /* Node indices                       */
            offset,
            status = FAILURE;    /* Error status                       */



        
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];




    /* Must check if DMode is compatible with underlying tree */
    if (tree_data->TreeType != TTYPE_FX2IR)
    {
        DR_Error("Inconsistent tree type for turbo pricing! "
                "(Hyb3_FwdTurbo_t)\n");
        goto RETURN;
    }
  


    /* If this is not an exercise date discount and return */
    if (!ExerFlag)                                              
    {

        if (Hyb3_Dev (FwdTurboPtr,  /* Disc expd value of a swap starting on the  */
                 t,            /* last exercise date                         */
                 T,
                 DCurve,
                 DISC_3D_CUPS,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
            
        }  /* if */
    
        status = SUCCESS;
        goto RETURN;
        
    }  /* if */


   
      
    /* "Basic" intrinsic value */
    for (i = Bottom1; i <= Top1; i ++)
    {
        for (j = Bottom2[i]; j <= Top2[i]; j++)
        {

            offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
            FwdTurbo = (double *)FwdTurboPtr + offset;
            UlTurbo  = (double *)UlTurboPtr + offset;

            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
            {
                FwdTurbo[k] = UlTurbo[k] - Strike;
            }  
        }
    }
   


  

    /* Calculate accrual if within accrual period */
    if (AccStartPaid <= CurrentDate)
    {


        if (DrDayCountFraction(AccStartPaid,
                               CurrentDate,                                 
                               FPmtDayCount,
                               &FDayCountAcc) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for dom pmt.(Hyb3_FwdTurbo_t)");
            goto RETURN;
        } 

        if (DrDayCountFraction(AccStartPaid,
                               CurrentDate,                                 
                               DPmtDayCount,
                               &DDayCountAcc) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for foreign pmt.(Hyb3_FwdTurbo_t)");
            goto RETURN;
        }  
        
        if (DrDayCountFraction(AccStartPaid,
                               CurrentDate,                                 
                               CapFloorDayCount,
                               &CapFloorDayCountAcc) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for cap/floor.(Hyb3_FwdTurbo_t)");
            goto RETURN;
        }  
 

        FloorAmt    = FloorRate * PrincForCapFloor * DCFracCapFloor;
        CapAmt      = CapRate   * PrincForCapFloor * DCFracCapFloor;
        FloorAmtAcc = FloorRate * PrincForCapFloor * CapFloorDayCountAcc;
        CapAmtAcc   = CapRate   * PrincForCapFloor * CapFloorDayCountAcc;

        FloorAmtStub = FloorAmt - FloorAmtAcc;
        CapAmtStub   = CapAmt   - CapAmtAcc;

       
        
        if (Arrears == 'Y')
        {
            
            /* Unless this is not a reset and the convention is */
            /* 'N', the recently paid coupon must be eliminated */
            /*   Assumption:  FRate and DRate  have been        */
            /*                DEV'd from the reset.             */
            if ( ResetFlag  ||  OptionStubConv != 'N') 
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    /* Prepare to address the 2-D slices */
                    offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                    DRate = (double *)DCouponRatePtr + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        /* Domestic */
                        if (DSimpOrComp ==  'C')
                        {
                            RateAux = pow(1.0 + DRate[j],DDCFracPaid)-1.0;
                        }  
                        else
                        {
                            RateAux =  DRate[j] * DDCFracPaid;
                        }
                        DPmtFull = RateAux * DPrincForCouponPaid;

                        /* Prepare to address the 3-D slices */
                        offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                        FRate = (double *)FXFCouponRatePtr + offset;
                        FwdTurbo = (double *)FwdTurboPtr + offset;

                        for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {

                            /* Foreign */
                            if (FSimpOrComp ==  'C')
                            {
                                RateAux = pow(1.0 + FRate[k],FDCFracPaid)- 1.0;
                            }
                            else
                            {
                                RateAux =  FRate[k] * FDCFracPaid;
                            }
                            FPmtFull = RateAux * FPrincForCouponPaid;

                            TurboPmtFull = MAXMIN(DPmtFull + FPmtFull, CapAmt
                                              , FloorAmt);
                            FwdTurbo[k] -= TurboPmtFull;

                        } /* For k */
                    }  /* For j */
                }   /* For i */
            }

            

            /* Now add stub accrued in case conv is B or S */
            if( (!ResetFlag) &&
                ((OptionStubConv == 'B') || (OptionStubConv == 'S')) )
            {    
                
                if (FSimpOrComp == 'C' || DSimpOrComp == 'C')
                {
                    DR_Error(" Not acceptable ... \n");
                    goto RETURN;
                }

                for (i = Bottom1; i <= Top1; i ++)
                {
                    /* Prepare to address the 2-D slices */
                    offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                    DRate = (double *)DCouponRatePtr + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {

                        /* Domestic */
                        RateAux =  DRate[j] * (DDCFracPaid - DDayCountAcc);
                        DPmtStub = RateAux * DPrincForCouponPaid;

                        /* Prepare to address the 3-D slices */
                        offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                        FRate = (double *)FXFCouponRatePtr + offset;
                        FwdTurbo = (double *)FwdTurboPtr + offset;


                        for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            /* Foreign */
                            RateAux =  FRate[k] 
                                      * (FDCFracPaid-FDayCountAcc);
                            FPmtStub = RateAux * FPrincForCouponPaid;

                            TurboPmtStub = MAXMIN(DPmtStub + FPmtStub,
                                                   CapAmtStub , FloorAmtStub);
                            FwdTurbo[k] += TurboPmtStub;
                            
                        }  /* For k  */
                    }   /* For j */
                }   /* For i */

            }  /* if then else */      
                        
            
        }
        else /*  Reset in advance */
        {

            if (ResetFlag) /* Nothing to do, next cp just made*/
            {
                    ;
            }    
            else if(FSimpOrComp == 'C' || DSimpOrComp == 'C')
            {
                DR_Error("Accrual calculations are not supported \n"
                        "when pmt is on a compounding rate."
                        "(Hyb3_FwdTurbo_t)\n");
                        goto RETURN;
            }
            else   
            {
                
                for (i = Bottom1; i <= Top1; i ++)
                {

                    /* Prepare to address the 2-D slices */
                    offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                    DRate = (double *)DCouponRatePtr + offset;
                    Zero = (double *)DZeroToPmtPtr + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        /* Domestic */
                        DPmtAcc  = DRate[j] * DDayCountAcc 
                                 * DPrincForCouponPaid;
                        DPmtFull = DRate[j] * DDCFracPaid 
                                 * DPrincForCouponPaid;
                        DPmtStub = DPmtFull - DPmtAcc;


                        /* Prepare to address the 3-D slices */
                        offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                        FRate = (double *)FXFCouponRatePtr + offset;
                        FwdTurbo = (double *)FwdTurboPtr + offset;


                        for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            /* Foreign */
                            FPmtAcc = FRate[k] * FDayCountAcc 
                                 * FPrincForCouponPaid;
                            FPmtFull = FRate[k] * FDCFracPaid 
                                 * FPrincForCouponPaid;
                            FPmtStub = FPmtFull - FPmtAcc;


                            TurboPmtAcc = MAXMIN(DPmtAcc + FPmtAcc, CapAmtAcc
                                             , FloorAmtAcc);
                            TurboPmtFull = MAXMIN(DPmtFull + FPmtFull, CapAmt
                                               , FloorAmt);
                            TurboPmtStub = MAXMIN(DPmtStub + FPmtStub, 
                                                   CapAmtStub , FloorAmtStub);


                            switch(OptionStubConv)
                            {
                                case 'B':
                                    FwdTurbo[k] += TurboPmtFull * Zero[j]
                                                         - TurboPmtAcc;
                                    break;

                                case 'S':
                                    FwdTurbo[k] += TurboPmtStub * Zero[j];
                                    break;

                                case 'N':
                                    FwdTurbo[k] += TurboPmtFull * Zero[j];
                                    break;
                            }

                        }  /* For k */
                    }  /* For j */
                }   /* For i */
            }   /* If not a reset date */    
               

        } /* Arrears or advance */
    }




    status = SUCCESS;

  RETURN:

    return (status);



}  /* Hyb3_FwdTurbo_t */




/*****  Hyb3_BTurbo_Price  ********************************************************/
/**                                                                             
         Price of a binary turbo swap including different amortisation         
         schedules for the domestic and foreign payments.                      
                                                                               
         Assumes dimensions as follows:                                        
             FX rates :              3-Dim                                     
             Domestic index:         2-Dim                                     
             Foreign index times FX: 3-Dim                                     
             Zero to payment:        2-Dim                                     
                                                                               
                                                                             */
int  Hyb3_BTurbo_Price(TSLICE      TurboPtr,         /**< (O) The binary turbo      */
                  /* Domestic part of formula */       
                  TSLICE      DCouponRatePtr,   /**< (I) Domestic coupon rate  */
                  double      DPrincForCoupon,  /**< (I) Domestic principal    */
                  double      DDCFrac,          /**< (I) Domestic day count fr */
                  char        DSimpOrComp,      /**< (I) Pmt simple or comp    */
                  /* Foreign part of formula  */
                  TSLICE      FXFCouponRatePtr, /**< (I) For'gn cp X FX rate   */
                  double      FPrincForCoupon,  /**< (I) Foreign principal     */
                  double      FDCFrac,          /**< (I) Foreign day count fr  */
                  char        FSimpOrComp,      /**< (I) Simple or compounding */
                  /* Common for all  */
                  double      FloorRate,        /**< (I) Floor level as rate   */
                  double      CapRate,          /**< (I) Cap level as rate     */
                  double      PrincForCapFloor, /**< (I) Princ for cap amount  */
                  double      DCFracCapFloor,   /**< (I) Dcf for cap amount    */

                  TSLICE      DZeroToPmtPtr,    /**< (I) Dom zero to turbo pmt */
                  int         ResetForAll,      /**< (I) Reset flag for turbo  */
                  char        ArrearsForAll,    /**< (I) Arrears flag for turbo*/

                  int         DPrincipalFlag,   /**< (I) Flag for dom principal*/
                  double      DPrincipal,       /**< (I) Domestic princ pmt    */
                  int         FPrincipalFlag,   /**< (I) Flag for foreign princ*/
                  double      FPrincipal,       /**< (I) Foreign princ pmt     */
                  TSLICE      FxSpotPtr,        /**< (I) FX rates at time point*/
                  /* Binary */
                  int         NbBarriers,       /**< (I) Nb FX barrier levels  */
                  double     *BarrierLevel,     /**< (I) Marking FX ranges     */
                  double     *BarrierWeight,    /**< (I) Weights for ranges    */
                  char        Smoothing,        /**< (I) TRUE if smoothing on  */

                  int         t,                /**< (I) Current time point    */
                  int         T,                /**< (I) Total nb of timepoints*/
                  int         DCurve,           /**< (I) Curve to discount on  */
                  HYB3_DEV_DATA   *dev_data,         /**< (I) Internal data         */
                  HYB3_TREE_DATA  *tree_data)        /**< (I) Tree data             */
{

        /* Locals for addressing the slices */
        double    *Turbo;
        double    *FRate;  
        double    *FXRate; 

        double    *Zero;   
        double    *DRate;  


        double     FloorAmt;
        double     CapAmt;

        double
                   DPmt,  /* Domestic payment */
                   FPmt,  /* Foreign payment  */
                   TurboPmt,
                   RateAux,
                   SmoothAux,
                   IndexStep;

        double     Barrier1 = 0.0;       /* Barrier levels applicable to    */
        double     Barrier2 = 0.0;       /* current FX level                */
        double     UpValue1, DownValue1; /* P/off above and below barrier 1 */
        double     UpValue2, DownValue2; /* P/off above and below barrier 2 */

        int        BoundLeftRight;       /* True when FX falls between barr */

        int
                Top1,      Bottom1,      /* Limits of the tree (1rst dim)   */
               *Top2,     *Bottom2,
              **Top3,    **Bottom3,
                                         /* Index used in search of barrier */
                BarrIdx,
                offset,
                i, j, k, l,             /* Node indices                     */
                status = FAILURE;       /* Error status = FAILURE initially */

        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];	                


    /* Discounted expected value function */        
    if (Hyb3_Dev(TurboPtr,         
            t,
            T,
            DCurve,
            DISC_3D_CUPS,
            dev_data,
            tree_data) != SUCCESS)
    {
        goto RETURN;
    }
            
  

        /* Manage principal payments for foreign and domestic */
        if (FPrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    Turbo = (double *)TurboPtr + offset;
                    FXRate = (double *)FxSpotPtr + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Turbo[k] += FPrincipal * FXRate[k];

                    }  /* for k */
                }
            }	
        }  /* if */

        if (DPrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    Turbo = (double *)TurboPtr + offset;
                    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Turbo[k] += DPrincipal;

                    }  /* for k */	
                }
            }
        } 




        /* If reset, add turbo payments */
        if (ResetForAll)
        {


            FloorAmt = FloorRate * PrincForCapFloor * DCFracCapFloor;
            CapAmt   = CapRate   * PrincForCapFloor * DCFracCapFloor;


            for (i = Bottom1; i <= Top1; i ++)
            {
                
                /* Prepare for accessing 2-D slices */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                DRate = (double *)DCouponRatePtr + offset;
                Zero  = (double *)DZeroToPmtPtr + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {

                    /* Prepare for accessing 3-D slices */
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    FXRate = (double *)FxSpotPtr + offset;
                    FRate  = (double *)FXFCouponRatePtr + offset;
                    Turbo  = (double *)TurboPtr + offset;

                    /* Domestic */
                    if (DSimpOrComp ==  'C')
                    {
                        RateAux = pow(1.0 + DRate[j] , DDCFrac) - 1.0;
                    }
                    else
                    {
                        RateAux =  DRate[j] * DDCFrac;
                    }
                    DPmt = RateAux * DPrincForCoupon;


                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {

                        /* Foreign */
                        if (FSimpOrComp ==  'C')
                        {
                            RateAux = pow(1.0 + FRate[k] , FDCFrac) - 1.0;
                        }
                        else
                        {
                            RateAux =  FRate[k] * FDCFrac;
                        }
                        FPmt = RateAux * FPrincForCoupon;

                        /* In principle, the pmt is dom + for */
                        TurboPmt = DPmt + FPmt;

                        /* Now deal with the binary clause */
                        /* First bracket the FX with the barrier levels */
                        BarrIdx = 0; 
                        BoundLeftRight = TRUE;
                        UpValue2 = 0.0; /* To avoid warnings */
                        for(l=0; l<NbBarriers; l++)
                        {
                            if (FXRate[k] < BarrierLevel[l])
                            {
                                break;
                            }
                            BarrIdx++;
                        }

                      
                        /* Now calculate pay off for smoothed and */
                        /* non-smoothed cases                     */
                        if (Smoothing == 'N')
                        {
                            /* The weights are multiplicative */
                            TurboPmt *= BarrierWeight[BarrIdx];           

                            /* Finally apply cap and floor */
                            TurboPmt = MAXMIN(TurboPmt,CapAmt,FloorAmt);

                        }
                        else
                        {
                            
                            /* Establish UpValues and DownValues  */
                            /* for smoothing of the payoff        */
                            if (BarrIdx == 0)                             
                            {                                             
                                BoundLeftRight = FALSE;
                                Barrier1 = BarrierLevel[0];                   
                                DownValue1 = BarrierWeight[0];            
                                UpValue1 = BarrierWeight[1];              
                            }                                             
                            else if (BarrIdx == NbBarriers)                 
                            {                                             
                                BoundLeftRight = FALSE;
                                Barrier1 = BarrierLevel[NbBarriers - 1];                   
                                DownValue1 = BarrierWeight[NbBarriers - 1]; 
                                UpValue1   = BarrierWeight[NbBarriers];     
                            }                                             
                            else                                          
                            {                                             
                                Barrier1 = BarrierLevel[BarrIdx - 1]; 
                                DownValue1 = BarrierWeight[BarrIdx -1];   
                                UpValue1   = BarrierWeight[BarrIdx];
                                Barrier2 = BarrierLevel[BarrIdx];      
                                DownValue2 = BarrierWeight[BarrIdx];      
                                UpValue2   = BarrierWeight[BarrIdx + 1];  
                            }                                             
                                                                         

                            IndexStep = Hyb3_GetIndexStep(FxSpotPtr,
                                                     3,
                                                     i, j, k,
                                                     t,
                                                     tree_data);

                            /* First barrier */
                            if (Smooth_Step(&(SmoothAux),
                                    MAXMIN(TurboPmt*UpValue1,  CapAmt,FloorAmt),
                                    MAXMIN(TurboPmt*DownValue1,CapAmt,FloorAmt),
                                    FXRate[k],
                                    Barrier1,
                                    IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }

                            
                            /* Second barrier, if required */
                            if (BoundLeftRight)
                            {
                                if (Smooth_Step(&(TurboPmt),
                                         MAXMIN(TurboPmt*UpValue2,CapAmt,FloorAmt),
                                         SmoothAux,
                                         FXRate[k],
                                         Barrier2,
                                         IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }
                            else
                            {
                                TurboPmt = SmoothAux;
                            }
                        }
                        
                        /* Add payment to total turbo */
                        if (ArrearsForAll == 'Y')
                        {
                            Turbo[k] += TurboPmt;
                        }
                        else
                        {
                            Turbo[k] += TurboPmt * Zero[j];
                        }


                    }  /* For k */

                } /* For j */

            } /* For i */
        
            
        }  /* if it is reset for all payments */

  

        status = SUCCESS;

      RETURN:

        return (status);



}  /* Hyb3_BTurbo_Price */
                  
 

/*****  Hyb3_FwdBTurbo_t  *********************************************************/
/**
        Calculates the  forward swap price on the lattice minus an offset
        level.  In other  words, this  function calculates  the intrinsic
        value of an option on the turbo fwd swap struck at the given off-
        set level.
  
        IMPORTANT ASSUMPTION:  This  forward calculator does  not support
        stubs. An  exercise must be an accrual start or else fall earlier
        than the absolute first accrual start of the turbo. 
  
        Assumptions are  made concerning the  dimensions of the variables
        being priced (see Hyb3_Turbo_Price above).
  
        In addition,  the  foreign coupon rate must  be pre-multiplied by
        the FX rate.   Both the foreign and the domestic coupon rates are
        assumed to be DEV'd outside the code in case the arrears turbo is
        being priced.
  
  
 */


 int  Hyb3_FwdBTurbo_t(TSLICE     FwdTurboPtr,    /**< (O) The forward         */
                  TSLICE     UlTurboPtr,          /**< (I) Underlying, no stubs*/

                  TSLICE     DCouponRatePtr,      /**< (I) Domestic rate       */
                  double     DPrincForCouponPaid, /**< (I) Domestic principal  */
                  double     DDCFracPaid,         /**< (I) Day count fraction  */
                  char       DSimpOrComp,         /**< (I) Simp or compounding */

                  TSLICE     FXFCouponRatePtr,    /**< (I) For rate X FX rate  */
                  double     FPrincForCouponPaid, /**< (I) Foreign principal   */
                  double     FDCFracPaid,         /**< (I) Foreign day count   */
                  char       FSimpOrComp,         /**< (I) Simp or compounding */
                  int        ResetFlag,           /**< (I) Turbo reset flag    */

                  double     FloorRate,           /**< (I) Floor level as rate */
                  double     CapRate,             /**< (I) Cap level as rate   */
                  double     PrincForCapFloor,    /**< (I) Princ for cap amount*/
                  double     DCFracCapFloor,      /**< (I) Dcc for cap amount  */
                
                  TSLICE     FxSpotPtr,           /**< (I) Current FX rates    */
                  int         NbBarriers,         /**< (I) Nb FX barr levels   */
                  double     *BarrierLevel,       /**< (I) Marking FX ranges   */
                  double     *BarrierWeight,      /**< (I) Weights for ranges  */
                  char       Smoothing,           /**< (I) Smothing on/off     */
                  
                  int        ExerFlag,            /**< (I) True if exercise    */
                  double     Strike,              /**< (I) Strike amount       */
                  char       Arrears,             /**< (I) Arrears flag        */

                  int        t,                   /**< (I) Curr timepoint      */
                  int        T,                   /**< (I) Total nb of t'points*/
                  int        DCurve,              /**< (I) Discount curve      */
                  HYB3_DEV_DATA   *dev_data,      /**< (I) Internal data       */
                  HYB3_TREE_DATA  *tree_data)     /**< (I) Tree data           */


{


    /* Locals for slice addressing */
    double     *FwdTurbo;
    double     *UlTurbo;   
    double     *DRate;     
    double     *FRate;
    double     *FXRate;


    double      RateAux;       /* Auxiliary for when rate is compounding */
    double      SmoothAux;     /* Auxiliary for the smoothing algorithm  */
    double      IndexStep;     /* Used in the smoothing algorithm        */

    double      TurboPmtFull;  /* Full turbo pmt just added              */
    double      FPmtFull;      /* Foreign portion of hte above           */ 
    double      DPmtFull;      /* Domestic portion of the above          */
    double      FloorAmt;
    double      CapAmt;


    int         BarrIdx;              /* Index for XF barrier look-up    */
    double      Barrier1 = 0.0;       /* Barrier levels applicable to    */
    double      Barrier2 = 0.0;       /* current FX level                */
    double      UpValue1, DownValue1; /* P/off above and below barrier 1 */
    double      UpValue2, DownValue2; /* P/off above and below barrier 2 */

    int         BoundLeftRight;       /* TRUE if FX level in between two */
                                      /* specified barrier values        */


    int         Top1,   Bottom1;   /* Tree limits (1rst dim)             */
    int        *Top2,  *Bottom2;   /* Tree limits (2nd dim)              */
    int       **Top3, **Bottom3;   /* Tree limits (3rd dim)              */
    int         i, j, k, l;        /* Node indices                       */
    int         offset;
    int         status = FAILURE;  /* Error status                       */



        
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];




    /* Must check if DMode is compatible with underlying tree */
    if (tree_data->TreeType != TTYPE_FX2IR)
    {
        DR_Error("Inconsistent tree type for turbo pricing! "
                "(Hyb3_FwdTurbo_t)\n");
        goto RETURN;
    }
  


    /* If this is not an exercise date discount and return */
    if (!ExerFlag)                                              
    {

        if (Hyb3_Dev (FwdTurboPtr,  /* Disc expd value of a swap starting on the  */
                 t,            /* last exercise date                         */
                 T,
                 DCurve,
                 DISC_3D_CUPS,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
            
        }  /* if */
    
        status = SUCCESS;
        goto RETURN;
        
    }  /* if */


   
      
    /* "Basic" intrinsic value */
    for (i = Bottom1; i <= Top1; i ++)
    {
        for (j = Bottom2[i]; j <= Top2[i]; j++)
        {

            offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
            FwdTurbo = (double *)FwdTurboPtr + offset;
            UlTurbo  = (double *)UlTurboPtr + offset;

            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
            {
                FwdTurbo[k] = UlTurbo[k] - Strike;
            }  
        }
    }
   

    /* IMPORTANT: Now we must trust the caller to have checked that */
    /*            exercises DO NOT result in stubs!                 */

  

    /* If reset  is in arrears and  has just occurred */
    /* then the recently added coupon must be removed */
    if ((Arrears == 'Y') && (ResetFlag))
    {

        FloorAmt    = FloorRate * PrincForCapFloor * DCFracCapFloor;
        CapAmt      = CapRate   * PrincForCapFloor * DCFracCapFloor;
  
        for (i = Bottom1; i <= Top1; i ++)
        {
            /* Prepare to address the 2-D slices */
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DRate = (double *)DCouponRatePtr + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                /* Domestic */
                if (DSimpOrComp ==  'C')
                {
                    RateAux = pow(1.0 + DRate[j],DDCFracPaid)-1.0;
                }  
                else
                {
                    RateAux =  DRate[j] * DDCFracPaid;
                }
                DPmtFull = RateAux * DPrincForCouponPaid;


                /* Prepare to address the 3-D slices */
                offset   = Hyb3_Node_Offset(3,i,j,t,tree_data);
                FXRate   = (double *)FxSpotPtr + offset;
                FRate    = (double *)FXFCouponRatePtr + offset;
                FwdTurbo = (double *)FwdTurboPtr + offset;


                for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                {

                    /* Foreign */
                    if (FSimpOrComp ==  'C')
                    {
                        RateAux = pow(1.0 + FRate[k],FDCFracPaid)- 1.0;
                    }
                    else
                    {
                        RateAux =  FRate[k] * FDCFracPaid;
                    }
                    FPmtFull = RateAux * FPrincForCouponPaid;


                    /* In principle, pmt is just the domestic plus the foreign */
                    TurboPmtFull = DPmtFull + FPmtFull;

                    /* Now deal with the binary clause */
                    /* First bracket the FX with the barrier levels */
                    BarrIdx = 0; 
                    BoundLeftRight = TRUE;
                    UpValue2 = 0.0; /* To avoid warnings */
                    for(l=0; l<NbBarriers; l++)
                    {
                        if (FXRate[k] < BarrierLevel[l])
                        {
                            break;
                        }
                        BarrIdx++;
                    }



                    /* Now calculate pay off for smoothed and */
                    /* non-smoothed cases                     */
                    if (Smoothing == 'N')
                    {

                        TurboPmtFull *= BarrierWeight[BarrIdx];

                        /* Finally apply cap and floor */
                        TurboPmtFull = MAXMIN(TurboPmtFull, CapAmt
                                          , FloorAmt);
                    }
                    else
                    {

                        /* Establish UpValues and DownValues  */
                        /* for smoothing of the payoff        */
                        if (BarrIdx == 0)                             
                        {                                             
                            BoundLeftRight = FALSE;
                            Barrier1 = BarrierLevel[0];                   
                            DownValue1 = BarrierWeight[0];            
                            UpValue1 = BarrierWeight[1];              
                        }                                             
                        else if (BarrIdx == NbBarriers)                 
                        {                                             
                            BoundLeftRight = FALSE;
                            Barrier1 = BarrierLevel[NbBarriers - 1];                   
                            DownValue1 = BarrierWeight[NbBarriers - 1]; 
                            UpValue1   = BarrierWeight[NbBarriers];     
                        }                                             
                        else                                          
                        {                                             
                            Barrier1 = BarrierLevel[BarrIdx - 1]; 
                            DownValue1 = BarrierWeight[BarrIdx -1];   
                            UpValue1   = BarrierWeight[BarrIdx];
                            Barrier2 = BarrierLevel[BarrIdx];      
                            DownValue2 = BarrierWeight[BarrIdx];      
                            UpValue2   = BarrierWeight[BarrIdx + 1];  
                        }                                             
                                     


                        IndexStep = Hyb3_GetIndexStep(FxSpotPtr,
                                                 3,
                                                 i, j, k,
                                                 t,
                                                 tree_data);

                        
                        if (Smooth_Step(&(SmoothAux),
                                MAXMIN(TurboPmtFull*UpValue1,  CapAmt,FloorAmt),
                                MAXMIN(TurboPmtFull*DownValue1,CapAmt,FloorAmt),
                                FXRate[k],
                                Barrier1,
                                IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        
                        /* Consider second barrier, if applicable */
                        if (BoundLeftRight)
                        {
                            if (Smooth_Step(&(TurboPmtFull),
                                    MAXMIN(TurboPmtFull*UpValue2,CapAmt,FloorAmt),
                                    SmoothAux,
                                    FXRate[k],
                                    Barrier2,
                                    IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                        else
                        {
                            TurboPmtFull = SmoothAux;
                        }
                    }
                    
                    FwdTurbo[k] -= TurboPmtFull;

                } /* For k */

            }  /* For j */

        }   /* For i */            
        
    }




    status = SUCCESS;

  RETURN:

    return (status);



}  /* Hyb3_FwdBTurbo_t */








/*****  Hyb3_Redeemer_Price  *****************************************************/
/**
 *      Calculates the price of a redeemer call spread   
 *
 *      Assumes dimensions as follows:                                      
 *           Redeemer:               3-Dim
 *           FX rates:               3-Dim                                   
 *           Dom zero to rdmr pmt:   2-Dim 
 *
 */
int   Hyb3_Redeemer_Price(
		     TSLICE     RedeemerPtr,      /**< (O) Underlying rdmr   */
                     TSLICE     FwdFXRdmrPtr,     /**< (I) Forward FX        */
                     TSLICE     DZeroToRdmrPmtPtr,/**< (I) Dom ZC to rdmr pmt*/

                     double     Notional,         /**< (I) Notional          */
                     double     LowStrike,        /**< (I) Rdmr low strike   */
                     double     HighStrike,       /**< (I) Rdmr high strike  */

                     int        t,                /**< (I) Current time prd  */
                     HYB3_TREE_DATA  *tree_data)  /**< (I) Tree data         */
{

    double
          *RedeemerPtrL,
          *FwdFXRdmrPtrL,
          *DZeroToRdmrPtrL;

    int
          Top1,     Bottom1,           /* Limits of the tree (1st dim) */
         *Top2,    *Bottom2,           /* Limits of the tree (2nd dim) */
        **Top3,   **Bottom3,           /* Limits of the tree (3rd dim) */

        offset,
        i, j, k;                       /* Node indices                 */
       

        
        Top1    = tree_data->Top1[t];	                
        Bottom1 = tree_data->Bottom1[t];
        Top2    = tree_data->Top2[t];	                
        Bottom2 = tree_data->Bottom2[t];
        Top3    = tree_data->Top3[t];	                
        Bottom3 = tree_data->Bottom3[t];	                

        
        for (i=Bottom1; i<=Top1; i++)
        {
            DZeroToRdmrPtrL = (double *)DZeroToRdmrPmtPtr +
                Hyb3_Node_Offset(2,i,0,t,tree_data);

            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {                
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                RedeemerPtrL  = (double *)RedeemerPtr + offset;
                FwdFXRdmrPtrL = (double *)FwdFXRdmrPtr + offset;
                
                for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                {
                    RedeemerPtrL[k] = Notional * (
                    MAX(FwdFXRdmrPtrL[k] - DZeroToRdmrPtrL[j]*HighStrike,0)
				   -MAX(DZeroToRdmrPtrL[j]*LowStrike  - FwdFXRdmrPtrL[k],0));                                                
                }
            }
        }
 
 
        return (SUCCESS);


}  /*  Hyb3_Redeemer_Price  */


/*****  Hyb3_BTurbo_Flows  ********************************************************/
/**                                                                            
         Price of a binary turbo swap including different amortisation         
         schedules for the domestic and foreign payments.                      
                                                                               
         Similar to Hyb3_BTurbo_Price, but here user needs to input 2 extra         
         Zeros to Pmt to discount dom/for principal payments                   
                                                                             */
int  Hyb3_BTurbo_Flows(
		  TSLICE      TurboPtr,         /**< (O) The binary turbo      */
                  /* Domestic part of formula */       
                  TSLICE      DCouponRatePtr,   /**< (I) Domestic coupon rate  */
                  double      DPrincForCoupon,  /**< (I) Domestic principal    */
                  double      DDCFrac,          /**< (I) Domestic day count fr */
                  /* Foreign part of formula  */
                  TSLICE      FXFCouponRatePtr, /**< (I) For'gn cp X FX rate   */
                  double      FPrincForCoupon,  /**< (I) Foreign principal     */
                  double      FDCFrac,          /**< (I) Foreign day count fr  */
                  /* Common for all  */
                  double      FloorRate,        /**< (I) Floor level as rate   */
                  double      CapRate,          /**< (I) Cap level as rate     */
                  double      PrincForCapFloor, /**< (I) Princ for cap amount  */
                  double      DCFracCapFloor,   /**< (I) Dcf for cap amount    */

                  TSLICE      DZeroToPmtPtr,    /**< (I) Dom zero to turbo pmt */
                  int         ResetForAll,      /**< (I) Reset flag for turbo  */
                  char        ArrearsForAll,    /**< (I) Arrears flag for turbo*/

                  int         DPrincipalFlag,   /**< (I) Flag for dom principal*/
                  double      DPrincipal,       /**< (I) Domestic princ pmt    */
                  TSLICE      DZeroToPrinPmtPtr, /**< (I) For zero to prin pmt */
                  int         FPrincipalFlag,   /**< (I) Flag for foreign princ*/
                  double      FPrincipal,       /**< (I) Foreign princ pmt     */
                  TSLICE      FZeroToPrinPmtPtr, /**< (I) Dom zero to prin pmt */
                  TSLICE      PrinFxPtr,        /**< (I) FX rates for foreign prin */
                  /* Binary */
                  TSLICE      BarrFxPtr,        /**< (I) FX rates for fx barrier */
                  int         NbBarriers,       /**< (I) Nb FX barrier levels  */
                  double     *BarrierLevel,     /**< (I) Marking FX ranges     */
                  double     *BarrierWeight,    /**< (I) Weights for ranges    */
                  char        Smoothing,        /**< (I) TRUE if smoothing on  */

                  int         t,                /**< (I) Current time point    */
                  int         T,                /**< (I) Total nb of timepoints*/
                  int         DCurve,           /**< (I) Curve to discount on  */
                  HYB3_DEV_DATA   *dev_data,         /**< (I) Internal data         */
                  HYB3_TREE_DATA  *tree_data)        /**< (I) Tree data             */
                  
                  
{


        /* Locals for addressing the slices */
        double    *Turbo;
        double    *FRate;  
        double    *BarrFXRate; 
        double    *PrinFXRate; 

        double    *Zero;   
        double    *DZero;   
        double    *FZero;   
        double    *DRate;  


        double     FloorAmt;
        double     CapAmt;

        double
                   DPmt,  /* Domestic payment */
                   FPmt,  /* Foreign payment  */
                   TurboPmt,
                   SmoothAux,
                   IndexStep;

        double     Barrier1 = 0.0;       /* Barrier levels applicable to    */
        double     Barrier2 = 0.0;       /* current FX level                */
        double     UpValue1, DownValue1; /* P/off above and below barrier 1 */
        double     UpValue2, DownValue2; /* P/off above and below barrier 2 */

        int        BoundLeftRight;       /* True when FX falls between barr */

        int
                Top1,      Bottom1,      /* Limits of the tree (1rst dim)   */
               *Top2,     *Bottom2,
              **Top3,    **Bottom3,
                                         /* Index used in search of barrier */
                BarrIdx,
                offset,
                i, j, k, l,             /* Node indices                     */
                status = FAILURE;       /* Error status = FAILURE initially */

        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];	                


    /* Discounted expected value function */        
    if (Hyb3_Dev(TurboPtr,         
            t,
            T,
            DCurve,
            DISC_3D_CUPS,
            dev_data,
            tree_data) != SUCCESS)
    {
        goto RETURN;
    }
            
  

        /* Manage principal payments for foreign and domestic */
        if (FPrincipalFlag)
        {
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            FZero = (double *)FZeroToPrinPmtPtr + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    Turbo = (double *)TurboPtr + offset;
                    PrinFXRate = (double *)PrinFxPtr + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Turbo[k] += FPrincipal * PrinFXRate[k] * FZero[i];

                    }  /* for k */
                }
            }	
        }  /* if */

        if (DPrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                /* Prepare for accessing 2-D slices */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                DZero  = (double *)DZeroToPrinPmtPtr + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    Turbo = (double *)TurboPtr + offset;
                    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Turbo[k] += DPrincipal * DZero[j];

                    }  /* for k */	
                }
            }
        } 




        /* If reset, add turbo payments */
        if (ResetForAll)
        {


            FloorAmt = FloorRate * PrincForCapFloor * DCFracCapFloor;
            CapAmt   = CapRate   * PrincForCapFloor * DCFracCapFloor;


            for (i = Bottom1; i <= Top1; i ++)
            {
                
                /* Prepare for accessing 2-D slices */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
                DRate = (double *)DCouponRatePtr + offset;
                Zero  = (double *)DZeroToPmtPtr + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {

                    /* Prepare for accessing 3-D slices */
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    BarrFXRate = (double *)BarrFxPtr + offset;
                    FRate  = (double *)FXFCouponRatePtr + offset;
                    Turbo  = (double *)TurboPtr + offset;

                    /* Domestic */
                    DPmt = DRate[j] * DDCFrac * DPrincForCoupon;


                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {

                        /* Foreign */
                        FPmt = FRate[k] * FDCFrac * FPrincForCoupon;

                        /* In principle, the pmt is dom + for */
                        TurboPmt = DPmt + FPmt;

                        /* Now deal with the binary clause */
                        /* First bracket the FX with the barrier levels */
                        BarrIdx = 0; 
                        BoundLeftRight = TRUE;
                        UpValue2 = 0.0; /* To avoid warnings */
                        for(l=0; l<NbBarriers; l++)
                        {
                            if (BarrFXRate[k] < BarrierLevel[l])
                            {
                                break;
                            }
                            BarrIdx++;
                        }

                      
                        /* Now calculate pay off for smoothed and */
                        /* non-smoothed cases                     */
                        if (Smoothing == 'N')
                        {
                            /* The weights are multiplicative */
                            TurboPmt *= BarrierWeight[BarrIdx];           

                            /* Finally apply cap and floor */
                            TurboPmt = MAXMIN(TurboPmt,CapAmt,FloorAmt);

                        }
                        else
                        {
                            
                            /* Establish UpValues and DownValues  */
                            /* for smoothing of the payoff        */
                            if (BarrIdx == 0)                             
                            {                                             
                                BoundLeftRight = FALSE;
                                Barrier1 = BarrierLevel[0];                   
                                DownValue1 = BarrierWeight[0];            
                                UpValue1 = BarrierWeight[1];              
                            }                                             
                            else if (BarrIdx == NbBarriers)                 
                            {                                             
                                BoundLeftRight = FALSE;
                                Barrier1 = BarrierLevel[NbBarriers - 1];                   
                                DownValue1 = BarrierWeight[NbBarriers - 1]; 
                                UpValue1   = BarrierWeight[NbBarriers];     
                            }                                             
                            else                                          
                            {                                             
                                Barrier1 = BarrierLevel[BarrIdx - 1]; 
                                DownValue1 = BarrierWeight[BarrIdx -1];   
                                UpValue1   = BarrierWeight[BarrIdx];
                                Barrier2 = BarrierLevel[BarrIdx];      
                                DownValue2 = BarrierWeight[BarrIdx];      
                                UpValue2   = BarrierWeight[BarrIdx + 1];  
                            }                                             
                                                                         

                            IndexStep = Hyb3_GetIndexStep(BarrFxPtr,
                                                     3,
                                                     i, j, k,
                                                     t,
                                                     tree_data);

                            /* First barrier */
                            if (Smooth_Step(&(SmoothAux),
                                    MAXMIN(TurboPmt*UpValue1,  CapAmt,FloorAmt),
                                    MAXMIN(TurboPmt*DownValue1,CapAmt,FloorAmt),
                                    BarrFXRate[k],
                                    Barrier1,
                                    IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }

                            
                            /* Second barrier, if required */
                            if (BoundLeftRight)
                            {
                                if (Smooth_Step(&(TurboPmt),
                                         MAXMIN(TurboPmt*UpValue2,CapAmt,FloorAmt),
                                         SmoothAux,
                                         BarrFXRate[k],
                                         Barrier2,
                                         IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }
                            else
                            {
                                TurboPmt = SmoothAux;
                            }
                        }
                        
                        /* Add payment to total turbo */
                        if (ArrearsForAll == 'Y')
                        {
                            Turbo[k] += TurboPmt;
                        }
                        else
                        {
                            Turbo[k] += TurboPmt * Zero[j];
                        }


                    }  /* For k */

                } /* For j */

            } /* For i */
        
            
        }  /* if it is reset for all payments */

  

        status = SUCCESS;

      RETURN:

        return (status);



}  /* Hyb3_BTurbo_Flows */


/*****  Hyb3_FwdBTurbo_Flows  *********************************************************/
/**
 *      Calculates the  forward swap price on the lattice minus an offset
 *      level.  In other  words, this  function calculates  the intrinsic
 *      value of an option on the turbo fwd swap struck at the given off-
 *      set level.
 *
 *      IMPORTANT ASSUMPTION:  This  forward calculator does  not support
 *      stubs. An  exercise must be an accrual start or else fall earlier
 *      than the absolute first accrual start of the turbo. 
 *
 *      Assumptions are  made concerning the  dimensions of the variables
 *      being priced (see Hyb3_Turbo_Price above).
 *
 *      In addition,  the  foreign coupon rate must  be pre-multiplied by
 *      the FX rate.   Both the foreign and the domestic coupon rates are
 *      assumed to be DEV'd outside the code in case the arrears turbo is
 *      being priced.
 *
 *      NOTE: this function is similar to Hyb3_FwdBTurbo_t, but it also remove
 *      dom and for principal pmts if it has been added on current date
 */


 int  Hyb3_FwdBTurbo_Flows(
		  TSLICE     FwdTurboPtr,         /**< (O) The forward     */
                  TSLICE     UlTurboPtr,          /**< (I) Underlying, no stubs*/

                  TSLICE     DCouponRatePtr,      /**< (I) Domestic rate       */
                  double     DPrincForCouponPaid, /**< (I) Domestic principal  */
                  double     DDCFracPaid,         /**< (I) Day count fraction  */

                  TSLICE     FXFCouponRatePtr,    /**< (I) For rate X FX rate  */
                  double     FPrincForCouponPaid, /**< (I) Foreign principal   */
                  double     FDCFracPaid,         /**< (I) Foreign day count   */
                  int        NeedToRemovePmtFlag, /**< (I) need to remove reset flag    */

                  double     FloorRate,           /**< (I) Floor level as rate */
                  double     CapRate,             /**< (I) Cap level as rate   */
                  double     PrincForCapFloor,    /**< (I) Princ for cap amount*/
                  double     DCFracCapFloor,      /**< (I) Dcc for cap amount  */

                  int        DPrincipalFlag,      /**< (I) Dom Prin Flag       */
                  double     DPrincipal,          /**< (I) Dom Prin pmt        */
                  int        FPrincipalFlag,      /**< (I) For Prin Flag       */
                  double     FPrincipal,          /**< (I) For Prin pmt        */
                
                  TSLICE     PrinFxPtr,           /**< (I) FX rates for foreign prin */
                  TSLICE     BarrFxPtr,           /**< (I) FX rates for barr   */

                  int         NbBarriers,         /**< (I) Nb FX barr levels   */
                  double     *BarrierLevel,       /**< (I) Marking FX ranges   */
                  double     *BarrierWeight,      /**< (I) Weights for ranges  */
                  char       Smoothing,           /**< (I) Smothing on/off     */
                  
                  int        ExerFlag,            /**< (I) True if exercise    */
                  double     Strike,              /**< (I) Strike amount       */

                  int        t,                   /**< (I) Curr timepoint      */
                  int        T,                   /**< (I) Total nb of t'points*/
                  int        DCurve,              /**< (I) Discount curve      */
                  HYB3_DEV_DATA   *dev_data,      /**< (I) Internal data       */
                  HYB3_TREE_DATA  *tree_data)     /**< (I) Tree data           */
{

    /* Locals for slice addressing */
    double     *FwdTurbo;
    double     *UlTurbo;   
    double     *DRate;     
    double     *FRate;
    double     *BarrFXRate;
    double     *PrinFXRate;


    double      SmoothAux;     /* Auxiliary for the smoothing algorithm  */
    double      IndexStep;     /* Used in the smoothing algorithm        */

    double      TurboPmtFull;  /* Full turbo pmt just added              */
    double      FPmtFull;      /* Foreign portion of hte above           */ 
    double      DPmtFull;      /* Domestic portion of the above          */
    double      FloorAmt;
    double      CapAmt;


    int         BarrIdx;              /* Index for XF barrier look-up    */
    double      Barrier1 = 0.0;       /* Barrier levels applicable to    */
    double      Barrier2 = 0.0;       /* current FX level                */
    double      UpValue1, DownValue1; /* P/off above and below barrier 1 */
    double      UpValue2, DownValue2; /* P/off above and below barrier 2 */

    int         BoundLeftRight;       /* TRUE if FX level in between two */
                                      /* specified barrier values        */


    int         Top1,   Bottom1;   /* Tree limits (1rst dim)             */
    int        *Top2,  *Bottom2;   /* Tree limits (2nd dim)              */
    int       **Top3, **Bottom3;   /* Tree limits (3rd dim)              */
    int         i, j, k, l;        /* Node indices                       */
    int         offset;
    int         status = FAILURE;  /* Error status                       */



        
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];




    /* Must check if DMode is compatible with underlying tree */
    if (tree_data->TreeType != TTYPE_FX2IR)
    {
        DR_Error("Inconsistent tree type for turbo pricing! "
                "(Hyb3_FwdTurbo_t)\n");
        goto RETURN;
    }
  


    /* If this is not an exercise date discount and return */
    if (!ExerFlag)                                              
    {

        if (Hyb3_Dev (FwdTurboPtr,  /* Disc expd value of a swap starting on the  */
                 t,            /* last exercise date                         */
                 T,
                 DCurve,
                 DISC_3D_CUPS,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
            
        }  /* if */
    
        status = SUCCESS;
        goto RETURN;
        
    }  /* if */


   
      
    /* "Basic" intrinsic value */
    for (i = Bottom1; i <= Top1; i ++)
    {
        for (j = Bottom2[i]; j <= Top2[i]; j++)
        {

            offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
            FwdTurbo = (double *)FwdTurboPtr + offset;
            UlTurbo  = (double *)UlTurboPtr + offset;

            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
            {
                FwdTurbo[k] = UlTurbo[k] - Strike;
            }  
        }
    }
   

    /* IMPORTANT: Now we must trust the caller to have checked that */
    /*            exercises DO NOT result in stubs!                 */

  

    /* If Need to remove pmt flag is ON */
    /* then the recently added coupon and principal pmts must be removed */
    if (NeedToRemovePmtFlag)
    {
        FloorAmt    = FloorRate * PrincForCapFloor * DCFracCapFloor;
        CapAmt      = CapRate   * PrincForCapFloor * DCFracCapFloor;
  
        for (i = Bottom1; i <= Top1; i ++)
        {
            /* Prepare to address the 2-D slices */
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DRate = (double *)DCouponRatePtr + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                /* Domestic */
                DPmtFull = DRate[j] * DDCFracPaid * DPrincForCouponPaid;


                /* Prepare to address the 3-D slices */
                offset   = Hyb3_Node_Offset(3,i,j,t,tree_data);
                BarrFXRate   = (double *)BarrFxPtr + offset;
                PrinFXRate   = (double *)PrinFxPtr + offset;
                FRate    = (double *)FXFCouponRatePtr + offset;
                FwdTurbo = (double *)FwdTurboPtr + offset;


                for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                {

                    /* Foreign */
                    FPmtFull = FRate[k] * FDCFracPaid * FPrincForCouponPaid;


                    /* In principle, pmt is just the domestic plus the foreign */
                    TurboPmtFull = DPmtFull + FPmtFull;

                    /* Now deal with the binary clause */
                    /* First bracket the FX with the barrier levels */
                    BarrIdx = 0; 
                    BoundLeftRight = TRUE;
                    UpValue2 = 0.0; /* To avoid warnings */
                    for(l=0; l<NbBarriers; l++)
                    {
                        if (BarrFXRate[k] < BarrierLevel[l])
                        {
                            break;
                        }
                        BarrIdx++;
                    }



                    /* Now calculate pay off for smoothed and */
                    /* non-smoothed cases                     */
                    if (Smoothing == 'N')
                    {

                        TurboPmtFull *= BarrierWeight[BarrIdx];

                        /* Finally apply cap and floor */
                        TurboPmtFull = MAXMIN(TurboPmtFull, CapAmt
                                          , FloorAmt);
                    }
                    else
                    {

                        /* Establish UpValues and DownValues  */
                        /* for smoothing of the payoff        */
                        if (BarrIdx == 0)                             
                        {                                             
                            BoundLeftRight = FALSE;
                            Barrier1 = BarrierLevel[0];                   
                            DownValue1 = BarrierWeight[0];            
                            UpValue1 = BarrierWeight[1];              
                        }                                             
                        else if (BarrIdx == NbBarriers)                 
                        {                                             
                            BoundLeftRight = FALSE;
                            Barrier1 = BarrierLevel[NbBarriers - 1];                   
                            DownValue1 = BarrierWeight[NbBarriers - 1]; 
                            UpValue1   = BarrierWeight[NbBarriers];     
                        }                                             
                        else                                          
                        {                                             
                            Barrier1 = BarrierLevel[BarrIdx - 1]; 
                            DownValue1 = BarrierWeight[BarrIdx -1];   
                            UpValue1   = BarrierWeight[BarrIdx];
                            Barrier2 = BarrierLevel[BarrIdx];      
                            DownValue2 = BarrierWeight[BarrIdx];      
                            UpValue2   = BarrierWeight[BarrIdx + 1];  
                        }                                             
                                     


                        IndexStep = Hyb3_GetIndexStep(BarrFxPtr,
                                                 3,
                                                 i, j, k,
                                                 t,
                                                 tree_data);

                        
                        if (Smooth_Step(&(SmoothAux),
                                MAXMIN(TurboPmtFull*UpValue1,  CapAmt,FloorAmt),
                                MAXMIN(TurboPmtFull*DownValue1,CapAmt,FloorAmt),
                                BarrFXRate[k],
                                Barrier1,
                                IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        
                        /* Consider second barrier, if applicable */
                        if (BoundLeftRight)
                        {
                            if (Smooth_Step(&(TurboPmtFull),
                                    MAXMIN(TurboPmtFull*UpValue2,CapAmt,FloorAmt),
                                    SmoothAux,
                                    BarrFXRate[k],
                                    Barrier2,
                                    IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                        else
                        {
                            TurboPmtFull = SmoothAux;
                        }
                    }
                    
                    FwdTurbo[k] -= TurboPmtFull;

                    if (DPrincipalFlag)
                    {
                        FwdTurbo[k] -= DPrincipal;
                    }

                    if (FPrincipalFlag)
                    {
                        FwdTurbo[k] -= FPrincipal * PrinFXRate[k];
                    }

                } /* For k */

            }  /* For j */

        }   /* For i */            
        
    }




    status = SUCCESS;

  RETURN:

    return (status);



}  /* Hyb3_FwdBTurbo_Flows */



/*****  Hyb3_BTurbo_Flows_Pmt  ****************************************************/
/**                                                                         
         Price of one binary turbo payment including amortisations             
                                                                               
         Similar to Hyb3_BTurbo_Flows, but no dev is done here. Also,               
         anything already stored in TurboPtr will be overwritten.              
                                                                             */
int  Hyb3_BTurbo_Flows_Pmt(
		  TSLICE      TurboPtr,         /**< (O) The binary turbo pmt  */
                  /* Domestic part of formula */       
                  TSLICE      DCouponRatePtr,   /**< (I) Domestic coupon rate  */
                  double      DPrincForCoupon,  /**< (I) Domestic principal    */
                  double      DDCFrac,          /**< (I) Domestic day count fr */
                  /* Foreign part of formula  */
                  TSLICE      FXFCouponRatePtr, /**< (I) For'gn cp X FX rate   */
                  double      FPrincForCoupon,  /**< (I) Foreign principal     */
                  double      FDCFrac,          /**< (I) Foreign day count fr  */
                  /* Common for all  */
                  double      FloorRate,        /**< (I) Floor level as rate   */
                  double      CapRate,          /**< (I) Cap level as rate     */
                  double      PrincForCapFloor, /**< (I) Princ for cap amount  */
                  double      DCFracCapFloor,   /**< (I) Dcf for cap amount    */

                  TSLICE      DZeroToPmtPtr,    /**< (I) Dom zero to turbo pmt */

                  double      DPrincipal,       /**< (I) Domestic princ pmt    */
                  TSLICE      DZeroToPrinPmtPtr, /**< (I) For zero to prin pmt */

                  double      FPrincipal,       /**< (I) Foreign princ pmt     */
                  TSLICE      FZeroToPrinPmtPtr, /**< (I) Dom zero to prin pmt */
                  TSLICE      PrinFxPtr,        /**< (I) FX rates for foreign prin */
                  /* Binary */
                  TSLICE      BarrFxPtr,        /**< (I) FX rates for fx barrier */
                  int         NbBarriers,       /**< (I) Nb FX barrier levels  */
                  double     *BarrierLevel,     /**< (I) Marking FX ranges     */
                  double     *BarrierWeight,    /**< (I) Weights for ranges    */
                  char        Smoothing,        /**< (I) TRUE if smoothing on  */

                  int         t,                /**< (I) Current time point    */
                  HYB3_TREE_DATA  *tree_data)        /**< (I) Tree data             */
                  
                  
{


        /* Locals for addressing the slices */
        double    *Turbo;
        double    *FRate;  
        double    *BarrFXRate; 
        double    *PrinFXRate; 

        double    *Zero;   
        double    *DZero;   
        double    *FZero;   
        double    *DRate;  


        double     FloorAmt;
        double     CapAmt;

        double
                   DPmt,  /* Domestic payment */
                   FPmt,  /* Foreign payment  */
                   TurboPmt,
                   SmoothAux,
                   IndexStep;

        double     Barrier1 = 0.0;       /* Barrier levels applicable to    */
        double     Barrier2 = 0.0;       /* current FX level                */
        double     UpValue1, DownValue1; /* P/off above and below barrier 1 */
        double     UpValue2, DownValue2; /* P/off above and below barrier 2 */

        int        BoundLeftRight;       /* True when FX falls between barr */

        int
                Top1,      Bottom1,      /* Limits of the tree (1rst dim)   */
               *Top2,     *Bottom2,
              **Top3,    **Bottom3,
                                         /* Index used in search of barrier */
                BarrIdx,
                offset,
                i, j, k, l,             /* Node indices                     */
                status = FAILURE;       /* Error status = FAILURE initially */

        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];	                



        /* Manage principal payments for foreign and domestic */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        FZero = (double *)FZeroToPrinPmtPtr + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            /* Prepare for accessing 2-D slices */
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DZero  = (double *)DZeroToPrinPmtPtr + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                Turbo = (double *)TurboPtr + offset;
                PrinFXRate = (double *)PrinFxPtr + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {
                    Turbo[k] = FPrincipal * PrinFXRate[k] * FZero[i]
                               + DPrincipal * DZero[j];

                }  /* for k */
            }
        }	


        /* add turbo payments */
        FloorAmt = FloorRate * PrincForCapFloor * DCFracCapFloor;
        CapAmt   = CapRate   * PrincForCapFloor * DCFracCapFloor;

        for (i = Bottom1; i <= Top1; i ++)
        {
            
            /* Prepare for accessing 2-D slices */
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            DRate = (double *)DCouponRatePtr + offset;
            Zero  = (double *)DZeroToPmtPtr + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {

                /* Prepare for accessing 3-D slices */
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                BarrFXRate = (double *)BarrFxPtr + offset;
                FRate  = (double *)FXFCouponRatePtr + offset;
                Turbo  = (double *)TurboPtr + offset;

                /* Domestic */
                DPmt = DRate[j] * DDCFrac * DPrincForCoupon;


                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {

                    /* Foreign */
                    FPmt = FRate[k] * FDCFrac * FPrincForCoupon;

                    /* In principle, the pmt is dom + for */
                    TurboPmt = DPmt + FPmt;

                    /* Now deal with the binary clause */
                    /* First bracket the FX with the barrier levels */
                    BarrIdx = 0; 
                    BoundLeftRight = TRUE;
                    UpValue2 = 0.0; /* To avoid warnings */
                    for(l=0; l<NbBarriers; l++)
                    {
                        if (BarrFXRate[k] < BarrierLevel[l])
                        {
                            break;
                        }
                        BarrIdx++;
                    }

                  
                    /* Now calculate pay off for smoothed and */
                    /* non-smoothed cases                     */
                    if (Smoothing == 'N')
                    {
                        /* The weights are multiplicative */
                        TurboPmt *= BarrierWeight[BarrIdx];           

                        /* Finally apply cap and floor */
                        TurboPmt = MAXMIN(TurboPmt,CapAmt,FloorAmt);

                    }
                    else
                    {
                        
                        /* Establish UpValues and DownValues  */
                        /* for smoothing of the payoff        */
                        if (BarrIdx == 0)                             
                        {                                             
                            BoundLeftRight = FALSE;
                            Barrier1 = BarrierLevel[0];                   
                            DownValue1 = BarrierWeight[0];            
                            UpValue1 = BarrierWeight[1];              
                        }                                             
                        else if (BarrIdx == NbBarriers)                 
                        {                                             
                            BoundLeftRight = FALSE;
                            Barrier1 = BarrierLevel[NbBarriers - 1];                   
                            DownValue1 = BarrierWeight[NbBarriers - 1]; 
                            UpValue1   = BarrierWeight[NbBarriers];     
                        }                                             
                        else                                          
                        {                                             
                            Barrier1 = BarrierLevel[BarrIdx - 1]; 
                            DownValue1 = BarrierWeight[BarrIdx -1];   
                            UpValue1   = BarrierWeight[BarrIdx];
                            Barrier2 = BarrierLevel[BarrIdx];      
                            DownValue2 = BarrierWeight[BarrIdx];      
                            UpValue2   = BarrierWeight[BarrIdx + 1];  
                        }                                             
                                                                     

                        IndexStep = Hyb3_GetIndexStep(BarrFxPtr,
                                                 3,
                                                 i, j, k,
                                                 t,
                                                 tree_data);

                        /* First barrier */
                        if (Smooth_Step(&(SmoothAux),
                                MAXMIN(TurboPmt*UpValue1,  CapAmt,FloorAmt),
                                MAXMIN(TurboPmt*DownValue1,CapAmt,FloorAmt),
                                BarrFXRate[k],
                                Barrier1,
                                IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }

                        
                        /* Second barrier, if required */
                        if (BoundLeftRight)
                        {
                            if (Smooth_Step(&(TurboPmt),
                                     MAXMIN(TurboPmt*UpValue2,CapAmt,FloorAmt),
                                     SmoothAux,
                                     BarrFXRate[k],
                                     Barrier2,
                                     IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                        else
                        {
                            TurboPmt = SmoothAux;
                        }
                    }
                    
                    /* Add payment to total turbo */
                    Turbo[k] += TurboPmt * Zero[j];

                }  /* For k */

            } /* For j */

        } /* For i */

  

        status = SUCCESS;

      RETURN:

        return (status);



}  /* Hyb3_BTurbo_Flows_Pmt */
 
