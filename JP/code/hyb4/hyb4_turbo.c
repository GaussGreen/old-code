/****************************************************************************/
/*      Calculation routines customised for a turbo swap                    */
/****************************************************************************/
/*      TURBO.C                                                             */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"
 
/*****  Hyb4_Turbo_Price  *********************************************************/
/**                                                                             
         Price of a turbo swap including different amortisation schedules      
         for the domestic and foreign payments.                                
                                                                               
         Assumes dimensions as follows:                                        
             FX rates :              3-Dim                                     
             Domestic index:         2-Dim                                     
             Foreign index times FX: 3-Dim                                     
             Zero to payment:        2-Dim                                     
                                                                             */
int   Hyb4_Turbo_Price(TSLICE      TurboPtr,         /**< (O) The turbo             */
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
                  HYB4_DEV_DATA   *dev_data,         /**< (I) Internal data         */
                  HYB4_TREE_DATA  *tree_data)        /**< (I) Tree data             */
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

        
        Top1    = tree_data->iMax[t];
        Top2    = tree_data->jMax[t];
        Top3    = tree_data->kMax[t];
        Bottom1 = tree_data->iMin[t];
        Bottom2 = tree_data->jMin[t];
        Bottom3 = tree_data->kMin[t];	                


    /* Discounted expected value function */        
    if (Hyb4_Dev(TurboPtr,         
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
                    offset = tree_data->NodeOffset2[t][i][j];
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
                    offset = tree_data->NodeOffset2[t][i][j];
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
                offset = tree_data->NodeOffset1[t][i];
                DRate = (double *)DCouponRatePtr + offset;
                Zero  = (double *)DZeroToPmtPtr + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {

                    /* Prepare for accessing 3-D slices */
                    offset = tree_data->NodeOffset2[t][i][j];
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



}  /* Hyb4_Turbo_Price */
 
int Hyb4_FwdTurbo_t(TSLICE     FwdTurboPtr,         /**< (O) The forward         */
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
                  HYB4_DEV_DATA  *dev_data,            /**< (I) Internal data       */
                  HYB4_TREE_DATA *tree_data)           /**< (I) Tree data           */
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



        
    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];

    /* VEZI: Must check if DMode is compatible with underlying tree */

    /* If this is not an exercise date discount and return */
    if (!ExerFlag)                                              
    {

        if (Hyb4_Dev (FwdTurboPtr,  /* Disc expd value of a swap starting on the  */
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

            offset = tree_data->NodeOffset2[t][i][j];
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
            DR_Error("Unable to calculate dcf for dom pmt.(Hyb4_FwdTurbo_t)");
            goto RETURN;
        } 

        if (DrDayCountFraction(AccStartPaid,
                               CurrentDate,                                 
                               DPmtDayCount,
                               &DDayCountAcc) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for foreign pmt.(Hyb4_FwdTurbo_t)");
            goto RETURN;
        }  
        
        if (DrDayCountFraction(AccStartPaid,
                               CurrentDate,                                 
                               CapFloorDayCount,
                               &CapFloorDayCountAcc) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for cap/floor.(Hyb4_FwdTurbo_t)");
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
                    offset = tree_data->NodeOffset1[t][i];
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
                        offset = tree_data->NodeOffset2[t][i][j];
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
                    offset = tree_data->NodeOffset1[t][i];
                    DRate = (double *)DCouponRatePtr + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {

                        /* Domestic */
                        RateAux =  DRate[j] * (DDCFracPaid - DDayCountAcc);
                        DPmtStub = RateAux * DPrincForCouponPaid;

                        /* Prepare to address the 3-D slices */
                        offset = tree_data->NodeOffset2[t][i][j];
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
                        "(Hyb4_FwdTurbo_t)\n");
                        goto RETURN;
            }
            else   
            {
                
                for (i = Bottom1; i <= Top1; i ++)
                {

                    /* Prepare to address the 2-D slices */
                    offset = tree_data->NodeOffset1[t][i];
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
                        offset = tree_data->NodeOffset2[t][i][j];
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



}  /* Hyb4_FwdTurbo_t */
