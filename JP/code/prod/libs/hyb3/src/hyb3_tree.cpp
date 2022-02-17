/****************************************************************************/
/*      Building of the tree: calculation of forward rates, spot volatility */
/*      stock forward volatility, interest rate drift at each time step and */
/*      tree limits.                                                        */
/****************************************************************************/
/*      TREE.c                                                              */
/****************************************************************************/

/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cupslib.h"

/*****  Hyb3_Build_Tree  **********************************************************/
/**
*       Main routine for the construction of the tree.
*/
int  Hyb3_Build_Tree 
        (int              CetOutputFlag,   /**< (I) On TRUE write CET to file  */
         T_CURVE          t_curve[2][3],   /**< (I) Zero curve data            */
         MKTVOL_DATA     *mktvol_data,     /**< (I) Market volatility data     */
         FX_DATA         *fx_data,         /**< (I) FX data (needed for CUPS)  */
         EQ_DATA         *eq_data,         /**< (I) Equity data (may be NULL)  */
         HYB3_TREE_DATA       *tree_data)  /**< (O) Tree data                  */

{

           /* The following six flags are  initialised to their most usual    */
           /* values and changed as required in the switch statement below.   */
    int    FxOn         = TRUE;     /* Flag to indicate whether FX relevant   */
    int    ThirdFXDimOn = FALSE;
    int    EquityOn   = FALSE;    /* Flag to specify the third asset          */
    int    NbIRatesOn = 2;        /* 1 for TTYPE_EQ1IR, 2 for all others      */
    double LRho[5];               /* Local var. for extension of correlations */ 
    int    RhoFxIrFacIndex[2];      /* Index of the 1st Fx/IR fact corr       */
                                    /* in Rho vector                          */ 


    int    FxMomentMatchOn;       /* AlexK: moment matching procedure       */
    
    int    EquityOnIr1 = FALSE;   /* Flag eq denominated on the 1st IR      */
    int    EquityOnIr2 = FALSE;   /* Flag eq denominated on the 2nd IR      */
    int    EquityComp  = FALSE;   /* Flag that specifies composite equity   */

    int    CvDiff[2];             /* Curve indices (three curves modelled)  */
    int    CvIdx1[2];             /* for each one fo the currencies  (dom   */
    int    CvIdx2[2];             /* and foreign), plus discount curve.     */
    int    CvDisc[2];   

    int    T  = tree_data->NbTP;
    int    xT = 0;

    int    i;                     /* Index                                  */
    int    j;                     /* Index                                  */
    int    status  = FAILURE;     /* Error status = FAILURE initially       */
  
    double  temp;

    /* error handling */
    char  ErrorMsg[MAXBUFF];
    char  CcyErrString[2][MAXBUFF];        /* defines the name for the tree dimensions */

    /* find the cutoff date index        */
    /* if not found, cutoff index = NbTP */

    while ((xT < T) &&
           (tree_data->TPDate[xT] != tree_data->xDate))
    {
        ++xT;
    }

    tree_data->xT = xT;

    /* disallow composite equity mode */

    if (tree_data->TreeType == TTYPE_EQC2IR)
    {
        DR_Error("Composite equity mode (EQC2IR) is disallowed.");
        goto RETURN;
    }

    /* if NbSigmaDBL has not been initialised yet, then set it to input Sigma*/
    if(tree_data->NbSigmaDBL < 0)
    {
        tree_data->NbSigmaDBL = (double)tree_data->NbSigmaMax;
    }

    /* initialise NbIRSigma to input NbSigma regardless of the dim of the tree  */
    /* The only case where NbIRSigma gets recomputed is if we are running a true*/
    /* FX tree                                                                  */
    tree_data->NbIRSigmaMax = tree_data->NbSigmaDBL;

    /* we apply a moment matching algorithm on the FX tree */
    FxMomentMatchOn  = tree_data->FxMomentMatching;

    /*  TreeType will  tell us if the tree is being  run in 2-D  mode  */
    /*  (i.e. two interest rates) or 3-D (i.e. plus equity or plus FX) */
    sprintf(CcyErrString[0] ,"%s","Foreign");    /* important: make sure that NbIRatesOn is <= 2 */
    sprintf(CcyErrString[1] ,"%s","Domestic");

    switch(tree_data->TreeType)
    {

        case TTYPE_1IR:
            NbIRatesOn = 1;
            FxOn = FALSE;
            sprintf(CcyErrString[0] ,"%s","Domestic");
            break;
        case TTYPE_1IR2F:
            NbIRatesOn  = 1;
            FxOn = FALSE;
            sprintf(CcyErrString[0] ,"%s","Domestic");
            break;
        case TTYPE_2IR:            
            break;
        case TTYPE_2IR2F1D:
            break;
        
        case TTYPE_FX2IR:
            ThirdFXDimOn = TRUE;
            break;
        
        case TTYPE_EQD2IR:
            EquityOn     = TRUE;
            EquityOnIr2  = TRUE;
            break;
        
        case TTYPE_EQF2IR:
            EquityOn     = TRUE;
            EquityOnIr1  = TRUE;
            break;
        
        case TTYPE_EQC2IR:
            EquityOn     = TRUE;
            EquityComp   = TRUE;
            break;
        
        case TTYPE_EQ1IR:
            EquityOn     = TRUE;
            NbIRatesOn   = 1;
            FxOn         = FALSE;
            EquityOnIr1  = TRUE;
            sprintf(CcyErrString[0] ,"%s","Domestic");
            break;
        
        case TTYPE_EQDFX2IR:
            NbIRatesOn   = 2;
            FxOn         = TRUE;
            ThirdFXDimOn = TRUE;
            EquityOn     = TRUE;
            EquityOnIr2  = TRUE;
            break;

        case TTYPE_EQFFX2IR:
            NbIRatesOn   = 2;
            FxOn         = TRUE;
            ThirdFXDimOn = TRUE;
            EquityOn     = TRUE;
            EquityOnIr1  = TRUE;
            break;
    }

    /* Check coherence of Nb of factor and tree geometry */
    switch (mktvol_data[0].NbFactor)
    {
        case 1:
            if ((tree_data->TreeType==TTYPE_1IR2F)||
                (tree_data->TreeType==TTYPE_2IR2F1D))
            {
                DR_Error("Nb of foreign IR factors incompatible with"
                         " tree geometry!(Hyb3_Build_Tree)\n");
                goto RETURN;
            }
            break;
        case 2:
            if ((tree_data->TreeType!=TTYPE_1IR2F)&&
                (tree_data->TreeType!=TTYPE_2IR2F1D))
            {
                DR_Error("Nb of IR factors incompatible with"
                         " tree geometry!(Hyb3_Build_Tree)\n");
                goto RETURN;
            }
            break;
        default:
            DR_Error("Nb of foreign IR factors must be 1 or 2!"
                     "(Hyb3_Build_Tree)\n");
            goto RETURN;
           
    }
    if ( (NbIRatesOn == 2)&&
        (mktvol_data[1].NbFactor != 1))
    {
        DR_Error("Only one domestic IR factor allowed!"
                 "(Hyb3_Build_Tree)\n");
        goto RETURN;
    }


    /* reconstruct the input flat files (useful for QLib debugging) */
    if (tree_data->debug)
    {
        if (ProduceFlatMarketFiles(t_curve, 
                                   mktvol_data, 
                                   fx_data,
                                   eq_data, 
                                   tree_data) != SUCCESS )
        {
            DR_Error("Hyb3_BuildTree:: failed to publish debug flat file information!");
            goto RETURN;
        }
    }


    /* Note that the swaption volatility is first bootstrapped */
    /* using NORMAL smile parameters so that the bootstrapped   */
    /* Fx volatility is independent of the smile and CET.      */
    /* We then perform CET and repeat the swaption volatility  */
    /* bootrapping, but this time using smile and CET          */  



    /* For all interest rates on, do the forwards and the vols */
    for (i = 0; i < NbIRatesOn; i++)   /* 0=Foreign, 1=Domestic */
    {

        CvDiff[i] = tree_data->CvDiff[i];
        CvIdx1[i] = tree_data->CvIdx1[i];
        CvIdx2[i] = tree_data->CvIdx2[i];
        CvDisc[i] = tree_data->CvDisc[i]; 

        for (j = 0; j < 3; j++)	  /* Index, cost of fund and credit curves */
        {
            /* Calculate forwards and zeros */
            if (Hyb3_Forward_Curve(tree_data->ZeroCoupon[i][j],
                              tree_data->ZeroRate[i][j],
                              tree_data->FwdRate[i][j],
                              &t_curve[i][j],
                              tree_data->TPDate,
                              tree_data->NbTP) == FAILURE)
            {
                sprintf (ErrorMsg, "Forward Curve Calculation failed for %s "
                            "dimension! (BuildTree)",CcyErrString[i]);
                DR_Error(ErrorMsg); 
                goto RETURN;
            }  /* if */ 
        }  /* for j */

        if (FxOn || EquityOn)
        {
        
            /* Bootstrap swaption volatility to get spot volatility curve */
            /* Use diffused t_curve and NORMAL smile parameters!          */
            if (Hyb3_SpotVol(mktvol_data[i].Aweight,
                        mktvol_data[i].BaseDate,
                        mktvol_data[i].NbVol,
                        mktvol_data[i].VolDate,
                        mktvol_data[i].Vol,
                        mktvol_data[i].VolUsed,
                        mktvol_data[i].Freq,
                        mktvol_data[i].DCC,
                        mktvol_data[i].SwapSt,
                        mktvol_data[i].SwapMat,
                        0.,                      /* "Normal" QLeft    */
                        0.,                      /* "Normal" QRight   */
                        0.,                      /* No forward shift  */
                        mktvol_data[i].NbFactor, /* Number of factors */
                        mktvol_data[i].Alpha,
                        mktvol_data[i].Beta,
                        mktvol_data[i].Rho,
                        mktvol_data[i].SkipFlag,
                        mktvol_data[i].CalibFlag,
                        &t_curve[i][CvDiff[i]]) != SUCCESS)
            {
                sprintf (ErrorMsg, "IR Spot Vol Function failed for %s "
                            "dimension! (BuildTree)",CcyErrString[i]);
                DR_Error(ErrorMsg); 
                goto RETURN;
            }

            /* Interpolate spot volatility curves */
            if (Hyb3_Interp_SpotVol (tree_data->IrAweight[i],
                                mktvol_data[i].BaseDate,
                                mktvol_data[i].NbVol,
                                mktvol_data[i].VolDate,
                                mktvol_data[i].Aweight,
                                mktvol_data[i].NbFactor,   
                                mktvol_data[i].Beta,
                                mktvol_data[i].CalibFlag,
                                tree_data->FwdRate[i][CvDiff[i]],
                                tree_data->Length,
                                tree_data->NbTP) == FAILURE)
            {
                sprintf (ErrorMsg, "Interpolatation of IR Spot Vol Curve failed for %s "
                            "dimension! (BuildTree)",CcyErrString[i]);
                DR_Error(ErrorMsg); 
                goto RETURN;
            }
            for (j = -1; j <= tree_data->NbTP; j++)
            {
                /*   Spot vol is stored in "Aweight" form, i.e spotvol=sigma*alpha(0)    */
                /*   for several factor IR, it is the Spot vol of the first factor       */
                                                               
                tree_data->SpotVol[i][j]=tree_data->IrAweight[i][0][j];
            }

        }/* if FxOn or EqOn */
    } /* for i */

    /* Interpolate all correlations alongside timeline in preparation for */
    /* CUPS drift and for volatility calculations (may be 1,3 or 5 curves)*/
    {
        long    l;

        if(FxOn)
        {
            RhoFxIrFacIndex[0]= mktvol_data[0].NbFactor * mktvol_data[1].NbFactor;
            RhoFxIrFacIndex[1]=RhoFxIrFacIndex[0]+mktvol_data[0].NbFactor;


            /* correlation between the two IR factors*/

            if( Hyb3_IR_IR_Corr (LRho,
                fx_data->Rho[0][0],
                mktvol_data[1].Rho,
                mktvol_data[0].Rho,
                mktvol_data[1].NbFactor,
                mktvol_data[0].NbFactor, 
                mktvol_data[1].Alpha,
                mktvol_data[0].Alpha,
                mktvol_data[1].Beta, 
                mktvol_data[0].Beta, 
                1,  /* 1 because corresl between exponential facs*/
                mktvol_data[1].CorrSwapSt,
                mktvol_data[1].CorrSwapMat,
                mktvol_data[0].CorrSwapSt, 
                mktvol_data[0].CorrSwapMat,
                mktvol_data[1].CorrSwapDCC,  
                mktvol_data[0].CorrSwapDCC,  
                mktvol_data[1].CorrSwapFreq, 
                mktvol_data[0].CorrSwapFreq, 
                fx_data->ValueDate, 
                fx_data->ValueDate,
                &t_curve[1][CvDiff[1]],
                &t_curve[0][CvDiff[0]]) == FAILURE)
            {
                DR_Error("IR_IR_Corr failed in Build_Tree");
                goto RETURN;
            }

            for(i = 0; i < NbIRatesOn; i++)
            {
                /* correlation fx vs IR factors */
            if( Hyb3_AssetIrCorr (&(LRho[RhoFxIrFacIndex[i]]), 
                fx_data->Rho[i+1][0],
                mktvol_data[i].Rho,
                mktvol_data[i].NbFactor,
                mktvol_data[i].Alpha,
                mktvol_data[i].Beta, 
                1,/* 1 because correl bewteen exponential facs*/
                mktvol_data[i].CorrSwapSt,
                mktvol_data[i].CorrSwapMat,
                mktvol_data[i].CorrSwapDCC,
                mktvol_data[i].CorrSwapFreq,
                fx_data->ValueDate,                                          
                &t_curve[i][CvDiff[i]]) == FAILURE)
            {
                DR_Error("Hyb3_AssetIrCorr failed for %s dimension!(BuildTree)",CcyErrString[i]);
                goto RETURN;
            }

            } /*for each IR*/
            


            for (l = -1; l <= tree_data->NbTP; l++)
            {
                for (i=0;i<(mktvol_data[0].NbFactor+1)*(mktvol_data[1].NbFactor+1)-1;i++)
                {
                tree_data->Rho[i][l] = LRho[i];
                }
            }
        }/* fx on */
        
        if (EquityOn)
        {   
            for (l = -1; l <= tree_data->NbTP; l++)
                tree_data->Rho[3][l] = eq_data->Rho[0][0];

            if(tree_data->TreeType != TTYPE_EQ1IR)
            {   
                for (l = -1; l <= tree_data->NbTP; l++)
                {
                    tree_data->Rho[4][l] = eq_data->Rho[1][0];
                    tree_data->Rho[5][l] = eq_data->Rho[2][0];
                }
            }
        }/* equity on*/
    }
    /* As long as there are two interest rates on, FX data is relevant.  */
    /* All the deterministic calculations (i.e. forwards, spot vols  and */
    /* total vols) related to the FX are performed regardless of whether */
    /* the third dimension is ON or is FX.  This is necessary because of */
    /* of the currency protected drift calculations.                     */

    if (FxOn)
    {   
        
        if(Hyb3_Tree_SmileParams(tree_data->SmileIndex,
                            tree_data->A1,
                            tree_data->A2,
                            tree_data->A3,
                            tree_data->A1C,
                            tree_data->A2C,
                            tree_data->A3C,
                            tree_data->TPDate,
                            tree_data->NbTP+1,
                            fx_data->a1,
                            fx_data->a2,
                            fx_data->a3,
                            fx_data->NbSmilePt,
                            fx_data->SmileDate, 
                            tree_data->tMin,
                            tree_data->tMax) == FAILURE)
        {
            DR_Error("unable to calculate required smile parameters in Build_tree\n");
            goto RETURN;
        }

        /* The last relevant tree time point in TPDate[NbTP], so that  */
        /* the last relevant vol point on the tree timeline is at      */
        /* TPDate[NbTP - 1]. Therefore to extended the spot vols up to */
        /* that point in ExtendSpotVol, we have to pass to             */
        /* Hyb3_Get_TreeFxSpotVol  (NbTP+1), as opposed to (NbTP)           */
   

        if(Hyb3_Get_TreeFxSpotVol(tree_data->SpotFxVol,
                             tree_data->TPDate,  
                             (tree_data->NbTP + 1), /* Cf comment above!*/
                             fx_data->FxVol,    
                             fx_data->VolDate,  
                             fx_data->NbVol,    
                             fx_data->InpSpotVol,
                             fx_data->InpSpotVolDate,
                             fx_data->NbInpSpotVol,
                             tree_data->SpotVol[1],
                             tree_data->SpotVol[0],
                             mktvol_data[1].Alpha,
                             mktvol_data[0].Alpha,
                             tree_data->FwdRate[1][CvDiff[1]],/*Diffuse curve*/
                             tree_data->FwdRate[0][CvDiff[0]],/*Diffuse curve*/
                             &(LRho[RhoFxIrFacIndex[1]]),     /*Cst corr FX/dom IR fact */
                             &(LRho[RhoFxIrFacIndex[0]]),     /*Cst corr FX/for IR fact */
                             &(LRho[0]),                      /*Cst corr dom IR /for IR fact */
                             mktvol_data[1].Rho,
                             mktvol_data[0].Rho,
                             fx_data->FxCutOffFlag,
                             fx_data->FxCutOffLast,
                             fx_data->FxCutOffLevel,
                             mktvol_data[1].Beta,
                             mktvol_data[0].Beta,
                             mktvol_data[1].NbFactor,                    /* Dom Nb Factor */
                             mktvol_data[0].NbFactor,                    /* For Nb Factor */
                             fx_data->ValueDate,
                             fx_data->FxBootStrapMode,

                             tree_data->A1C,
                             tree_data->A2C,
                             tree_data->A3C) == FAILURE)
        {   
            DR_Error("Hyb3_Get_TreeFxSpotVol has failed in Hyb3_Build_Tree\n");
            goto RETURN;
        }
 


        if (Hyb3_FX_Vol(tree_data->FxVol,
                   tree_data->SpotFxVol,
                   mktvol_data[1].NbFactor,
                   mktvol_data[0].NbFactor,
                   &(LRho[RhoFxIrFacIndex[1]]),
                   &(LRho[RhoFxIrFacIndex[0]]),
                   &(LRho[0]),
                   mktvol_data[1].Rho,
                   mktvol_data[0].Rho,
                   tree_data->SpotVol,
                   tree_data->FwdRate[1][CvDiff[1]],
                   tree_data->FwdRate[0][CvDiff[0]],
                   (mktvol_data[1].Alpha),
                   (mktvol_data[0].Alpha),
                   tree_data->TPDate,       /* &(fx_data->VolDate[1]),*/          
                   tree_data->TPDate,       /* &(fx_data->VolDate[1]),*/
                   tree_data->A1C,
                   tree_data->A2C,
                   tree_data->A3C,
                   (mktvol_data[1].Beta),                   
                   (mktvol_data[0].Beta),
                   tree_data->TPDate,          
                   (tree_data->NbTP+1)) == FAILURE)
        {
           goto RETURN;
        }


        if (Hyb3_Forward_FX(tree_data->FwdFx,
                       tree_data->NbTP,
                       tree_data->ZeroCoupon,
                       tree_data->CvDisc,     /* Use Discount curve*/
                       fx_data->Spot) == FAILURE)
        {
            goto RETURN;
        }



        if(Hyb3_FXMidNode(tree_data->FxMidNode,    /* (O) Centre of the tree in FX dim */
                     tree_data->NbTP,         /* (I) Total nb of nodes in timeline*/
                     tree_data->Length,       /* (I) Time step in years           */
                     tree_data->SpotFxVol,    /* (I) FX spotvol alongside timeline*/
                     tree_data->SmileIndex,
                     tree_data->FxVol,        /* new parameter  needed for estimated delta K drift adj */
                     tree_data->NbSigmaMax,
                     tree_data->A1,
                     tree_data->A2,
                     tree_data->A3) == FAILURE)
        {
            goto RETURN;
        }
        
        if(ThirdFXDimOn)
        {
            if(Hyb3_Find_IR_StdDev(&(tree_data->NbIRSigmaMax),
                              tree_data->NbSigmaDBL,
                              mktvol_data[1].Beta[0],
                              mktvol_data[0].Beta[0],
                              tree_data->SpotVol,
                              fx_data->Rho[0][0],
                              tree_data->NbTP,
                              tree_data->FwdRate,
                              tree_data->FxVol,
                              tree_data->Length,
                              tree_data->TPDate,
                              mktvol_data[1].QLeft,
                              mktvol_data[1].QRight,
                              mktvol_data[1].FwdShift,
                              mktvol_data[0].QLeft,
                              mktvol_data[0].QRight,
                              mktvol_data[0].FwdShift
                              ) == FAILURE) goto RETURN;

            /* check Fwd FX */
            tree_data->CalcCheckSlices = TRUE;
        }
    } /* Fx is on */


    /* There are three types of equity trees: domestic (TTYPE_EQD2IR),    */
    /* foreign cups (TTYPE_EQF2IR) and foreign (TTYPE_EQC2IR).            */
    /* The domestic case is straightforward. In the cups case we only     */
    /* have to add a drift and we do that by adjusting the center of the  */
    /* equity tree. The foreign case is more complex: we first bootstrapp */
    /* the equity vol and calculate the forward stock normally. But then  */
    /* we create an artificial domestic asset by multiplying the forward  */
    /* stock by the fx forward, replace the volatility by the composite   */
    /* volatility and adjust the correlations accordingly.                */

    if (EquityOn)
    {
        if(Hyb3_Tree_SmileParams(tree_data->SmileIndex,
                            tree_data->A1,
                            tree_data->A2,
                            tree_data->A3,
                            tree_data->A1C,
                            tree_data->A2C,
                            tree_data->A3C,
                            tree_data->TPDate,
                            tree_data->NbTP+1,
                            eq_data->a1,
                            eq_data->a2,
                            eq_data->a3,
                            eq_data->NbSmilePt,
                            eq_data->SmileDate, 
                            tree_data->tMin,
                            tree_data->tMax) == FAILURE)
        {
            DR_Error("unable to calculate required smile parameters in Build_tree\n");
            goto RETURN;
        }

        /* Calculate forward  stock and settlement dates */
        /* (Discount curves passed for the ccy of stock) */

        if (Hyb3_Forward_Stock (tree_data->FwdEq,
                           tree_data->NodeSettleDate,
                           tree_data->NodeSettleTime,
                           eq_data->ValueDate,
                           eq_data->Spot,
                           eq_data->NbFwd,
                           eq_data->Fwd,
                           eq_data->FwdDate,
                           eq_data->FwdType,
                           eq_data->NbSettle,
                           eq_data->LastTrading,
                           eq_data->SettleDate,
                           eq_data->SettleType,
                           eq_data->NbBorrow,
                           eq_data->BorrowDate,
                           eq_data->Borrow,
                           EquityOnIr2 ?
                             &t_curve[1][CvDisc[1]] :
                             &t_curve[0][CvDisc[0]],
                           tree_data->TPDate,
                           tree_data->Length,
                           EquityOnIr2 ?
                             tree_data->ZeroCoupon[1][CvDisc[1]] :
                             tree_data->ZeroCoupon[0][CvDisc[0]],
                           tree_data->NbTP) == FAILURE)
        {
            goto RETURN; 
        }

        if(Hyb3_Get_TreeEqSpotVol(tree_data->SpotEqVol,
                             tree_data->FwdEq,
                             tree_data->TPDate,
                             tree_data->NbTP+1,     
                             eq_data->Vol,
                             eq_data->VolDate,
                             eq_data->NbVol,  
                             eq_data->InpSpotVol,      /*(I)Input Spot Eq vols          */
                             eq_data->InpSpotVolDate,  /*(I)Input SpotVol Dates         */
                             eq_data->NbInpSpotVol,    /*(I)Nb of Input spot vols       */
                             EquityOnIr2 ?
                             tree_data->SpotVol[1] :
                             tree_data->SpotVol[0],    /*(I)dom ir (tree) Aweight       */
                             EquityOnIr2 ?
                             mktvol_data[1].Alpha :
                             mktvol_data[0].Alpha,     /*(I)dom factor weights          */
                             EquityOnIr2 ?
                             tree_data->FwdRate[1][CvDiff[1]] :/*Diffuse curve*/
                             tree_data->FwdRate[0][CvDiff[0]], /*Diffuse curve*/
                             EquityOnIr2 ?
                             tree_data->Rho[4] :
                             tree_data->Rho[3],
                             eq_data->EqCutOffFlag,    /*(I)                            */
                             eq_data->EqCutOffLast,    /*(I)                            */
                             eq_data->EqCutOffLevel,   /*(I)                            */
                             EquityOnIr2 ?
                             mktvol_data[1].Beta :
                             mktvol_data[0].Beta,
                             1,                        /*(I)                            */
                             tree_data->TPDate[0],     /*(I)                            */
                             eq_data->EqBootStrapMode,
                             tree_data->A1C,
                             tree_data->A2C,
                             tree_data->A3C) == FAILURE)
        {
            goto RETURN;
        }

        if(Hyb3_FXMidNode(tree_data->EqMidNode,         /* (O) Centre of the tree in Eq dim  */
                          tree_data->NbTP,              /* (I) Total nb of nodes in timeline */
                          tree_data->Length,            /* (I) Time step in years            */
                          tree_data->SpotEqVol,         /* (I) Eq spotvol alongside timeline */
                          tree_data->SmileIndex,
                          tree_data->FxVol,
                          tree_data->NbSigmaMax,
                          tree_data->A1,
                          tree_data->A2,
                          tree_data->A3) == FAILURE)
        {
            goto RETURN;
        }

        /* Adjust forward stock,volatility and correlations  */
        /* of composite equity                               */

        /* $$$ Must address this in presence of eq smile, if */
        /* we ever want to extend the smile formulation to   */
        /* the 3F cups mode.                                 */
        if (EquityComp == TRUE)
        {
            if (Hyb3_Composite_Eq(tree_data->FwdEq,
                             tree_data->EqMidNode,
                             tree_data->SpotEqVol,
                             tree_data->Rho,
                             tree_data->NodeSettleDate,
                             &t_curve[0][CvDisc[0]],

                             &t_curve[1][CvDisc[1]],
                             tree_data->FwdFx,
                             tree_data->SpotFxVol,
                             tree_data->Length,
                             tree_data->NbTP) == FAILURE)
            {
                goto RETURN;
            }  /* if */
        }  /* if */

        /* Calculate equity option volatility at each node in the tree */
        /* Now the composite equity behaves like a domestic one        */

        if(Hyb3_Eq_Vol( tree_data->EqVol,  
                        tree_data->SpotEqVol,
                        EquityOnIr1?
                            tree_data->Rho[3] :
                            tree_data->Rho[4],
                        EquityOnIr1?
                            tree_data->SpotVol[0] :
                            tree_data->SpotVol[1],
                        EquityOnIr1?
                            tree_data->FwdRate[0][CvDisc[0]] :
                            tree_data->FwdRate[1][CvDisc[1]],
                        EquityOnIr2 ?
                                  mktvol_data[1].Alpha :
                                  mktvol_data[0].Alpha,    /*(I)dom factor weights          */
                        tree_data->TPDate,                 /* &(fx_data->VolDate[1]),*/          
                        tree_data->A1C,
                        tree_data->A2C,
                        tree_data->A3C,
                        EquityOnIr1?
                        (mktvol_data[0].Beta) :                   
                        (mktvol_data[1].Beta),
                        tree_data->TPDate,          
                        (tree_data->NbTP+1)))
        {
            goto RETURN;
        }  /* if */
    }

    /*  CET: Called to enhance Vladimir's approximation */
    /*   (output to CET.prn)  */
    for (i = 0; i < NbIRatesOn; i++)
    {  
        /* if we are running an FX tree with CET, we need to run CET with the    */
        /* (reduced) NbIRSigma                                                   */
        /* Given that we systematically initialise IRSigma to SigmaDBL at        */
        /* the start of Build Tree, we need to trick it by temporarily replacing */
        /* SigmaDBL by IRSigma                                                   */

        temp                  = tree_data->NbSigmaDBL ;
        tree_data->NbSigmaDBL = tree_data->NbIRSigmaMax;

        if (Hyb3_Cet_Main (CetOutputFlag,
                      i,
                      t_curve,
                      mktvol_data,
                      tree_data) == FAILURE)
        {
            goto RETURN;
        } 

        /* restore the value of SigmaDBL*/
        tree_data->NbSigmaDBL = temp;

    }


    for (i = 0; i < NbIRatesOn; i++)
    {
        /* Bootstrap CET'd swaption volatilities to get spot */
        /* volatility curve. Use diffused t_curve !          */
        if (Hyb3_SpotVol(mktvol_data[i].Aweight,
                    mktvol_data[i].BaseDate,
                    mktvol_data[i].NbVol,
                    mktvol_data[i].VolDate,
                    mktvol_data[i].Vol,
                    mktvol_data[i].VolUsed,
                    mktvol_data[i].Freq,
                    mktvol_data[i].DCC,
                    mktvol_data[i].SwapSt,
                    mktvol_data[i].SwapMat,
                    mktvol_data[i].QLeft,
                    mktvol_data[i].QRight,
                    mktvol_data[i].FwdShift,
                    mktvol_data[i].NbFactor,
                    mktvol_data[i].Alpha,
                    mktvol_data[i].Beta,
                    mktvol_data[i].Rho,
                    mktvol_data[i].SkipFlag,
                    mktvol_data[i].CalibFlag,
                    &t_curve[i][CvDiff[i]]) != SUCCESS)
        {
            goto RETURN;
        }

        /* Interpolate spot volatility curves */
        if (Hyb3_Interp_SpotVol (tree_data->IrAweight[i],
                            mktvol_data[i].BaseDate,
                            mktvol_data[i].NbVol,
                            mktvol_data[i].VolDate,
                            mktvol_data[i].Aweight,
                            mktvol_data[i].NbFactor,
                            mktvol_data[i].Beta,
                            mktvol_data[i].CalibFlag,
                            tree_data->FwdRate[i][CvDiff[i]],
                            tree_data->Length,
                            tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }
        for (j = -1; j <= tree_data->NbTP; j++)
        {
            /*   Spot vol is stored in "Aweight" form, i.e spotvol=sigma*alpha(0)    */
                                                         
            tree_data->SpotVol[i][j]=tree_data->IrAweight[i][0][j]; 
        }

        /*NO NEED TO EXTRAPOLATE CORRELATION AGAIN FOR SEVERAL FACTOR IR: INDEPENDANT OF SCALING OF AWEIGHT*/
    }


    /* Weights */
    if (Hyb3_Tree_Weights(tree_data->Aweight,
                          tree_data->Rho,
                          tree_data->TreeType,
                          tree_data->NbTP,
                          tree_data->IrAweight,
                          tree_data->SpotFxVol,
                          tree_data->SpotEqVol) == FAILURE)
    {
        goto RETURN;
    }

    for (i=0; i < NbIRatesOn; i++)
    {

        /* Provisional tree limits just for purposes of drift calibration    */
        if (Hyb3_Tree_Limits_IR2(&(tree_data->Width[0]), /* These widths are */
                            &(tree_data->HalfWidth[0]),  /* recal'd later    */
                            tree_data->Top1,             /* (O) 1D Upper lim */
                            tree_data->Bottom1,          /* (O) 1D Lower lim */
                            tree_data->Top2,             /* (O) 2D Upper lim */
                            tree_data->Bottom2,          /* (O) 2D Lower lim */
                            tree_data->IrAweight[i],
                            mktvol_data[i].NbFactor,
                            tree_data->Length,
                            tree_data->LengthJ,
                            tree_data->NbTP,
                            (tree_data->NbIRSigmaMax), 
                            mktvol_data[i].Beta) == FAILURE)
        {
            goto RETURN;
        }

        if (Hyb3_Tree_Limits_IR2(&(tree_data->Width[0]), /* These widths are */
                            &(tree_data->HalfWidth[0]),  /* recal'd later    */
                            tree_data->OutTop1,          /* (O) 1D Upper lim */
                            tree_data->OutBottom1,       /* (O) 1D Lower lim */
                            tree_data->OutTop2,          /* (O) 2D Upper lim */
                            tree_data->OutBottom2,       /* (O) 2D Lower lim */
                            tree_data->IrAweight[i],
                            mktvol_data[i].NbFactor,
                            tree_data->Length,
                            tree_data->LengthJ,
                            tree_data->NbTP,
                            (tree_data->NbIRSigmaMax),
                            mktvol_data[i].Beta) == FAILURE)
        {
            goto RETURN;
        }

        /* Get one factor IR grid drift */
        if (mktvol_data[i].NbFactor == 1)
        {
            /* Get one factor IR grid drift */
            if (Hyb3_Find_Drift_1D(tree_data->IrZCenter[i],
                              tree_data->ZeroCoupon[i][CvDiff[i]],
                              tree_data->FwdRate[i][CvDiff[i]],
                              mktvol_data[i].QLeft,
                              mktvol_data[i].QRight,
                              mktvol_data[i].FwdShift,
                              mktvol_data[i].Beta[0],
                              tree_data->SpotVol[i],
                              tree_data->Top1,
                              tree_data->Bottom1,
                              tree_data->OutTop1,   
                              tree_data->OutBottom1,
                              tree_data->NbTP,
                              tree_data->Length,
                              tree_data->LengthJ,
                              tree_data) == FAILURE)
           {
                goto RETURN;
           }
        }
        else if (mktvol_data[i].NbFactor==2)
        {
            if (Hyb3_Find_Drift_2D(tree_data->IrZCenter[i],
                              tree_data->ZeroCoupon[i][CvDiff[i]],
                              tree_data->FwdRate[i][CvDiff[i]],
                              mktvol_data[i].QLeft,
                              mktvol_data[i].QRight,
                              mktvol_data[i].FwdShift,
                              mktvol_data[i].Beta,
                              tree_data->IrAweight[i],
                              tree_data->Top1,
                              tree_data->Bottom1,
                              tree_data->OutTop1,   
                              tree_data->OutBottom1,
                              tree_data->Top2,
                              tree_data->Bottom2,
                              tree_data->OutTop2,   
                              tree_data->OutBottom2,
                              tree_data->NbTP,
                              tree_data->Length,
                              tree_data->LengthJ,
                              tree_data) == FAILURE)
           {
                goto RETURN;
           }
        }

        /* Free Provisional limit if not CET*/
        if ((tree_data->TreeType != TTYPE_1IR)&&
        (tree_data->TreeType != TTYPE_1IR2F))
        {
        if (Hyb3_Free_TreeLimits_IR2(tree_data->Top1,             /* (O) 1D Upper lim */
                            tree_data->Bottom1,          /* (O) 1D Lower lim */
                            tree_data->Top2,             /* (O) 2D Upper lim */
                            tree_data->Bottom2,          /* (O) 2D Lower lim */
                            mktvol_data[i].NbFactor,
                            tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }
        if (Hyb3_Free_TreeLimits_IR2(tree_data->OutTop1,          /* (O) 1D Upper lim */
                                     tree_data->OutBottom1,       /* (O) 1D Lower lim */
                                     tree_data->OutTop2,          /* (O) 2D Upper lim */
                                     tree_data->OutBottom2,       /* (O) 2D Lower lim */
                                     mktvol_data[i].NbFactor,
                                     tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }
        }
    
    }


    /* Calculate the cups adjustment */
    if (NbIRatesOn == 2)
    {
        /* CAREFUL: This cups adjustment does not take into account FX smile */
        /* it is still based on the Log-Normal FX + Normal IR approximation  */
           if (Hyb3_CUP_Drift(tree_data->DriftCUPS,
                    tree_data->IrAweight[0],
                    mktvol_data[0].NbFactor,
                    tree_data->SpotFxVol,
                    &(tree_data->Rho[RhoFxIrFacIndex[0]]),/* FX vs IRFor Fac */
                    tree_data->Length,
                    tree_data->NbTP) == FAILURE)
        {
            goto RETURN; 
        }
    }

    /* Add cups drift to foreign equity */
    if ((EquityOnIr1 == TRUE) && (NbIRatesOn == 2))
    {
        if (Hyb3_CUPSEqMidNode(tree_data->EqMidNode,
                          tree_data->SpotEqVol,
                          tree_data->SpotFxVol,
                          tree_data->Rho[5], 
                          tree_data->Length,
                          tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }
    }

    /* Final tree limits, if not running in CET mode (i.e. one */
    /* interest rate only).                   */
    if ((tree_data->TreeType != TTYPE_1IR)&&
        (tree_data->TreeType != TTYPE_1IR2F))
    {

        /* Tree_limits for inner ellipse  */
        if (Hyb3_Tree_Limits (tree_data->Width, 
                         tree_data->HalfWidth,
                         tree_data->xWidth, 
                         tree_data->xHalfWidth,
                         tree_data->Top1,
                         tree_data->Bottom1,
                         tree_data->Top2,
                         tree_data->Bottom2,
                         tree_data->Top3,
                         tree_data->Bottom3,
                         tree_data->Top4,
                         tree_data->Bottom4,
                         tree_data->NbSigmaDBL,
                         tree_data->NbIRSigmaMax,
                         tree_data->NbTP,
                         tree_data->xT,
                         tree_data->TreeType,
                         tree_data->FwdRate,
                         tree_data->CvDisc,
                         tree_data->IrAweight,               /* IR spotvols */
                         tree_data->SpotFxVol,  /* Spotvol  of */
                         tree_data->SpotEqVol,  /* third asset */
                         tree_data->Rho,
                         mktvol_data[0].Beta,
                         (NbIRatesOn == 1) ?
                           NULL :
                           mktvol_data[1].Beta, /* Unused if TTYPE_EQ1IR */
                         tree_data->Length,
                         tree_data->LengthJ) == FAILURE)
        {
            goto RETURN;
        }

        if ((tree_data->TreeType != TTYPE_EQDFX2IR) &&
            (tree_data->TreeType != TTYPE_EQFFX2IR))
        {
            /* Tree_limits for outer ellipse */
            if (Hyb3_Tree_Limits (tree_data->Width, 
                         tree_data->HalfWidth,
                         tree_data->xWidth, 
                         tree_data->xHalfWidth,
                         tree_data->OutTop1,
                         tree_data->OutBottom1,
                         tree_data->OutTop2,
                         tree_data->OutBottom2,
                         tree_data->OutTop3,
                         tree_data->OutBottom3,
                         tree_data->OutTop4,
                         tree_data->OutBottom4,
                         tree_data->NbSigmaDBL,
                         tree_data->NbIRSigmaMax,
                         tree_data->NbTP,
                         tree_data->xT,
                         tree_data->TreeType,
                         tree_data->FwdRate,
                         tree_data->CvDisc,
                         tree_data->IrAweight,
                         tree_data->SpotFxVol,
                         tree_data->SpotEqVol,
                         tree_data->Rho,
                         mktvol_data[0].Beta,
                         (NbIRatesOn == 1) ?
                           NULL :
                           mktvol_data[1].Beta, /* Unused if TTYPE_EQ1IR */
                         tree_data->Length,
                         tree_data->LengthJ) == FAILURE)
            {               
                goto RETURN;
            }
       
        
            /* Tree_limits for inner ellipse  */
            if (Hyb3_Tree_Limits (tree_data->InnerTreeWidth, 
                             tree_data->InnerTreeHalfWidth,
                             tree_data->xWidth, 
                             tree_data->xHalfWidth,
                             tree_data->InnerTreeTop1,
                             tree_data->InnerTreeBottom1,
                             tree_data->InnerTreeTop2,
                             tree_data->InnerTreeBottom2,
                             tree_data->InnerTreeTop3,
                             tree_data->InnerTreeBottom3,
                             tree_data->InnerTreeTop4,
                             tree_data->InnerTreeBottom4,
                             tree_data->NbIRSigmaMax, 
                             tree_data->NbIRSigmaMax - INNER_TREE_SIZE,
                             tree_data->NbTP,
                             tree_data->xT,
                             tree_data->TreeType,
                             tree_data->FwdRate,
                             tree_data->CvDisc,
                             tree_data->IrAweight,             /* IR spotvols */
                             tree_data->SpotFxVol,  /* Spotvol  of */
                             tree_data->SpotEqVol,  /* third asset */
                             tree_data->Rho,
                             mktvol_data[0].Beta,
                             (NbIRatesOn == 1) ?
                               NULL :
                               mktvol_data[1].Beta, /* Unused if TTYPE_EQ1IR */
                             tree_data->Length,
                             tree_data->LengthJ) == FAILURE)
            {
                goto RETURN;
            }
        }

        /* Last Step: correct the CUPS drift  */
        if ( (tree_data-> TreeType == TTYPE_FX2IR) && (FxMomentMatchOn)  )
        {

            /* calculate the adjustment factor first */
            if ( Hyb3_CUPSAdjust_FXSmile( tree_data, mktvol_data, fx_data ) == FAILURE)
            {
               goto RETURN;
            }
        }

    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Build_Tree */


/***************************************************************************************/
/*

 ***************************************************************************************/
int Hyb3_Find_IR_StdDev(
                   double   *IRSigmaMax,    /**< (O) allowable Nb of std devs on IR dims    */
                   double   NbFxSigma,      /**< (I) input FX std devs                      */
                   double   BetaD,          /**< (I) dom mr                                 */
                   double   BetaF,          /**< (I) for mr                                 */
                   double   *SpotVol[2],    /**< (I) IR spot vols                           */
                   double   Rho,            /**< (I) dom/for ir correl                      */
                   int      NbTP,           /**< (I) Nb TP in tree                          */
                   double   *FwdRate[2][3], /**< (I) Fwd Rates                              */
                   double   *CompositeVol,  /**< (I) Fx Vol                                 */
                   double   *Length,        /**< (I)                                        */
                   long     *TPDate,        /**< (I)                                        */
                   double   QLeftD,         /**< (I)                                        */
                   double   QRightD,        /**< (I)                                        */
                   double   FwdShiftD,      /**< (I)                                        */
                   double   QLeftF,         /**< (I)                                        */                                      
                   double   QRightF,        /**< (I)                                        */
                   double   FwdShiftF)      /**< (I)                                        */

{
    int     status = FAILURE;
    long    StartDate = 0L;
    long    EndDate   = 0L;
    int     NbEvDates;
    long    *EvDates = NULL;

    double  NextLength;
    double  Var1L = 0.0;
    double  Var2L = 0.0;    
    double  CovL  = 0.0;
    double  Vol1,Vol2,CRho, Cov,Delta_t;

    double  *Var1   = NULL,
            *Var2   = NULL,
            *Cov12    = NULL,
            *TotalT = NULL;

    int     i,t;

    double  TimeToNextEvDate = 0.0;
    int     ThisOffset,NextOffset;
    
    double  T ;
   
    double  ZCenterD,ZCenterF;
    int     CurrOffs;
    double  *MaxSigma = NULL;
    double  Vol2TimesRho,Vol2TimesSqrtRho;
   
    double  FXstddev,loc, MaxAllowed, MinAllowed;
    double  FwdFor,FwdDom;
    double  QM;
    double  St;

    
    double  X;
    int     jEps1, kX;
    
    double  NormalF,NormalD;
    double  Eps1;
    double  RateF;
    double  MinDrift,MaxDrift;

    double  MinStdev;
    double  DelX,DelEps;
    int     NbStepX, NbStepEps;
    int     HasExceeded;

    char    message[MAXBUFF];

    int     *EvDatesOffset = NULL;
    double  *TimeToEvDates = NULL;

    double  Ratio;

    /* checks */
    


    NbEvDates       = 0L;

    if(DrDateFwdAny(TPDate[0],
                    1,   
                    'A', /*  1 year  */
                    'F', /* Forward*/
                    &StartDate) == FAILURE) goto RETURN;

    EndDate = TPDate[NbTP];

    if (StartDate >= EndDate)/* time line is < 1 year then no need to do any search*/
    {   
        MinStdev   = NbFxSigma;
       *IRSigmaMax = IR_MINSIZE + (MinStdev - IR_MINSIZE)/3.0;

        if((*IRSigmaMax < IR_MINSIZE) || (*IRSigmaMax < INNER_TREE_SIZE))
        {   
            sprintf(message,"IR Stdev too small: %lf\n",*IRSigmaMax);
            DR_Error(message);
            return(FAILURE);
        }
        status = SUCCESS;
        return(status);
    }
    else
    {
        CurrOffs = -1L;  
        if(DateListFromFreq(StartDate,
                            EndDate,
                            'Q',    /* Frequency*/
                            'F',    /* Stub     */
                            &NbEvDates,
                            &EvDates) == FAILURE) goto RETURN;

        

        EvDatesOffset = (int*)    DR_Array(INT   ,0,NbEvDates-1L);
        TimeToEvDates = (double*) DR_Array(DOUBLE,0,NbEvDates-1L);
        
        if((EvDatesOffset == NULL) ||
           (TimeToEvDates == NULL))
            goto RETURN;

        for(i = 0; i<NbEvDates;i++)
        {
            CurrOffs = GetDLOffset(NbTP+1,
                                   TPDate,
                                   EvDates[i],
                                   CbkLOWER);
            
            if(CurrOffs < 0L) goto RETURN;
            EvDatesOffset[i]    = CurrOffs;
            TimeToEvDates[i]    = Daysact(TPDate[0], EvDates[i])/365.0;
        }/* for i */
    }

    
    if(NbEvDates < 2L) goto RETURN;
    MaxSigma = (double*) DR_Array(DOUBLE,0,NbEvDates -2L);
    if(MaxSigma == NULL) goto RETURN;


    Var1    = (double*) DR_Array(DOUBLE,0,NbTP);
    Var2    = (double*) DR_Array(DOUBLE,0,NbTP);
    Cov12   = (double*) DR_Array(DOUBLE,0,NbTP);
    TotalT  = (double*) DR_Array(DOUBLE,0,NbTP);

    Var1L = Var2L  = 0.0;
    CovL           = 0.0;
    T              = 0.0;


    for(t = 0; t<= NbTP;t++)
    {
        Var1[t]   = Var1L;
        Var2[t]   = Var2L;
        Cov12[t]  = CovL;
        TotalT[t] = T;

        NextLength = Length[t];
        T += NextLength;


        Var1L *= exp (- 2. * BetaF * NextLength);
        Var1L += SQUARE(SpotVol[0][t])* NextLength * Hyb3_ExpDecay(2.*BetaF,NextLength);
        

        Var2L *= exp (- 2. * BetaD * NextLength);
        Var2L += SQUARE(SpotVol[1][t]) * NextLength * Hyb3_ExpDecay(2.*BetaD,NextLength);
        

        CovL *= exp (- (BetaD + BetaF) * NextLength);
        CovL += Rho * SpotVol[0][t] * SpotVol[1][t] * NextLength * 
                Hyb3_ExpDecay(BetaD + BetaF, NextLength);
    }


    /* StepSize for X */
    DelX  = 0.05;
    NbStepX = (int) ((NbFxSigma)/DelX);

    /* Step Size for Eps1 */
    DelEps = 0.05;

    for(i = 0; i <= NbEvDates - 2L;i++)
    {
        ThisOffset       = EvDatesOffset[i];
        NextOffset       = EvDatesOffset[i+1];
        TimeToNextEvDate = Daysact(EvDates[i],EvDates[i+1])/365.;

        Delta_t           = Daysact(TPDate[ThisOffset],EvDates[i])/365.;

        /* Adjust variances if necessary*/
        if(Delta_t > TINY)
        {
            Vol1 = Var1[ThisOffset] * exp (- 2. * BetaF * Delta_t);
            Vol1 += SQUARE(SpotVol[0][ThisOffset])* Delta_t * Hyb3_ExpDecay(2.*BetaF,Delta_t);

            Vol2 = Var2[ThisOffset] * exp (- 2. * BetaD * Delta_t);
            Vol2 += SQUARE(SpotVol[1][ThisOffset])* Delta_t * Hyb3_ExpDecay(2.*BetaD,Delta_t);

            Cov = Cov12[ThisOffset] * exp (- (BetaD + BetaF) * Delta_t);
            Cov += Rho * SpotVol[0][ThisOffset] * SpotVol[1][ThisOffset] * Delta_t * 
                   Hyb3_ExpDecay(BetaD + BetaF, Delta_t);
        }
        else
        {
            Vol1 = Var1[ThisOffset];
            Vol2 = Var2[ThisOffset];
            Cov  = Cov12[ThisOffset];
        }


        Vol1 = sqrt(Vol1);
        Vol2 = sqrt(Vol2);
        CRho = Cov/ Vol1 / Vol2;
                         
        if ((CRho - 1 > TINY) || (CRho + 1 < TINY))
        {
            DR_Error("CRho1 is not a correl number ( >1 or <-1)\n");
            goto RETURN;
        }

        Vol2TimesRho     = Vol2 * CRho;
        Vol2TimesSqrtRho = Vol2 * sqrt(1-CRho * CRho);

        /* Composite Vol is strictly speaking not the right compvol to use*/
        FXstddev     = CompositeVol[NextOffset] * sqrt(TimeToEvDates[i+1]);
        loc          = STDDEV_FACTOR * NbFxSigma * FXstddev;
        
        /* approximate the Fwd rates */
        FwdFor     = FwdRate[0][1][ThisOffset] * TimeToNextEvDate/Length[ThisOffset];
        FwdDom     = FwdRate[1][1][ThisOffset] * TimeToNextEvDate/Length[ThisOffset];

        /* Foreign Calibration constant */
        QM  = (QLeftF + QRightF) / 2;
        St  = Vol1 * (1 + FwdShiftF) / (1 + QM * FwdShiftF);

        if (ConvexityC_BS2Q (St, QLeftF, QRightF, FwdShiftF, &ZCenterF) == FAILURE)
        {
            DR_Error ("Could not find calib const for Foreign Rate !");
            goto RETURN;
        }
        /* Domestic Calibration constant */
        QM  = (QLeftD + QRightD) / 2;
        St  = Vol2 * (1 + FwdShiftD) / (1 + QM * FwdShiftD);

        if (ConvexityC_BS2Q (St, QLeftD, QRightD, FwdShiftD, &ZCenterD) == FAILURE)
        {
            DR_Error ("Could not find calib const for Dom Rate !");
            goto RETURN;
        }


        if(Mapping_2Q(&Ratio,
                      QLeftD,
                      QRightD,
                      FwdShiftD,
                      ZCenterD) ==  FAILURE) goto RETURN;

        Ratio *= FwdDom/(1+FwdShiftD);
        Ratio += 1.0;

        if(Mapping_2Q(&RateF,
                      QLeftF,
                      QRightF,
                      FwdShiftF,
                      ZCenterF) ==  FAILURE) goto RETURN;

        RateF *= FwdFor/(1+FwdShiftF);
        RateF += 1.0;

        if(fabs(RateF) < TINY) goto RETURN;

        Ratio /= RateF;
        MaxAllowed   = Ratio * exp( loc );
        MinAllowed   = Ratio * exp(-loc );


        X = NbFxSigma;
        for(kX = 0 ; kX < NbStepX ; kX++)
        {   
            /* this should not happen but check just to be on the safe side */
            if(X < 0) goto RETURN;

            NbStepEps = (int)((2*X)/DelEps); 
            Eps1     = X;

            HasExceeded = FALSE;
            for(jEps1 = 0; jEps1 <= NbStepEps; jEps1++)
            {   
                NormalD = ZCenterD + Vol2TimesRho* Eps1 
                          + Vol2TimesSqrtRho * sqrt(MAX(X*X - Eps1 * Eps1,0.0));

                if(Mapping_2Q(&MaxDrift,
                              QLeftD,
                              QRightD,
                              FwdShiftD,
                              NormalD) ==  FAILURE) goto RETURN;
                MaxDrift *= FwdDom/(1+FwdShiftD);
                MaxDrift += 1.0;

                /* We should really be adding a CUPS adj , we are not*/
                NormalF = ZCenterF + Vol1 * Eps1;
                if(Mapping_2Q(&RateF,
                              QLeftF,
                              QRightF,
                              FwdShiftF,
                              NormalF) ==  FAILURE) goto RETURN;

                RateF *= FwdFor/(1+FwdShiftF);
                RateF += 1.0;

                if(fabs(RateF) < TINY) goto RETURN;
                MaxDrift /= RateF;

                NormalD = ZCenterD + Vol2TimesRho* Eps1 - 
                          Vol2TimesSqrtRho * sqrt(MAX(X*X - Eps1 * Eps1,0.0));

                if(Mapping_2Q(&MinDrift,
                              QLeftD,
                              QRightD,
                              FwdShiftD,
                              NormalD) ==  FAILURE) goto RETURN;

                MinDrift *= FwdDom/(1+FwdShiftD);
                MinDrift += 1.0;
                MinDrift /= RateF;

                HasExceeded = ((MaxDrift > MaxAllowed) || (MinDrift < MinAllowed));

                if(HasExceeded) break;
                
                Eps1 -= DelEps;
            }/* for jEps */

            if(!HasExceeded) break;
            X -= DelX;

        }/* for kX */

        MaxSigma[i] = X;
    }/* for i */

    MinStdev = MaxSigma[0];

    for(i=1;i<= NbEvDates-2L;i++)
    {
        MinStdev = MIN(MaxSigma[i], MinStdev);
    }


    *IRSigmaMax = IR_MINSIZE + (MinStdev - IR_MINSIZE)/3.0;

    if((*IRSigmaMax < IR_MINSIZE) || (*IRSigmaMax < INNER_TREE_SIZE))
    {   
        sprintf(message,"IR Stdev too small: %lf\n",*IRSigmaMax);
        DR_Error(message);
        goto RETURN;
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("Find IR_Stddev failed\n");
    }
    Free_DR_Array(Var1,DOUBLE,0,NbTP);
    Free_DR_Array(Var2,DOUBLE,0,NbTP);
    Free_DR_Array(TotalT,DOUBLE,0,NbTP);
    Free_DR_Array(Cov12,DOUBLE,0,NbTP);
    Free_DR_Array(MaxSigma,DOUBLE,0,NbEvDates -2L);
    Free_DR_Array(EvDatesOffset, INT,0,NbEvDates-1L);
    Free_DR_Array(TimeToEvDates, DOUBLE,0,NbEvDates-1L);
    Free_DR_Array(EvDates, LONG,0,NbEvDates-1L);

    return(status);

}

/*****  Hyb3_Tree_Limits  ********************************************************/
/**
        Calculates tree limits: we calculate  the total vol of the forwards
        (up to three dims) and cut the tree according to a parabola in each
        dimension. Also allocate on the go memory for the tree limits.
 
        Assumes that the TreeWeights have been computed prior to calling.
 
*/
int     Hyb3_Tree_Limits (
                  int       *Width,          /**< (O) Min nb of nodes each dim */
                  int       *HalfWidth,      /**< (O) Max nb of nodes each dim */
                  int       *xWidth,         /**< (O) Min nb of nodes each dim */
                  int       *xHalfWidth,     /**< (O) Max nb of nodes each dim */
                  int       *Top1,           /**< (O) Upper limits of 1D tree  */
                  int       *Bottom1,        /**< (O) Lower limits             */
                  int       **Top2,          /**< (O) Upper limits of 2D tree  */
                  int       **Bottom2,       /**< (O) Lower limits             */
                  int       ***Top3,         /**< (O) Upper limits of 3D tree  */
                  int       ***Bottom3,      /**< (O) Lower limits             */
                  int       ****Top4,        /**< (O) Upper limits of 4D tree  */
                  int       ****Bottom4,     /**< (O) Lower limits             */
                  double    NbSigmaDBL,      /**< (I) Nb of std dev used to cut*/
                  double    NbIRSigmaMax,
                  int       NbTP,            /**< (I) Total number of nodes    */
                  int       xT,
                  int       TreeType,        /**< (I) Nb and type of factors   */
                  double   *FwdRate[2][3],   /**< (I) Use COF fwd rate at nodes*/
                  int      *CvDiff,
                  double   *IrAweight[2][3], /**< (I) IR spotvols in timeline  */
                  double   *SpotVolFx,       /**< (I) Spotvols of non-IR asset */
                  double   *SpotVolEq,       /**< (I) Spotvols of non-IR asset */
                  double   *Rho[6],          /**< (I) Correlations             */
                  double   *Beta1,           /**< (I) Mean reversion of IR's   */
                  double   *Beta2,           /**< (I) Mean reversion of IR's   */
                  double   *Length,          /**< (I) Length of time steps     */
                  double   *LengthJ          /**< (I) Time step for jump size  */
                     )
{
        double
                *SpotVolAsset1 = NULL,
                *SpotVolAsset2 = NULL,
                *SpotVolAsset3 = NULL,
                *SpotVolAsset4 = NULL,
                *SpotVol[2]    = {NULL, NULL},
                *Rho1 = NULL,                 /* IR1 vs IR2    (for convenience)  */
                *Rho2 = NULL,                 /* IR1 vs Asset3 (for convenience)  */
                *Rho3 = NULL,                 /* IR2 vs Asset3 (for convenience)  */
                *Rho4 = NULL,
                *Rho5 = NULL,
                *Rho6 = NULL,
                 LBeta[3],                    /* Local mean reversion             */
                *AIr1 = NULL,                 /* Bond vol integral for foreign    */
                *AIr2 = NULL,                 /* Bond vol integral for domestic   */
                CurrLengthJ,                  /* Length of curr t step for jump   */
                NextLengthJ,                  /* Lenght of next t step for jump   */
                NextLength,                   /* Lenght of next time step         */
                Var1,  Var2,  Var3,  Var4,    /* Variance of factors              */
                Cov1,  Cov2,  Cov3,  Cov4,    /* Covariance of factors            */
                Vol1,  Vol2,  Vol3,  Vol4,    /* Volatility of factors            */
                CRho1, CRho2, CRho3, CRho4,   /* Composite correlations           */
                SpotCov1,
                SpotCov2,                     /* Cov term for spot processes      */
                SpotCov3,                     /* Cov term for spot processes      */
                SpotCov4,                     /* Cov term for spot processes      */
                SpotVar2,
                SpotVar3,                     /* Variance term for spot process   */
                SpotVar4,                     /* Variance term for spot process   */
                I[16],
                Delta,                    /*                                  */
            
                det,
                dt,
                T,                            /* Total time for integration       */
                du;                           /* Scaled delta time                */
            

        int
                Centre, h, h2, h3, h4,
                t, i, j, k, u, v, z,       /* Counters                         */

                ThirdDimOn  = FALSE,          /* Flag to indicate that tree is 3-D*/
                FourthDimOn = FALSE,          /* Flag to indicate that tree is 4-D*/
                TotDim      = 0,              /* Total nb of tree dims (2 or 3)   */
                EquityOn    = FALSE,          /* True if 3rd dim is eq (not FX)   */
                EquityOnIr1 = FALSE,          /* True if eq denomination is FOR   */
                ThirdDimIR  = FALSE,          /* True if 3rd dim is IR            */

                status = FAILURE;             /* Error status = FAILURE initially */
        
        double  DomA, ForA;

        double  xx[4];
        double  y[4];

        double  J[4][4];

        double  r[4][4];
        double  R[4][4];
   
        double  C[4][4];

        double  n[4][4]; 

        double  N[4][4];
        double  G[4][4];

        double  a[4][4];

        double  vol[4];
        double  Vol[4];

        double  C20[MAXNBDATE];
        double  C21[MAXNBDATE];
        double  C22[MAXNBDATE];
        double  C30[MAXNBDATE];
        double  C31[MAXNBDATE];
        double  C32[MAXNBDATE];
        double  C33[MAXNBDATE];


        for (u=0; u<4; ++u)
        {
            vol[u] = 0;
            Vol[u] = 0;

            xx[u] = 0;
            y[u]  = 0;

            for (v=0; v<4; ++v)
            {
                J[u][v] = 0;

                r[u][v] = 0;
                R[u][v] = 0;
            
                C[u][v] = 0;

                n[u][v] = 0; 

                N[u][v] = 0;
                G[u][v] = 0;

                a[u][v] = 0;
            }
        }

        /* Initialise accumulators used for variances etc.     */
        Var1  =  Var2 =  Var3 =  Var4 = 0.;  /* Total variances        */
        Vol1  =  Vol2 =  Vol3 =  Vol4 = 0.;  /* Total volatilities     */
        Cov1  =  Cov2 =  Cov3 =  Cov4 = 0.;  /* Covariances            */
        CRho1 = CRho2 = CRho3 = CRho4 = 0.;  /* Composite correlations */
        SpotCov1 = 0.0;
        SpotCov2 = 0.0; /* 'Basic' term for IR1 vs Asset3 cov  */
        SpotCov3 = 0.0; /* 'Basic' term for IR2 vs Asset3 cov  */
        SpotCov4 = 0.0; /* 'Basic' term for IR2 vs Asset3 cov  */
        SpotVar2 = 0.0;
        SpotVar3 = 0.0; /* 'Basic' term for variance of Asset3 */
        SpotVar4 = 0.0; /* 'Basic' term for variance of Asset3 */

        /*Default case */
        SpotVol[0] = IrAweight[0][0];
        SpotVol[1] = IrAweight[1][0];
        LBeta[0] = Beta1[0];                    /* Foreign mean reversion    */
        LBeta[1] = (Beta2==NULL) ? 0: Beta2[0]; /* Domestic mean reversion   */
        LBeta[2] = 0;                           /* No mean rev. on third dim */

        SpotVolAsset1 = SpotVol[0];

        /*  TreeType will tell us if the tree is being run in 2-D mode (i.e. */
        /*  2 IR's) or 3-D (i.e. plus eq or FX) and therefore correlations   */
        if (TreeType == TTYPE_EQ1IR)
        {
            SpotVolAsset2 = SpotVolEq;
            ThirdDimOn = FALSE;
            FourthDimOn = FALSE;
            TotDim     = 2;
            Rho1       = Rho[3];
        }
        else if (TreeType == TTYPE_2IR)
        {
            SpotVolAsset2 = SpotVol[1];
            ThirdDimOn = FALSE;
            FourthDimOn = FALSE;
            TotDim     = 2;
            Rho1       = Rho[0];
        }
        else if (TreeType == TTYPE_2IR2F1D)
        {
            /* Due to Aweight extrapolation, effective correlation
             * may be different from input factor correlation
             */
            SpotVolAsset2 = (double*)DR_Array(DOUBLE, -1, NbTP);
            Rho1          = (double*)DR_Array(DOUBLE, -1, NbTP);
            for (t = -1; t <= NbTP; t++)
            {
                SpotVolAsset2[t] = sqrt(SQUARE(IrAweight[0][1][t])+
                                           SQUARE(IrAweight[0][2][t]));
                Rho1[t]          = IrAweight[0][1][t]/SpotVolAsset2[t]; 
            }
            SpotVolAsset3 = IrAweight[1][0];
            ThirdDimOn = TRUE;
            ThirdDimIR = TRUE;
            TotDim     =3;
            EquityOn = FALSE;
            Rho2     = Rho[0]; /*IrD vs 1st factor IrF */
            Rho3     = Rho[1]; /*IrD vs 2nd factor IrF */
            LBeta[1] = Beta1[1];
            LBeta[2] = Beta2[0];
        }
        else if (TreeType == TTYPE_FX2IR)
        {
            SpotVolAsset2 = SpotVol[1];
            SpotVolAsset3 = SpotVolFx;
            ThirdDimOn = TRUE;
            FourthDimOn = FALSE;
            TotDim     = 3;
            Rho1       = Rho[0];
            EquityOn = FALSE;
            Rho2 = Rho[1]; /* IrF vs FX */
            Rho3 = Rho[2]; /* IrD vs FX */
        }
        else if (TreeType == TTYPE_EQD2IR)
        {
            SpotVolAsset2 = SpotVol[1];
            SpotVolAsset3 = SpotVolEq;
            ThirdDimOn = TRUE;
            FourthDimOn = FALSE;
            TotDim     = 3;
            Rho1       = Rho[0];
            EquityOn = TRUE;
            Rho2 = Rho[3]; /* IrF vs Eq */
            Rho3 = Rho[4]; /* IrD vs Eq */
            EquityOnIr1 = FALSE;
        }
        else if (TreeType == TTYPE_EQF2IR)
        {
            SpotVolAsset2 = SpotVol[1];
            SpotVolAsset3 = SpotVolEq;
            ThirdDimOn = TRUE;
            FourthDimOn = FALSE;
            TotDim     = 3;
            Rho1       = Rho[0];
            EquityOn = TRUE;
            Rho2 = Rho[3]; /* IrF vs Eq */
            Rho3 = Rho[4]; /* IrD vs Eq */
            EquityOnIr1 = TRUE;
        }
        else if (TreeType == TTYPE_EQDFX2IR)
        {
            SpotVolAsset2 = SpotVol[1];
            SpotVolAsset3 = SpotVolFx;
            SpotVolAsset4 = SpotVolEq;
            TotDim        = 4;
            ThirdDimOn    = TRUE;
            FourthDimOn   = TRUE;

            EquityOn      = TRUE;
            EquityOnIr1   = FALSE;

            Rho1 = Rho[0];
            Rho2 = Rho[1];
            Rho3 = Rho[2];
            Rho4 = Rho[3];
            Rho5 = Rho[4];
            Rho6 = Rho[5];
        }
        else if (TreeType == TTYPE_EQFFX2IR)
        {
            SpotVolAsset2 = SpotVol[1];
            SpotVolAsset3 = SpotVolFx;
            SpotVolAsset4 = SpotVolEq;
            TotDim        = 4;
            ThirdDimOn    = TRUE;
            FourthDimOn   = TRUE;
        
            EquityOn      = TRUE;
            EquityOnIr1   = TRUE;

            Rho1 = Rho[0];
            Rho2 = Rho[1];
            Rho3 = Rho[2];
            Rho4 = Rho[3];
            Rho5 = Rho[4];
            Rho6 = Rho[5];
        }
        
        if (FourthDimOn) /* do fx & eq dimensions in the new way */
        {
            double  value;

            double  Time    [MAXNBDATE];
            double  One     [MAXNBDATE];

            double  FwIrFor [MAXNBDATE];
            double  FwIrDom [MAXNBDATE];

            double  MrIrFor = Beta1[0];
            double  MrIrDom = Beta2[0];

            double *SvIrFor = SpotVolAsset1;
            double *SvIrDom = SpotVolAsset2;

            double *SvFx    = SpotVolAsset3;
            double *SvEq    = SpotVolAsset4;

            double *CrIrForIrDom = Rho1;
            double *CrIrForFx    = Rho2;
            double *CrIrDomFx    = Rho3;
            double *CrIrForEq    = Rho4;
            double *CrIrDomEq    = Rho5;
            double *CrFxEq       = Rho6;

            for (t=0; t<MAXNBDATE; ++t)
            {
                One [t] = 1;
            }

            value = 0;

            for (t=0; t<=NbTP; ++t)
            {
                Time[t] = value;

                dt = Length[t];

                FwIrFor [t] = FwdRate[0][CvDiff[0]][t] / dt;
                FwIrDom [t] = FwdRate[1][CvDiff[1]][t] / dt;

                value += dt;
            }

            TreeCovariance_IrFx (C20,
                                 NbTP,
                                 Time,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrFor,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvFx,
                                 One,
                                 CrIrForIrDom,
                                 CrIrForFx);

            TreeCovariance_IrFx (C21,
                                 NbTP,
                                 Time,
                                 MrIrDom,
                                 SvIrDom,
                                 FwIrFor,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvFx,
                                 CrIrForIrDom,
                                 One,
                                 CrIrDomFx);

            TreeCovariance_FxFx (C22,
                                 NbTP,
                                 Time,
                                 FwIrFor,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 FwIrFor,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvFx,
                                 SvFx,
                                 One,
                                 CrIrForIrDom,
                                 CrIrForFx,
                                 CrIrForIrDom,
                                 One,
                                 CrIrDomFx,
                                 CrIrForFx,
                                 CrIrDomFx,
                                 One);

            if (TreeType == TTYPE_EQDFX2IR)
            {
            TreeCovariance_IrEq (C30,
                                 NbTP,
                                 Time,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvEq,
                                 CrIrForIrDom,
                                 CrIrForEq);

            TreeCovariance_IrEq (C31,
                                 NbTP,
                                 Time,
                                 MrIrDom,
                                 SvIrDom,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvEq,
                                 One,
                                 CrIrDomEq);

            TreeCovariance_FxEq (C32,
                                 NbTP,
                                 Time,
                                 FwIrFor,
                                 MrIrFor,
                                 SvIrFor,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvFx,
                                 SvEq,
                                 CrIrForIrDom,
                                 CrIrForEq,
                                 One,
                                 CrIrDomEq,
                                 CrIrDomFx,
                                 CrFxEq);

            TreeCovariance_EqEq (C33,
                                 NbTP,
                                 Time,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 FwIrDom,
                                 MrIrDom,
                                 SvIrDom,
                                 SvEq,
                                 SvEq,
                                 One,
                                 CrIrDomEq,
                                 CrIrDomEq,
                                 One);
        }
            else /* TreeType == TTYPE_EQFFX2IR */
            {
                TreeCovariance_IrEq (C30,
                                     NbTP,
                                     Time,
                                     MrIrFor,
                                     SvIrFor,
                                     FwIrFor,
                                     MrIrFor,
                                     SvIrFor,
                                     SvEq,
                                     One,
                                     CrIrForEq);

                TreeCovariance_IrEq (C31,
                                     NbTP,
                                     Time,
                                     MrIrDom,
                                     SvIrDom,
                                     FwIrFor,
                                     MrIrFor,
                                     SvIrFor,
                                     SvEq,
                                     CrIrForIrDom,
                                     CrIrDomEq);

                TreeCovariance_FxEq (C32,
                                     NbTP,
                                     Time,
                                     FwIrFor,
                                     MrIrFor,
                                     SvIrFor,
                                     FwIrDom,
                                     MrIrDom,
                                     SvIrDom,
                                     FwIrFor,
                                     MrIrFor,
                                     SvIrFor,
                                     SvFx,
                                     SvEq,
                                     One,
                                     CrIrForEq,
                                     CrIrForIrDom,
                                     CrIrDomEq,
                                     CrIrForFx,
                                     CrFxEq);

                TreeCovariance_EqEq (C33,
                                     NbTP,
                                     Time,
                                     FwIrFor,
                                     MrIrFor,
                                     SvIrFor,
                                     FwIrFor,
                                     MrIrFor,
                                     SvIrFor,
                                     SvEq,
                                     SvEq,
                                     One,
                                     CrIrForEq,
                                     CrIrForEq,
                                     One);
            }
        }

        /* Cover first period as a special case (5 nodes regardless) */

        for (i = 0; i < TotDim; i++)
        {
            Width[i]     = 5;
            HalfWidth[i] = 2;

            R[i][i] = 1;
                
        }  /* for i */

        /* Dimensions 1 and 2 are always there */  
        Bottom1[0] = -2;
        Top1[0]    =  2;
        Bottom2[0] = (int *)DR_Array(INT, -2, 2);
        Top2[0]    = (int *)DR_Array(INT, -2, 2);
        if (Top2[0] == NULL || Bottom2[0] == NULL)
        {
            goto RETURN;
        }

        /* Bond vol integrals may be needed in two instances */
        if (TreeType == TTYPE_EQ1IR)
        {
            AIr1 = (double *)DR_Array(DOUBLE, 0, NbTP);
            if (AIr1 == NULL)
            {
                goto RETURN;
            }

        }
        if (ThirdDimOn)
        {    
            Bottom3[0] = (int **)DR_Array(INT_PTR, -2, 2);
            Top3[0]    = (int **)DR_Array(INT_PTR, -2, 2);
            if (Top3[0] == NULL || Bottom3[0] == NULL)
            {                                         
                goto RETURN;                          
            }                                        

            /* Allocate also the Af and Ad integrals */
            AIr2 = (double *)DR_Array(DOUBLE, 0, NbTP);
            AIr1 = (double *)DR_Array(DOUBLE, 0, NbTP);
            if (AIr1 == NULL || AIr2 == NULL)
            {                                         
                goto RETURN;                          
            }  

        }
        if (FourthDimOn)
        {    
            Bottom4[0] = (int ***)DR_Array(INT_D_PTR, -2, 2);
            Top4[0]    = (int ***)DR_Array(INT_D_PTR, -2, 2);
            if (Top4[0] == NULL || Bottom4[0] == NULL)
            {                                         
                goto RETURN;                          
            }                                        
        }

        for (i = -2; i <= 2; i++)
        {
            Bottom2[0][i] = -2;
            Top2[0][i]    =  2;
            if (ThirdDimOn)
            {    
                Bottom3[0][i] = (int *)DR_Array(INT, -2, 2);
                Top3[0][i]    = (int *)DR_Array(INT, -2, 2);
                if (Bottom3[0][i]==NULL || Top3[0][i]==NULL)
                {
                    goto RETURN;
                }

                if (FourthDimOn)
                {
                    Bottom4[0][i] = (int **)DR_Array(INT_PTR, -2, 2);
                    Top4[0][i]    = (int **)DR_Array(INT_PTR, -2, 2);
                    if (Bottom4[0][i]==NULL || Top4[0][i]==NULL)
                    {
                        goto RETURN;
                    }
                }

                for (j = -2; j <= 2; j++)
                {
                    Bottom3[0][i][j] = -2;
                    Top3[0][i][j]    =  2;

                    if (FourthDimOn)
                    {
                        Bottom4[0][i][j] = (int *)DR_Array(INT, -2, 2);
                        Top4[0][i][j]    = (int *)DR_Array(INT, -2, 2);
                        if (Bottom4[0][i][j]==NULL || Top4[0][i][j]==NULL)
                        {
                            goto RETURN;
                        }

                        for (k = -2; k <= 2; k++)
                        {
                            Bottom4[0][i][j][k] = -2;
                            Top4[0][i][j][k]    =  2;

                        }

                    } /* for k */

                }  /* for j */

            }  /* if */

        }  /* for i */


        /*  TREE LIMITS  */
        /*  Walk forward on the tree, calculating the tree limits and */
        /*  updating the total variances and composite correlations.  */
        /*                                                            */
        /*  Here we try to unify the approach  for all tree types  by */
        /*  using the variables SpotVolAsset2 and SpotVolAsset3.      */
        

                
        for (t = 0; t <= NbTP; t++)
        {   
            if (TreeType == TTYPE_EQ1IR)
            {
                vol[0] = SpotVol[0][t-1];
                vol[1] = SpotVolEq [t-1];

                r[1][0] = Rho1[t-1]; /* VEZI */
            }
            else if (TreeType == TTYPE_2IR)
            {
                vol[0] = SpotVol[0][t-1];
                vol[1] = SpotVol[1][t-1];

                r[1][0] = Rho1[t-1];
            }
            else if (TreeType == TTYPE_FX2IR)
            {
                vol[0] = SpotVol[0][t-1];
                vol[1] = SpotVol[1][t-1];
                vol[2] = SpotVolFx [t-1];

                r[1][0] = Rho1[t-1];
                r[2][0] = Rho2[t-1];
                r[2][1] = Rho3[t-1];
            }
            else if (TreeType == TTYPE_EQD2IR) /* VEZI */
            {
                vol[0] = SpotVol[0][t-1];
                vol[1] = SpotVol[1][t-1];
                vol[2] = SpotVolEq [t-1];

                r[1][0] = Rho1[t-1];
                r[2][0] = Rho2[t-1];
                r[2][1] = Rho3[t-1];
            }
            else if (TreeType == TTYPE_EQF2IR)
            {
                vol[0] = SpotVol[0][t-1];
                vol[1] = SpotVol[1][t-1];
                vol[2] = SpotVolEq [t-1];

                r[1][0] = Rho1[t-1];
                r[2][0] = Rho2[t-1];
                r[2][1] = Rho3[t-1];
            }
            else if ((TreeType == TTYPE_EQDFX2IR) ||
                     (TreeType == TTYPE_EQFFX2IR))
            {
                vol[0] = SpotVol[0][t-1];
                vol[1] = SpotVol[1][t-1];
                vol[2] = SpotVolFx [t-1];
                vol[3] = SpotVolEq [t-1];

                r[1][0] = Rho1[t-1];
                r[2][0] = Rho2[t-1];
                r[2][1] = Rho3[t-1];
                r[3][0] = Rho4[t-1];
                r[3][1] = Rho5[t-1];
                r[3][2] = Rho6[t-1];
            }
            else if (TreeType == TTYPE_2IR2F1D)
            {
                vol[0] = IrAweight[0][0][t-1];
                vol[1] = sqrt(SQUARE(IrAweight[0][1][t-1])+
                                           SQUARE(IrAweight[0][2][t-1]));
                vol[2] = IrAweight[1][0][t-1];

                /* Due to Aweight extrapolation, effective correlation
                 * may be different from input factor correlation
                 */
                r[1][0] = IrAweight[0][1][t-1]/vol[1]; 
                r[2][0] = Rho2[t-1];
                r[2][1] = Rho3[t-1];
            }

        
            /*   Jump coefficients. */
            CurrLengthJ = LengthJ[t-1];
            NextLengthJ = LengthJ[t];
            
            NextLength  = Length[t];
            du     = sqrt (JUMPCOEFF * CurrLengthJ);



               
            /* Limits of 1D tree */
            if (t > 0)
            {   
                /* cholesky decomposition for the increment */
                if(Hyb3_Triangulation4D(n, TotDim, r) != SUCCESS)
                {
                    DR_Error("Bad correlation structure (Hyb3_Tree_Limits)\n");
                    goto RETURN;
                }
       
                /*spatial discretization*/
                for (u=0; u < TotDim; ++u)
                {
                    for (v=0; v <= u; ++v)
                    {
                        J[u][v] = vol[u] * n[u][v] * du;
                    }
                }
                /* cholesky decomposition for the terminal asset distribution */
                if(Hyb3_Triangulation4D(N, TotDim, R) != SUCCESS)
                {
                    DR_Error("Bad terminal correlation structure (Hyb3_Tree_Limits)\n");
                    goto RETURN;
                }


                /* inverse of N */
                G[0][0] = 1;
                G[1][1] = 1 / N[1][1];
                G[1][0] = - G[1][1] * N[1][0] / N[0][0];
                if (ThirdDimOn)       
                {
                    G[2][2] = 1 / N[2][2];
                    G[2][1] = - G[2][2] * N[2][1] / N[1][1];
                    G[2][0] = - (G[2][1] * N[1][0] + G[2][2] * N[2][0]) / N[0][0];
                    if (FourthDimOn)  
                    {
                        G[3][3] = 1 / N[3][3];
                        G[3][2] = - G[3][3] * N[3][2] / N[2][2];
                        G[3][1] = - (G[3][2] * N[2][1] + G[3][3] * N[3][1]) / N[1][1];
                        G[3][0] = - (G[3][1] * N[1][0] + G[3][2] * N[2][0] + G[3][3] * N[3][0]) / N[0][0];
                    }
                }


                /* slopes */
                for (z=0; z<TotDim; ++z) 
                {
                    for (k=0; k<z; ++k)
                    {
                        a[z][k] = 0;
                        for (v=k; v<z; ++v)
                        {
                             a[z][k] += G[z][v] * J[v][k] / Vol[v];
                        }
                        a[z][k] *= N[z][z] * Vol[z];
                        a[z][k] += J[z][k];
                        a[z][k] /= - J[z][z];
                    }
                }

 

                /* Note the fabs () as they are no guarantee Jump[] are positive */
                /* and make sure we leave at least 5 nodes in the tree           */
                h = (int) ceil (NbIRSigmaMax * Vol[0] / fabs (J[0][0])); 
                h = MAX (h,  2);
                Centre = 0;
 
                /*  Minimum half size of the ellipse: slope + 1 node. This is to
                 *   avoid case(1) where there is no square of 9 nodes to branch to.
                 *   (case(2) is correct):
                 *                                       x
                 *                                   x   x
                 *               x               O   O   O
                 *           x   x               O   O   O
                 *       x   x   x   =======>    O   O   O
                 *       x   x                   x   x
                 *       x                       x
                 *
                 *           (1)                     (2)
                 */
                h2 = (int)ceil(fabs(a[1][0])) + 2;
                                             
                Bottom1[t] = Centre - h;
                Top1[t]    = Centre + h;
                
                /* Record biggest limits of the tree for memory allocation */
                HalfWidth[0] = MAX (HalfWidth[0], h);
                Width[0]     = 2 * HalfWidth[0] + 1;
                
            
                /* Allocate for the 2D limits and the 3D and 4D matrices of limits */
                Bottom2[t] = (int *)DR_Array(INT, Bottom1[t], Top1[t]);
                Top2[t]    = (int *)DR_Array(INT, Bottom1[t], Top1[t]);
                if (Bottom2[t] == NULL || Top2[t] == NULL)
                {
                    goto RETURN;
                }

                if (ThirdDimOn)
                {    
                    Bottom3[t] = (int **)DR_Array(INT_PTR, Bottom1[t], Top1[t]);
                    Top3[t]    = (int **)DR_Array(INT_PTR, Bottom1[t], Top1[t]);
                    if (Bottom3[t] == NULL || Top3[t] == NULL)  
                    {    
                        goto RETURN;
                    }

                    if (FourthDimOn && (t <= xT))
                    {
                        Bottom4[t] = (int ***)DR_Array(INT_D_PTR, Bottom1[t], Top1[t]);
                        Top4[t]    = (int ***)DR_Array(INT_D_PTR, Bottom1[t], Top1[t]);
                        if (Bottom4[t] == NULL || Top4[t] == NULL)  
                        {    
                            goto RETURN;
                        }
                    }
                }   
                for (i = Bottom1[t]; i <= Top1[t]; i++)
                {
                    Centre = NEAR_INT(i * a[1][0]);

                    xx[0] = i * J[0][0]/ Vol[0];
                    
                    y[0]  = G[0][0] * xx[0];

                    Delta  = (SQUARE(NbIRSigmaMax) - SQUARE(xx[0]));

                    /* Due to rounding, xx[0] may have been higher than NbSigmaDBL */
                    Delta = MAX(TINY, Delta); 

                    Delta = sqrt (Delta);
                    Delta *= N[1][1] / (J[1][1] / Vol[1]);

                    h = (int) ceil(Delta);
                    h = MAX(h, h2);


                    
                    Bottom2[t][i] = Centre - h;
                    Top2[t][i]    = Centre + h;

                    HalfWidth[1] = MAX (HalfWidth[1], h);
                    Width[1] = 2 * HalfWidth[1] + 1;
                    

                    /* Limits of 3D tree.They depend on first two dimensions */
                    if (ThirdDimOn)
                    {    

                        /* Memory allocation */
                        Bottom3[t][i] = (int *)DR_Array(INT, Bottom2[t][i], Top2[t][i]);
                        Top3[t][i]    = (int *)DR_Array(INT, Bottom2[t][i], Top2[t][i]);
                        if (Bottom3[t][i] == NULL || Top3[t][i] == NULL)
                        {
                            goto RETURN;
                        }

                        if (FourthDimOn && (t <= xT))
                        {
                            Bottom4[t][i] = (int **)DR_Array(INT_PTR, Bottom2[t][i], Top2[t][i]);
                            Top4[t][i]    = (int **)DR_Array(INT_PTR, Bottom2[t][i], Top2[t][i]);
                            if (Bottom4[t][i] == NULL || Top4[t][i] == NULL)
                            {
                                goto RETURN;
                            }

                        }


                        

                        h3 = (int) ceil(fabs(a[2][0]) + fabs(a[2][1])) + 2;

                        for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
                        {
                            Centre = NEAR_INT(i * a[2][0] + j * a[2][1]);

                            
                            xx[1] = (i * J[1][0] + j * J[1][1]) / Vol[1];
                             y[1] = G[1][0] * xx[0] + G[1][1] * xx[1];

                            Delta = SQUARE(NbSigmaDBL) - SQUARE(y[0]) - SQUARE(y[1]);

                            Delta = MAX(TINY, Delta); 

                            Delta = sqrt(Delta);

                            Delta *= N[2][2] / (J[2][2] / Vol[2]);

                            h = (int)ceil(Delta);
                            h = MAX(h,h3);

                            Bottom3[t][i][j] = Centre - h;
                            Top3[t][i][j]    = Centre + h;

                            HalfWidth[2] = MAX(HalfWidth[2], h);
                            Width[2] = 2 * HalfWidth[2] + 1;
                            


                            if (FourthDimOn)
                            {   
                                            
                                h4 = (int) ceil(fabs(a[3][0]) + fabs(a[3][1]) + fabs(a[3][2])) + 2;

                            /* Memory allocation */ 
                            if (t <= xT)
                            {
                                Bottom4[t][i][j] = (int *)DR_Array(INT, Bottom3[t][i][j], Top3[t][i][j]);
                                Top4[t][i][j]    = (int *)DR_Array(INT, Bottom3[t][i][j], Top3[t][i][j]);
                                if (Bottom4[t][i][j] == NULL || Top4[t][i][j] == NULL)
                                {
                                    goto RETURN;
                                }
                            }

                                for (k = Bottom3[t][i][j]; k <= Top3[t][i][j]; k++)
                                {
                                    Centre = NEAR_INT(i * a[3][0] + j * a[3][1] + k * a[3][2]);

                                    xx[2] = (i * J[2][0] + j * J[2][1] + k * J[2][2]) / Vol[2];

                                    y[2] = G[2][0] * xx[0] + G[2][1] * xx[1] + G[2][2] * xx[2];

                                    Delta = SQUARE(NbSigmaDBL) - SQUARE(y[0]) - SQUARE(y[1]) - SQUARE(y[2]);

                                    Delta = MAX(TINY, Delta); 
    
                                    Delta = sqrt(Delta);

                                    Delta *= N[3][3] / (J[3][3] / Vol[3]);

                                    h = (int)ceil(Delta);
                                    h = MAX(h,h4);

                                if (t <= xT)
                                {
                                    Bottom4 [t][i][j][k] = Centre - h;
                                    Top4    [t][i][j][k] = Centre + h;
                                }

                                    HalfWidth [3] = MAX(HalfWidth[3], h);
                                    Width     [3] = 2 * HalfWidth[3] + 1;

                                } /* for k */
                            }  /* if 4th dim*/
                        }  /* for j */
                    }  /* if 3rd dim*/
                }  /* for i */
            } /* if t>0 */
            
             
    

            /* CORRELATION MATRIX */
            /* Compute correlation matrix at next time step. */


            /*  1 - First dimension */

            Var1 *= exp (- 2. * LBeta[0] * NextLength);
            Var1 += SQUARE(SpotVol[0][t])*NextLength * Hyb3_ExpDecay(2.*LBeta[0],NextLength);
            Vol1 = sqrt (Var1);

            C[0][0] = Var1;

            Vol[0]  = Vol1;

            /*  2 - Second dimension (may be IR, EQ, or second fgn factor ) */
            if (TreeType == TTYPE_EQ1IR)
            {
                /* In this case, 2nd dim is equity */

                /* Equity fwd is composite asset */
                SpotVar2 += SQUARE(SpotVolAsset2[t]) * NextLength;

                SpotCov1 *= exp (- LBeta[0] * Length[t]);
                SpotCov1 += Rho1[t] * SpotVol[0][t] * SpotVolAsset2[t] * NextLength 
                            * Hyb3_ExpDecay (LBeta[0], NextLength);

                /* Refresh accumulators */
                for (k = 0; k <= 3; k++)
                    I[k] = 0.;

                
                for (j = t, T = 0.; j >=0; j--)
                {    

                    /* Locals just to improve performance */
                    double   deltaT       = Length[j];
                    double   expBeta1     = exp(-LBeta[0]*deltaT);
                    double   sqrSpotVolF  = SQUARE(SpotVol[0][j]);
                    double   ABeta1       = Hyb3_ExpDecay(LBeta[0],deltaT);


                    /* Calculate the A integrals (the A's correct from vol */
                    /* of the spot rate to vol of bond of finite maturity) */    
                    AIr1[j] += FwdRate[0][CvDiff[0]][t]
                             * Hyb3_ExpDecay(LBeta[0],NextLength) * exp(-LBeta[0]*T);
                    
                    /* Variance terms */                          
                    I[0] += sqrSpotVolF * AIr1[j] * AIr1[j] * deltaT;
                    I[1] += 2.*Rho1[j]*SpotVolAsset2[j]*SpotVol[0][j]*AIr1[j]*deltaT;

                    /* Covariance terms */
                    I[2] *= expBeta1;
                    I[2] += sqrSpotVolF * AIr1[j] * ABeta1 * deltaT;
                                                
                    T += deltaT;
                        

                }  /* for j */

                    
                Var2 = SpotVar2 + I[0] + I[1];
                Cov1 = SpotCov1 + I[2];    

                Vol2  = sqrt (Var2);
                CRho1 = Cov1 / Vol1 / Vol2;
            }
            else
            {

                /* In all other cases, 2nd dim is the domestic IR 
                 * or fgn second factor 
                 */

                Var2 *= exp (- 2. * LBeta[1] * NextLength);
                Var2 += SQUARE(SpotVolAsset2[t]) * NextLength 
                      * Hyb3_ExpDecay(2.*LBeta[1],NextLength);
                Vol2 = sqrt (Var2);

                Cov1 *= exp (- (LBeta[0] + LBeta[1]) * NextLength);
                Cov1 +=   Rho1[t] * SpotVol[0][t] * SpotVolAsset2[t] 
                    * NextLength * Hyb3_ExpDecay(LBeta[0] + LBeta[1], NextLength);

                CRho1 = Cov1 / Vol1 / Vol2;
                if ((CRho1 - 1 > TINY) || (CRho1 + 1 < TINY))
                {
                    DR_Error("CRho1 is not a correl number (>1 or < -1)\n");
                    goto RETURN;
                }
            }

            C[1][0] = Cov1;
            C[1][1] = Var2;

            Vol[1]  = Vol2;

            R[1][0] = CRho1;

            /*  3 - Third dimension (may be FX, EQ or domestic IR) */
            if (ThirdDimOn && !FourthDimOn)
            {    
                
                /* Terms in var and covariance which are same for FX & EQ */
                SpotVar3 *= exp (- 2. * LBeta[2] * NextLength); 
                SpotVar3 += SQUARE(SpotVolAsset3[t]) * NextLength
                           * Hyb3_ExpDecay(2.*LBeta[2],NextLength);

                SpotCov2 *= exp (- (LBeta[0] + LBeta[2]) * NextLength);
                SpotCov2 += Rho2[t] * SpotVol[0][t] * SpotVolAsset3[t] 
                            * NextLength * Hyb3_ExpDecay(LBeta[0] + LBeta[2], NextLength);

                SpotCov3 *= exp (- (LBeta[1] + LBeta[2]) * NextLength);
                SpotCov3 += Rho3[t] * SpotVolAsset2[t] * SpotVolAsset3[t] 
                            * NextLength * Hyb3_ExpDecay(LBeta[1] + LBeta[2], NextLength);
                
                if (ThirdDimIR == FALSE ) /* should be FX or EQ then*/
                {

                    /* Refresh accumulators */
                    for (k = 0; k <= 8; k++)
                        I[k] = 0.;

                    DomA = 0.0;
                    ForA = 0.0;

                
                    for (j = t, T = 0.; j >=0; j--)
                    {    

                        /* Locals just to improve performance */
                        double   deltaT = Length[j];
                        double   sqrSpotVolF  = SQUARE(SpotVol[0][j]);
                        double   sqrSpotVolD  = SQUARE(SpotVolAsset2[j]);
                        double   ABeta1  = Hyb3_ExpDecay(LBeta[0],deltaT);
                        double   ABeta2  = Hyb3_ExpDecay(LBeta[1],deltaT);
                    

                    
                        DomA *= exp(-LBeta[1] * Length[j]);
                        DomA += FwdRate[1][CvDiff[1]][j] * Hyb3_ExpDecay(LBeta[1],Length[j]);

                        ForA *= exp(-LBeta[0] * Length[j]);
                        ForA += FwdRate[0][CvDiff[0]][j] * Hyb3_ExpDecay(LBeta[0],Length[j]);


                        /* Variance terms */  
                    
                        I[0] += sqrSpotVolF * ForA * ForA * Length[j];
                        I[1] += sqrSpotVolD * DomA * DomA * Length[j];
                    
                        I[2] -= 2.*Rho1[j] * SpotVol[0][j]    * SpotVolAsset2[j]
                                 * DomA * ForA * Length[j];

                        I[3] -= 2.*Rho2[j] * SpotVolAsset3[j] * SpotVol[0][j]
                                 * ForA * Length[j];

                        I[4] += 2.*Rho3[j] * SpotVolAsset3[j] * SpotVolAsset2[j]
                                 * DomA * Length[j];

                        /* Covariance terms */
                    
                        I[5] += sqrSpotVolF * ForA * ABeta1 * Length[j] * exp(-LBeta[0] * T);
                    
                        I[6] += sqrSpotVolD * DomA * ABeta2 * Length[j] * exp(-LBeta[1] * T);

                   
                        I[7] += SpotVol[0][j] * SpotVolAsset2[j] * DomA * ABeta1 * Length[j]
                                * Rho1[j] * exp(-LBeta[0] * T);

                        I[8] += SpotVol[0][j] * SpotVolAsset2[j] * ForA * ABeta2 * Length[j]
                                * Rho1[j] * exp(-LBeta[1] * T);
                            
                        T += deltaT;
                        
                    }  /* for j */
                
                
                
                    if (EquityOn)
                    {
                        if (EquityOnIr1)
                        {
                            Var3 = SpotVar3 + I[0] - I[3];
                            Cov2 = SpotCov2 + I[5];
                            Cov3 = SpotCov3 + I[7];
                        }
                        else
                        {
                            Var3 = SpotVar3 + I[1] + I[4];
                            Cov2 = SpotCov2 + I[8];
                            Cov3 = SpotCov3 + I[6];
                        }
                    }
                    else  /* FX is on */
                    {
                        Var3 = SpotVar3 + I[0] + I[1] + I[2] + I[3] + I[4];
                        Cov2 = SpotCov2 - I[5] + I[7];
                        Cov3 = SpotCov3 + I[6] - I[8];
                    }
                }
                else
                {
                    Var3 = SpotVar3;
                    Cov2 = SpotCov2;
                    Cov3 = SpotCov3;

                }


                Vol3  = sqrt (Var3);
                CRho2 = Cov2 / Vol1 / Vol3;
                CRho3 = Cov3 / Vol2 / Vol3;
                if((CRho2 - 1.0 > TINY) || (CRho2 + 1.0 < TINY))
                {
                    DR_Error("CRho2 is ill-defined\n");
                    goto RETURN;
                }

                if((CRho3 - 1.0 > TINY) || (CRho3 + 1.0 < TINY))
                {
                    DR_Error("CRho3 is ill-defined\n");
                    goto RETURN;
                }

                if((1- CRho1*CRho1 - CRho2*CRho2 - CRho3*CRho3 + 2*CRho1*CRho2*CRho3) < TINY)
                {
                    DR_Error("Resulting composite correl mtx is not positive\n");
                    goto RETURN;
                }

                C[2][0] = Cov2;
                C[2][1] = Cov3;
                C[2][2] = Var3;

                Vol[2]  = Vol3;

                R[2][0] = CRho2;
                R[2][1] = CRho3;

            }  /* if 3rd dim ON and not 4th dim ON */

            if (FourthDimOn) /* do dimensions 3 and 4 in the new way */
            {
                i = t+1;

                if (i <= NbTP)
                {
                    C[2][0] = C20[i];
                    C[2][1] = C21[i];
                    C[2][2] = C22[i];

                    C[3][0] = C30[i];
                    C[3][1] = C31[i];
                    C[3][2] = C32[i];
                    C[3][3] = C33[i];

                    if (C[2][2] < 0)
                    {
                        DR_Error("Negative Fx Variance.\n");
                        goto RETURN;
                    }

                    if (C[3][3] < 0)
                    {
                        DR_Error("Negative Eq Variance.\n");
                        goto RETURN;
                    }

                    Vol[2] = sqrt(C[2][2]);
                    Vol[3] = sqrt(C[3][3]);

                    R[2][0] = C[2][0] / Vol[2] / Vol[0];
                    R[2][1] = C[2][1] / Vol[2] / Vol[1];
                
                    R[3][0] = C[3][0] / Vol[3] / Vol[0];
                    R[3][1] = C[3][1] / Vol[3] / Vol[1];
                    R[3][2] = C[3][2] / Vol[3] / Vol[2];

                    if((R[2][0] - 1.0 > TINY) || (R[2][0] + 1.0 < TINY))
                    {
                        DR_Error("Rho20 is ill-defined\n");
                        goto RETURN;
                    }

                    if((R[2][1] - 1.0 > TINY) || (R[2][1] + 1.0 < TINY))
                    {
                        DR_Error("Rho21 is ill-defined\n");
                        goto RETURN;
                    }

                    if((R[3][0] - 1.0 > TINY) || (R[3][0] + 1.0 < TINY))
                    {
                        DR_Error("Rho30 is ill-defined\n");
                        goto RETURN;
                    }

                    if((R[3][1] - 1.0 > TINY) || (R[3][1] + 1.0 < TINY))
                    {
                        DR_Error("Rho31 is ill-defined\n");
                        goto RETURN;
                    }

                    if((R[3][2] - 1.0 > TINY) || (R[3][2] + 1.0 < TINY))
                    {
                        DR_Error("Rho32 is ill-defined\n");
                        goto RETURN;
                    }

                    det = 1 -       R[1][0] * R[1][0]
                            -       R[2][0] * R[2][0]
                            -       R[2][1] * R[2][1]
                            -       R[3][0] * R[3][0]
                            -       R[3][1] * R[3][1]
                            -       R[3][2] * R[3][2]

                            + 2 * ( R[2][1] * R[3][1] * R[3][2]
                            +       R[2][0] * R[3][0] * R[3][2]
                            +       R[1][0] * R[3][0] * R[3][1]
                            +       R[1][0] * R[2][0] * R[2][1] )

                            +       R[1][0] * R[1][0] * R[3][2] * R[3][2]
                            +       R[2][0] * R[2][0] * R[3][1] * R[3][1]
                            +       R[2][1] * R[2][1] * R[3][0] * R[3][0]

                            - 2 * ( R[1][0] * R[2][1] * R[3][2] * R[3][0]
                            +       R[1][0] * R[2][0] * R[3][1] * R[3][2]
                            +       R[2][0] * R[2][1] * R[3][0] * R[3][1] );

                    if (det < TINY)
                    {
                        DR_Error("Resulting composite correl mtx is not positive\n");
                        goto RETURN;
                    }
                }
            }

            if (t <= xT)
            {
                xHalfWidth [0] = HalfWidth [0];
                xWidth     [0] = Width     [0];

                if (TotDim > 1)
                {
                    xHalfWidth [1] = HalfWidth [1];
                    xWidth     [1] = Width     [1];

                    if (TotDim > 2)
                    {
                        xHalfWidth [2] = HalfWidth [2];
                        xWidth     [2] = Width     [2];

                        if (TotDim > 3)
                        {
                            xHalfWidth [3] = HalfWidth [3];
                            xWidth     [3] = Width     [3];
                        }
                    }
                }
            }
        }

        /* Finally repeat the limit values for the node NbTP + 1 */

        Bottom1[NbTP+1] = Bottom1[NbTP];
        Top1[NbTP+1]    = Top1[NbTP];

        Bottom2[NbTP+1] = (int *)DR_Array(INT, Bottom1[NbTP+1], Top1[NbTP+1]);
        Top2[NbTP+1]    = (int *)DR_Array(INT, Bottom1[NbTP+1], Top1[NbTP+1]);
        if (Bottom2[NbTP+1] == NULL || Top2[NbTP+1] == NULL)
        {
            goto RETURN;
        }

        if (ThirdDimOn)
        {    
            Bottom3[NbTP+1] = (int **)DR_Array(INT_PTR, Bottom1[NbTP+1], Top1[NbTP+1]);
            Top3[NbTP+1]    = (int **)DR_Array(INT_PTR, Bottom1[NbTP+1], Top1[NbTP+1]);
            if (Bottom3[NbTP+1] == NULL || Top3[NbTP+1] == NULL)  
            {                                           
                goto RETURN;                            
            }                                           
        }  

        if (FourthDimOn && (NbTP+1 <= xT))
        {    
            Bottom4[NbTP+1] = (int ***)DR_Array(INT_D_PTR, Bottom1[NbTP+1], Top1[NbTP+1]);
            Top4[NbTP+1]    = (int ***)DR_Array(INT_D_PTR, Bottom1[NbTP+1], Top1[NbTP+1]);
            if (Bottom4[NbTP+1] == NULL || Top4[NbTP+1] == NULL)  
            {                                           
                goto RETURN;                            
            }                                           
        }

        for (i = Bottom1[NbTP+1]; i <= Top1[NbTP+1]; i++)
        {
            Bottom2[NbTP+1][i] = Bottom2[NbTP][i];
            Top2[NbTP+1][i]    = Top2[NbTP][i];

            if (ThirdDimOn)
            {    
                Bottom3[NbTP+1][i] = (int *)DR_Array(INT, Bottom2[NbTP+1][i], 
                                                Top2[NbTP+1][i]);
                Top3[NbTP+1][i]    = (int *)DR_Array(INT, Bottom2[NbTP+1][i], 
                                                Top2[NbTP+1][i]);
                if (Bottom3[NbTP+1][i]==NULL || Top3[NbTP+1][i]==NULL)
                {
                    goto RETURN;
                }


                for (j = Bottom2[NbTP+1][i]; j <= Top2[NbTP+1][i]; j++)
                {
                    Bottom3[NbTP+1][i][j] = Bottom3[NbTP][i][j];
                    Top3[NbTP+1][i][j]    = Top3[NbTP][i][j];
                }
            }

            if (FourthDimOn && (NbTP+1 <= xT))
            {    
                Bottom4[NbTP+1][i] = (int **)DR_Array(INT_PTR, Bottom2[NbTP+1][i], Top2[NbTP+1][i]);
                Top4   [NbTP+1][i] = (int **)DR_Array(INT_PTR, Bottom2[NbTP+1][i], Top2[NbTP+1][i]);

                if (Bottom4 [NbTP+1][i]==NULL ||
                    Top4    [NbTP+1][i]==NULL)
                {
                    goto RETURN;
                }

                for (j = Bottom2[NbTP+1][i]; j <= Top2[NbTP+1][i]; j++)
                {
                    Bottom4[NbTP+1][i][j] = (int *)DR_Array(INT, Bottom3[NbTP+1][i][j], Top3[NbTP+1][i][j]);
                    Top4   [NbTP+1][i][j] = (int *)DR_Array(INT, Bottom3[NbTP+1][i][j], Top3[NbTP+1][i][j]);

                    if (Bottom4 [NbTP+1][i][j] == NULL ||
                        Top4    [NbTP+1][i][j] == NULL)
                    {
                        goto RETURN;
                    }

                    for (k = Bottom3[NbTP+1][i][j]; k <= Top3[NbTP+1][i][j]; k++)
                    {
                        Bottom4 [NbTP+1][i][j][k] = Bottom4 [NbTP][i][j][k];
                        Top4    [NbTP+1][i][j][k] = Top4    [NbTP][i][j][k];
                    }
                }
            }
        }


        
        status = SUCCESS;

      RETURN:

 
        Free_DR_Array(AIr1, DOUBLE, 0, NbTP);
        Free_DR_Array(AIr2, DOUBLE, 0, NbTP);

        if (TreeType == TTYPE_2IR2F1D)
        {
            Free_DR_Array(SpotVolAsset2, DOUBLE, -1, NbTP);
            Free_DR_Array(Rho1, DOUBLE, -1, NbTP);
        }



        return (status);

}  /* Hyb3_Tree_Limits */



/*****  Hyb3_Tree_Limits_IR1  *****************************************************/
/**
 	Calculates the limits of one factor interest rate tree.
 
 
*/
int     Hyb3_Tree_Limits_IR1(
            int     *Width,                 /**< (O) Maximum node index        */
            int     *HalfWidth,             /**< (O) Minimum node index        */
            int     *Top,                   /**< (O) Upper limit of the tree   */
            int     *Bottom,                /**< (O) Lower limit of the tree   */
            double  *SpotVol,               /**< (I) Inst vol of short rate    */
            double  *Length,                /**< (I) Length of each time step  */
            double  *LengthJ,               /**< (I) Time step for jump size   */
            int     NbTP,                   /**< (I) Total number of steps     */
            double  NbSigmaDBL,             /**< (I) Nb of std devs for trim   */
            double  Beta)                   /**< (I) Mean reversion coeff      */
{
    double

        BaseVol2,             /* Square of base volatility at current node  */
        BaseVol,              /* Base volatility at current node            */
        Jump,                 /* Jump size at current period                */
        du;                   /* Jump coefficient                           */

    int 
              
        CurrentMax,           /* Current maximum of nodes in each dimension */
        i,                    /*                                            */
        status = FAILURE;     /* Error status = FAILURE initially           */





        BaseVol  = BaseVol2 = 0.;
        *Width = 5;
        *HalfWidth = 2;
                
        for (i = 0; i <= NbTP; i++)
        {        
            du = sqrt (JUMPCOEFF * LengthJ[i-1]);
            Jump = SpotVol[i-1] * du;

            
            /* CurrentMax is such that CurrentMax*Jump = NbSigmaDBL*BaseVol */
            CurrentMax = (int) ceil (NbSigmaDBL * BaseVol / Jump); 
            CurrentMax = MAX (CurrentMax, 2);
        
            /* The bottom and the top of the tree can't be more */
            /* than CurrentMax away from the center.            */
            Bottom[i] = -CurrentMax;     
            Top[i]    =  CurrentMax;    

            /* The base volatility is the integral   */ 
            /* of the spot volatility (cf basevol.c) */                                         
            BaseVol2 *= exp (- 2. * Beta * Length[i]);                                 
            BaseVol2 += SpotVol[i]*SpotVol[i]*Length[i]* Hyb3_ExpDecay(2.*Beta,Length[i]);
            BaseVol   = sqrt (BaseVol2);

            /* Record the largest node number for memory allocation */                                                                                    
            *HalfWidth = MAX (*HalfWidth,  CurrentMax);     
            *Width = 2 * (*HalfWidth) + 1;  

        }  /* for i */

        Bottom[NbTP+1] = Bottom[NbTP];
        Top[NbTP+1]    = Top[NbTP];

        
        status = SUCCESS;


        return (status);



}  /* Hyb3_Tree_Limits_IR1 */
/*****  Hyb3_Tree_Limits_IR2  *****************************************************/
/**
     Calculates the tree limits for an interest rate with up to two factor
 
 
*/
int     Hyb3_Tree_Limits_IR2(
            int     *Width,                 /**< (O) Maximum node index         */
            int     *HalfWidth,             /**< (O) Minimum node index         */
            int     *Top1,                  /**< (O) Upper limit of the 1D tree */
            int     *Bottom1,               /**< (O) Lower limit of the 1D tree */
            int     **Top2,                 /**< (O) Upper limit of the 2D tree */
            int     **Bottom2,              /**< (O) Lower limit of the 2D tree */
            double  **Aweight,              /**< (I) Aweight of short rate      */
            int     NbFactor,               /**< (I) Number of Factors          */
            double  *Length,                /**< (I) Length of each time step   */
            double  *LengthJ,               /**< (I) Time step for jump size    */
            int     NbTP,                   /**< (I) Total number of steps      */
            double  NbSigmaDBL,             /**< (I) Nb of std devs for trim    */
            double  *Beta)                  /**< (I) Mean reversion coeff       */
{
    double

        Var1, Var2,          /* Variance of factors at current node        */
        Vol1, Vol2,          /* Volatility of factors                      */
        Cov1,                /* Covariance of factors                      */
        Rho1,                /* Correlation of factors                     */
        a1,                  /* Ellipsoid axis                             */
        Jump[3],             /* Jump size at current period                */
        du,                  /* Jump coefficient                           */
        Delta, x1;               /* Intermediate results                        */

    int 
              
        CurrentMax,          /* Current maximum of nodes in each dimension */
        Center,              /* Center of the ellipse                      */
        h2,                  /* Minimum width of ellipsoid                 */
        i,t,                 /*                                            */
        status = FAILURE;    /* Error status = FAILURE initially           */





    Vol1  = Var1 = 0.;
    Vol2  = Var2 = 0.;
    Cov1  = Rho1 = 0 ;

    /* 
     *  Start the tree with 5 nodes in each dimension.
     */

    for (i = 0; i < 2; i++)
    {
        Width[i]     = 5;
        HalfWidth[i] = 2;
    }

    Bottom1[0] = -2;
    Top1[0]    =  2;

    if (NbFactor >= 2)
    {  
        Bottom2[0] = (int *) DR_Array (INT, -2, 2);
        Top2[0]    = (int *) DR_Array (INT, -2, 2);

        if (  (Bottom2[0] == NULL)
            || (Top2[0]    == NULL))
        {
            goto RETURN;
        }


        for (i = -2; i <= 2; i++)
        {
            Bottom2[0][i] = -2;
            Top2[0][i]    =  2;

        }  /* for i */
    }  /* if */


    for (t = 0; t <= NbTP; t++)
    {        
        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        Jump[0] = Aweight[0][t-1] * du;
        if (NbFactor >= 2)
        {  
            Jump[1] = Aweight[1][t-1] * du;
            Jump[2] = Aweight[2][t-1] * du;
        }


        if (t>0)
        {
            /* CurrentMax is such that CurrentMax*Jump = NbSigmaDBL*Vol1 */
            CurrentMax = (int) ceil (NbSigmaDBL * Vol1 / Jump[0]); 
            CurrentMax = MAX (CurrentMax, 2);

            /* The bottom and the top of the tree can't be more */
            /* than CurrentMax away from the center.            */
            Bottom1[t] = -CurrentMax;     
            Top1[t]    =  CurrentMax;    



            /* Record the largest node number for memory allocation */
            HalfWidth[0] = MAX (HalfWidth[0],  CurrentMax);     
            Width[0] = 2 * (HalfWidth[0]) + 1;  

            /*
             *    Limits of the 2D tree. They depend on first dimension.
             */

            if (NbFactor >= 2)
            {

                /*
                 *  Allocate memory for the 2D limits and the 3D matrices of limits.
                 */

                Bottom2[t] = (int *) DR_Array (INT, Bottom1[t], Top1[t]);
                Top2[t]    = (int *) DR_Array (INT, Bottom1[t], Top1[t]);

                if (  (Bottom2[t] == NULL)
                    || (Top2[t]    == NULL))
                {
                    goto RETURN;
                }
                /* Slope of 2D ellipse axis */
                a1 = (Rho1 * Vol2 / Vol1 * Jump[0] - Jump[1]) / fabs (Jump[2]);

                /*  Minimum half size of the ellipse: slope + 1 node. This is to
                 *   avoid case(1) where there is no square of 9 nodes to branch to.
                 *   (case(2) is correct):
                 *                                       x
                 *                                   x   x
                 *               x               O   O   O
                 *           x   x               O   O   O
                 *       x   x   x   =======>    O   O   O
                 *       x   x                   x   x
                 *       x                       x
                 *
                 *           (1)                     (2)
                 */

                h2 = (int) ceil (fabs(a1)) + 2;

                for (i = Bottom1[t]; i <= Top1[t]; i++)
                {
                    x1 = i * Jump[0] / Vol1;

                    Delta = MAX (TINY, (1. - Rho1*Rho1) * (NbSigmaDBL*NbSigmaDBL - x1*x1) * SQUARE (Vol2 / Jump[2]));
                    Delta = sqrt (Delta);

                    CurrentMax = (int) ceil (Delta);
                    CurrentMax = MAX (CurrentMax, h2);

                    Center = NEAR_INT (i * a1);

                    Bottom2[t][i] = Center - CurrentMax;
                    Top2[t][i]    = Center + CurrentMax;

                    HalfWidth[1] = MAX (HalfWidth[1], CurrentMax);
                    Width[1]     = 2 * HalfWidth[1] + 1;
                } /* for i */
            } /*if 2 Factor */
        } /*if t>0*/

        /*
         *  Compute correlation matrix at next time step.
         */

        /* The base volatility is the integral   */ 
        /* of the spot volatility (cf basevol.c) */                                         
        Var1 *= exp (- 2. * Beta[0] * Length[t]);                                 
        Var1 += Aweight[0][t]*Aweight[0][t]*Length[t]* Hyb3_ExpDecay(2.* Beta[0],Length[t]);
        Vol1   = sqrt (Var1);


        if (NbFactor >= 2)
        {  
            Var2 *= exp (- 2. * Beta[1] * Length[t]);
            Var2 += Aweight[1][t] * Aweight[1][t] * Length[t] * Hyb3_ExpDecay (2. * Beta[1], Length[t]);
            Var2 += Aweight[2][t] * Aweight[2][t] * Length[t] * Hyb3_ExpDecay (2. * Beta[1], Length[t]);

            Vol2 = sqrt (Var2);

            Cov1 *= exp (- (Beta[0] + Beta[1]) * Length[t]);
            Cov1 += Aweight[0][t] * Aweight[1][t] * Length[t] * Hyb3_ExpDecay (Beta[0] + Beta[1], Length[t]);

            Rho1 = Cov1 / Vol1 / Vol2;
        } /*if 2 factor */



    }  /* for t */

    Bottom1[NbTP+1] = Bottom1[NbTP];
    Top1[NbTP+1]    = Top1[NbTP];
    
    if(NbFactor >= 2)
    {
        Bottom2[NbTP+1] = (int *)DR_Array(INT, Bottom1[NbTP+1], Top1[NbTP+1]);
        Top2[NbTP+1]    = (int *)DR_Array(INT, Bottom1[NbTP+1], Top1[NbTP+1]);
        if (Bottom2[NbTP+1] == NULL || Top2[NbTP+1] == NULL)
        {
            goto RETURN;
        }

        for (i = Bottom1[NbTP+1]; i <= Top1[NbTP+1]; i++)
        {
            Bottom2[NbTP+1][i] = Bottom2[NbTP][i];
            Top2[NbTP+1][i]    = Top2[NbTP][i];
        }
    }




    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Hyb3_Tree_Limits_IR2: could not allocate memory for the tree!");
    }


    return (status);



}  /* Hyb3_Tree_Limits_IR2 */

/*****  Hyb3_Tree_FreeLimits_IR2  *****************************************************/
/**
     Free the memory allocated for the tree limits
 
 
*/
int     Hyb3_Free_TreeLimits_IR2(
            int     *Top1,                  /**< (O) Upper limit of the 1D tree */
            int     *Bottom1,               /**< (O) Lower limit of the 1D tree */
            int     **Top2,                 /**< (O) Upper limit of the 2D tree */
            int     **Bottom2,              /**< (O) Lower limit of the 2D tree */
            int     NbFactor,               /**< (I) Number of Factors          */
            int     NbTP)                   /**< (I) Total number of steps      */
{


    int t;



    if (NbFactor >=2)
    {

        if ((Top1 != NULL) && (Bottom1 != NULL) &&
            (Top2 != NULL) && (Bottom2 != NULL))
        {
            for (t = 0; t <= NbTP+1; t++)
            {    
                if (Top2[t] != NULL)
                {
                    Free_DR_Array(Top2[t], INT, Bottom1[t], Top1[t]);
                }

                if (Bottom2[t] != NULL)
                {
                    Free_DR_Array(Bottom2[t], INT, Bottom1[t], Top1[t]);
                }

            }  /* for t */

        }  /* if */
    }/* if */

 

    return (SUCCESS);



}  /* Hyb3_Free_TreeLimits_IR2 */


/*****  Hyb3_Tree_Weights  ********************************************************/
/**
     Calculates the tree diffusion weights (i.e. orhtogonalisation factors 
     multiplied by  the spot volatilities and the time steps).
  
 */
int  Hyb3_Tree_Weights
         (double     **Aweights,     /**< (O) Max six weights per time step    */
          double      *Rho[6],       /**< (I) Correlations at each time step   */
          int          TreeType,     /**< (I) 2 or 3-D, Equity or FX           */
          int          NbTP,         /**< (I) Total number of nodes in tree    */
          double      *IrAweight[2][3], /**< (I) Spot vols of interest rate dims  */
          double      *SpotVolFx,
          double      *SpotVolEq) /**< (I) Spot vol of non-IR asset (EQ/FX) */
                  
{
    double

       *SpotVol1 = NULL,
       *SpotVol2 = NULL,
       *SpotVol3 = NULL,
       *SpotVol4 = NULL,
       *Rho1 = NULL, /* Local var used for corr of 1st vs. 2nd assets */
       *Rho2 = NULL, /* Local var used for corr of 1st vs. 3rd assets */
       *Rho3 = NULL, /* Local var used for corr of 2nd vs. 3rd assets */
       *Rho4 = NULL,
       *Rho5 = NULL,
       *Rho6 = NULL;
   double 
        Vol2ndFact, x;

    double n[4][4];
    double r[4][4];
    double vol[4];

    int 
        t, u, v, w,
        NbDims = 0,
        status = FAILURE;     /* Error status = FAILURE initially */

  

        switch (TreeType)
        {

            case TTYPE_1IR:   /* Used for CET runs */

                NbDims = 1;
                SpotVol1 = IrAweight[0][0];
                break;

            case TTYPE_2IR:

                NbDims = 2;
                SpotVol1 = IrAweight[0][0];
                SpotVol2 = IrAweight[1][0];
                Rho1 = Rho[0];
                break;

            case TTYPE_EQ1IR:

                NbDims = 2;
                SpotVol1 = IrAweight[0][0];
                SpotVol2 = SpotVolEq;
                Rho1 = Rho[3]; /* Corr IR1 vs. EQ */
                break;

            case TTYPE_FX2IR:

                NbDims = 3;
                SpotVol1 = IrAweight[0][0];
                SpotVol2 = IrAweight[1][0];
                SpotVol3 = SpotVolFx;
                Rho1 = Rho[0];
                Rho2 = Rho[1];
                Rho3 = Rho[2];
                break;

            case TTYPE_2IR2F1D:

                NbDims=3;
                SpotVol3= IrAweight[1][0];
                Rho2=Rho[0];
                Rho3=Rho[1];
                break;

            case TTYPE_1IR2F:

                NbDims=2;
                break;

            case TTYPE_EQD2IR:
            case TTYPE_EQF2IR:
            case TTYPE_EQC2IR:

                NbDims = 3;
                SpotVol1 = IrAweight[0][0];
                SpotVol2 = IrAweight[1][0];
                SpotVol3 = SpotVolEq;
                Rho1 = Rho[0];
                Rho2 = Rho[3];
                Rho3 = Rho[4];
                break;

            case TTYPE_EQDFX2IR:
            case TTYPE_EQFFX2IR:

                NbDims = 4;
                SpotVol1 = IrAweight[0][0];
                SpotVol2 = IrAweight[1][0];
                SpotVol3 = SpotVolFx;
                SpotVol4 = SpotVolEq;
                Rho1 = Rho[0];
                Rho2 = Rho[1];
                Rho3 = Rho[2];
                Rho4 = Rho[3];
                Rho5 = Rho[4];
                Rho6 = Rho[5];
                break;

            default:

                DR_Error("Hyb3_Tree_Weights: Unknown tree type !\n");
                goto RETURN;

        } /* End switch */

        switch(TreeType)
        {
        case TTYPE_2IR2F1D:     /*In these case, Aweight for the 2 factors   */
        case TTYPE_1IR2F:       /* foreign IR are given                      */
                    /* Weights */       
            for (t = -1; t<= NbTP; t++)
            {        

                /* replace by t-1 */
                Aweights[0][t] = IrAweight[0][0][t];          
            
                /* Second dimension ON */
                if (NbDims > 1)
                {
                    /* Aweights have already been computed */
                    Aweights[1][t] = IrAweight[0][1][t];    
                    Aweights[2][t] = IrAweight[0][2][t];

                
                    /* Third dimension ON */
                    if (NbDims > 2)
                    {
                        Aweights[3][t] = Rho2[t];
                        Vol2ndFact=sqrt(Aweights[1][t]*Aweights[1][t]+Aweights[2][t]*Aweights[2][t]);


                        Aweights[4][t] = (Rho3[t]*Vol2ndFact - Aweights[1][t]*Rho2[t]) 
                                      / Aweights[2][t];
                        x=  1.0 - Aweights[3][t]*Aweights[3][t]
                                -Aweights[4][t]*Aweights[4][t];
                        if (x<TINY)
                        {
                           goto RETURN;
                        }
                        Aweights[5][t] = sqrt(x);
                    
                        Aweights[3][t] *= SpotVol3[t];   
                        Aweights[4][t] *= SpotVol3[t];   
                        Aweights[5][t] *= SpotVol3[t];
                    }
                }
             
            }  /* for i */
            break;

        default:
        /* Weights */       
            for (t = -1; t<= NbTP; t++)
            {       
                /* replace by t-1 */
                Aweights[0][t] = 1.;          
                Aweights[0][t] *= SpotVol1[t] ;  

                /* Second dimension ON */
                if (NbDims > 1)
                {
                    /* Rho1 and sqrt(1 - Rho1^2) are orthogonalisation factors */
                    Aweights[1][t] = Rho1[t];    
                    Aweights[2][t] = sqrt(1.0 - Rho1[t]*Rho1[t]);
            
                    Aweights[1][t] *= SpotVol2[t] ;   
                    Aweights[2][t] *= SpotVol2[t] ; 
                
                    /* Third dimension ON */
                    if (NbDims > 2)
                    {
                        Aweights[3][t] = Rho2[t];

                        Aweights[4][t] = (Rho3[t] - Rho1[t]*Rho2[t]) 
                                      / sqrt(1.0 - Rho1[t]*Rho1[t]);       
                        Aweights[5][t] = sqrt(  1.0 - Rho1[t]*Rho1[t] - Rho2[t]*Rho2[t]
                                       - Rho3[t]*Rho3[t]+2*Rho1[t]*Rho2[t]*Rho3[t])
                                      /sqrt(1.0 - Rho1[t]*Rho1[t]);
                        /* NB: weight[5] should be above a minimum level, but */
                        /* this should have been checked at input level.      */
                        Aweights[3][t] *= SpotVol3[t];   
                        Aweights[4][t] *= SpotVol3[t];   
                        Aweights[5][t] *= SpotVol3[t];

                        /* Fourth dimension ON */
                        if (NbDims > 3)
                        {
                            r[1][0] = Rho1[t];
                            r[2][0] = Rho2[t];
                            r[2][1] = Rho3[t];
                            r[3][0] = Rho4[t];
                            r[3][1] = Rho5[t];
                            r[3][2] = Rho6[t];

                            vol[0] = SpotVol1[t];
                            vol[1] = SpotVol2[t];
                            vol[2] = SpotVol3[t];
                            vol[3] = SpotVol4[t];

                            n[0][0] = 1;
    
                            n[1][0] = r[1][0] / n[0][0];

                            n[1][1] = sqrt(1 - SQUARE(n[1][0]));

                            n[2][0] = r[2][0] / n[0][0];

                            n[2][1] = (r[2][1] - n[2][0] * n[1][0]) / n[1][1];

                            n[2][2] = sqrt(1 - SQUARE(n[2][0]) - SQUARE(n[2][1]));

                            n[3][0] = r[3][0] / n[0][0];

                            n[3][1] = (r[3][1] - n[3][0] * n[1][0]) / n[1][1];

                            n[3][2] = (r[3][2] - n[3][0] * n[2][0] - n[3][1] * n[2][1]) / n[2][2];

                            n[3][3] = sqrt(1 - SQUARE(n[3][0]) - SQUARE(n[3][1]) - SQUARE(n[3][2]));

                            w = 0;

                            for (u=0; u<4; ++u)
                            {
                                for (v=0; v<=u; ++v)
                                {
                                    Aweights[w][t] = vol[u] * n[u][v];

                                    ++w;
                                }
                            }
                        }
                    }
                }
            }
        } /*End switch */
            

            



    status = SUCCESS;

  RETURN:

    return (status);



}  /* Hyb3_Tree_Weights */



/*****  Hyb3_Tree_Correlations  ****************************************************/
/**
 *   Interpolates the input correlation curves to obtain correlations along-
 *   side the  tree timeline .
 */
int  Hyb3_Tree_Correlations
        (double  **Rho,              /**< (O) Correlations at each time step   */
         int       NbTP,             /**< (I) Total number of nodes in tree    */
         long     *TPDate,           /**< (I) Dates in tree                    */
         int       TreeType,         /**< (I) Whether 3D, EqDom, EqFor etc.    */
         int       NbCorrFx,         /**< (I) Nb of pts in FX corr curves      */
         long      CorrDateFx[MAXNBDATE+1],   /**< (I) Benchmark FX corr dates */
         double    RhoBenchFx[3][MAXNBDATE+1],/**< (I) Corr on benchmark dates */
         int       NbCorrEq,                  /**< (I) Nb of pts in corr curves*/
         long      CorrDateEq[MAXNBDATE],     /**< (I) Benchmark eq corr dates */
         double    RhoBenchEq[3][MAXNBDATE+1])/**< (I) Corr on benchmark dates */
{
    double

        Rho1=9.9,  /*  Local var used for corr of IR1 vs. IR2     */
        Rho2=9.9,  /*  Local var used for corr of IR1 vs. FX      */
        Rho3=9.9,  /*  Local var used for corr of IR2 vs. FX      */
        Rho4=9.9,  /*  Local var used for corr of IR1 vs. Equity  */
        Rho5=9.9,  /*  Local var used for corr of IR2 vs. Equity  */
        Rho6=9.9;  /*  Local var used for corr of FX  vs. Equity  */

        /* NB: In the multi-currency case, Ir1 is the foreign */
        /*     currency, Ir2 is the domestic currency.        */

        /* NB2: In the single currency case,the only relevant */
        /*      correlation will be Rho4.                     */


    long
        CorrDate1, 
        CorrDate2;

    int 
        i, j, 
        i1 = 0, 
        i2 = 0,
        status = FAILURE;     /* Error status = FAILURE initially  */

       
        /* There are two very distinct cases to consider: */
        /* single-currency and multi-currency.            */

        /* 1 - SINGLE CURRENCY  (TTYPE_1IR)  */
        if (TreeType == TTYPE_1IR)
        {
            /* One factor, one IR, nothing to be done. */
            /* Must be running in CET mode.            */
            return(SUCCESS);

        }
        /* 2 - SINGLE CURRENCY WITH EQUITY (TTYPE_EQ1IR)  */
        else if (TreeType == TTYPE_EQ1IR)
        {

            for (i = 1, i1 = 0; i <= NbCorrEq; i++) 
            {        
                CorrDate1 = CorrDateEq[i-1];
                CorrDate2 = CorrDateEq[i];
                
                i2 = NbTP;
                while ((TPDate[i2] > CorrDate2) && (i2 > 0))
                    i2--;        	                       
                
                if (i2 == i1 - 1)	  
                    continue;
                
                for (j = i1; j <= i2; j++) 
                {        
                    dlinterp (TPDate[j], &Rho4, 
                              CorrDate1, CorrDate2, 
                              RhoBenchEq[0][i-1], RhoBenchEq[0][i]);
                    
                    Rho[3][j] = Rho4;
                        
                }  /* for j */
                            
                i1 = i2 + 1;
            
            }  /* for i */
            
            
            for (i = i2 + 1; i <= NbTP; i++)  
            {         
                Rho[3][i] = Rho4;   
                        
            }
            
            /* We need one extra node to build the lattice */
            Rho[3][-1] = Rho[3][0]; 

        }
        else
        /* 3 - MULTI-CURRENCY (all other types) */
        {


            /* First three correlation curves will be FX-related */
            
            /* Interpolate 3 correlation curves in the FX data */
            for (i = 1, i1 = 0; i <= NbCorrFx; i++) 
            {        
                CorrDate1 = CorrDateFx[i-1];
                CorrDate2 = CorrDateFx[i];
                
                /* Node index i2 falling just before the current vol point  */
                i2 = NbTP;
                while ((TPDate[i2] > CorrDate2) && (i2 > 0)) 
                    i2--;        	
                
                /* Ignore this point if it is falling onto the previous one */
                if (i2 == i1 - 1)
                    continue;
                    
                for (j = i1; j <= i2; j++) /* Note last node i2 is included */
                {        
                    dlinterp (TPDate[j],  &Rho1, 
                              CorrDate1,    CorrDate2, 
                              RhoBenchFx[0][i-1], RhoBenchFx[0][i]);
            
                    dlinterp (TPDate[j],  &Rho2, 
                              CorrDate1,    CorrDate2, 
                              RhoBenchFx[1][i-1], RhoBenchFx[1][i]);
            
                    dlinterp (TPDate[j],  &Rho3, 
                              CorrDate1,    CorrDate2, 
                              RhoBenchFx[2][i-1], RhoBenchFx[2][i]);
                
                    Rho[0][j] = Rho1;
                    Rho[1][j] = Rho2;
                    Rho[2][j] = Rho3;
                            
                }  /* for j */
                                
                i1 = i2 + 1;
            
            }  /* for i */
            
            /* If tree extend beyond last vol point fill correlations with  */
            /* the last value found (correlations will be const thereafter) */ 
            for (i = i2 + 1; i <= NbTP; i++)     
            {                               
                    Rho[0][i] = Rho1;    
                    Rho[1][i] = Rho2;
                    Rho[2][i] = Rho3;
                            
            }  /* for i */
            
            /* Need one extra node to build the lattice */
            for (j = 0; j < 3; j++)
            {
                Rho[j][-1] = Rho[j][0];
            }
            
            
            
            /* THREE MORE CORRELATION CURVES IF EQUITY ON */
            
            /* If equity ON, interpolate 3 correlation curves in the eq data */
            if (   (TreeType == TTYPE_EQD2IR)
                 ||(TreeType == TTYPE_EQF2IR)
                 ||(TreeType == TTYPE_EQC2IR)
                 ||(TreeType == TTYPE_EQDFX2IR)
                 ||(TreeType == TTYPE_EQFFX2IR))
            {
                for (i = 1, i1 = 0; i <= NbCorrEq; i++) 
                {        
                    CorrDate1 = CorrDateEq[i-1];
                    CorrDate2 = CorrDateEq[i];
                    
                    i2 = NbTP;
                    while ((TPDate[i2] > CorrDate2) && (i2 > 0))
                        i2--;        	                       
                    
                    if (i2 == i1 - 1)	  
                        continue;
                    
                    for (j = i1; j <= i2; j++) 
                    {        
                        dlinterp (TPDate[j], &Rho4, 
                                  CorrDate1, CorrDate2, 
                                  RhoBenchEq[0][i-1], RhoBenchEq[0][i]);
                        dlinterp (TPDate[j], &Rho5, 
                                  CorrDate1, CorrDate2, 
                                  RhoBenchEq[1][i-1], RhoBenchEq[1][i]);
                        dlinterp (TPDate[j], &Rho6, 
                                  CorrDate1, CorrDate2, 
                                  RhoBenchEq[2][i-1], RhoBenchEq[2][i]);
                        
                        Rho[3][j] = Rho4;
                        Rho[4][j] = Rho5;
                        Rho[5][j] = Rho6;
                            
                    }  /* for j */
                                
                    i1 = i2 + 1;
            
                }  /* for i */
            
            
                for (i = i2 + 1; i <= NbTP; i++)  
                {         
                    Rho[3][i] = Rho4;
                    Rho[4][i] = Rho5;                   
                    Rho[5][i] = Rho6;        
                            
                }  /* for i */
            
                /* We need one extra node to build the lattice */
                for (j = 3; j<6 ; j++)
                    Rho[j][-1] = Rho[j][0];    
            
            } /* if Equity ON */

        } /* If then else on single/multi ccy */



    status = SUCCESS;

    return (status);



}  /* Hyb3_Tree_Correlations */

/**** Hyb3_Transform_SmileParams ********************************************
 * 
 * The function does:
 *
 * a) converts Skew/Crv/Level to internal A1, A2, A3 in the tanh/sech formula
 * b) check validity (smile functions is well defined and positive)
 *
 *****************************************************************************/

int Hyb3_Transform_SmileParams(double Sk,
                               double Crv,
                               double Level,
                               double *aa,
                               double *bb,
                               double *cc,
                               char   *Error)
{
    int status = FAILURE;

    double a = -9999.0;
    double b = -9999.0;
    double c =  9999.0;

    if (Level < 0.0)
    {
        strcpy(Error, "Smile level is negative at higher end");
        goto RETURN;
    }

    if (IS_EQUAL(Level, 0.0)) /* flat local vol */ 
    {
        if (!(IS_EQUAL(Crv, 0.0)) || 
            !(IS_EQUAL(Sk,  0.0)))
        {
            { 
                strcpy(Error, "If a3 is zero, a1 and a2 must be zero too");
                goto RETURN;
            }
        }

        a = 0.0;
        b = 0.0;
        c = 1.0;
    }
    else /* Level > 0 */
    {
        /* no curvature */
        if (IS_EQUAL(Crv, 0.0)) 
        {
            /* no skew */
            if (IS_EQUAL(Sk, 0.0))
            {
                a = 0.0;
                b = 0.0;
                c = 1.0;
            }
            else /* skew */
            {
                c = fabs(Sk) / Level;	
                a = Sk / c;
                b = 0.0;
            }
        }
        else /* curvature, pos or neg */
        { 
            b =  pow(- fabs(Sk) + sqrt(Sk * Sk + 4.0 * fabs(Crv) * Level), 2.0) /
                 (4.0 * Crv); 
            c =  sqrt(Crv / b);
            a =  Sk / c;
        }
    }

    if ((1 - a + b <= 0) ||
        (1 + a + b <= 0))
    {
        strcpy(Error, "Smile function is negative at one or two ends");
        goto RETURN;
    }

    /* check for roots of the smile function, as a function of exp(cx) */
    /* any real roots must be negative                                 */

    if ((a*a) - (2*b) - 1 >= 0)
    {
        if (b + sqrt((a*a) - (2*b) - 1) > 0)
        {
            strcpy(Error, "Smile function has positive roots");
            goto RETURN;
        }
    }

    status = SUCCESS;

  RETURN:

    *aa = a;
    *bb = b;
    *cc = c;

    return status;
}

/**** Hyb3_Tree_SmileParams **************************************************
 * 
 * The function does:
 *
 * a) converts Skew/Crv/Level to internal A1, A2, A3 in the tanh/sech formula
 * b) extend parameters to every timepoint in the tree timeline
 * c) populates tMin and tMax which are then used in Hyb3_FillGrid
 *
 * NOTE1: for every t, the same smile parameters are used from tMin[t] to tMax[t]
 * NOTE2: Smile Idx is Forward Looking: SmileIdx[t] points
 *        to the smile parameters that apply to the period [t,t+1]
 * NOTE3: We currently ignore smile dates which are denser than
 * the tree time-line.
 *
 *****************************************************************************/

int Hyb3_Tree_SmileParams(int *SmileIdx, /**< (O)*/
                          double *A1,    /**< (O)*/
                          double *A2,    /**< (O)*/
                          double *A3,    /**< (O)*/
                          double *A1C,   /**< (O)*/
                          double *A2C,   /**< (O)*/
                          double *A3C,   /**< (O)*/
                          long   *TPDate,
                          int    NbTP,
                          double a1[MAXNBDATE],
                          double a2[MAXNBDATE],
                          double a3[MAXNBDATE],
                          int    NbSmilePt,
                          long   SmileDate[MAXNBDATE],
                          long   tMin[MAXNBDATE],
                          long   tMax[MAXNBDATE])
{
    int status = FAILURE;
    int i,j;
    long    OffsetBefore,OffsetSofar;
    double Sk, Crv, Level, a, b, c;
    long    Curr_tMin, CurrIdx;

    char Error[MAXBUFF] = "";

    if(NbSmilePt> MAXNBDATE) goto RETURN;

    for(i = 0;i < NbSmilePt; i++)
    {   
        Sk    = a1[i];
        Crv   = a2[i];
        Level = a3[i];

        if (Hyb3_Transform_SmileParams(Sk,
                                       Crv,
                                       Level,
                                       &a,
                                       &b,
                                       &c,
                                       Error) == FAILURE)
        {
            DR_Error("Bad smile params at line %d.\n"
                     "%s.\n",
                     i+1,
                     Error);

            goto RETURN;
        }
        
        A1[i] = a;
        A2[i] = b;
        A3[i] = c;
    }
    

    /* fill extra point to populate Fx @ TP = 0*/
    A1[-1] = A1[0];
    A2[-1] = A2[0];
    A3[-1] = A3[0];

    /* fill in SmileIdx array */
    OffsetBefore = OffsetSofar = 0L;
    for(i = 0; i < NbSmilePt; i++)
    {
        OffsetSofar = GetDLOffset(NbTP,
                                  TPDate,
                                  SmileDate[i],
                                  CbkLOWER);

        if(OffsetSofar < 0L) goto RETURN;
        if(OffsetSofar == OffsetBefore) continue;

        /* remark: <= to ensure filling of last data point */
        for(j = OffsetBefore; j <= OffsetSofar; j++)
        {
            SmileIdx[j] = i;
        }
        OffsetBefore = OffsetSofar;
    }

    /* AK: original version: if(OffsetSofar < NbTP -1) */
    /* not sure why last point (NbTP) was not filled   */
    if(OffsetSofar < NbTP-1 )
    {
        for(j=OffsetSofar ; j < NbTP ; j++)
        {
            SmileIdx[j] = NbSmilePt-1;
        }
    }


    /* fill in tMin, tMax arrays */
    CurrIdx = SmileIdx[0];
    Curr_tMin = 0;
    for(j = 0; j < NbTP; j++)
    {
        if (CurrIdx != SmileIdx[j])
        {
            for (i=Curr_tMin; i<j; i++)
            {
                tMin[i] = Curr_tMin;
                tMax[i] = j-1;
            }
            CurrIdx = SmileIdx[j];
            Curr_tMin = j;
        }
    }

    /* AK: Fill up to NbTP effectively: algorithm might have to be cleaned  */
    /* at some stage (to avoid confusion on length of tree NbTP or NbTP+1)  */
    for (i=Curr_tMin; i<=j; i++)
    {
        tMin[i] = Curr_tMin;
        tMax[i] = j-1;
    }

    SmileIdx[-1] = SmileIdx[0];

    for(i = 0; i <= NbTP; i++)
    {
       A1C[i] = A1[SmileIdx[i]];
       A2C[i] = A2[SmileIdx[i]];
       A3C[i] = A3[SmileIdx[i]];
    }

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Hyb3_Tree_SmileParams: failed.");
    }
    return(status);

}


