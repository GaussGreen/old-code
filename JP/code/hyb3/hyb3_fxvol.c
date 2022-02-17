/****************************************************************************/
/*      Calculation of fx spot volatility by bootstrapping benchmark fx     */
/*      option volatilities.                                                */
/****************************************************************************/
/*      FXVOL.C                                                             */
/****************************************************************************/

/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

#undef   MAXNBTP
#define  MAXNBTP            2048  /* Max nb of time points in array passed to FX vol funcs -  */


/** * NOTE: Composite vols functions in this file *****************************
 *
 *  Top level function that compute composite fx vol for a given list of expiry
 *  ---------------------------------------------------------------------------
 *  Hyb3_Get_FxVol     (calls v1)  (stickyturbo uses this -- needs to change [20-Jan-04]) (LOGNORMAL FX)
 *  Hyb3_Get_FxVol2    (calls v1)
 *
 *  Top level func that calc comp. fx vol for every pt in the tree timeline
 *  -----------------------------------------------------------------------
 *  Hyb3_FX_Vol (calls v2)
 *
 *	Top level func that computes the spot FX vol at every pt in the timeline
 *  ------------------------------------------------------------------------
 *  Hyb3_Get_FxSpotVol (calls v3)  (stickyturbo uses this -- needs to change [20-Jan-04]) (LOGNORMAL FX)
 *  Hyb3_Get_TreeFxSpotVol (calls v4)  -- used in Hyb3_Build_Tree
 *
 *  (v1) Hyb3_MultiFac_FxVol       - low level with no FX smile (i.e. lognormal)
 *  (v2) Hyb3_MultiFac_FxVol2      - low level with FX smile
 *  (V3) Hyb3_MultiFac_Spot_FxVol  - low level with no FX smile (i.e. lognormal)
 *  (V4) Hyb3_MultiFac_Spot_FxVol2 - low level with FX smile
 *
 *  products should stop calling functions that use (v1) or (v3) (e.g. stickyturbo [20-Jan-04])
 *
 */

/** * NOTE: Foward FX functions in this file *****************************
 *
 *  Hyb3_Forward_FX - computes FFX for every pt in the tree-timeline
 *  Hyb3_FX_Fwd     - computes FFX for one given maturity
 *  Hyb3_ForwardFX  - more clever FFX function for a given maturity
 *
 *  NOTE: 20-Jan-04: Hyb3_FX_Fwd should be replaced by Hyb3_ForwardFX. 
 *        However, stickyturbo uses Hyb3_FX_Fwd it.
 *        
 */



/** **  Hyb3_FX_Vol  ***********************************************************************
*
*  Computes the composite FX Vols at every time point of the tree timeline.
*  
*  Effectively, it calls the Hyb3_MultiFac_FxVol2 which 
*  a) uses the analytical approximation based on Henrik's paper
*  b) is the exact same function as in srm3
*
*
*******************************************************************************************/


int Hyb3_FX_Vol(double   *FxVolCurve,        /**< (O) Fx volatility at each time step */
           double   *SpotFxVol,         /**< (I) Inst fx vol at each time step   */
           int      DomNbFac,           /**< (I) Dom Nb of Factor   */
           int      ForNbFac,           /**< (I) For Nb of Factor   */
           double  *RhoFxDomFac,      /**<(I) Corr. between fx and Dom factors     */
           double  *RhoFxForFac,      /**<(I) Corr between fx and For factors      */
           double  *RhoDomFacForFac,  /**<(I) Corr between Dom and For factors     */
           double  *RhoDomFac,        /**<(I) Corr between Dom factors             */
           double  *RhoForFac,        /**<(I) Corr between For fcators             */
           double   *SpotVol[2],        /**< (I) Inst vol of short rate at node  */
           double   *FwdRateDom,        /**< (I) Pass Forward Domestic           */
           double   *FwdRateFor,        /**< (I) Pass Forward Foreign            */
           double   DomAlpha[3],        /**< (I) dom factor weights              */ 
           double   ForAlpha[3],        /**< (I) for factor weights              */
           long     *VolDate,           /**< (I) fwd fx maturity                 */
           long     *VolIntegrationDate,/**< (I) expiry of option                */
           double   *A1C,               /**< (I) FX smile param, [0,NbTP-2]      */
           double   *A2C,               /**< (I) FX smile param, [0,NbTP-2]      */
           double   *A3C,               /**< (I) FX smile param, [0,NbTP-2]      */
           double   BetaDom[3],         /**< (I) Mean reversion coefficient (dom)*/
           double   BetaFor[3],         /**< (I) Mean reversion coefficient (for)*/
           long     *TPDate,            /**< (I) date of each pt.                */
           int      NbTP)               /**< (I) Total nb of nodes or time steps */
{

    int
        k,
        status = FAILURE;        /* Error status = FAILURE initially */

    double *SDomSpotVol = NULL;
    double *SForSpotVol = NULL;

    SDomSpotVol    =   (double *)DR_Array(DOUBLE, -1, NbTP+1);
    SForSpotVol    =   (double *)DR_Array(DOUBLE, -1, NbTP+1);

    if ((SDomSpotVol == NULL) ||
        (SForSpotVol == NULL))
    {
        DR_Error("Unable to allocate memory.");
        goto RETURN;
    }

    /*Change spot vol format to match srm3 definitions (division by alpha factors)*/
    for (k = 0 ; k < NbTP; k++)
    {
        SDomSpotVol[k] = SpotVol[1][k] / DomAlpha[0];
        SForSpotVol[k] = SpotVol[0][k] / ForAlpha[0];
    }

    if(Hyb3_MultiFac_FxVol2( FxVolCurve,
                        SDomSpotVol,       
                        SForSpotVol,     
                        DomAlpha,        
                        ForAlpha,     
                        BetaDom,         
                        BetaFor,        
                        FwdRateDom,
                        FwdRateFor,
                        DomNbFac,                 
                        ForNbFac,                
                        RhoFxDomFac,          
                        RhoFxForFac,          
                        RhoDomFacForFac,        
                        RhoDomFac,              
                        RhoForFac,            
                        SpotFxVol,      
                        A1C,
                        A2C,
                        A3C,
                        VolDate,       
                        VolIntegrationDate,        
                        (NbTP+1),                    
                        TPDate,          
                        (NbTP+1)) != SUCCESS){

        goto RETURN;
    }


    /*Add -1 point*/
    FxVolCurve[-1] = FxVolCurve[0];


    status = SUCCESS;

    RETURN:
    if(SDomSpotVol!=NULL){
        Free_DR_Array(SDomSpotVol, DOUBLE, -1, NbTP+1);
    }
    if(SForSpotVol!=NULL){
        Free_DR_Array(SForSpotVol, DOUBLE, -1, NbTP+1);
    }

    if (status == FAILURE)
    {
        DR_Error("Hyb3_FX_Vol: failed");
    }
    return (status);
}





int Hyb3_For_FX_Vol(
                    double   *FxVolCurve,        /**< (O) Fx volatility at each time step */
                    double   *SpotFxVol,         /**< (I) Inst fx vol at each time step   */
                    double   *RhoCurve[6],       /**< (I) Correlations at each time step  */
                    double   *SpotVol[2],        /**< (I) Inst vol of short rate at node  */
                    double   *FwdRateDom,        /**< (I) Pass Forward Domestic           */
                    double   *FwdRateFor,        /**< (I) Pass Forward Foreign            */
                    double   DomAlpha[3],        /**< (I) dom factor weights              */ 
                    double   ForAlpha[3],        /**< (I) for factor weights              */
                    long     VolDate,            /**< (I) fwd fx maturity                 */
                    long     ReviseDate,         /**< (I) expiry of option                */
                    double   *A1C,               /**< (I) FX smile param, [0,NbTP-2]      */
                    double   *A2C,               /**< (I) FX smile param, [0,NbTP-2]      */
                    double   *A3C,               /**< (I) FX smile param, [0,NbTP-2]      */
                    double   BetaDom[3],         /**< (I) Mean reversion coefficient (dom)*/
                    double   BetaFor[3],         /**< (I) Mean reversion coefficient (for)*/
                    long     *TPDate,            /**< (I) date of each pt.                */
                    int      NbTP)               /**< (I) Total nb of nodes or time steps */
{

    int
        k,
        status = FAILURE;        /* Error status = FAILURE initially */

    double *SDomSpotVol = NULL;
    double *SForSpotVol = NULL;

    SDomSpotVol    =   (double *)DR_Array(DOUBLE, -1, NbTP+1);
    SForSpotVol    =   (double *)DR_Array(DOUBLE, -1, NbTP+1);

    if ((SDomSpotVol == NULL) ||
        (SForSpotVol == NULL))
    {
        DR_Error("Unable to allocate memory.");
        goto RETURN;
    }

    /*Change spot vol format to match srm3 definitions (division by alpha factors)*/
    for (k = 0 ; k < NbTP; k++)
    {
        SDomSpotVol[k] = SpotVol[1][k] / DomAlpha[0];
        SForSpotVol[k] = SpotVol[0][k] / ForAlpha[0];
    }

    if(Hyb3_MultiFac_FxVol2(
                        FxVolCurve,
                        SDomSpotVol,
                        SForSpotVol,
                        DomAlpha,
                        ForAlpha,
                        BetaDom,
                        BetaFor,
                        FwdRateDom,
                        FwdRateFor,
                        1,
                        1,     
                        RhoCurve[2],
                        RhoCurve[1],
                        RhoCurve[0],
                        NULL,
                        NULL,
                        SpotFxVol,
                        A1C,
                        A2C,
                        A3C,
                        &VolDate,
                        &ReviseDate,
                        1,
                        TPDate,
                        (NbTP+1))!=SUCCESS)
    {
        goto RETURN;
    }


    /*Add -1 point*/
    /* FxVolCurve[-1] = FxVolCurve[0]; */ /* should be commented out ?!? */


    status = SUCCESS;

    RETURN:
    if(SDomSpotVol!=NULL){
        Free_DR_Array(SDomSpotVol, DOUBLE, -1, NbTP+1);
    }
    if(SForSpotVol!=NULL){
        Free_DR_Array(SForSpotVol, DOUBLE, -1, NbTP+1);
    }

    if (status == FAILURE)
    {
        DR_Error("Hyb3_FX_Vol: failed");
    }
    return (status);
}




 
/*****  Hyb3_Forward_FX  **********************************************************/
/**
 *	 Calculates the fx forward at every timepoint in the tree timeline
 *
 */
int   Hyb3_Forward_FX(double  *FxFwd,        /**< (O) FX forward at each node      */                 
                 int      NbTP,         /**< (I) Total nb of nodes in timeline*/                 
                 double  *ZCoupon[2][3],/**< (I) Zero cp price at each node   */
                 int     *ZDisc,        /**< (I) Zero to use for discounting  */
                 double   SpotFx)       /**< (I) Spot fx                      */
{
    int
        i,
        status = FAILURE;           /* Error status = FAILURE initially     */

  
    /* FX forwards are calculated using the COF zero curve */
    for (i = 0; i <= NbTP; i++)
    {
        FxFwd[i]     = SpotFx * ZCoupon[0][ZDisc[0]][i] / ZCoupon[1][ZDisc[1]][i];
                
    }  /* for i */

    /* Finally calculate the FX forward for node NbTP + 1 */
    FxFwd[NbTP+1] = SpotFx * ZCoupon[0][ZDisc[0]][NbTP+1] 
                         / ZCoupon[1][ZDisc[1]][NbTP+1];
          
    status = SUCCESS;
    return (status);

}  /* Hyb3_Forward_FX */

/** Hyb3_FXMidNode************************************************************
 *
 *  Computes the mid-node for obtaining the centre of the tree at each timepoint
 *  This mid-node is added to the FwdFX in lattice
 *
 *************************************************************************/
int   Hyb3_FXMidNode(double  *MidNode,      /**< (O) Centre of the tree in FX dim */
                int      NbTP,         /**< (I) Total nb of nodes in timeline*/
                double  *Length,       /**< (I) Time step in years           */
                double  *SpotFxVol,    /**< (I) FX spotvol alongside timeline*/
                int     *SmileIdx,
                double  *FxVol, 
                double  NbSigmaMax,
                double  *A1,
                double  *A2,
                double  *A3)               
{
    int
        i,Idx,
        status = FAILURE;           /* Error status = FAILURE initially     */

    double
        Vol2 = 0.0;                 /* Square of vol at the current node    */
    
    double  gdashL;

    
    for (i = 0; i <= NbTP; i++)
    {
        MidNode[i] = -0.5* Vol2 ;
        Idx        = SmileIdx[i];


        if(Hyb3_GDashfunc(&(gdashL),
                     1.0,           /* evaluate function at x = 1.0 */
                     A1[Idx],
                     A2[Idx],
                     A3[Idx]) == FAILURE) goto RETURN;
                     
        Vol2 += SpotFxVol[i] * SpotFxVol[i] * Length[i] * gdashL;
        
    }  /* for i */

    
    MidNode[NbTP+1] = -0.5 * Vol2;

        
    status = SUCCESS;
RETURN:
    return (status);

}  /* Hyb3_FXMidNode */

/*************** Hyb3_FX_Fwd  ************************************/
/** calculates the forward fx rate for a specified Maturity    
  
  
 *************************************************************/
int Hyb3_FX_Fwd(    double      *FwdFx,     /**< (O) */
               double      SpotFx,     /**< (I) */
               T_CURVE     *DomCurve,  /**< (I) */
               T_CURVE     *ForCurve,  /**< (I) */
               long        ValueDate,
               long        Maturity)   /**< (I) */
                
{

    double  DomDisc,ForDisc; /* disc. fact. in dom and for maturing  */
                            /* on input Maturity Date               */

    double  ZeroRateL;      /* needed for GetZeroPriceRate function, no other use*/

    int     status = FAILURE;

    if ( FwdFx == NULL || DomCurve == NULL
         || ForCurve == NULL)
    {
        DR_Error ("Memory Alloc failed in Hyb3_FX_Fwd\n");
        goto RETURN;
    }

    /* get DomDisc */
    if (GetZeroPriceRate(&ZeroRateL,
                         &DomDisc,
                         Maturity,
                         DomCurve) != SUCCESS)
    {
        goto RETURN;
    }

    /* ForDisc */
    if (GetZeroPriceRate(&ZeroRateL,
                         &ForDisc,
                         Maturity,
                         ForCurve) != SUCCESS)
    {
        goto RETURN;
    }

    *FwdFx = SpotFx * ForDisc / DomDisc;

    status = SUCCESS;
RETURN:
    return(status);
}/* Hyb3_FX_Fwd */


/*****  Hyb3_FX_Fwd_Ratio  ****************************************************/
/**
 *       Utility routine to calculate the ratio of two deterministic fwd fx
 *       
 */
int Hyb3_FX_Fwd_Ratio (double      *Ratio,        /**< (O) FX1 / FX2 */
                  T_CURVE     *DomCurve,      /**< (I) */
                  T_CURVE     *ForCurve,      /**< (I) */
                  long        Maturity1,
                  long        Maturity2)       /**< (I) */
                
{

    double  DomDisc1, DomDisc2, ForDisc1, ForDisc2; 

    double  ZeroRateL;      /* needed for GetZeroPriceRate function, no other use*/

    int     status = FAILURE;

    if ( Ratio == NULL || DomCurve == NULL
         || ForCurve == NULL)
    {
        DR_Error ("Memory Alloc failed in Hyb3_FX_Fwd\n");
        goto RETURN;
    }

    /* get DomDisc */
    if (GetZeroPriceRate(&ZeroRateL,
                  &DomDisc1,
                  Maturity1,
                  DomCurve) != SUCCESS)
    {
        goto RETURN;
    }

    if (GetZeroPriceRate(&ZeroRateL,
                  &DomDisc2,
                  Maturity2,
                  DomCurve) != SUCCESS)
    {
        goto RETURN;
    }

    /* ForDisc */
    if (GetZeroPriceRate(&ZeroRateL,
                  &ForDisc1,
                  Maturity1,
                  ForCurve) != SUCCESS)
    {
        goto RETURN;
    }
    
    if (GetZeroPriceRate(&ZeroRateL,
                  &ForDisc2,
                  Maturity2,
                  ForCurve) != SUCCESS)
    {
        goto RETURN;
    }

    *Ratio = (ForDisc1 / DomDisc1) / (ForDisc2 / DomDisc2) ;

    status = SUCCESS;
RETURN:
    return(status);
}/* Hyb3_FX_Fwd */


/********* Hyb3_MultiFac_FxVol **********************************/
/** 1-Assumes and tests whether Voldate is on the input time line.
    2-Assumes TPDate[0] = FxBaseDate
    3- the correlation inputs are with respect to the exponential factors
  
 *************************************************************/

int Hyb3_MultiFac_FxVol(
        double  *CompositeVol,  /**< (O) Composite FX Vol for FX Expiries*/
        double  *DomSpotVol,    /**< (I)Instant.dom int.rate vol*/
        double  *ForSpotVol,    /**< (I)Instant.for.int rate vol*/
        double  *DomAlpha,      /**< (I)*/
        double  *ForAlpha,      /**< (I)*/
        double  *DomBeta,       /**< (I) domestic mean reversions*/
        double  *ForBeta,       /**< (I) foreign mean reversions*/
        double  *DomFwdRate,    /**< (I)dom rates used in Christian's A factors*/
        double  *ForFwdRate,    /**< (I)for rates used in Christian's A factors*/
        int     DomNbFac,       /**< (I)Nb factors for Dom yield curve*/
        int     ForNbFac,       /**< (I)Nb factors for For yield curve*/
        double  *RhoFxDom,      /**< (I)Corr. between fx and Dom factors*/
        double  *RhoFxFor,      /**< (I)Corr between fx and For factors*/
        double  *RhoDomFor,     /**< (I)Corr between Dom and For factors*/
        double  *RhoDom,        /**< (I)Corr between Dom factors  */
        double  *RhoFor,        /**< (I)Corr between For fcators*/
        double  *SpotFxVol,     /**< (I)Inst fx vol           */
        long    *VolExpDate,    /**< (I)Expiries of CompositeVol */
        int     NbVol,          /**< (I)Nb of output CompositeVols*/
        long    *TPDate,        /**< (I)Dates on time line for vol integration*/
        int     NbTP)           /**< (I)Nb of time points*/
{

    int
        NbNoFxInt = 0,    /*Nb of integrals not containing spot fx vol*/
        NbWithFxInt = 0,  /*Nb of integrals due to covar fx/Dom,fx/For*/
        status=FAILURE;

    int
        i,k,m,l,index,
        CurrentNbInt,   /* used to index integrlas not involving Fx*/
        found;      /* tests whether the Vol expiry date is on the time line*/


    double
            *NoFxInt = NULL,   /*Integrals not containing spot fxvol*/
            *WithFxInt = NULL, /*Integrals due to covar Fx/Dom and Fx/For*/
            *DomA = NULL,
            *ForA = NULL,
            *Length = NULL;
    double
            VarFx;  /* used to calulate the integral due to Spot Fx var.*/
    char
            ErrorMsg[MAXBUFF];
    




    /* Quick checks first */
    if (DomNbFac != 1 && DomNbFac != 2 && DomNbFac != 3)
     
    {
        DR_Error("Invalid Number of domestic yield curve factors "
                "in Hyb3_MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if (ForNbFac != 1 && ForNbFac != 2 && ForNbFac != 3)
     
    {
        DR_Error("Invalid Number of foreign yield curve factors "
                "in Hyb3_MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }



    if( NbTP < 2)
    {
        DR_Error ("Invalid number of time points in Hyb3_MultiFac_FxVol: "
                    "should be at least 2\n");
        return (status);
    }

    if( NbVol < 1)
    {
        DR_Error ("Invalid number of Vol points in MultiFac_fxVol: "
                    "should be at least 1");
        return (status);
    }



    if (CompositeVol == NULL ||
        DomSpotVol == NULL ||
        ForSpotVol == NULL ||
        DomAlpha == NULL ||
        ForAlpha == NULL ||
        DomBeta == NULL ||
        ForBeta == NULL ||
        DomFwdRate == NULL ||
        ForFwdRate == NULL ||
        RhoFxDom == NULL ||
        RhoFxFor == NULL ||
        RhoDomFor == NULL ||
        RhoDom == NULL ||
        RhoFor == NULL ||
        SpotFxVol == NULL ||
        VolExpDate == NULL ||
        TPDate == NULL)

    {
        DR_Error ("invalid pointer inputs to MultiFx_FxVol\n");
        return (status);
    }


    NbNoFxInt = (DomNbFac + ForNbFac) + 
                (DomNbFac + ForNbFac - 1) * (DomNbFac + ForNbFac) / 2;

    NbWithFxInt = (DomNbFac + ForNbFac);


    NoFxInt = (double*) DR_Array (DOUBLE,0,NbNoFxInt - 1);
    WithFxInt = (double*) DR_Array (DOUBLE,0,NbWithFxInt -1);
    DomA = (double*) DR_Array (DOUBLE,0,DomNbFac - 1);
    ForA = (double*) DR_Array (DOUBLE,0,ForNbFac - 1);
    Length=(double*) DR_Array (DOUBLE,0,NbTP - 2);


    if( NoFxInt == NULL ||
        WithFxInt == NULL ||
        DomA == NULL ||
        ForA == NULL ||
        Length == NULL)

    {        
        goto RETURN;
    }

    for (i = 0; i < NbTP - 1; i++)
    {   
         if ( TPDate[i+1] <= TPDate[i])
        {
            sprintf("ErrorMsg","TPDate[%d] and TPDate[%d] in input time line"
                                "are not in a strictly"
                                " ascending order\n",i,i+1);
            DR_Error ("ErrorMsg");
            return (status);
        }
        Length[i] = Daysact(TPDate[i],TPDate[i+1]) / 365.;
    }

    for (k = 0; k < NbVol; k++)
    {
        
        /* first make sure VolExpdate[k] is on time line*/

        found = 0;

        for (index = 0; index < NbTP; index++)
        {   
            found = (TPDate[index] == VolExpDate[k]);

            if(found)
            {                   
                break;
            }
        }/*  for index */

        if (!found)
        {
            sprintf("ErrorMsg", "Vol Exp date #%ld is not"
                                " on the input time line\n",VolExpDate[k]);
            DR_Error("ErrorMsg");
            goto RETURN;
        }
        

        if (VolExpDate[k] == TPDate[0])
        {
            CompositeVol[k] = 0.0;
            continue;
        }
            

        /* Some initialisations */

       for(i = 0; i < DomNbFac; i++)
        {
           DomA[i] = 0.;

        }/* for i*/
        
        for(i = 0; i < ForNbFac; i++)
        {
            ForA[i] = 0.;

        }/* for i*/

        for(i = 0; i < NbNoFxInt; i++)
        {
            NoFxInt[i] = 0.;
        }

        for(i = 0; i < NbWithFxInt; i++)
        {
            WithFxInt[i] = 0.;
        }
        
        VarFx = 0.;

        /* initialisations are now complete*/
        
        for(i = index - 1; i >= 0; i--)
        {
            CurrentNbInt = -1;
            
            /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-DomBeta[l] * Length[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(DomBeta[l],Length[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0; l < ForNbFac; l++)
            {
                ForA[l] *= exp(-ForBeta[l] * Length[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(ForBeta[l],Length[i]);

            }/* for l */

            
            /* variance of domestic factors*/
            
            for (l = 0; l < DomNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += DomAlpha[l] * DomAlpha[l] * 
                                    DomSpotVol[i] * DomSpotVol[i]
                                    * DomA[l] * DomA[l] * Length[i];
            }
            
            /* covariance between domestic factors*/

            if (DomNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* DomAlpha[0] * DomAlpha[1] * 
                                             DomA[0] * DomA[1]
                                            * DomSpotVol[i] * DomSpotVol[i] *
                                            RhoDom[0]*Length[i];

            }/* end of case: DomNbFac > 1 */
            
            if (DomNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * DomAlpha[0] * DomAlpha[2] 
                                        * DomA[0]*DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] 
                                        * RhoDom[1] * Length[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * DomAlpha[1] * DomAlpha[2] 
                                        * DomA[1] * DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] 
                                        * RhoDom[2] * Length[i];

            }/* end of case DomNbFac > 2 */

            /* variance of foreign factors*/
            
            for (l = 0; l < ForNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += ForAlpha[l] * ForAlpha[l] 
                                        * ForSpotVol[i]*ForSpotVol[i]
                                        * ForA[l] * ForA[l]*Length[i];

            } /* for l */

            /* covariance between foreign factors*/

            if (ForNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* ForAlpha[0] * ForAlpha[1]
                                        * ForA[0] * ForA[1]
                                        * ForSpotVol[i] * ForSpotVol[i] 
                                        * RhoFor[0]*Length[i];
            }
            
            if (ForNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * ForAlpha[0] * ForAlpha[2] 
                                        * ForA[0]*ForA[2]
                                        * ForSpotVol[i] * ForSpotVol[i] 
                                        * RhoFor[1] * Length[i];
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * ForAlpha[1] * ForAlpha[2] 
                                        * ForA[1] * ForA[2]
                                        * ForSpotVol[i] * ForSpotVol[i] 
                                        * RhoFor[2] * Length[i];
            }
            
            /* covariance between domestic and foreign factors */

            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {
                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt] -= 2.0 * DomAlpha[l] * DomA[l] *
                                            DomSpotVol[i] *
                                            ForAlpha[m] * ForA[m] * 
                                            ForSpotVol[i] * 
                                            RhoDomFor[m + l * ForNbFac] 
                                            * Length[i];
                }/* for m */

            }/* for l */
            
            /*covariance between Fx and domestic*/

            for (l = 0; l < DomNbFac; l++)
            {
                WithFxInt[l] += 2.0 * SpotFxVol[i]* DomSpotVol[i]
                                * DomAlpha[l] * DomA[l] * 
                                RhoFxDom[l]*Length[i];

            }/* for l*/


            /*covariance between fx and for*/

            for (l = 0; l < ForNbFac; l++)
            {
                WithFxInt[l+DomNbFac] -= 2.0 * SpotFxVol[i]* ForSpotVol[i]
                                        * ForAlpha[l] * ForA[l] 
                                        * RhoFxFor[l]*Length[i];
            }/* for l*/

            VarFx += SpotFxVol[i] * SpotFxVol[i] * Length[i];

        }   /* for i*/

        
        CompositeVol[k]=0.;

        for(i = 0; i < NbNoFxInt; i++)
        {
            CompositeVol[k] += NoFxInt[i];
        }

        for(i = 0; i < NbWithFxInt; i++)
        {
            CompositeVol[k] += WithFxInt[i];
        }

        CompositeVol[k] += VarFx;

        if( CompositeVol[k] < TINY)
        {
            sprintf(ErrorMsg, "Problem in calulating Forward FX Vol at"
                    "Expiry #%ld\n",VolExpDate[k]);

            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* Recall that TPDate[0]=FxBaseDate and that*/
        /* TPDate[0] < VolExpDate[k] because of earlier check
           at beginning of kVol loop, so that we can divide!*/
        
        CompositeVol[k] /= (Daysact(TPDate[0],VolExpDate[k]) / 365.);
        CompositeVol[k] = sqrt (CompositeVol[k]);
                    
    }/* for k*/

    status = SUCCESS;

RETURN:
    Free_DR_Array (NoFxInt,DOUBLE,0,NbNoFxInt - 1);
    Free_DR_Array (WithFxInt,DOUBLE,0,NbWithFxInt - 1);
    Free_DR_Array (DomA,DOUBLE,0,DomNbFac - 1);
    Free_DR_Array (ForA,DOUBLE,0,ForNbFac - 1);
    Free_DR_Array (Length,DOUBLE,0,NbTP - 2);

    return (status);


}/* Multifac_FxVol */





/********* Hyb3_For_MultiFac_FxVol **********************************/
/** 1-Assumes and tests whether Voldate is on the input time line.
    2-Assumes TPDate[0] = FxBaseDate
    3- the correlation inputs are with respect to the exponential factors
    Calculates the Implied volatility of the forward fx (as defined by the expiry date), 
    Vol is integrated between 0 and Revise date.
 *************************************************************/

int Hyb3_For_MultiFac_FxVol(
        double  *CompositeVol,  /*< *(O) Composite FX Vol for FX Expiries*/
        double  *DomSpotVol,    /*< *(I)Instant.dom int.rate vol*/
        double  *ForSpotVol,    /*< *(I)Instant.for.int rate vol*/
        double  *DomAlpha,      /*< *(I)*/
        double  *ForAlpha,      /*< *(I)*/
        double  *DomBeta,       /*< *(I) domestic mean reversions*/
        double  *ForBeta,       /*< *(I) foreign mean reversions*/
        double  *DomFwdRate,    /*< *(I)dom rates used in Christian's A factors*/
        double  *ForFwdRate,    /*< *(I)for rates used in Christian's A factors*/
        int     DomNbFac,       /*< *(I)Nb factors for Dom yield curve*/
        int     ForNbFac,       /*< *(I)Nb factors for For yield curve*/
        double  *RhoFxDom,      /*< *(I)Corr. between fx and Dom factors*/
        double  *RhoFxFor,      /*< *(I)Corr between fx and For factors*/
        double  *RhoDomFor,     /*< *(I)Corr between Dom and For factors*/
        double  *RhoDom,        /*< *(I)Corr between Dom factors  */
        double  *RhoFor,        /*< *(I)Corr between For fcators*/
        double  *SpotFxVol,     /*< *(I)Inst fx vol           */
        long    VolExpDate,     /*< *(I)Expiries of CompositeVol */
        long    ReviseDate,     /*< *(I)Integration of the spot between 0 and Revise date */
        long    *TPDate,        /*< *(I)Dates on time line for vol integration*/
        int     NbTP)           /*< *(I)Nb of time points*/
                

{

    int
        NbNoFxInt = 0,    /*Nb of integrals not containing spot fx vol*/
        NbWithFxInt = 0,  /*Nb of integrals due to covar fx/Dom,fx/For*/
        status=FAILURE;

    int
        i,m,l,index1, index2,
        CurrentNbInt,   /* used to index integrlas not involving Fx*/
        found1, found2;      /* tests whether the Vol expiry date is on the time line*/


    double
            *NoFxInt = NULL,   /*Integrals not containing spot fxvol*/
            *WithFxInt = NULL, /*Integrals due to covar Fx/Dom and Fx/For*/
            *DomA = NULL,
            *ForA = NULL,
            *Length = NULL;
    double
            VarFx;  /* used to calulate the integral due to Spot Fx var.*/
    char
            ErrorMsg[MAXBUFF];
    

    /* Quick checks first */
    if (DomNbFac != 1 && DomNbFac != 2 && DomNbFac != 3)
     
    {
        DR_Error("Invalid Number of domestic yield curve factors "
                "in Hyb3_MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if (ForNbFac != 1 && ForNbFac != 2 && ForNbFac != 3)
     
    {
        DR_Error("Invalid Number of foreign yield curve factors "
                "in Hyb3_MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    /*Check that revise date is before expiry date */
    if(ReviseDate>VolExpDate)
        goto RETURN;

    if( NbTP < 2)
    {
        DR_Error ("Invalid number of time points in Hyb3_MultiFac_FxVol: "
                    "should be at least 2\n");
        return (status);
    }


/*
    if (
        DomSpotVol == NULL ||
        ForSpotVol == NULL ||
        DomAlpha == NULL ||
        ForAlpha == NULL ||
        DomBeta == NULL ||
        ForBeta == NULL ||
        DomFwdRate == NULL ||
        ForFwdRate == NULL ||
        RhoFxDom == NULL ||
        RhoFxFor == NULL ||
        RhoDomFor == NULL ||
        RhoDom == NULL ||
        RhoFor == NULL ||
        SpotFxVol == NULL ||
        TPDate == NULL)

    {
        DR_Error ("invalid pointer inputs to MultiFx_FxVol\n");
        return (status);
    }

*/
    NbNoFxInt = (DomNbFac + ForNbFac) + 
                (DomNbFac + ForNbFac - 1) * (DomNbFac + ForNbFac) / 2;

    NbWithFxInt = (DomNbFac + ForNbFac);


    NoFxInt = (double*) DR_Array (DOUBLE,0,NbNoFxInt - 1);
    WithFxInt = (double*) DR_Array (DOUBLE,0,NbWithFxInt -1);
    DomA = (double*) DR_Array (DOUBLE,0,DomNbFac - 1);
    ForA = (double*) DR_Array (DOUBLE,0,ForNbFac - 1);
    Length=(double*) DR_Array (DOUBLE,0,NbTP - 2);


    if( NoFxInt == NULL ||
        WithFxInt == NULL ||
        DomA == NULL ||
        ForA == NULL ||
        Length == NULL)

    {        
        goto RETURN;
    }

    for (i = 0; i < NbTP - 1; i++)
    {   
         if ( TPDate[i+1] <= TPDate[i])
        {
            sprintf("ErrorMsg","TPDate[%d] and TPDate[%d] in input time line"
                                "are not in a strictly"
                                " ascending order\n",i,i+1);
            DR_Error ("ErrorMsg");
            return (status);
        }
        Length[i] = Daysact(TPDate[i],TPDate[i+1]) / 365.;
    }

    
        
        /* first make sure VolExpdate and Revise date are on time line*/

        found1 = 0;
        found2 = 0;

        for (index1 = 0; index1 < NbTP; index1++)
        {   
            found1 = (TPDate[index1] == ReviseDate);

            if(found1)
            {                   
                break;
            }
        }/*  for index1 */

        for (index2 = 0; index2 < NbTP; index2++)
        {   
            found2 = (TPDate[index2] == VolExpDate);

            if(found2)
            {                   
                break;
            }
        }/*  for index1 */

        if ((!found1)||(!found2))
        {
            sprintf("ErrorMsg", "Vol Exp date #%ld is not"
                                " on the input time line\n",(VolExpDate));
            DR_Error("ErrorMsg");
            goto RETURN;
        }
        

        if ((VolExpDate == TPDate[0])||(ReviseDate == TPDate[0]))
        {
            (*CompositeVol) = 0.0;
        }
            

        /* Some initialisations */

       for(i = 0; i < DomNbFac; i++)
        {
           DomA[i] = 0.;

        }/* for i*/
        
        for(i = 0; i < ForNbFac; i++)
        {
            ForA[i] = 0.;

        }/* for i*/

        for(i = 0; i < NbNoFxInt; i++)
        {
            NoFxInt[i] = 0.;
        }

        for(i = 0; i < NbWithFxInt; i++)
        {
            WithFxInt[i] = 0.;
        }
        
        VarFx = 0.;

        /* initialisations are now complete*/
        /* Calculate DomA(TRevise, VolExpDate) and ForA(TRevise, VolExpDate) factors */
        
        for(i = index2 - 1; i >= index1; i--)
        {
            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-DomBeta[l] * Length[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(DomBeta[l],Length[i]);
            }/* for l */

            for (l = 0; l < ForNbFac; l++)
            {
                ForA[l] *= exp(-ForBeta[l] * Length[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(ForBeta[l],Length[i]);

            }/* for l */
        }


        for(i = index1 - 1; i >= 0; i--)
        {
            CurrentNbInt = -1;
            
            /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-DomBeta[l] * Length[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(DomBeta[l],Length[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0; l < ForNbFac; l++)
            {
                ForA[l] *= exp(-ForBeta[l] * Length[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(ForBeta[l],Length[i]);

            }/* for l */

            
            /* variance of domestic factors*/
            
            for (l = 0; l < DomNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += DomAlpha[l] * DomAlpha[l] * 
                                    DomSpotVol[i] * DomSpotVol[i]
                                    * DomA[l] * DomA[l] * Length[i];
            }
            
            /* covariance between domestic factors*/

            if (DomNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* DomAlpha[0] * DomAlpha[1] * 
                                             DomA[0] * DomA[1]
                                            * DomSpotVol[i] * DomSpotVol[i] *
                                            RhoDom[0]*Length[i];

            }/* end of case: DomNbFac > 1 */
            
            if (DomNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * DomAlpha[0] * DomAlpha[2] 
                                        * DomA[0]*DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] 
                                        * RhoDom[1] * Length[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * DomAlpha[1] * DomAlpha[2] 
                                        * DomA[1] * DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] 
                                        * RhoDom[2] * Length[i];

            }/* end of case DomNbFac > 2 */

            /* variance of foreign factors*/
            
            for (l = 0; l < ForNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += ForAlpha[l] * ForAlpha[l] 
                                        * ForSpotVol[i]*ForSpotVol[i]
                                        * ForA[l] * ForA[l]*Length[i];

            } /* for l */

            /* covariance between foreign factors*/

            if (ForNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* ForAlpha[0] * ForAlpha[1]
                                        * ForA[0] * ForA[1]
                                        * ForSpotVol[i] * ForSpotVol[i] 
                                        * RhoFor[0]*Length[i];
            }
            
            if (ForNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * ForAlpha[0] * ForAlpha[2] 
                                        * ForA[0]*ForA[2]
                                        * ForSpotVol[i] * ForSpotVol[i] 
                                        * RhoFor[1] * Length[i];
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * ForAlpha[1] * ForAlpha[2] 
                                        * ForA[1] * ForA[2]
                                        * ForSpotVol[i] * ForSpotVol[i] 
                                        * RhoFor[2] * Length[i];
            }
            
            /* covariance between domestic and foreign factors */

            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {
                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt] -= 2.0 * DomAlpha[l] * DomA[l] *
                                            DomSpotVol[i] *
                                            ForAlpha[m] * ForA[m] * 
                                            ForSpotVol[i] * 
                                            RhoDomFor[m + l * ForNbFac] 
                                            * Length[i];
                }/* for m */

            }/* for l */
            
            /*covariance between Fx and domestic*/

            for (l = 0; l < DomNbFac; l++)
            {
                WithFxInt[l] += 2.0 * SpotFxVol[i]* DomSpotVol[i]
                                * DomAlpha[l] * DomA[l] * 
                                RhoFxDom[l]*Length[i];

            }/* for l*/


            /*covariance between fx and for*/

            for (l = 0; l < ForNbFac; l++)
            {
                WithFxInt[l+DomNbFac] -= 2.0 * SpotFxVol[i]* ForSpotVol[i]
                                        * ForAlpha[l] * ForA[l] 
                                        * RhoFxFor[l]*Length[i];
            }/* for l*/

            VarFx += SpotFxVol[i] * SpotFxVol[i] * Length[i];

        }   /* for i*/

        
        (*CompositeVol)=0.;

        for(i = 0; i < NbNoFxInt; i++)
        {
            (*CompositeVol) += NoFxInt[i];
        }

        for(i = 0; i < NbWithFxInt; i++)
        {
            (*CompositeVol) += WithFxInt[i];
        }

        (*CompositeVol) += VarFx;

        if( (*CompositeVol) < TINY)
        {
            sprintf(ErrorMsg, "Problem in calulating Forward FX Vol at"
                    "Expiry #%ld\n",VolExpDate);

            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* Recall that TPDate[0]=FxBaseDate and that*/
        /* TPDate[0] < VolExpDate[k] because of earlier check
           at beginning of kVol loop, so that we can divide!*/
        
        (*CompositeVol) /= (Daysact(TPDate[0],ReviseDate) / 365.);
        (*CompositeVol) = sqrt (*CompositeVol);
                    

    status = SUCCESS;

RETURN:
    Free_DR_Array (NoFxInt,DOUBLE,0,NbNoFxInt - 1);
    Free_DR_Array (WithFxInt,DOUBLE,0,NbWithFxInt - 1);
    Free_DR_Array (DomA,DOUBLE,0,DomNbFac - 1);
    Free_DR_Array (ForA,DOUBLE,0,ForNbFac - 1);
    Free_DR_Array (Length,DOUBLE,0,NbTP - 2);

    return (status);
}/* FOr_Multifac_FxVol */






/************ Hyb3_Get_FxVol ******************************/
/**Outputs composite fx vols for the specified input
   VolExpDates.
   This function bootstraps IR vols using NORMAL smile
   parameters.
  
   ASSUMES:
    1-The input SpotFxVol array to be of size 
      (fx_data.NbVol + fx_data.NbInpSpotVol);
    
    2-SpotFxVols are BACKWARD looking
  
    3-The IMPLICIT expiries/dates along side SpotFxVols are:
  
       -fx_data.VolDate (i.e input composite expiry dates) for the first NbVol
        elements.
       - fx_data.InpSpotVolDate (i.e input spotvol dates) for the remaining 
        elements
      
 ******************************************************/


int Hyb3_Get_FxVol (
            T_CURVE      *DomCurve,         /**< (I)*/
            T_CURVE      *ForCurve,         /**< (I)*/
            MKTVOL_DATA  *DomMktVol,        /**< (I)*/
            MKTVOL_DATA  *ForMktVol,        /**< (I)*/
            FX_DATA      *fx_data,          /**< (I)*/
          
            /* Ppy for creating time line */
            int         Ppy,                /**< (I)*/

            /* Fx details */
            double      *SpotFxVol,        /**< (I)*/
            int         NbFxVol,           /**< (I)Nb of Vols we are outputting*/
            long        *VolExpDates,      /**< (I)composite vols expiries*/
            double      *CompositeVol)     /**< (O)*/

{

    int     
            i,k,status = FAILURE;
    double
            
            *Discrete_SpotVol = NULL;
    long    
            *TPDate = NULL,
            *FxSpotVolDates = NULL;
    double  
            *DomFwd = NULL,
            *ForFwd = NULL;

   
    double  *DomSpotVol = NULL,
            *ForSpotVol = NULL;
    double  
            alpha;  /* this is Alpha[0] */
    int     
            NbTP = 0,
            TotNbFxSpotVols = 0;
            

    double 
           test_mtx [3][3]; /* triangualtion matrices*/
    double 
            test_rho[3];   /* to test validity of input corr.structure*/  
           
    /* Some basic checks first */
    

    
    test_rho[0] = fx_data->Rho[0][0];
    test_rho[1] = fx_data->Rho[1][0];
    test_rho[2] = fx_data->Rho[2][0];
    
    if (Triangulation (test_mtx,3,test_rho) == FAILURE)
    {
        DR_Error ("Invalid input correlation structure"
                    "between Fx rate and swap rates in Hyb3_Get_FxVol\n");
        return (status);
    }

    if (DomCurve == NULL ||
        ForCurve == NULL ||
        DomMktVol == NULL ||
        ForMktVol == NULL ||
        fx_data == NULL ||
        SpotFxVol == NULL ||
        VolExpDates == NULL ||
        CompositeVol == NULL)
    {
        DR_Error("Invalid pointer inputs to Hyb3_Get_FxVol\n");
        return (status);
    }
    
    if((fx_data->NbVol < 0) || (fx_data->NbInpSpotVol < 0))
    {
        DR_Error("Nb input vols must be >= 0\n");
        goto RETURN;
    }

    if (fx_data->NbInpSpotVol > 0)
    {
        if(fx_data->InpSpotVolDate[0] <= fx_data->ValueDate)
        {
            DR_Error("First spotvol date must be > fx value date\n");
            goto RETURN;
        }
    }
    if ((fx_data->NbVol > 0) && (fx_data->NbInpSpotVol > 0))
    {
        if(fx_data->VolDate[fx_data->NbVol] >= fx_data->InpSpotVolDate[0])
        {
            DR_Error("First internla composite voldate must be "
                        " stricylt before first spotvol date\n");
            goto RETURN;
        }
    }
    
    if (fx_data->NbVol > 0)
    {
        if (fx_data->VolDate[0] != fx_data->Today) goto RETURN;
    }


    TotNbFxSpotVols = fx_data->NbVol + fx_data->NbInpSpotVol;
    if (TotNbFxSpotVols < 1)
    {
        DR_Error("Total number of spot fxvols is <=0\n");
        goto RETURN;
    }

    FxSpotVolDates = (long*) DR_Array(LONG,0,TotNbFxSpotVols - 1);
    if (FxSpotVolDates == NULL) goto RETURN;
    
    for (k = 0 ; k < fx_data->NbVol ; k++)
    {
        FxSpotVolDates[k] = fx_data->VolDate[k+1];
    }
    
    for (k = fx_data->NbVol ; k < TotNbFxSpotVols; k++)
    {
        FxSpotVolDates[k] = fx_data->InpSpotVolDate[k - fx_data->NbVol];
    }

    if (FxSpotVolDates[0] <= fx_data->Today) goto RETURN;
        
    if (TotNbFxSpotVols > 1)
    {   
        long    l;

        for(l = 0 ; l < TotNbFxSpotVols - 1; l++)
        {
            if (FxSpotVolDates[l] >= FxSpotVolDates[l+1]) goto RETURN;
        }
    }

    /* Create time line */
    if (MergeDateLists(TotNbFxSpotVols,
                       FxSpotVolDates,
                       &(NbTP),
                       &(TPDate)) == FAILURE)
    {
        DR_Error("MergeDateLists failed for TotNbFxVols in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    if (MergeDateLists(
                        DomMktVol->NbVol,
                        DomMktVol->VolDate,
                        &(NbTP),
                        &(TPDate)) == FAILURE)
    {
        DR_Error("MergeDateLists failed in Hyb3_Get_FxVol for Dom vols\n");
        goto RETURN;
    }

    if (MergeDateLists(
                        ForMktVol->NbVol,
                        ForMktVol->VolDate,
                        &NbTP,
                        &TPDate) == FAILURE)
    {
        DR_Error("MergeDateLists failed in Hyb3_Get_FxVol for For vols\n");
        goto RETURN;
    }
        
    if (MergeDateLists(
                        NbFxVol,
                        VolExpDates,
                        &NbTP,
                        &TPDate) == FAILURE)
    {
        DR_Error("MergeDateLists failed in Hyb3_Get_FxVol for input fx vols\n");
        goto RETURN;
    }

    /* Add Fx base to the TPDate list */
    if (MergeDateLists(
                        1,
                        &(fx_data->Today),
                        &NbTP,
                        &TPDate) == FAILURE)
    {
        DR_Error("MergedDateLists failed in Hyb3_Get_FxVol for fx base date\n");
        goto RETURN;
    }
       
    /* recall that DrExtendAmerSch will remove all dates STRICTLY
      before fx_data->Today*/
    if (DrExtendAmerSch(fx_data->Today,
                        Ppy,
                        &NbTP,
                        &TPDate,
                        NULL,
                        'L',
                        NULL,
                        NULL,
                        NULL,
                        NULL,
                        NULL) == FAILURE)
    {
        DR_Error ("DrExtendAmerSch failed in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    if (TPDate[0] != fx_data->Today)
    {
        DR_Error ("TPDate[0] is not fx_data.Today (Hyb3_Get_FxVol)\n");
        goto RETURN;
    }

    /* bootstrap IR spot vol using NORMAL parameters */
     if (Hyb3_SpotVol(   DomMktVol->Aweight,
                    DomMktVol->BaseDate,
                    DomMktVol->NbVol,
                    DomMktVol->VolDate,
                    DomMktVol->Vol,
                    DomMktVol->VolUsed,
                    DomMktVol->Freq,
                    DomMktVol->DCC,
                    DomMktVol->SwapSt,
                    DomMktVol->SwapMat,
                    0., /* "Normal" QLeft    */                         
                    0., /* "Normal" QRight   */                   
                    0., /* No forward shift  */                    
                    1,  /* 1 Factor          */
                    DomMktVol->Alpha,
                    DomMktVol->Beta,
                    DomMktVol->Rho,
                    DomMktVol->SkipFlag,
                    DomMktVol->CalibFlag,
                    DomCurve) != SUCCESS)
    {
        DR_Error("Hyb3_SpotVol for domestic failed in Get_SpotFxVol\n");
        goto RETURN;
    }

    if (Hyb3_SpotVol(    ForMktVol->Aweight,
                    ForMktVol->BaseDate,
                    ForMktVol->NbVol,
                    ForMktVol->VolDate,
                    ForMktVol->Vol,
                    ForMktVol->VolUsed,
                    ForMktVol->Freq,
                    ForMktVol->DCC,
                    ForMktVol->SwapSt,
                    ForMktVol->SwapMat,
                    0.,  /* "Normal" QLeft    */                  
                    0.,  /* "Normal" QRight   */                  
                    0.,  /* No forward shift  */                  
                    1,   /* 1 Factor          */
                    ForMktVol->Alpha,
                    ForMktVol->Beta,
                    ForMktVol->Rho,
                    ForMktVol->SkipFlag,
                    ForMktVol->CalibFlag,
                    ForCurve) != SUCCESS)
    {
        DR_Error("Hyb3_SpotVol for foreign failed in Get_SpotFxVol\n");
        goto RETURN;
    }
        
    /* Discretise spot vol onto time line */
    DomSpotVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    
    if (DomSpotVol == NULL)
    {
        DR_Error ("Invalid memory allocation for domestic"
                    "vol in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* The spot Vol needed for the rest of the function
     * is the Spot Vol of the EXPONENTIAL factors,as opposed to     
     * the orthogonal ones produced by triangulation,
     * therefore use only Aweight[0]                    */
    
    if (Hyb3_DiscretiseInput(DomSpotVol,
                        DomMktVol->Aweight[0],
                        DomMktVol->CalibFlag,
                        DomMktVol->NbVol,
                        DomMktVol->VolDate,
                        NbTP,
                        TPDate) == FAILURE)
    {
        DR_Error("Hyb3_DiscretiseInput failed for domestic"
                    "Spot vol in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* spot vol input to Hyb3_MultiFac_FxVol should be "weight free",
     * so rescale.                                           */
    
    alpha = DomMktVol->Alpha[0];

    /* Hyb3_Param_Check ensures that Alpha[0] non zero*/
    for (i = 0; i < NbTP - 1;i++)
    {
        DomSpotVol[i] /= alpha; 
    }

    ForSpotVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    
    if (ForSpotVol == NULL)
    {
        DR_Error ("Invalid memory allocation for foreign vol"
                    "in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
   
    if (Hyb3_DiscretiseInput(ForSpotVol,
                        ForMktVol->Aweight[0],
                        ForMktVol->CalibFlag,
                        ForMktVol->NbVol,
                        ForMktVol->VolDate,
                        NbTP,
                        TPDate) == FAILURE)
    {
        DR_Error("Hyb3_DiscretiseInput failed for foreign"
                    "Spot vol in Get_Spot_fxVol\n");
        goto RETURN;
    }
                        
    alpha = ForMktVol->Alpha[0];
    for (i = 0; i < NbTP - 1; i++)
    {
        ForSpotVol[i] /= alpha;
    }
    
    Discrete_SpotVol = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    if (Discrete_SpotVol == NULL)
    {
        DR_Error("memory allocation failure in Hyb3_Get_FxVol\n");
        goto RETURN;
    }

    
    if (Hyb3_DiscretiseInput(Discrete_SpotVol,
                        SpotFxVol,
                        TRUE,   /*CalibFlag*/
                        TotNbFxSpotVols,
                        FxSpotVolDates,
                        NbTP,
                        TPDate) == FAILURE)
    {
        DR_Error("Hyb3_DiscretiseInput failed for fx"
                    "Spot vol in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* interpolate relevant forward rates*/
    DomFwd = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    ForFwd = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    
    if (DomFwd == NULL || ForFwd == NULL)
    {
        DR_Error ("Failed to allocate memory for forward rates"
                    "in Hyb3_Get_FxVol\n");
        goto RETURN;
    }

    if (Hyb3_SimpleFwdCurve(DomFwd,
                        DomCurve,
                        TPDate,
                        NbTP) == FAILURE)
    {
        DR_Error ("Hyb3_SimpleFwdCurve failed for domestic"
                    "curve in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
  

    if (Hyb3_SimpleFwdCurve( ForFwd,
                        ForCurve,
                        TPDate,
                        NbTP) == FAILURE)
    {
        DR_Error ("Hyb3_SimpleFwdCurve failed for foreign"
                    "curve in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* we are now ready to call MultiFac_FXVol*/

    if (Hyb3_MultiFac_FxVol(
                        CompositeVol,
                        DomSpotVol,
                        ForSpotVol,
                        DomMktVol->Alpha,
                        ForMktVol->Alpha,
                        DomMktVol->Beta,
                        ForMktVol->Beta,
                        DomFwd,
                        ForFwd,
                        1,/* 1 factor*/
                        1,/* 1 factor*/
                        &(fx_data->Rho[2][0]),/* Rho Fx/Dom*/
                        &(fx_data->Rho[1][0]),/* Rho Fx/For*/
                        &(fx_data->Rho[0][0]),/* Rho Dom/For */
                        DomMktVol->Rho,
                        ForMktVol->Rho,
                        Discrete_SpotVol,
                        VolExpDates,
                        NbFxVol,
                        TPDate,
                        NbTP) == FAILURE)
    {
        DR_Error ("Hyb3_MultiFac_FxVol failed in get_FxVol\n");
        goto RETURN;
    }

    status = SUCCESS;
RETURN: 
    if (NbTP >= 1)
    {
        Free_DR_Array (TPDate,LONG,0,NbTP - 1);
    }
    if (TotNbFxSpotVols >= 1)
    {
        Free_DR_Array (FxSpotVolDates,LONG, 0, TotNbFxSpotVols - 1);
    }

    if (NbTP >= 2)
    {
        Free_DR_Array (DomSpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (ForSpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (Discrete_SpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (DomFwd,DOUBLE,0,NbTP - 2);
        Free_DR_Array (ForFwd,DOUBLE,0,NbTP - 2);
    }
    
    return (status);

}/* Hyb3_Get_FxVol*/


/*********** Hyb3_DiscretiseInput ********************************************/
/**
  
    given NbTP time points, outputs NbTP-1 Vol values.
  
  
    NOTE:1- skipps VolDates < TPDate[0]
         2- tests that each VolDate >=TPDate[0] coincides with one TPDate 
         3- assumes VolDate are in a STRICTLY ascending order.
  
 ***************************************************************************/



int Hyb3_DiscretiseInput(
        double  *DiscretisedVol,        /**< (O)                               */
        double  *Vol,                   /**< (I) Vol to discretise onto TPdates*/
        int     CalibFlag,              /**< (I)                               */
        int     NbVol,                  /**< (I) Nb of input Vol points        */
        long    *VolDate,               /**< (I) Vol dates                     */
        int     NbTP,                   /**< (I) Total Nb of TP                */
        long    *TPDate)                /**< (I) TP dates                      */
{
    int i,j,k,status = FAILURE;
    
    int found;

    int TPSofar,TPBefore;


    if( DiscretisedVol == NULL  ||
        Vol == NULL             ||
        VolDate == NULL         ||
        TPDate == NULL)
    {
        DR_Error ("Invalid pointer inputs to Hyb3_DiscretiseInput\n");
        return (status);
    }
   

    if (NbTP < 1)
    {
        DR_Error("Invalid number of input time points in Hyb3_DiscretiseInput: "
                    "should be at least 1\n");
        
        return(status);
    }
    

    if (CalibFlag  == FALSE)
    {
        for (j = 0; j < (NbTP - 1); j++)
        {
            DiscretisedVol[j] = Vol[0];
        }

        return (SUCCESS);
    }
   

    TPSofar = TPBefore = 0;

    for (i = 0; i < NbVol; i++)
    {   
        /* if VolDate < TPDate[0], skip i */
        if (VolDate[i] < TPDate[0])
            continue;

        /* first make sure Vol is on time line*/

        found = 0;
        for ( ; TPSofar < NbTP; TPSofar++)
        {
            found = (TPDate[TPSofar] == VolDate[i]);

            if (found)
            {
                break;
            }

        }/* for TPSofar */

        if (!found)
        {
            sprintf("ErrorMsg","VolDate %ld is not on the time "
                                "line (Hyb3_DiscretiseInput)\n",
                                VolDate[i]);
            DR_Error("ErrorMsg");
            goto RETURN;
        }
        
        for (k = TPBefore; k < TPSofar;k++)
        {
            DiscretisedVol[k] = Vol[i];

        } /* for k */

        TPBefore = TPSofar;

    } /* for i*/

    if (TPSofar < NbTP - 1)
    {
        /* extend with last vol level*/

        for (j = TPSofar; j < NbTP - 1; j++)
        {
            DiscretisedVol[j] = Vol[NbVol - 1];

        }/* for j*/
    }

    status = SUCCESS;


RETURN:
    return(status);

}/*Discretise Input*/





/***** Hyb3_Get_TreeFxSpotVol *****************************************************/
/** CAREFUL: The Input IR vols are assumed to be in an "Aweight" form 
   i.e "alpha*spotvol". 
   
   1-Assumes that TPDate ranges from TPDate[-1] to TPDate[NbTP - 1]
     The last meaningful TreeSpotVol is TreeSpotVol[NbTP - 2], 
     The first meaningful TreeSpotVol is TreeSpotVol[0].
     TreeSpotVol is then extended by:
     TreeSpotVol[NbTP - 1] ,which is set to TreeSpotVol[NbTP - 2],and
     TreeSpotVol[-1], which is set to TreeSpotVol[0];
  
   2-The array of InpCompVol is offset by one, i.e, the first
     effective expiry is InpVompVol[1], the last one is InpCompVol[NbCompVol]
  
   3-The array of InpSpotVol is NOT offset by one, i.e, 
     the first spot vol is InpSpotVol[0], and the last one is 
     InpSpotVol[NbInpSpotVol - 1]  
  
 *****************************************************************************/

int Hyb3_Get_TreeFxSpotVol(double    *TreeSpotVol, /**< (O)Bootst'd spotvol on tree    */
                      long      *TPDate,      /**< (I)tree time pts               */
                      long      NbTP,         /**< (I)Last TP is TPDate[NbTP-1]   */
                      double    *InpCompVol,  /**< (I)input composite vols        */
                      long      *InpCompVDate,/**< (I)option expiry dates         */
                      int       NbCompVol,    /**< (I)Number of composite vols    */
                      double    *InpSpotVol,  /**< (I)Input Spot Fx vols          */
                      long      *InpSpotVDate,/**< (I)Input Hyb3_SpotVol Dates         */
                      int       NbInpSpotVol, /**< (I)Nb of Input spot vols       */
                      double    *DomSpotVol,  /**< (I)dom ir (tree) Aweight       */
                      double    *ForSpotVol,  /**< (I)for ir (tree) Aweight       */
                      double    *DomAlpha,    /**< (I)dom factor weights          */ 
                      double    *ForAlpha,    /**< (I)for factor weights          */
                      double    *DomFwdRate,  /**< (I)dom (tree) fwd rates        */
                      double    *ForFwdRate,  /**< (I)for (tree) fwd rates        */
                      double    *RhoFxDom,    /**< (I)                            */
                      double    *RhoFxFor,    /**< (I)                            */
                      double    *RhoDomFor,   /**< (I)                            */
                      double    *RhoDom,      /**< (I)                            */
                      double    *RhoFor,      /**< (I)                            */
                      int       FxCutOffFlag, /**< (I)                            */
                      int       FxCutOffLast, /**< (I)                            */
                      double    FxCutOffLevel,/**< (I)                            */
                      double    *DomBeta,     /**< (I)                            */
                      double    *ForBeta,     /**< (I)                            */
                      int       DomNbFac,     /**< (I)                            */
                      int       ForNbFac,     /**< (I)                            */
                      long      Today,        /**< (I)                            */
                      int       FxBootStrapMode,
                      double    *A_Fx,        /**<  (I) FX smile param, [0,NbTP-2]*/
                      double    *B_Fx,        /**<  (I) FX smile param, [0,NbTP-2]*/
                      double    *C_Fx)        /**<  (I)                           */
{


    int status = FAILURE;
    long    PseudoCompVDate[MAXNBDATE+1];/*Stores interpolated expiry dates  */
    double  PseudoCompVol[MAXNBDATE+1];  /*stores interpolated composite vols*/
    long    *RawSpotVDate   = NULL;      /*spot vol dates before extension   */
    double  *RawSpotVol     = NULL;      /*spotvols bef. extension on TPDate */
    int     NbRawSpotVol    = 0;     /*total Nb of spot vols bef. extension  */
    int     NbPseudoCompVol = 0;     /*Nb of valid interpolated compositevols*/


    int     i,k;
    long    TPSofar, TPBefore;
    double  Vol;
    double  *LDomVol = NULL;        /* dom Aweight/alpha*/  
    double  *LForVol = NULL;        /* for Aweight/alpha*/
    
    /* basic checks */
    if (TreeSpotVol == NULL ||
        TPDate      == NULL ||
        DomSpotVol  == NULL ||
        ForSpotVol  == NULL ||
        DomAlpha    == NULL ||
        ForAlpha    == NULL ||
        DomBeta     == NULL ||
        ForBeta     == NULL ||
        DomFwdRate  == NULL ||
        ForFwdRate  == NULL ||
        RhoFxDom    == NULL ||
        RhoFxFor    == NULL ||
        RhoDomFor   == NULL ||
        RhoDom      == NULL ||
        RhoFor      == NULL)
    {
        DR_Error("invalid input pointers to Hyb3_Get_TreeFxSpotVol\n");
        goto RETURN;
    }
    
    
    if (NbTP < 2) goto RETURN;

    /* Do the nil calibration case and exit*/
    if (FxBootStrapMode == FX_CONSTANT_SPOT_VOL)
    {
        for (k = -1 ; k < NbTP ; k++)
        {
            TreeSpotVol[k] = FxCutOffLevel;
        }
        return (SUCCESS);
    }

    if ((NbInpSpotVol < 0) || (NbCompVol < 0)) goto RETURN;
    if ((NbCompVol > 0) && ((InpCompVol == NULL) || (InpCompVDate == NULL)))
    {        
        goto RETURN;
    }


    if ((NbInpSpotVol > 0) && ((InpSpotVol == NULL) || (InpSpotVDate == NULL)))
    {
        goto RETURN;
    }
    
    /* Check that composite vol dates have been read in using Hyb3_Fx_Input_W*/
    /* which offsets vols and dates by 1 and sets VolDate[0] to Today   */
    if (NbCompVol > 0)
    { 
        if(InpCompVDate[0] != Today)
        {
            DR_Error("First Composite Vol Date is not Today\n");
            goto RETURN;
        }
        if(!IS_EQUAL(InpCompVol[0],InpCompVol[1])) goto RETURN;
    }
    
    
    if (TPDate[0] != Today)
    {
        DR_Error("TPDate[0] is not Today \n");
        goto RETURN;
    }
        
    /* Prepare Ir vols */
    if (DomAlpha[0] < TINY || ForAlpha[0] < TINY) goto RETURN;

    LDomVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    LForVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);

    if (LDomVol == NULL || LForVol == NULL) goto RETURN;

    for (k = 0 ; k < NbTP - 1 ; k++)
    {
        LDomVol[k] = DomSpotVol[k] / DomAlpha[0];
        LForVol[k] = ForSpotVol[k] / ForAlpha[0];
    }


    /* Find Pseudo composite vols and expiries */
    /*************************************************/
    /* Input composite vols are offset by 1, whereas */
    /* Pseudocomposite vols are offset by 0          */
    /*************************************************/

    TPSofar = TPBefore = 0;
    NbPseudoCompVol    = 0;

    for (k = 1; k <= NbCompVol ; k++) /*input composite vols are offset by 1*/
    {
        /* Look for the tree time  point  TPSoFar  falling  */ 
        /* immediately before the current volatility input  */         
        TPSofar = NbTP - 1;
        while ((TPDate[TPSofar] > InpCompVDate[k]) && (TPSofar > 0))  
                    TPSofar--;
        
        /* Ignore this point if it is falling onto the previous one */        
        if (TPSofar == TPBefore)
            continue;

        /* Interpolate fx volatility at current node as the input */
        /* volatility points may not be falling on a node.        */    
        dlinterp (TPDate[TPSofar],     &Vol,                     
                  InpCompVDate[k-1] ,  InpCompVDate[k],
                  InpCompVol[k-1]   ,  InpCompVol[k]);
        
        PseudoCompVDate[NbPseudoCompVol] = TPDate[TPSofar];
        PseudoCompVol[NbPseudoCompVol]   = Vol;

        NbPseudoCompVol++;
        TPBefore = TPSofar;
    }

    /* quick checks*/    
    if ((NbPseudoCompVol > NbCompVol) || (NbPseudoCompVol < 0)) goto RETURN;    
    if ((NbInpSpotVol > 0) && (NbPseudoCompVol > 0))
    {
        if(PseudoCompVDate[NbPseudoCompVol-1] >= InpSpotVDate[0]) goto RETURN;
    }
    
    if (NbPseudoCompVol > 0)
    {
        if(PseudoCompVDate[0] <= Today) goto RETURN;
    }

    if (NbPseudoCompVol > 1)
    {
        for (i = 0 ; i < NbPseudoCompVol - 1; i++)
        {
            if (PseudoCompVDate[i] >= PseudoCompVDate[i+1]) goto RETURN;
        }
    }

    NbRawSpotVol = NbPseudoCompVol + NbInpSpotVol;
    if (NbRawSpotVol < 1) 
    {
        DR_Error("No vols to bootstrap (no spot vols and "
                "no composite vols left after interpolation!)\n");
        goto RETURN;
    }
    
    RawSpotVDate    = (long*)   DR_Array(LONG,  0, NbRawSpotVol - 1);
    RawSpotVol      = (double*) DR_Array(DOUBLE,0, NbRawSpotVol - 1);

    if ((RawSpotVDate == NULL) || (RawSpotVol == NULL)) goto RETURN;


    /* fill in vol dates and spot vols */
    for (k = 0 ; k < NbPseudoCompVol; k++)
    {
        RawSpotVDate[k] = PseudoCompVDate[k];
    }
  
    if (NbPseudoCompVol > 0)
    {   
        
        if (Hyb3_MultiFac_Spot_FxVol2(RawSpotVol,
                                LDomVol,
                                LForVol,
                                DomAlpha,
                                ForAlpha,
                                DomFwdRate,
                                ForFwdRate,
                                RhoFxDom,
                                RhoFxFor,
                                RhoDomFor,
                                RhoDom,
                                RhoFor,
                                NbPseudoCompVol,
                                PseudoCompVDate,
                                PseudoCompVDate,
                                PseudoCompVol,
                                FxCutOffFlag,
                                FxCutOffLast,
                                FxCutOffLevel,
                                DomBeta,
                                ForBeta,
                                DomNbFac,
                                ForNbFac,
                                A_Fx,
                                B_Fx,
                                C_Fx,
                                TPDate,
                                NbTP) == FAILURE) goto RETURN;        

    }/* if there are composite vols to bootstrap*/
    
    /* Append inputspot vol dates and vols*/
    for (k = NbPseudoCompVol ; k < NbRawSpotVol ; k++)
    {
        RawSpotVDate[k] = InpSpotVDate[k - NbPseudoCompVol];
        RawSpotVol[k]   = InpSpotVol[k - NbPseudoCompVol];
    }
    

    if(ExtendSpotVol(TreeSpotVol,
                     NbRawSpotVol,
                     RawSpotVDate,
                     RawSpotVol,
                     NbTP,
                     TPDate) == FAILURE) goto RETURN;


    /* extend by one */
    TreeSpotVol[-1] = TreeSpotVol[0];
    status = SUCCESS;

RETURN:
    
    if (NbRawSpotVol >= 1)
    {
        Free_DR_Array(RawSpotVDate,LONG,0,NbRawSpotVol - 1);
        Free_DR_Array(RawSpotVol, DOUBLE,0,NbRawSpotVol - 1);
    }
    if (NbTP >= 2)
    {
        Free_DR_Array(LDomVol, DOUBLE, 0, NbTP - 2);
        Free_DR_Array(LForVol, DOUBLE, 0, NbTP - 2);
    }
    
    return(status);
}



/*****  Hyb3_ForwardFX ************************************************************/
/**
 *      Calculate a deterministic forward FX level as at the FX value date.
 *      The maturity of the forward FX is denoted by T
 *      The function accounts for the different value dates in the 2 discount
 *      zero curves
 */
int    Hyb3_ForwardFX(double  *FwdFX,   /**<  (O) fwd spread value           */
                 long     T,            /**< (I) fwd fx maturity             */
                 double   SpotFX,       /**< (I) Spot level at FX Value Date */
                 long     FXValueDate,  /**< (I) in YYYYMMDD                 */

                 T_CURVE const* fcrv,   /** foreign zero curve               */
                 T_CURVE const* dcrv)   /** domestic zero curve              */
{
    int      status = FAILURE;

    double   ForZero;   /* foreign fwd zero from FXValueDate to T  */
    double   DomZero;   /* domestic fwd zero from FXValueDate to T */
    double   ZeroPriceL;

    /* basic checks */

    if (FwdFX    == NULL) 
        goto RETURN;

    /* compute the foreign fwd zero */

    ForZero    = GetZeroPrice(T, fcrv);
    ZeroPriceL = GetZeroPrice(FXValueDate, fcrv);

    if (ForZero < TINY || ZeroPriceL < TINY)
        goto RETURN;

    if (FXValueDate >= fcrv->ValueDate)
        ForZero /= ZeroPriceL;

    DomZero    = GetZeroPrice(T, dcrv);
    ZeroPriceL = GetZeroPrice(FXValueDate, dcrv);

    if (DomZero < TINY || ZeroPriceL < TINY)
        goto RETURN;

    if (FXValueDate >= dcrv->ValueDate)
        DomZero /= ZeroPriceL;

    *FwdFX = SpotFX * ForZero / DomZero;

    status = SUCCESS;

RETURN:

    return (status);

} /* Hyb3_ForwardFX */


/************ Hyb3_Get_FxVol2 ******************************/
/**Outputs composite fx vols for the specified input
   VolExpDates.
   This function bootstraps IR vols using NORMAL smile
   parameters.
  
    THIS FUNCTION IS DIFFERENT FROM Hyb3_Get_FxVol because:
  
    1-The input SpotFxVol array can be of any size
    
    2-SpotFxVols are FORWARD looking
  
    3-The dates along side SpotFxVols are specified as an input array
  
*/

int Hyb3_Get_FxVol2 (
            T_CURVE      *DomCurve,         /**<(I)*/
            T_CURVE      *ForCurve,         /**<(I)*/
            MKTVOL_DATA  *DomMktVol,        /**<(I)*/
            MKTVOL_DATA  *ForMktVol,        /**<(I)*/
            FX_DATA      *fx_data,          /**<(I)*/
          
            /* Ppy for creating time line */
            int         Ppy,                /**<(I)*/

            /* Fx details */
            long        NbFxSpotVols,
            long       *FxSpotVolDates,
            double     *FxSpotVol,

            int         NbFxVol,           /**<(I)Nb of Vols we are outputting*/
            long        *VolExpDates,      /**<(I)composite vols expiries*/
            double      *CompositeVol)     /**<(O)*/

{

    int     
            i,k,status = FAILURE;
    double
           
            *Discrete_SpotVol = NULL;
    long    
            *TPDate = NULL;
    double  
            *DomFwd = NULL,
            *ForFwd = NULL;

    double  *DomSpotVol = NULL,
            *ForSpotVol = NULL;

    double  
            alpha;  /* this is Alpha[0] */

    int     
            NbTP = 0;
            

    double 
           test_mtx [3][3]; /* triangualtion matrices*/
    double 
            test_rho[3];   /* to test validity of input corr.structure*/  
           
    /* Some basic checks first */
    

    
    test_rho[0] = fx_data->Rho[0][0];
    test_rho[1] = fx_data->Rho[1][0];
    test_rho[2] = fx_data->Rho[2][0];
    
    if (Triangulation (test_mtx,3,test_rho) == FAILURE)
    {
        DR_Error ("Invalid input correlation structure"
                    "between Fx rate and swap rates in Hyb3_Get_FxVol\n");
        return (status);
    }

    if (DomCurve == NULL ||
        ForCurve == NULL ||
        DomMktVol == NULL ||
        ForMktVol == NULL ||
        fx_data == NULL ||
        FxSpotVol == NULL ||
        FxSpotVolDates == NULL ||
        VolExpDates == NULL ||
        CompositeVol == NULL)
    {
        DR_Error("Invalid pointer inputs to Hyb3_Get_FxVol\n");
        return (status);
    }
    
    /* Create time line */
    if (MergeDateLists(NbFxSpotVols,
                       FxSpotVolDates,
                       &(NbTP),
                       &(TPDate)) == FAILURE)
    {
        DR_Error("MergeDateLists failed for TotNbFxVols in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    if (MergeDateLists(
                        DomMktVol->NbVol,
                        DomMktVol->VolDate,
                        &(NbTP),
                        &(TPDate)) == FAILURE)
    {
        DR_Error("MergeDateLists failed in Hyb3_Get_FxVol for Dom vols\n");
        goto RETURN;
    }

    if (MergeDateLists(
                        ForMktVol->NbVol,
                        ForMktVol->VolDate,
                        &NbTP,
                        &TPDate) == FAILURE)
    {
        DR_Error("MergeDateLists failed in Hyb3_Get_FxVol for For vols\n");
        goto RETURN;
    }
        
    if (MergeDateLists(
                        NbFxVol,
                        VolExpDates,
                        &NbTP,
                        &TPDate) == FAILURE)
    {
        DR_Error("MergeDateLists failed in Hyb3_Get_FxVol for input fx vols\n");
        goto RETURN;
    }

    /* Add Fx base to the TPDate list */
    if (MergeDateLists(
                        1,
                        &(fx_data->Today),
                        &NbTP,
                        &TPDate) == FAILURE)
    {
        DR_Error("MergedDateLists failed in Hyb3_Get_FxVol for fx base date\n");
        goto RETURN;
    }
       
    /* recall that DrExtendAmerSch will remove all dates STRICTLY
     *before fx_data->Today*/
    if (DrExtendAmerSch(fx_data->Today,
                        Ppy,
                        &NbTP,
                        &TPDate,
                        NULL,
                        'L',
                        NULL,
                        NULL,
                        NULL,
                        NULL,
                        NULL) == FAILURE)
    {
        DR_Error ("DrExtendAmerSch failed in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    if (TPDate[0] != fx_data->Today)
    {
        DR_Error ("TPDate[0] is not fx_data.Today (Hyb3_Get_FxVol)\n");
        goto RETURN;
    }

    /* bootstrap IR spot vol using NORMAL parameters */
     if (Hyb3_SpotVol(   DomMktVol->Aweight,
                    DomMktVol->BaseDate,
                    DomMktVol->NbVol,
                    DomMktVol->VolDate,
                    DomMktVol->Vol,
                    DomMktVol->VolUsed,
                    DomMktVol->Freq,
                    DomMktVol->DCC,
                    DomMktVol->SwapSt,
                    DomMktVol->SwapMat,
                    0., /* "Normal" QLeft    */                         
                    0., /* "Normal" QRight   */                   
                    0., /* No forward shift  */                    
                    1,  /* 1 Factor          */
                    DomMktVol->Alpha,
                    DomMktVol->Beta,
                    DomMktVol->Rho,
                    DomMktVol->SkipFlag,
                    DomMktVol->CalibFlag,
                    DomCurve) != SUCCESS)
    {
        DR_Error("Hyb3_SpotVol for domestic failed in Get_SpotFxVol\n");
        goto RETURN;
    }

    if (Hyb3_SpotVol(    ForMktVol->Aweight,
                    ForMktVol->BaseDate,
                    ForMktVol->NbVol,
                    ForMktVol->VolDate,
                    ForMktVol->Vol,
                    ForMktVol->VolUsed,
                    ForMktVol->Freq,
                    ForMktVol->DCC,
                    ForMktVol->SwapSt,
                    ForMktVol->SwapMat,
                    0.,  /* "Normal" QLeft    */                  
                    0.,  /* "Normal" QRight   */                  
                    0.,  /* No forward shift  */                  
                    1,  /* 1 Factor          */
                    ForMktVol->Alpha,
                    ForMktVol->Beta,
                    ForMktVol->Rho,
                    ForMktVol->SkipFlag,
                    ForMktVol->CalibFlag,
                    ForCurve) != SUCCESS)
    {
        DR_Error("Hyb3_SpotVol for foreign failed in Get_SpotFxVol\n");
        goto RETURN;
    }
        
    /* Discretise spot vol onto time line */
    DomSpotVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    
    if (DomSpotVol == NULL)
    {
        DR_Error ("Invalid memory allocation for domestic"
                    "vol in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* The spot Vol needed for the rest of the function
     * is the Spot Vol of the EXPONENTIAL factors,as opposed to     
     * the orthogonal ones produced by triangulation,
     * therefore use only Aweight[0]                    */
    
    if (Hyb3_DiscretiseInput(DomSpotVol,
                        DomMktVol->Aweight[0],
                        DomMktVol->CalibFlag,
                        DomMktVol->NbVol,
                        DomMktVol->VolDate,
                        NbTP,
                        TPDate) == FAILURE)
    {
        DR_Error("Hyb3_DiscretiseInput failed for domestic"
                    "Spot vol in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* spot vol input to Hyb3_MultiFac_FxVol should be "weight free",
     * so rescale.                                           */
    
    alpha = DomMktVol->Alpha[0];

    /* Hyb3_Param_Check ensures that Alpha[0] non zero*/
    for (i = 0; i < NbTP - 1;i++)
    {
        DomSpotVol[i] /= alpha; 
    }

    ForSpotVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    
    if (ForSpotVol == NULL)
    {
        DR_Error ("Invalid memory allocation for foreign vol"
                    "in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
   
    if (Hyb3_DiscretiseInput(ForSpotVol,
                        ForMktVol->Aweight[0],
                        ForMktVol->CalibFlag,
                        ForMktVol->NbVol,
                        ForMktVol->VolDate,
                        NbTP,
                        TPDate) == FAILURE)
    {
        DR_Error("Hyb3_DiscretiseInput failed for foreign"
                    "Spot vol in Get_Spot_fxVol\n");
        goto RETURN;
    }
                        
    alpha = ForMktVol->Alpha[0];
    for (i = 0; i < NbTP - 1; i++)
    {
        ForSpotVol[i] /= alpha;
    }
    
    Discrete_SpotVol = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    if (Discrete_SpotVol == NULL)
    {
        DR_Error("memory allocation failure in Hyb3_Get_FxVol\n");
        goto RETURN;
    }

    
/*    if (Hyb3_DiscretiseInput(Discrete_SpotVol,
                        SpotFxVol,
                        TRUE,   
                        TotNbFxSpotVols,
                        FxSpotVolDates,
                        NbTP,
                        TPDate) == FAILURE)
    {
        DR_Error("Hyb3_DiscretiseInput failed for fx"
                    "Spot vol in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
*/    
    
    for (i = 0, k = 1 ; i <= NbTP - 2; i++)
    {
        if (k <= (NbFxSpotVols - 1))
        {
            if (TPDate[i] >= FxSpotVolDates[k]) k++;
            Discrete_SpotVol[i] = FxSpotVol[k - 1];
        }
        else
            Discrete_SpotVol[i] = FxSpotVol[k - 1];
    }/* for i */
    
    
    /* interpolate relevant forward rates*/
    DomFwd = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    ForFwd = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    
    if (DomFwd == NULL || ForFwd == NULL)
    {
        DR_Error ("Failed to allocate memory for forward rates"
                    "in Hyb3_Get_FxVol\n");
        goto RETURN;
    }

    if (Hyb3_SimpleFwdCurve( DomFwd,
                        DomCurve,
                        TPDate,
                        NbTP) == FAILURE)
    {
        DR_Error ("Hyb3_SimpleFwdCurve failed for domestic"
                    "curve in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
  

    if (Hyb3_SimpleFwdCurve( ForFwd,
                        ForCurve,
                        TPDate,
                        NbTP) == FAILURE)
    {
        DR_Error ("Hyb3_SimpleFwdCurve failed for foreign"
                    "curve in Hyb3_Get_FxVol\n");
        goto RETURN;
    }
    
    /* we are now ready to call MultiFac_FXVol*/

    if (Hyb3_MultiFac_FxVol(
                        CompositeVol,
                        DomSpotVol,
                        ForSpotVol,
                        DomMktVol->Alpha,
                        ForMktVol->Alpha,
                        DomMktVol->Beta,
                        ForMktVol->Beta,
                        DomFwd,
                        ForFwd,
                        1,/* 1 factor*/
                        1,/* 1 factor*/
                        &(fx_data->Rho[2][0]),/* Rho Fx/Dom*/
                        &(fx_data->Rho[1][0]),/* Rho Fx/For*/
                        &(fx_data->Rho[0][0]),/* Rho Dom/For */
                        DomMktVol->Rho,
                        ForMktVol->Rho,
                        Discrete_SpotVol,
                        VolExpDates,
                        NbFxVol,
                        TPDate,
                        NbTP) == FAILURE)
    {
        DR_Error ("Hyb3_MultiFac_FxVol failed in get_FxVol\n");
        goto RETURN;
    }

    status = SUCCESS;
RETURN: 
    if (NbTP >= 1)
    {
        Free_DR_Array (TPDate,LONG,0,NbTP - 1);
    }

    if (NbTP >= 2)
    {
        Free_DR_Array (DomSpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (ForSpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (Discrete_SpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (DomFwd,DOUBLE,0,NbTP - 2);
        Free_DR_Array (ForFwd,DOUBLE,0,NbTP - 2);
    }
    
    return (status);

}/* Hyb3_Get_FxVol2*/

/**********  Hyb3_CalcVolA_xx  ************************************************************************/
/**
*  Returns volA_xx, as defined in the appendix of "Long-dated FX smile model", HR, 15 May 2003
*  
*  Output: array of doubles (volA_xx[0,n],...,volA_xx[n-1,n])
*          where volA_xx[i,j] is the time T_i curvature for an ATMF option with expiry T_j      
*
*  Notes:
*
*  1) Implements Eq. 20 with the following approximation (in latex terms)
* 
*                     [ \tau \Sigma_{A,x}^2 ]_{\tau} \approx \Sigma_{A,x}^2 
*
*  2) if "vol" is fwd fx vol - maturity t[n] - then "loc_vol_x" and "loc_vol_xx" need to be the
*     (normalised) slope and curvature of the local fwd fx vol rather than of the local spot fx 
*     vol  
*
***************************************************************************************************/

int Hyb3_CalcVolA_xx (
                      double  *volA_xx,       /**< (O) indexed [0,n-1]                                 */
                      long    n,              /**< (I) index of last time point accessed, t[n] etc     */
                      double  *vol,           /**< (I) spot fx vol, [0,n-1]                            */
                      double  *loc_vol_x,     /**< (I) slope of local vol wrt log-money, [0,n]         */
                      double  *loc_vol_xx,    /**< (I) curvature of local vol wrt log-money, [0,n]     */
                      double  *t,             /**< (I) time grid, Time[0] = 0, [0,n]                   */
                      double  *dt)            /**< (I) intervals in grid, [0,n-1]                      */ 
{

       double   s1    = 0.0, 
                s2    = 0.0, 
                s3    = 0.0, 
                g,
                u1, u2, 
                tau_0, tau_1;

       int      i, status = FAILURE;

       if ( 
           (!volA_xx)    ||
           (!vol)        ||
           (!loc_vol_x)  ||
           (!loc_vol_xx) ||
           (!t)          ||
           (!dt)        
          )
       {
           DR_Error("CalcVolA_xx: NULL pointer passed\n");
           goto RETURN;
       }
       
       if (n <= 0) 
       {
           DR_Error("CalcVolA_xx: number of time points must be positive\n");
           goto RETURN;
        }
       
       if (n >= MAXNBTP) 
       {
           DR_Error("CalcVolA_xx: too many time points\n");
           goto RETURN;
       }
       
       /* explicit Euler integration - requires dense time line */
       for (i = n-1; i >= 0; i--)
       {

           if ( (vol[i] < TINY) ) 
           {
               DR_Error("CalcVolA_xx: vol too small\n");
               goto RETURN;
           }

           tau_0 = t[n] - t[i];         
           tau_1 = t[n] - t[i+1];
   
           s1   += vol[i] * vol[i] * dt[i];
           /* u1 = Sigma_A(t[i],t[n])   */
           u1    = sqrt(s1 / tau_0);    

           s2   += 0.5 * vol[i] * vol[i] * loc_vol_x[i] * u1 * u1  
               * (tau_0 * tau_0 - tau_1 * tau_1);
           /* u2 = Sigma_A_x(t[i],t[n]) */	
           u2    = s2 * u1 / (s1 * s1);
 
           g     = - u2 * u2 +
                   ( loc_vol_x[i] * loc_vol_x[i] + loc_vol_xx[i] 
                     - 4.0 * loc_vol_x[i] * u2 / u1
                     + 3.0 * u2 * u2 / (u1 * u1) ) 
                     * vol[i] * vol[i]; 
           /* s3 uses approx: time deriv of u2 is negligible */                     
           s3   += (tau_0 * tau_0 * tau_0 - tau_1 * tau_1 * tau_1)
                    * u1 * u1 * u1 * u1 
                    * g / 3.0;
           /* volA_xx[t[i],t[n]] */
           volA_xx[i] = s3 * u1 / (s1 * s1 * s1);

       } /* i */
       
       status = SUCCESS;


RETURN: 

       return(status);

}


/**********  Hyb3_CalcPerturb  ************************************************************************/
/**
*  Returns integral of vol^2(s) * volA(s,T_n) * volA_xx(s,T_n) / vol^2(T_{n-1}) over (t,T_{n-1}) 
*  output[i] contains the value of integral over [T_i,T_{n-1}]
*
*  Notes:
*
*  1) uses two approximations, 
* 
*                     [ \tau \Sigma_{A,x}^2 ]_{\tau} \approx \Sigma_{A,x}^2        
*                     sigma \approx sigma(t)  instead of sigma(t,T)
*
*  2) t[n] used  
*
***************************************************************************************************/

int Hyb3_CalcPerturb ( 
             double  *output,        /**< (O) indexed [0,n-1]                                  */
             long    n,              /**< (I) index of last time point accessed, t[n] etc      */
             double  *vol,           /**< (I) spot fx vol, [0,n-1]                             */
             double  *loc_vol_x,     /**< (I) slope of local vol wrt log-money, [0,n]          */
             double  *loc_vol_xx,    /**< (I) curvature of local vol wrt log-money, [0,n]      */
             double  *t,             /**< (I) time grid, Time[0] = 0, [0,n]                    */
             double  *dt)            /**< (I) intervals in grid, [0,n-1]                       */ 
{

       double   s1    = 0.0, 
                s2    = 0.0, 
                s3    = 0.0,
                s4    = 0.0,
                s5    = 0.0,
                g,
                volA_xx = 0.0,
                u1, u2, 
                tau_0, tau_1;

       int      i, status = FAILURE;

       if ( 
           (!output)     ||
           (!vol)        ||
           (!loc_vol_x)  ||
           (!loc_vol_xx) ||
           (!t)          ||
           (!dt)        
           )
       {
           DR_Error("CalcPerturb: NULL pointer passed\n");
           goto RETURN;
       }
       
       if (n <= 0) 
       {
       DR_Error("CalcPerturb: number of time points must be positive\n");
           goto RETURN;
       }
       
       if (n >= MAXNBTP) 
       {
           DR_Error("CalcPerturb: too many time points\n");
           goto RETURN;
       }
       
       /* explicit Euler integration, so requires dense time line     */
       for (i = n-1; i >= 0; i--)
       {
           if ( (vol[i] < TINY)) 
           {
               DR_Error("CalcPerturb: vol too small\n");
               goto RETURN;
           }
           
           tau_0 = t[n] - t[i];         
           tau_1 = t[n] - t[i+1];
   
           s1   += vol[i] * vol[i] * dt[i];
           /* u1 = Sigma_A(t[i],t[n])   */
           u1    = sqrt(s1 / tau_0);      

           s2   += 0.5 * vol[i] * vol[i] * loc_vol_x[i] * u1 * u1  
                   * (tau_0 * tau_0 - tau_1 * tau_1);
           /* u2 = Sigma_A_x(t[i],t[n]) */	
           u2    =  s2 * u1 / (s1 * s1);

           g     = - u2 * u2 +
                   (loc_vol_x[i] * loc_vol_x[i] + loc_vol_xx[i] 
                     - 4.0 * loc_vol_x[i] * u2 / u1
                     + 3.0 * u2 * u2 / (u1 * u1) ) 
                    * vol[i] * vol[i]; 
           /* s3 uses approx: time deriv of u2 is negligible */                     
           s3   += (tau_0 * tau_0 * tau_0 - tau_1 * tau_1 * tau_1)
                     * u1 * u1 * u1 * u1 
                 * g / 3.0;
           /* volA_xx[t[i],t[n]] */
           volA_xx = s3 * u1 / (s1 * s1 * s1);

           /* calculate smile pertubation                                              */
           /* \int_t^T vol^2 d [ tau * volA * volA_xx ] / dT du / vol[i]^2 at T = t[i] */ 

           s4   += (tau_0 * tau_0 - tau_1 * tau_1) * u1 * u1 * g / 2.0;

           s5   +=  vol[i] * vol[i] * ( - 2.0 * volA_xx / u1
                                        + 2.0 * s4 / (s1 * s1) ) * dt[i];

           output[i] = s5;

       } /* i */

       status = SUCCESS;


RETURN: 

        return(status);

}


/********** MultiFac_Spot_FXVol2 *********************************************/
/**Outputs as many spot fx vols as there are Vol Expiry dates (i.e NbVol).
   The correlation factors are with respect to the exponential factors
  
   NOTES:
  
   1) Assumes that TPDate[0]=fx_data.Today.
   2) Assumes and checks that each VolDate and VolIntegrationDate is on the time 
   line.
   3) VolDate and VolIntegrationDate must be entered in a strictly ascending 
   order.
   
 ********************************************************************************/


int Hyb3_MultiFac_Spot_FxVol2(
        double  *SpotFxVol,           /**< (O) Instant.fx vol                */
        double  *SpotVolDom,          /**< (I) Instant dom. int. rate vol    */
        double  *SpotVolFor,          /**< (I) Instant for. int.rate vol     */
        double  *AlphaDom,            /**< (I) dom. relative factor weights  */
        double  *AlphaFor,            /**< (I) for. relative factor weights  */
        double  *DomFwdRate,          /**< (I) dom. FwdRate @ each time pt   */
        double  *ForFwdRate,          /**< (I) for. FwdRate @ each time pt   */
        double  *RhoFxDomFac,         /**< (I) correl fx/dom. curve          */
        double  *RhoFxForFac,         /**< (I) correl fx/for. curve          */
        double  *RhoDomFacForFac,     /**< (I) correl dom./for. curves       */
        double  *RhoDomFac,           /**< (I) correl between dom. factors   */
        double  *RhoForFac,           /**< (I) correl between for. factors   */
        int     NbImpVol,             /**< (I) Nb of fx implied vol points   */
        long    *VolDate,             /**< (I) fwd fx maturity               */
        long    *VolIntegrationDate,  /**< (I) expiry of option              */
        double  *FxVol,               /**< (I) fx vol for option which       */
                                      /**<    expirs on VolIntegrationDate, */
                                      /**<     the underlying fwd matures on */
                                      /**<     VolDate                       */
        int     FxCutOffFlag,         /**< (I) True if cut off allowed       */
        int     FxCutOffLast,         /**< (I) True if c/off at last level   */
        double  FxCutOffLevel,        /**< (I) User-given or from cups.h     */
        double  *BetaDom,             /**< (I) Domestic mean reversion       */
        double  *BetaFor,             /**< (I) Foreign mean reversion        */
        int     DomNbFac,             /**< (I) Nb of fac. for dom. curve     */
        int     ForNbFac,             /**< (I) Nb of fac. for for.curve      */
        double  *A_Fx,                /**< (I) FX smile param, [0,NbTP-2]    */
        double  *B_Fx,                /**< (I) FX smile param, [0,NbTP-2]    */
        double  *C_Fx,                /**< (I) FX smile param, [0,NbTP-2]    */
        long    *TPDate,              /**< (I) date of each pt.              */
        int     NbTP)                 /**< (I) Total Nb of time steps        */
{

    double  Poly[3],                     /* coeffs of second order polynomial               */
            *DomA             = NULL,    /* A factors for dom.                              */
            *ForA             = NULL,    /* A factors for for.                              */
            *NoFxInt          = NULL,    /* integrals not containing the FX rate            */
            *WithFxInt        = NULL,    /* integrals due to covar. FX with dom/for         */
            *DeltaTime        = NULL,
            *Time             = NULL,   
            *a                = NULL,    /* short for A_Fx                                  */
            *b                = NULL,    /* short for B_Fx                                  */
            *c                = NULL,    /* short for C_Fx                                  */
            *volA             = NULL,    /* lowest order approximation                      */
            *loc_vol_x        = NULL,    /* local vol slope NORMALISED by FX spot vol       */
            *loc_vol_xx       = NULL,    /* local vol curvature NORMALISED by FX spot vol   */
            *SpotFxVolSoFar   = NULL,    /* values of spot FX Vol @ each TP                 */
            *FxVolExt         = NULL,    /* FxVol extended to TPDate timeline               */
            *FxVolUnSmile     = NULL,    /* Imp Vol, "FxVol" minus 1st order smile contrib  */
            *Perturb          = NULL;      

    double  FxVolPrev,
            epsilon, 
            tau_0, t_low, delta_t, slope, /* times used for interpolation of implied vols   */ 
            temp1, temp2,              
            Sigma = 0.0,                  /* Fx Spot Vol in current bucket                  */
            discriminant;                 /* discriminant of sec. order polynomial          */

    int     kVol,status = FAILURE,
            TPIntBefore,           /* integration time index of prev.Vol point */
            TPIntSoFar,            /* integration time index of curr.Vol point */
            TPMatSoFar,            /* time index of curr. fx mat date          */
            TPFirstImpVol,         /* time index for first implied vol         */
            TPLastImpVol,          /* time index for last implied vol          */
            NextTP,                /* temp variable used for TPDate offset     */ 
            VolTooLow,             /* flag for spot vol criterion              */
            i,k,l,m,n,
            CurrentNbInt,   /*used to index integrals not involving Fx rate */
            NbNoFxInt,      /*Total Nb of integrals not containing  FX rate */
            NbWithFxInt;    /*Total Nb of int. due to covar. Fx with For/Dom*/

    char    ErrorMsg[MAXBUFF];  /*Error message to be sent to DR_Error      */


    /* Quick checks first */

    if (DomNbFac != 1 && DomNbFac != 2 && DomNbFac != 3)   
    {
        DR_Error("Invalid Number of domestic yield curve factors "
                 "in Hyb3_MultiFac_Spot_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if (ForNbFac != 1 && ForNbFac != 2 && ForNbFac != 3)   
    {
        DR_Error("Invalid Number of foreign yield curve factors "
                "in Hyb3_MultiFac_Spot_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }

    if (NbTP < 2)
    {
        DR_Error ("Invalid number of time points in Hyb3_MultiFac_Spot_FxVol2: "
                    "should be at least 2\n");
        return (status);
    }

    if (NbImpVol < 1)
    {
        DR_Error ("Invalid number of Vol points in Hyb3_MultiFac_Spot_FxVol2: "
                    "should be at least 1");
        return (status);
    }
    
    if (RhoDomFac == NULL && DomNbFac > 1)
    {
        DR_Error("Dom factor correlation array cannot be NULL"
                " when more than 1 factor (Hyb3_MultiFac_Spot_FxVol2)\n");
        return(status);
    }
    if (RhoForFac == NULL && ForNbFac > 1)
    {
        DR_Error("Foreign factor correlation array cannot be NULL"
                " when more than 1 factor (Hyb3_MultiFac_Spot_FxVol2)\n");
        return(status);
    }

    if( SpotVolDom          == NULL ||
        SpotVolFor          == NULL ||
        AlphaDom            == NULL ||
        AlphaFor            == NULL ||
        DomFwdRate          == NULL ||
        ForFwdRate          == NULL ||
        RhoFxDomFac         == NULL ||
        RhoFxForFac         == NULL ||
        RhoDomFacForFac     == NULL ||
        VolDate             == NULL ||
        VolIntegrationDate  == NULL ||
        FxVol               == NULL ||
        BetaDom             == NULL ||
        BetaFor             == NULL ||
        TPDate              == NULL)
    {   
        DR_Error("Invalid pointer inputs to Hyb3_MultiFac_Spot_FxVol2\n");
        return (status);
    }
    
    NbNoFxInt = (DomNbFac + ForNbFac) + 
                (DomNbFac + ForNbFac - 1) * (DomNbFac + ForNbFac) / 2;

    NbWithFxInt = (DomNbFac + ForNbFac);


    NoFxInt           = (double*) DR_Array (DOUBLE,0,NbNoFxInt - 1);
    WithFxInt         = (double*) DR_Array (DOUBLE,0,NbWithFxInt - 1);
    DomA              = (double*) DR_Array (DOUBLE,0,DomNbFac - 1);
    ForA              = (double*) DR_Array (DOUBLE,0,ForNbFac - 1);
    DeltaTime         = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    Time              = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    a                 = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    b                 = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    c                 = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA              = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_x         = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_xx        = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    SpotFxVolSoFar    = (double*) DR_Array (DOUBLE,0,NbTP - 1); 
    FxVolExt          = (double*) DR_Array (DOUBLE,0,NbTP - 1); 
    FxVolUnSmile      = (double*) DR_Array (DOUBLE,0,NbTP - 1); 
    Perturb           = (double*) DR_Array (DOUBLE,0,NbTP - 1); 

    if( NoFxInt         == NULL ||
        WithFxInt       == NULL ||
        DomA            == NULL ||
        ForA            == NULL ||
        DeltaTime       == NULL ||
        Time            == NULL ||
        a               == NULL ||
        b               == NULL ||
        c               == NULL ||
        volA            == NULL ||
        loc_vol_x       == NULL ||
        loc_vol_xx      == NULL ||
        SpotFxVolSoFar  == NULL ||
        FxVolExt        == NULL ||
        FxVolUnSmile    == NULL ||
        Perturb         == NULL)
    {
        DR_Error("Memory allocation failed in Hyb3_MultiFac_Spot_FxVol2\n");
        status = SUCCESS;
        goto RETURN;
    }

    for (i = 0; i < NbTP - 1;i++)
    {   
    
        if ( TPDate[i+1] <= TPDate[i])
        {
            sprintf(ErrorMsg,"TPDate[%d] and TPDate[%d] in input time line"
                      "are not in a strictly ascending order\n",i,i+1);
            DR_Error (ErrorMsg);
            return (status);
        }
    
        /* time step between time points, and corresponding times in years */
        DeltaTime[i]  =    Daysact(TPDate[i],TPDate[i+1]) / 365.;
        Time[i]       =    Daysact(TPDate[0],TPDate[i]) / 365.; 
        a[i]          =    A_Fx[i];
        b[i]          =    B_Fx[i];
        c[i]          =    C_Fx[i];

        /* atm slope wrt log-fwd-money */  
        loc_vol_x[i]    =   c[i] * a[i];  
        /* atm curvature wrt log-fwd-money */  
        loc_vol_xx[i]   =   c[i] * c[i] * b[i];

    }

    /* we have already checked that NbTp >= 2 */
    Time[NbTP-1] = Daysact(TPDate[0],TPDate[NbTP-1]) / 365.; 

    /* find the first time point on or after the FIRST imp vol date */

    TPFirstImpVol = GetDLOffset(NbTP,
                                TPDate,
                                VolIntegrationDate[0],
                                CbkHIGHER);
    
    if (TPFirstImpVol < 0) 
    {
        sprintf(ErrorMsg, "Vol Integration date %ld is beyond the"
            "input time line (Hyb3_MultiFac_Spot_FxVol2)\n",VolIntegrationDate[0]);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    TPLastImpVol = GetDLOffset(NbTP,
                               TPDate,
                               VolIntegrationDate[NbImpVol-1],
                               CbkHIGHER);

    if (TPLastImpVol < 0)
    {
        sprintf(ErrorMsg, "Vol Integration date %ld is beyond the"
            "input time line (Hyb3_MultiFac_Spot_FxVol2)\n",VolIntegrationDate[NbImpVol-1]);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /************************************************************************/
    /*                                                                      */
    /*              remove 1st order FX smile contribution                  */
    /*                                                                      */
    /************************************************************************/

    /* first estimate zero expiry limit of implied vol = FxVolUnSmile[0]    */ 
    /* the estimate is obtained from Eq. 8 of "Long-dated FX smile model",  */
    /* by assuming FX spot vol constant and approximating RHS               */

    /* time to first option expiry */
    tau_0      =  Daysact(TPDate[0],VolIntegrationDate[0]) / 365.;

    /* approximation assumes const FX spot vol until first expiry           */
    epsilon = 2. * tau_0 * FxVol[0] * FxVol[0] * ( loc_vol_xx[0] 
                            - 0.5 * loc_vol_x[0] * loc_vol_x[0] ) / 3.0;

    if (epsilon < -0.5) /* sqr root of negative number */
    {
        sprintf(ErrorMsg, "Can't calculate initial FX spot vol (MultiFac_Spot_FxVol2)\n");
        DR_Error(ErrorMsg);
        goto RETURN;
    }  
    else if (fabs(epsilon) < TINY) 
    {
        FxVolUnSmile[0] = FxVol[0];
    }
    else 
    {
        FxVolUnSmile[0] = FxVol[0] * 
                          sqrt( (- 1. + sqrt(1 + 2. * epsilon) ) / epsilon);
    }

                                   /* extend the implied vols to full time line         */    
    t_low     = 0.0;               /* initialise as today                               */
    delta_t   = tau_0;             /* first step size equals first expiry               */
    FxVolPrev = FxVolUnSmile[0];   /* guess for implied vol extrapolated to zero expiry */
    slope     = FxVol[0] * FxVol[0] - FxVolUnSmile[0] * FxVolUnSmile[0];

    for (i = 0, k = 0 ; i <= NbTP - 1;i++)
    {   
        if (k <= NbImpVol - 1)
        {
            if (TPDate[i] >=  VolIntegrationDate[k])       
            {
                if (k < NbImpVol - 1)
                {
                    FxVolPrev  =  FxVol[k];
                    slope      =  FxVol[k+1] * FxVol[k+1] - FxVol[k] * FxVol[k] ;
                    t_low     +=  delta_t;
                    /* time step for next update   */
                    delta_t    =  Daysact(VolIntegrationDate[k],VolIntegrationDate[k+1]) / 365.;
                    k++;
                }
                else 
                {
                    /*  flat after last implied vol */
                    slope      = 0.0;
                    FxVolPrev  = FxVol[NbImpVol-1];
                    k++; /* i.e., k = NbImpVol */
                }
            }
            /* interpolate linearly on squares of vols - in the short end, in particular,  */ 
            /* this is more accurate than interpolation on the vols                        */
            FxVolExt[i] =  sqrt( FxVolPrev * FxVolPrev  + 
                slope * (Time[i] - t_low) / delta_t );
            /* instead of FxVolExt[i] /= ( IS_EQUAL(Time[i],0.0) ? 1.0 : sqrt(Time[i] ) );? */
            
        }
        else /* flat after last implied vol    */
        {         
            FxVolExt[i] =  FxVol[NbImpVol - 1];
        }
    }

    /* initialise dummy for integration */
    temp1   = FxVolUnSmile[0] * FxVolUnSmile[0] * DeltaTime[0];
    volA[1] = FxVolUnSmile[0];

    /* only volA is used outside this loop */
    for (i = 1; i < TPLastImpVol; i++)
    {   

        if (Hyb3_CalcPerturb( 
                         Perturb,                 /* (O) [0,i]...[i-1,i]                       */
                         i,                       /* (I) number of time points                 */
                         FxVolUnSmile,            /* (I) implied fx vol                        */
                         loc_vol_x,               /* (I) slope of local vol wrt log-money      */
                         loc_vol_xx,              /* (I) curvature of local vol wrt log-money  */
                         Time,                    /* (I) time grid, Time[0] = 0                */
                         DeltaTime)               /* (I) intervals in grid                     */
                         == FAILURE) goto RETURN;

        /* FxVol[i] etc is the implied vol at Time[i] - not Time[i+1] */

        temp2   =    ( Time[i+1] * FxVolExt[i+1] * FxVolExt[i+1] - 
                       Time[i] * FxVolExt[i] * FxVolExt[i] ) / 
                     ( 1.0 + Perturb[0] );

        if (temp2 < 0.0)
        {
            DR_Error("Hyb3_MultiFac_Spot_FxVol2: negative fx vol in calibrating interpolated implied vol at %ld\n", TPDate[i+1]);
            goto RETURN;
        }
        else
        {
            FxVolUnSmile[i] = sqrt(temp2/DeltaTime[i]);
        }
        temp1   +=  FxVolUnSmile[i] * FxVolUnSmile[i] * DeltaTime[i];

        /* Time[i] > 0 because i > 0 and TPDate is strictly increasing  */
        /* volA[i] refers to Time[i]                                    */
        volA[i+1]  =  sqrt(temp1/Time[i+1]); 

    }

    /************* 1st order FX smile contribution removed ******************/


    TPIntBefore = TPIntSoFar = 0;
    TPMatSoFar = 0;

    for (kVol = 0; kVol < NbImpVol; kVol++)
    {   

        /* check that IntegrationDate and VolDate are in the right order*/
        if (VolIntegrationDate[kVol] > VolDate[kVol])
        {
            sprintf(ErrorMsg, "Vol Integration date %ld is after"
                    "corresponding vol mat date (Hyb3_MultiFac_Spot_FxVol2)\n",
                    VolIntegrationDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
        if ( VolIntegrationDate[kVol] < TPDate[0])
        {
            DR_Error("VolIntegrationdate cannot be before TPDate[0].\n");
            goto RETURN;
        }

        /* check that VolIntegrationDate is on input time line         */

        /* - GetDLOffset searches for the new TPIntSoFar               */
        /*   in the range {TPDate[TPIntSoFar],TPDate[NbTP-1]}          */
        /* - since function returns offset relative to TPIntSoFar,     */
        /*   we add it to old value of TPIntSoFar                      */
        NextTP = GetDLOffset(NbTP - TPIntSoFar,
                             TPDate + TPIntSoFar,
                             VolIntegrationDate[kVol],
                             CbkEXACT);

        TPIntSoFar += NextTP;

        if (NextTP < 0)
        {
            sprintf(ErrorMsg, "Vol Integration date %ld is not on the"
                    "input time line (or dates are not in ascending order)"
                    "(Hyb3_MultiFac_Spot_FxVol2)\n",VolIntegrationDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
        /* check that Voldate is on input time line                    */

        /* - GetDLOffset searches for the new TPMatSoFar               */
        /*   in the range {TPDate[TPIntSoFar],TPDate[NbTP-1]}          */
        /* - since function returns offset relative to TPIntSoFar,     */
        /*   we add it to old value of TPIntSoFar                      */
        NextTP = GetDLOffset(NbTP - TPIntSoFar,
                             TPDate + TPIntSoFar,
                             VolDate[kVol],
                             CbkEXACT);

        TPMatSoFar = TPIntSoFar + NextTP;

        if (NextTP < 0)
        {
            sprintf(ErrorMsg, "Vol date %ld is not on the"
                    "input time line (Hyb3_MultiFac_Spot_FxVol2)\n",VolDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* the following should not normally happen as fx   */ 
        /* dates are input using New_Fx_Input_W,            */
        /* which checks for strictly ascending order and    */
        /* eliminates dates before FxBaseDate=TPDate[0]     */

        if (TPIntSoFar == TPIntBefore)
        {
          continue;
        }
        
        /* first,some initialisations*/
        
        for (i = 0; i < DomNbFac;i++)
        {
            DomA[i] = 0.;

        } /* for i*/
        
        for (i = 0; i < ForNbFac;i++)
        {
            ForA[i] = 0.;

        } /* for i*/

        for (i = 0; i < NbNoFxInt; i++)
        {
            NoFxInt[i] = 0.;

        } /* for i*/
        
        for (i = 0; i < NbWithFxInt; i++)
        {
            WithFxInt[i] = 0.;

        } /* for i*/

        for (i = 0; i < 3; i++)
        {
            Poly[i] = 0.;

        } /* for i*/
        
        /* prepare A-factors    */

        for (i = TPMatSoFar - 1; i>= TPIntSoFar; i--)
        {
            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-BetaDom[l] * DeltaTime[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0;l < ForNbFac;l++)
            {
                ForA[l] *= exp(-BetaFor[l] * DeltaTime[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(BetaFor[l],DeltaTime[i]);

            }/* for l */


        }/* for i */
        
        /* we can start integrating now */
        for (i = TPIntSoFar - 1; i >= TPIntBefore; i--)
        {
            CurrentNbInt = -1;
            
            /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-BetaDom[l] * DeltaTime[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0;l < ForNbFac;l++)
            {
                ForA[l] *= exp(-BetaFor[l] * DeltaTime[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(BetaFor[l],DeltaTime[i]);

            }/* for l */

            
            /* variance of domestic factors*/
            
            for (l = 0;l < DomNbFac;l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += AlphaDom[l] * AlphaDom[l] * 
                                         SpotVolDom[i] * SpotVolDom[i]
                                         * DomA[l] * DomA[l] * DeltaTime[i];
            }
            
            /* covariance between domestic factors*/

            if (DomNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* AlphaDom[0] * AlphaDom[1] * 
                                         DomA[0] * DomA[1]
                                         * SpotVolDom[i] * SpotVolDom[i] *
                                         RhoDomFac[0]*DeltaTime[i];

            }/* end of case: DomNbFac > 1 */
            
            if (DomNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * AlphaDom[0] * AlphaDom[2] *
                                         DomA[0] * DomA[2]
                                         * SpotVolDom[i] * SpotVolDom[i] * 
                                         RhoDomFac[1] * DeltaTime[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * AlphaDom[1] * AlphaDom[2] * 
                                         DomA[1] * DomA[2]
                                         * SpotVolDom[i] * SpotVolDom[i] * 
                                         RhoDomFac[2] * DeltaTime[i];

            }/* end of case DomNbFac > 2 */

            /* variance of foreign factors*/
            
            for (l = 0; l < ForNbFac;l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += AlphaFor[l] * AlphaFor[l] * 
                                         SpotVolFor[i] * SpotVolFor[i]
                                         * ForA[l] * ForA[l] * DeltaTime[i];

            } /* for l */

            /* covariance between foreign factors*/

            if (ForNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* AlphaFor[0] * AlphaFor[1] *
                                         ForA[0] * ForA[1]
                                         * SpotVolFor[i] * SpotVolFor[i] 
                                         * RhoForFac[0] * DeltaTime[i];
            }
            
            if (ForNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* AlphaFor[0] * AlphaFor[2] * 
                                         ForA[0]*ForA[2]
                                         * SpotVolFor[i] * SpotVolFor[i] *
                                         RhoForFac[1] * DeltaTime[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * AlphaFor[1] * AlphaFor[2] * 
                                         ForA[1] * ForA[2]
                                         * SpotVolFor[i] * SpotVolFor[i] * 
                                         RhoForFac[2] * DeltaTime[i];
            }
            
            /* covariance between domestic and foreign factors */

            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {
                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt] -= 
                                2.0 * AlphaDom[l] * DomA[l] * SpotVolDom[i] *
                                AlphaFor[m] * ForA[m] * SpotVolFor[i] * 
                                RhoDomFacForFac[m + l * ForNbFac]
                                * DeltaTime[i];
                }/* for m */

            }/* for l */
            
            /*covariance between Fx and domestic*/

            for (l = 0; l < DomNbFac; l++)
            {
                WithFxInt[l] += 2.0 * SpotVolDom[i] * AlphaDom[l] * 
                                DomA[l] * RhoFxDomFac[l] * DeltaTime[i];

            }/* for l*/


            /*covariance between fx and for*/

            for (l = 0; l < ForNbFac; l++)
            {
                WithFxInt[l+DomNbFac] -= 2.0 * SpotVolFor[i] * 
                                         AlphaFor[l] * ForA[l] * 
                                         RhoFxForFac[l] * DeltaTime[i];

            }/* for l*/

        }   /* for i*/
        
        for (i = TPIntBefore - 1; i >= 0; i--)
        {               
            CurrentNbInt = -1;

          /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-BetaDom[l] * DeltaTime[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0; l < ForNbFac; l++)
            {
                ForA[l] *= exp(-BetaFor[l] * DeltaTime[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(BetaFor[l],DeltaTime[i]);

            }/* for l */  
                        
            /* variance of each domestic factor*/
            
            for (l = 0; l < DomNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += AlphaDom[l] * AlphaDom[l] * 
                                         SpotVolDom[i] * SpotVolDom[i]
                                         * DomA[l] * DomA[l] * DeltaTime[i];
            }
            
            /* covariance between domestic factors*/

            if (DomNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* AlphaDom[0] * AlphaDom[1] * 
                                         DomA[0] * DomA[1]
                                         * SpotVolDom[i] * SpotVolDom[i] * 
                                         RhoDomFac[0]*DeltaTime[i];
            }
            
            if (DomNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * AlphaDom[0] * AlphaDom[2] * 
                                         DomA[0]*DomA[2]
                                         * SpotVolDom[i] * SpotVolDom[i] * 
                                         RhoDomFac[1] * DeltaTime[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * AlphaDom[1] * AlphaDom[2] * 
                                         DomA[1] * DomA[2]
                                         * SpotVolDom[i] * SpotVolDom[i] * 
                                         RhoDomFac[2] * DeltaTime[i];

            }

            /* variance of foreign factors*/
            
            for (l = 0; l < ForNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += AlphaFor[l] * AlphaFor[l] * 
                                         SpotVolFor[i]*SpotVolFor[i]
                                         * ForA[l] * ForA[l]*DeltaTime[i];

            } /* for l */
            
            /* covariance between foreign factors*/

            if (ForNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* AlphaFor[0] * AlphaFor[1] * 
                                         ForA[0] * ForA[1]
                                         * SpotVolFor[i] * SpotVolFor[i] *
                                         RhoForFac[0]*DeltaTime[i];

            }
            
            if (ForNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * AlphaFor[0] * AlphaFor[2] *
                                         ForA[0]*ForA[2]
                                         * SpotVolFor[i] * SpotVolFor[i] * 
                                         RhoForFac[1] * DeltaTime[i];
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * AlphaFor[1] * AlphaFor[2] * 
                                         ForA[1] * ForA[2]
                                         * SpotVolFor[i] * SpotVolFor[i] *
                                         RhoForFac[2] * DeltaTime[i];
            }
           
            /* covariance between domestic and foreign factors */

            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {
                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt] -= 2.0 * AlphaDom[l] * DomA[l] * 
                                             SpotVolDom[i] * SpotVolFor[i] *
                                             AlphaFor[m] * ForA[m] *
                                             RhoDomFacForFac[m + l * ForNbFac] * DeltaTime[i];
                }/* for m */

            }/* for l */
                        
            /* constant coefficients due to covar fx/dom*/
            
            for (l = 0; l < DomNbFac; l++)
            {
                Poly[0] += 2.0 * SpotVolDom[i] * AlphaDom[l] * DomA[l] *
                                                 SpotFxVolSoFar[i] 
                                               * RhoFxDomFac[l] * DeltaTime[i];
            }

            /*constant coeffs due to covar fx/for*/

            for (l = 0; l < ForNbFac; l++)
            {
                Poly[0] -= 2.0 * SpotVolFor[i] * AlphaFor[l] * 
                                                 ForA[l] * SpotFxVolSoFar[i]
                                               * RhoFxForFac[l] * DeltaTime[i];
            }

            /* constant coeffs due to fx variance*/

            Poly[0] += SpotFxVolSoFar[i] * SpotFxVolSoFar[i] * DeltaTime[i];

        }/* for i */

        for (l = 0; l < NbNoFxInt; l++)
        {
            Poly[0] += NoFxInt[l];
        }

        /* 1) we have assumed TPDate[0] = FxBaseDate                         */
        /* 2) volA is defined on the "finer" time line  {TPDate}             */
        /*    so we use the index TPIntSoFar, for which                      */  
        /*    TPdate[TPIntSoFar] = VolDate[kVol]                             */
        
        Poly[0] -= volA[TPIntSoFar] * volA[TPIntSoFar] * 
                            Daysact(TPDate[0],TPDate[TPIntSoFar]) / 365.;

        for (l = 0; l < NbWithFxInt; l++)
        {
            Poly[1] += WithFxInt[l];
        }
        
        /* notice that TPDate[TPIntBefore] < VolIntegrationDate[kVol]        */
        /* i.e (Poly[2] non zero)                                            */
        /* because of earlier check on TPIntBefore==TPIntSoFar               */

        Poly[2] = (double) Daysact(TPDate[TPIntBefore],VolIntegrationDate[kVol])
                           / 365.;

        discriminant = Poly[1] * Poly[1] - 4.* Poly[2] * Poly[0];

        VolTooLow = (discriminant < TINY);
  
        if (discriminant >= TINY)
        {       
            /* using the larger root ensures a continuous spot fx vol        */ 
            Sigma  = (-Poly[1] + sqrt(discriminant)) / 2. / Poly[2];
                            
            /* Evaluate whether spotvol satisfies appropriate criterion */
            if((FxCutOffLast) || (!FxCutOffFlag))
            {
                VolTooLow  = (FxVol[kVol] > FxCutOffLevel * Sigma);
            }
            else
            {
                VolTooLow = (Sigma < FxCutOffLevel);
            }
        }
        /* now implement cutoff actions as the case may be */

        if (VolTooLow)
        {
            if (FxCutOffFlag)
            {
                if (FxCutOffLast)
                {
                    if(kVol == 0)
                    {
                        DR_Error("Unable to bootstrap at least one FX\n"
                                "vol point and therefore unable to\n"
                                "cut off at last spot vol level.\n");

                        goto RETURN;
                    }
                    else
                    {   
                        for (n = kVol; n < NbImpVol; n++)
                        {
                            SpotFxVol[n] = SpotFxVol[kVol-1];
                        }
                    }
                }/* end of FxCutOffLast=true */

                else
                {                           
                    for(n = kVol; n < NbImpVol; n++)
                    {                            
                        SpotFxVol[n] = FxCutOffLevel;
                    }
                }

                break;

            }/* end of FxCutOffFlag */

            else
            {
                sprintf(ErrorMsg,"Problem in calculating fx vol at %ld:less"
                    " than %4.2lf%% \n(level determined by ratio to fwd vol)"
                    ".\n(MultiFac_FX_Spot_Vol)\n",
                    VolDate[kVol],100.*FxVol[kVol]/FxCutOffLevel);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        } /* end of VolTooLow */

        SpotFxVol[kVol] = Sigma;

        for (i = TPIntBefore; i < TPIntSoFar; i++)
        {
            SpotFxVolSoFar[i] = Sigma;

        }/* for i*/

        TPIntBefore  = TPIntSoFar;

    
    }/*for kVol*/
    
    status = SUCCESS;

RETURN:

    Free_DR_Array (NoFxInt,          DOUBLE,0,NbNoFxInt - 1);
    Free_DR_Array (WithFxInt,        DOUBLE,0,NbWithFxInt - 1);
    Free_DR_Array (DomA,             DOUBLE,0,DomNbFac - 1);
    Free_DR_Array (ForA,             DOUBLE,0,ForNbFac - 1);
    Free_DR_Array (DeltaTime,        DOUBLE,0,NbTP - 1);
    Free_DR_Array (Time,             DOUBLE,0,NbTP - 1);
    Free_DR_Array (a,                DOUBLE,0,NbTP - 1);     
    Free_DR_Array (b,                DOUBLE,0,NbTP - 1);    
    Free_DR_Array (c,                DOUBLE,0,NbTP - 1);   
    Free_DR_Array (volA,             DOUBLE,0,NbTP - 1);
    Free_DR_Array (loc_vol_x,        DOUBLE,0,NbTP - 1);   
    Free_DR_Array (loc_vol_xx,       DOUBLE,0,NbTP - 1);  
    Free_DR_Array (SpotFxVolSoFar,   DOUBLE,0,NbTP - 1);
    Free_DR_Array (FxVolExt,         DOUBLE,0,NbTP - 1);
    Free_DR_Array (FxVolUnSmile,     DOUBLE,0,NbTP - 1);   
    Free_DR_Array (Perturb,          DOUBLE,0,NbTP - 1);  

    return(status);

}/* Hyb3_MultiFac_Spot_FxVol2*/


/*************** MultiFac_FXVol2 *********************************************/
/**Outputs as many implied fx vols as there are Vol Expiry dates (i.e NbVol).
   The correlation factors are with respect to the exponential factors
  
   NOTES:
   1-Assumes that TPDate[0]=fx_data.Today.
   2-Assumes and checks that each VolDate and FxRstDate is on the time 
   line.
   3-VolDate and FxMatDate must be entered in a strictly ascending 
   order.
   
 *****************************************************************************/

int Hyb3_MultiFac_FxVol2(
                double  *ImpVol,           /**<(O) Implied FX Vol for FX Expiries       */
                double  *SpotVolDom,       /**<(I) Instant.dom int.rate vol             */
                double  *SpotVolFor,       /**<(I) Instant.for.int rate vol             */
                double  *AlphaDom,         /**<(I) factor weights                       */
                double  *AlphaFor,         /**<(I) factor weights                       */
                double  *BetaDom,          /**<(I) domestic mean reversions             */
                double  *BetaFor,          /**<(I) foreign mean reversions              */
                double  *DomFwdRate,       /**<(I) dom rates used in Christian's A facs */
                double  *ForFwdRate,       /**<(I) for rates used in Christian's A facs */
                int      DomNbFac,         /**<(I) Nb factors for Dom yield curve       */
                int      ForNbFac,         /**<(I) Nb factors for For yield curve       */
                double  *RhoFxDomFac,      /**<(I) Corr. between fx and Dom factors     */
                double  *RhoFxForFac,      /**<(I) Corr between fx and For factors      */
                double  *RhoDomFacForFac,  /**<(I) Corr between Dom and For factors     */
                double  *RhoDomFac,        /**<(I) Corr between Dom factors             */
                double  *RhoForFac,        /**<(I) Corr between For fcators             */
                double  *SpotFxVol,        /**<(I) Inst fx vol                          */
                double  *A_Fx,
                double  *B_Fx,
                double  *C_Fx,
                long    *FxMatDate,        /**<(I) maturity date                        */
                long    *FxRstDate,        /**<(I) fx reset date                        */
                int      NbVol,            /**<(I) Nb of output ImpVols                 */
                long    *TPDate,           /**<(I) Dates on time line for vol integrat. */
                int      NbTP)             /**<(I) Nb of time points                    */
{

    int           i, k, m, l, status=FAILURE;

    int
                  NbNoFxInt,           /* Nb of integrals not containing spot fx vol  */
                  NbWithFxInt,         /* Nb of integrals due to covar fx/Dom,fx/For  */
                  IndexRst,IndexMat,   /* FX reset and maturity indices on timeline   */
                  CurrentNbInt;        /* used to index integrlas not involving Fx    */

    double
                  *NoFxInt        = NULL,    /* Integrals not containing spot fxvol         */
                  *WithFxInt      = NULL,    /* Integrals due to covar Fx/Dom and Fx/For    */
                  *DomA           = NULL,
                  *ForA           = NULL,
                  *DeltaTime      = NULL,
                  *Time           = NULL,
                  *a              = NULL,    /* short for A_Fx                        */
                  *b              = NULL,    /* short for B_Fx                        */
                  *c              = NULL,    /* short for C_Fx                        */
                  *volA           = NULL,    /* lowest order approximation            */
                  *volA2          = NULL,    /* square of lowest order approx         */
                  *volA_x         = NULL,    /* slope wrt log moneyness               */
                  *volA_xx        = NULL,    /* curvature wrt log moneyness           */
                  *loc_vol_x      = NULL,    /* local vol slope                       */
                  *loc_vol_xx     = NULL;    /* local vol curvature                   */ 

    double        
                  varFX,
                  volB_sum = 0.;
            
    char
                  ErrorMsg[MAXBUFF];
    

    /* Quick checks first */
    if (DomNbFac != 1 && DomNbFac != 2 && DomNbFac != 3)
     
    {
        DR_Error("Invalid Number of domestic yield curve factors "
                 "in Hyb3_MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if (ForNbFac != 1 && ForNbFac != 2 && ForNbFac != 3)
     
    {
        DR_Error("Invalid Number of foreign yield curve factors "
                 "in Hyb3_MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if( NbTP < 2)
    {
        DR_Error ("Invalid number of time points in Hyb3_MultiFac_FxVol: "
                  "should be at least 2\n");
        return (status);
    }

    if( NbVol < 1)
    {
        DR_Error ("Invalid number of Vol points in MultiFac_fxVol: "
                  "should be at least 1");
        return (status);
    }

    if (ImpVol             == NULL ||
        SpotVolDom         == NULL ||
        SpotVolFor         == NULL ||
        AlphaDom           == NULL ||
        AlphaFor           == NULL ||
        BetaDom            == NULL ||
        BetaFor            == NULL ||
        DomFwdRate         == NULL ||
        ForFwdRate         == NULL ||
        RhoFxDomFac        == NULL ||
        RhoFxForFac        == NULL ||
        RhoDomFacForFac    == NULL ||
        SpotFxVol          == NULL ||
        FxMatDate          == NULL ||
        FxRstDate          == NULL ||
        TPDate             == NULL)

    {
        DR_Error ("invalid pointer inputs to MultiFx_FxVol\n");
        return (status);
    }

    if(DomNbFac > 1 && RhoDomFac == NULL)
    {
        DR_Error("cannot have more than 1 factor"
                 " and a NULL factor correlation inputs\n");
        return (status);
    }

    if(ForNbFac > 1 && RhoForFac == NULL)
    {
        DR_Error("cannot have more than 1 factor"
                 " and a NULL factor correlation inputs\n");
        return (status);
    }

    NbNoFxInt = (DomNbFac + ForNbFac) + 
                (DomNbFac + ForNbFac - 1) * (DomNbFac + ForNbFac) / 2;

    NbWithFxInt = (DomNbFac + ForNbFac);


    NoFxInt         = (double*) DR_Array (DOUBLE,0,NbNoFxInt - 1);
    WithFxInt       = (double*) DR_Array (DOUBLE,0,NbWithFxInt -1);
    DomA            = (double*) DR_Array (DOUBLE,0,DomNbFac - 1);
    ForA            = (double*) DR_Array (DOUBLE,0,ForNbFac - 1);
    DeltaTime       = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    Time            = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    a               = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    b               = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    c               = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA            = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA2           = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA_x          = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA_xx         = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_x       = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_xx      = (double*) DR_Array (DOUBLE,0,NbTP - 1);

    if( NoFxInt      == NULL ||
        WithFxInt    == NULL ||
        DomA         == NULL ||
        ForA         == NULL ||
        DeltaTime    == NULL ||
        Time         == NULL ||
        a            == NULL ||
        b            == NULL ||
        c            == NULL ||
        volA         == NULL ||
        volA2        == NULL ||
        volA_x       == NULL ||
        volA_xx      == NULL ||
        loc_vol_x    == NULL ||
        loc_vol_xx   == NULL)
    {
        goto RETURN;
    }

    for (i = 0; i < NbTP - 1; i++)
    {   

        if (TPDate[i+1] <= TPDate[i])
        {

            sprintf(ErrorMsg,"TPDate[%d] and TPDate[%d] in input time line"
                             "are not in a strictly"
                             " ascending order\n",i,i+1);
            DR_Error (ErrorMsg);
            return (status);

        }

        /* time step between time points, and corresponding times in years */
        DeltaTime[i] = Daysact(TPDate[i],TPDate[i+1]) / 365.;
        Time[i] = Daysact(TPDate[0],TPDate[i]) / 365.; 

        /* change fx smile arrays to short hand */
        a[i]    =    A_Fx[i];
        b[i]    =    B_Fx[i];
        c[i]    =    C_Fx[i];

        /* atm slope wrt log-fwd-money */  
        loc_vol_x[i]    =   c[i] * a[i];  
        /* atm curvature wrt log-fwd-money */  
        loc_vol_xx[i]   =   c[i] * c[i] * b[i];

    }

    Time[NbTP-1] =  Daysact(TPDate[0],TPDate[NbTP-1]) / 365.; 

    IndexMat = 0;  /* TPDate[IndexMat] = FxMatDate[k] */

    for (k = 0; k < NbVol; k++)
    {

        /* check that FxRstDate[i] <= FxMatDate[i] */
        if (FxRstDate[k] > FxMatDate[k])
        {

            DR_Error ("FX reset date is after maturity date\n");
                return(status);

        }

        /* find maturity date on time line */
        IndexMat = GetDLOffset(NbTP,
                               TPDate,
                               FxMatDate[k],
                               CbkEXACT);		

        if (IndexMat < 0)
        {

            sprintf(ErrorMsg, "Option expiry date #%ld is not"
                              " on the input time line\n",FxMatDate[k]);
            DR_Error(ErrorMsg);
            goto RETURN;

        }
        
        /* reset date = base date */
        if (FxRstDate[k] == TPDate[0])
        {

            ImpVol[k] = 0.0;
            continue; /* next k */

        }

        /* find FxRstDate[k] on time line*/        
        IndexRst = GetDLOffset(NbTP,
                               TPDate,
                               FxRstDate[k],
                               CbkEXACT);	

        if (IndexRst < 0)
        {

            sprintf(ErrorMsg, "Option expiry date #%ld is not"
                              " on the input time line\n",
                              FxRstDate[k]);
            DR_Error(ErrorMsg);
            goto RETURN;

        }

        /*   now TPDate[IndexRst] == FxRstDate[k] */
        
        if (IndexRst > IndexMat )    
        {

            DR_Error("IndexRst is strictly greater than"
                     "IndexMat\n");
            goto RETURN;

        }

        /* Some initialisations */
        varFX = 0.;

        for(l = 0; l < DomNbFac; l++)
        {

           DomA[l]      =   0.;

        }/* for i*/
        
        for(l = 0; l < ForNbFac; l++)
        {

            ForA[l]      =   0.;

        }/* for i*/

        for(l = 0; l < NbNoFxInt; l++)
        {

            NoFxInt[l]   =   0.;

        }

        for(l = 0; l < NbWithFxInt; l++)
        {
        
            WithFxInt[l] =   0.;

        }
        /******************************************/
        /*      initialisations are complete      */
        /******************************************/  
        
        for(i = IndexMat - 1; i >= IndexRst; i--)
        {
                        
            /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {

                DomA[l]   *=   exp(-BetaDom[l] * DeltaTime[i]);
                DomA[l]   +=   DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0; l < ForNbFac; l++)
            {
                ForA[l]   *=   exp(-BetaFor[l] * DeltaTime[i]);
                ForA[l]   +=   ForFwdRate[i] * Hyb3_ExpDecay(BetaFor[l],DeltaTime[i]);

            }/* for l */

        }/* for i */

        for (i = IndexRst - 1; i >= 0;i --)
        {

            CurrentNbInt     = 0;

           /* domestic A's*/
            for (l = 0; l < DomNbFac; l++)
            {

                DomA[l]   *=   exp(-BetaDom[l] * DeltaTime[i+1]);
                DomA[l]   +=   DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            
            /* foreign A's */
            for (l = 0; l < ForNbFac; l++)
            {

                ForA[l]   *=   exp(-BetaFor[l] * DeltaTime[i+1]);
                ForA[l]   +=   ForFwdRate[i] * Hyb3_ExpDecay(BetaFor[l],DeltaTime[i]);

            }/* for l */
            

            /* variance of domestic factors */
            for (l = 0; l < DomNbFac; l++)
            {   

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]    +=    AlphaDom[l] * AlphaDom[l] * 
                                                 SpotVolDom[i] * SpotVolDom[i]
                                                 * DomA[l] * DomA[l] * DeltaTime[i];

            }
            
            /* covariance between domestic factors*/
            if (DomNbFac > 1)
            {

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]    +=    2. * AlphaDom[0] * AlphaDom[1] * 
                                                 DomA[0] * DomA[1]
                                                 * SpotVolDom[i] * SpotVolDom[i] *
                                                 RhoDomFac[0]*DeltaTime[i];

            }/* end of case: DomNbFac > 1 */
            
            if (DomNbFac > 2)
            {   

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]   +=     2. * AlphaDom[0] * AlphaDom[2] 
                                                 * DomA[0]*DomA[2]
                                                 * SpotVolDom[i] * SpotVolDom[i] 
                                                 * RhoDomFac[1] * DeltaTime[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]   +=     2. * AlphaDom[1] * AlphaDom[2] 
                                                 * DomA[1] * DomA[2]
                                                 * SpotVolDom[i] * SpotVolDom[i] 
                                                 * RhoDomFac[2] * DeltaTime[i];

            }/* end of case DomNbFac > 2 */

            /* variance of foreign factors*/
            for (l = 0; l < ForNbFac; l++)
            {   

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]   +=     AlphaFor[l] * AlphaFor[l] 
                                                 * SpotVolFor[i]*SpotVolFor[i]
                                                 * ForA[l] * ForA[l]*DeltaTime[i];

            } /* for l */

            /* covariance between foreign factors*/
            if (ForNbFac > 1)
            {

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]   +=     2. * AlphaFor[0] * AlphaFor[1]
                                                 * ForA[0] * ForA[1]
                                                 * SpotVolFor[i] * SpotVolFor[i] 
                                                 * RhoForFac[0]*DeltaTime[i];

            }
            
            if (ForNbFac > 2)
            {   

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]   +=     2. * AlphaFor[0] * AlphaFor[2] 
                                                 * ForA[0] * ForA[2]
                                                 * SpotVolFor[i] * SpotVolFor[i] 
                                                 * RhoForFac[1] * DeltaTime[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]   +=     2. * AlphaFor[1] * AlphaFor[2] 
                                                 * ForA[1] * ForA[2]
                                                 * SpotVolFor[i] * SpotVolFor[i] 
                                                 * RhoForFac[2] * DeltaTime[i];

            }
            
            /* covariance between domestic and foreign factors */
            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {

                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt-1] -=   2. * AlphaDom[l] * DomA[l] *
                                                  SpotVolDom[i] *
                                                  AlphaFor[m] * ForA[m] * 
                                                  SpotVolFor[i] * 
                                                  RhoDomFacForFac[m + l * ForNbFac] 
                                                  * DeltaTime[i];

                }/* for m */

            }/* for l */
            
            /* covariance between Fx and domestic */
            for (l = 0; l < DomNbFac; l++)
            {

                WithFxInt[l]             +=   2. * SpotFxVol[i] * SpotVolDom[i]
                                              * AlphaDom[l] * DomA[l] * 
                                              RhoFxDomFac[l] * DeltaTime[i];

            }/* for l*/

            /* covariance between fx and for */
            for (l = 0; l < ForNbFac; l++)
            {

                WithFxInt[l+DomNbFac]    -=   2. * SpotFxVol[i] * SpotVolFor[i]
                                              * AlphaFor[l] * ForA[l] 
                                              * RhoFxForFac[l] * DeltaTime[i];

            }/* for l */

            /* integral of "composite" vol square over (t_i,T) */
            varFX      +=  SpotFxVol[i] * SpotFxVol[i] * DeltaTime[i];  
            volA2[i]    =  varFX; 
            
            for(l = 0; l < NbNoFxInt; l++)
            {
                volA2[i]    +=  NoFxInt[l];
            }

            for(l = 0; l < NbWithFxInt; l++)
            {
                volA2[i]    +=  WithFxInt[l];
            }

        } /* i */
        /* finished IR contribution */ 

        /* Note that volA[i] = \Sigma_A(t_i,T), where T = Time[IndexRst] */
        for (i = IndexRst - 1; i >= 0; i--)
        {
            volA[i]      = sqrt( volA2[i] / (Time[IndexRst] - Time[i]) );
        } 

        /* set up array of curvatures */
        if (Hyb3_CalcVolA_xx( 
                        volA_xx,      /* (O)                                      */
                        IndexRst,     /* (I) number of time points                */
                        SpotFxVol,    /* (I) spot fx vol                          */
                        loc_vol_x,    /* (I) slope of local vol wrt log-money     */
                        loc_vol_xx,   /* (I) curvature of local vol wrt log-money */
                        Time,         /* (I) time grid, evaluation time = 0       */
                        DeltaTime)    /* (I) intervals in grid                    */
                        == FAILURE) goto RETURN;

        volB_sum = 0.;
        for (i = IndexRst - 1; i >= 0; i--)
        {
            /* integrate result */
            /* using volA[i-1] because volA[i] corresponds to Time[i+1] */

            /* the RHS equals T - t integrated over (t_i,t_{i+1}) */
            volB_sum +=  SpotFxVol[i] * SpotFxVol[i] * volA[i] * volA_xx[i] 
                                      * 0.5 * ( - Time[i+1] * Time[i+1] + Time[i] * Time[i] 
                                      + 2.0 * Time[IndexRst] * (Time[i+1] - Time[i]) );
        }   /* for i */
        
        /* Since the implied vol approximation is unsafe for large pertubations, an error   */
        /* is returned in those cases - the user should then lower implied vol curvature.   */                                                                  
        /* Having the cutoff at 100% for vol^2 terms is somewhat arbitary                   */ 
        if( fabs(volB_sum) > 3.0 * volA2[0])
        {
            sprintf(ErrorMsg, "Volatility or its curvature is too large, bootstrapping" 
                              "fails at expiry #%ld\n",FxMatDate[k]);

            DR_Error(ErrorMsg);
            goto RETURN;
        }


        ImpVol[k] = sqrt( (volA2[0] + volB_sum) / (Time[IndexRst] - Time[0]) ); 
        
        if( ImpVol[k] < TINY)
        {
            sprintf(ErrorMsg, "Problem in calculating Implied FX Vol at"
                              "Expiry #%ld\n",FxMatDate[k]);

            DR_Error(ErrorMsg);
            goto RETURN;
        }

    }/* for k*/

    status = SUCCESS;

RETURN:

    Free_DR_Array (NoFxInt,     DOUBLE,0,NbNoFxInt - 1);
    Free_DR_Array (WithFxInt,   DOUBLE,0,NbWithFxInt - 1);
    Free_DR_Array (DomA,        DOUBLE,0,DomNbFac - 1);
    Free_DR_Array (ForA,        DOUBLE,0,ForNbFac - 1);
    Free_DR_Array (DeltaTime,   DOUBLE,0,NbTP - 1);
    Free_DR_Array (Time,        DOUBLE,0,NbTP - 1);
    Free_DR_Array (a,           DOUBLE,0,NbTP - 1);
    Free_DR_Array (b,           DOUBLE,0,NbTP - 1);
    Free_DR_Array (c,           DOUBLE,0,NbTP - 1);
    Free_DR_Array (volA,        DOUBLE,0,NbTP - 1);
    Free_DR_Array (volA2,       DOUBLE,0,NbTP - 1);
    Free_DR_Array (volA_x,      DOUBLE,0,NbTP - 1);
    Free_DR_Array (volA_xx,     DOUBLE,0,NbTP - 1);
    Free_DR_Array (loc_vol_x,   DOUBLE,0,NbTP - 1);
    Free_DR_Array (loc_vol_xx,  DOUBLE,0,NbTP - 1);

    return (status);

}/* Multifac_FxVol2*/


/*************** Hyb3_Get_FxSpotVol ****************************************/
/**This function:
    1- bootstraps IR vols assuming NORMAL smile parameters
    2- bootstraps FX composite vols
    3- appends the Input Spot vols carried in fx_data structure
    4- Handles the FX_CONSTANT_SPOT_VOL bootstrap mode directly
       
  
   On OUTPUT: SpotFxVol contains (fx_data.NbVol + fx_data.NbInpSpotVol)
              BACKWARD looking spot vols.
              The (implicit) expiries along side these spot vols are as follows
              The first fx_data.NbVol SpotFxVols correspond to 
              composite fx_data.VolDates, the next block corresponds
              to fx_data.InpSpotVolDates.
 **********************************************************************/
int Hyb3_Get_FxSpotVol(
            T_CURVE         *DomCurve,       /**<(I)*/
            T_CURVE         *ForCurve,       /**<(I)*/
            MKTVOL_DATA     *DomMktVol,      /**<(I)*/
            MKTVOL_DATA     *ForMktVol,      /**<(I)*/
            FX_DATA         *fx_data,        /**<(I)*/
                       
            /* Ppy for creating time schedule*/
            int             Ppy,
            double          *SpotFxVol)     /*(O)*/
{

    int     i,status = FAILURE;
    double
            *DomSpotVol = NULL,
            *ForSpotVol = NULL;
    long    
            *TPDate = NULL;
    double  
            *DomFwd = NULL,
            *ForFwd = NULL;

    double  alpha;  /* this is Alpha[0] */
    int     
            NbTP = 0;
    double 
           test_mtx [3][3]; /* triangulation matrices*/

    double test_rho[3];     /* used to test validity of correlation inputs*/
   
   
    
    /* Some basic checks first */
   
    test_rho[0] = fx_data->Rho[0][0];
    test_rho[1] = fx_data->Rho[1][0];
    test_rho[2] = fx_data->Rho[2][0];
    
    if (Triangulation (test_mtx,3,test_rho) == FAILURE)
    {
        DR_Error ("Invalid input correlation structure"
                    "between Fx rate and swap rates in Get_SpotFxVol\n");
        return (status);
    }
      

    if (DomCurve == NULL ||
        ForCurve == NULL ||
        DomMktVol == NULL ||
        ForMktVol == NULL ||
        fx_data == NULL ||        
        SpotFxVol == NULL )
    {
        DR_Error("Invalid pointer inputs to Get_SpotFxVol\n");
        return (status);
    }
    
    if ((fx_data->NbVol < 0) || (fx_data->NbInpSpotVol < 0))
    {
        DR_Error("number of vols must be >= 0 \n");
        goto RETURN;
    }
    if (fx_data->NbVol + fx_data->NbInpSpotVol <= 0)
    {
        DR_Error("invalid total number of VolDates (composite +spot) "
                 "in Hyb3_Get_FxSpotVol\n");

        goto RETURN;
    }

    /* deal with nil calibration first */
    if (fx_data->FxBootStrapMode == FX_CONSTANT_SPOT_VOL)
    {
        long l;
        for (l = 0; l < (fx_data->NbVol + fx_data->NbInpSpotVol); l++)
        {
            SpotFxVol[l] = fx_data->FxCutOffLevel;
        }
        return (SUCCESS);
    }
    
    if((fx_data->NbVol > 0) && (fx_data->NbInpSpotVol > 0))
    {
        if(fx_data->VolDate[fx_data->NbVol] >= fx_data->InpSpotVolDate[0])
        {
            DR_Error("First input spotvol date must be > last "
                     " composite voldate in Hyb3_Get_FxSpotVol\n");
            goto RETURN;
        }
    }
    
    if (fx_data->NbInpSpotVol > 0)
    {
        if(fx_data->InpSpotVolDate[0] <= fx_data->ValueDate)
        {
            DR_Error("First fx spot vol date must be > fx value date \n");
            goto RETURN;
        }
    }
    /*********************************************************/
    /* if there are any composite vols to be bootstraped, do */
    /* all the necessary work.                               */
    /*********************************************************/
    if (fx_data->NbVol > 0)
    {
        /* Create time line */

        /* first call to MergeDateLists with TPDate = NULL and NbTP = 0 */
        if (MergeDateLists(
                            DomMktVol->NbVol,
                            DomMktVol->VolDate,
                            &(NbTP),
                            &(TPDate)) == FAILURE)  
        {
            DR_Error("MergeDateLists failed in Get_SpotFxVol for " 
                        "DomVol dates\n");
            goto RETURN;
        }

        if (MergeDateLists(
                            ForMktVol->NbVol,
                            ForMktVol->VolDate,
                            &(NbTP),
                            &(TPDate)) == FAILURE)
        {
            DR_Error("MergeDateLists failed in Get_SpotFxVol for "
                     "ForVol dates\n");
            goto RETURN;
        }
        /* check first fx VolDate*/
        if (fx_data->VolDate[0] != fx_data->Today)
        {
            DR_Error("the first Fx VolDate should be fx_data.Today\n");
            goto RETURN;
        }
        if (MergeDateLists((fx_data->NbVol) + 1,
                            fx_data->VolDate,
                            &NbTP,
                            &TPDate) == FAILURE)
        {
            DR_Error ("MergeDateLists failed in Get_SpotFxVol for FxVol dates\n");
            goto RETURN;
        }


        /* Add Fx base to the TPDate list */
        if (MergeDateLists(
                            1,
                            &(fx_data->Today),
                            &NbTP,
                            &TPDate) == FAILURE)
        {
            DR_Error("MergedDateLists failed in Get_SpotfxVol for fx base date\n");
            goto RETURN;
        }
    
        /*note: DrExtendAmerSch removes all the dates Strictly before Base date*/
        if (DrExtendAmerSch(fx_data->Today,
                            Ppy,
                            &NbTP,
                            &TPDate,
                            NULL,
                            'L',
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            NULL) == FAILURE)
        {
            DR_Error ("DrExtendAmerSch failed in Get_SpotFxVol\n");
            goto RETURN;
        }
    
        if (TPDate[0] != fx_data->Today)
        {
            DR_Error ("TPDate[0] is not fx_data.Today (Hyb3_Get_FxSpotVol)\n");
            goto RETURN;
        }
        /* generate domestic and foreign spot vols */
        /* using NORMAL smile parameters           */
        if (Hyb3_SpotVol(    DomMktVol->Aweight,
                        DomMktVol->BaseDate,
                        DomMktVol->NbVol,
                        DomMktVol->VolDate,
                        DomMktVol->Vol,
                        DomMktVol->VolUsed,
                        DomMktVol->Freq,
                        DomMktVol->DCC,
                        DomMktVol->SwapSt,
                        DomMktVol->SwapMat,
                        0., /* "Normal" QLeft    */                         
                        0., /* "Normal" QRight   */                   
                        0., /* No forward shift  */                    
                        1,  /* 1 Factor          */
                        DomMktVol->Alpha,
                        DomMktVol->Beta,
                        DomMktVol->Rho,
                        DomMktVol->SkipFlag,
                        DomMktVol->CalibFlag,
                        DomCurve) != SUCCESS)
        {
            DR_Error("Hyb3_SpotVol for domestic failed in Get_SpotFxVol\n");
            goto RETURN;
        }

        if (Hyb3_SpotVol(    ForMktVol->Aweight,
                        ForMktVol->BaseDate,
                        ForMktVol->NbVol,
                        ForMktVol->VolDate,
                        ForMktVol->Vol,
                        ForMktVol->VolUsed,
                        ForMktVol->Freq,
                        ForMktVol->DCC,
                        ForMktVol->SwapSt,
                        ForMktVol->SwapMat,
                        0.,  /* "Normal" QLeft    */                  
                        0.,  /* "Normal" QRight   */                  
                        0.,  /* No forward shift  */                  
                        1,  /* 1 Factor          */
                        ForMktVol->Alpha,
                        ForMktVol->Beta,
                        ForMktVol->Rho,
                        ForMktVol->SkipFlag,
                        ForMktVol->CalibFlag,
                        ForCurve) != SUCCESS)
        {
            DR_Error("Hyb3_SpotVol for foreign failed in Get_SpotFxVol\n");
            goto RETURN;
        }
    
        /* Discretise spot vol onto time line */
        DomSpotVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    
        if (DomSpotVol == NULL)
        {
            DR_Error ("Invalid memory allocation for domestic"
                      "vol in Get_SpotFxVol\n");
            goto RETURN;
        }
    
        ForSpotVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);
    
        if (ForSpotVol == NULL)
        {
            DR_Error ("Invalid memory allocation for foreign vol"
                        "in Get_SpotFxVol\n");
            goto RETURN;
        }
    
        /* The spot Vol needed for the rest of the function
         * is the Spot Vol of the EXPONENTIAL factors,as opposed to the 
         * orthogonal the orthogonal ones produced by triangulation,
         * therefore use only Aweight[0]                    */

        if (Hyb3_DiscretiseInput(DomSpotVol,
                            DomMktVol->Aweight[0],
                            DomMktVol->CalibFlag,
                            DomMktVol->NbVol,
                            DomMktVol->VolDate,
                            NbTP,
                            TPDate) == FAILURE)
        {
            DR_Error("Hyb3_DiscretiseInput failed for domestic"
                      "Spot vol in Get_Spot_fxVol\n");
            goto RETURN;
        }
    
        /* spot vol input to Hyb3_MultiFac_FxVol should be "weight free".
         * Since Aweight[0]=spotvol*Alpha[0], we have to rescale */

        alpha = DomMktVol->Alpha[0];

        /* Hyb3_Param_Check, ensures that Alpha[0] non zero*/

        for (i = 0; i < NbTP - 1;i++)
        {
            DomSpotVol[i] /= alpha; 
        }

        if (Hyb3_DiscretiseInput(ForSpotVol,
                            ForMktVol->Aweight[0],
                            ForMktVol->CalibFlag,
                            ForMktVol->NbVol,
                            ForMktVol->VolDate,
                            NbTP,
                            TPDate) == FAILURE)
        {
            DR_Error("Hyb3_DiscretiseInput failed for foreign"
                     "Spot vol in Get_Spot_fxVol\n");
            goto RETURN;
        }
                        
        alpha = ForMktVol->Alpha[0];
        for (i = 0; i < NbTP - 1; i++)
        {
            ForSpotVol[i] /= alpha;
        }

        /* obtain appropriate forward rates*/
        DomFwd = (double*) DR_Array (DOUBLE,0,NbTP - 2);
        ForFwd = (double*) DR_Array (DOUBLE,0,NbTP - 2);
    
        if (DomFwd == NULL || ForFwd == NULL)
        {
            DR_Error ("Failed to allocate memory in Get_SpotFxVol\n");
            goto RETURN;
        }

        if (Hyb3_SimpleFwdCurve( DomFwd,
                            DomCurve,
                            TPDate,
                            NbTP) == FAILURE)
        {
            DR_Error ("Hyb3_SimpleFwdCurve failed for domestic"
                     "curve in Get_SpotFxVol\n");
            goto RETURN;
        }
    
        if (Hyb3_SimpleFwdCurve( ForFwd,
                            ForCurve,
                            TPDate,
                            NbTP) == FAILURE)
        {
            DR_Error ("Hyb3_SimpleFwdCurve failed for foreign"
                        "curve in Get_SpotFxVol\n");
            goto RETURN;
        }
        /***********************************************************/
        /* we can now bootstrapp fxvol!                            */
        /* CAREFUL: the function requires the input VolDates to be */
        /* STRICTLY GREATER than TPDate[0] (= fx_data.Today),      */
        /* so we need to skip                                      */
        /* the first VolDate, and decrease the Nb of vols by 1     */
        /***********************************************************/

        if (Hyb3_MultiFac_Spot_FxVol(SpotFxVol,
                                DomSpotVol,
                                ForSpotVol,
                                DomMktVol->Alpha,
                                ForMktVol->Alpha,
                                DomFwd,
                                ForFwd,
                                &(fx_data->Rho[2][0]),/* Rho fx/dom*/
                                &(fx_data->Rho[1][0]),/* Rho Fx/For*/
                                &(fx_data->Rho[0][0]),/* Rho Dom/For*/
                                DomMktVol->Rho,
                                ForMktVol->Rho,
                                fx_data->NbVol,
                                &(fx_data->VolDate[1]),
                                &(fx_data->VolDate[1]),
                                &(fx_data->FxVol[1]),
                                fx_data->FxCutOffFlag,
                                fx_data->FxCutOffLast,
                                fx_data->FxCutOffLevel,
                                DomMktVol->Beta,
                                ForMktVol->Beta,
                                1,/* 1 factor for dom*/
                                1,/* 1 factor for For*/
                                TPDate,
                                NbTP) == FAILURE)
        {
            DR_Error ("Hyb3_MultiFac_Spot_FxVol failed \n");
            goto RETURN;
        }
    }/* if fx_data->NbVol > 0 */
    
    /* Append any input spot vols*/
    if(fx_data->NbInpSpotVol > 0)
    {
        long    l;
        long    TotNbVol;
        
        TotNbVol = fx_data->NbVol + fx_data->NbInpSpotVol;

        for (l = fx_data->NbVol ; l < TotNbVol; l++)
        {
            SpotFxVol[l] = fx_data->InpSpotVol[l - fx_data->NbVol];
        }
    }


    status = SUCCESS;
    
RETURN:
    if (NbTP > 0)
    {
        Free_DR_Array (TPDate,LONG,0,NbTP - 1);    
    }

    if (NbTP >= 2)
    {
        Free_DR_Array (DomSpotVol,DOUBLE,0,NbTP - 2);
        Free_DR_Array (ForSpotVol, DOUBLE,0,NbTP - 2);
        Free_DR_Array (DomFwd,DOUBLE,0,NbTP - 2);
        Free_DR_Array (ForFwd,DOUBLE,0,NbTP - 2);
    }
    
    return (status);
    
}/* Hyb3_Get_FxSpotVol */


/********** MultiFac_Spot_FXVol *********************************************/
/**Outputs as many spot fx vols as there are Vol Expiry dates (i.e NbVol).
   The correlation factors are with respect to the exponential factors
  
   NOTES:
   1-Assumes that TPDate[0]=fx_data.Today.
   2-Assumes and checks that each VolDate and VolIntegrationDate is on the time 
   line.
   3-VolDate and VolIntegrationDate must be entered in a strictly ascending 
   order.
   
 *****************************************************************************/

int Hyb3_MultiFac_Spot_FxVol(
        double  *SpotFxVol,         /* (O) Instant.fx vol                */
        double  *DomSpotVol,        /* (I) Instant dom. int. rate vol    */
        double  *ForSpotVol,        /* (I) Instant for. int.rate vol     */
        double  *DomAlpha,          /* (I) dom. relative factor weights  */
        double  *ForAlpha,          /* (I) for. relative factor weights  */
        double  *DomFwdRate,        /* (I) dom. FwdRate @ each time pt   */
        double  *ForFwdRate,        /* (I) for. FwdRate @ each time pt   */
        double  *RhoFxDom,          /* (I) correl fx/dom. curve          */
        double  *RhoFxFor,          /* (I) correl fx/for. curve          */
        double  *RhoDomFor,         /* (I) correl dom./for. curves       */
        double  *RhoDom,            /* (I) correl between dom. factors   */
        double  *RhoFor,            /* (I) correl between for. factors   */
        int     NbVol,              /* (I) Nb of fx vol points           */
        long    *VolDate,           /* (I) fwd fx maturity               */
        long    *VolIntegrationDate,/* (I) expiry of option              */
        double  *FxVol,             /* (I) fx vol for option which       */
                                    /*     expirs on VolIntegrationDate, */
                                    /*     the underlying fwd matures on */
                                    /*     VolDate                       */
        int     FxCutOffFlag,       /* (I) True if cut off allowed       */
        int     FxCutOffLast,       /* (I) True if c/off at last level   */
        double  FxCutOffLevel,      /* (I) User-given or from cups.h     */
        double  *DomBeta,           /* (I) Domestic mean reversion       */
        double  *ForBeta,           /* (I) Foreign mean reversion        */
        int     DomNbFac,           /* (I) Nb of fac. for dom. curve     */
        int     ForNbFac,           /* (I) Nb of facc for for.curve      */
        long    *TPDate,            /* (I) date of each pt.              */
        int     NbTP)               /* (I) Total Nb of time steps        */
{
    double  a[3],               /* coeffs of second order polynomial        */
            *DomA = NULL,       /* A factors for dom.                       */
            *ForA = NULL,       /* A factors for for.                       */
            *NoFxInt = NULL,    /* integrals not containing the FX rate     */
            *WithFxInt = NULL,  /* integrals due to covar. FX with dom/for  */
            *Length = NULL;     /* length of timesteps                      */
        
    double  *TPSpotFxVol = NULL,    /* values of spot FX @ each TP          */
            Sigma = 0.0,            /* Fx Spot Vol in current bucket        */
            discriminant;           /* discriminant of sec. order polynomial*/
            
    int     kVol,status=FAILURE,
            TPIntBefore,           /*integration time index of prev.Vol point*/
            TPIntSofar,            /*integration time index of curr.Vol point*/
            TPMatSofar,            /*time index of curr. fx mat date         */
            VolTooLow,             /*Flag for spot vol criterion             */
            i,l,m,n,
            CurrentNbInt,   /*used to index integrals not involving Fx rate */
            NbNoFxInt,      /*Total Nb of integrals not containing  FX rate */
            NbWithFxInt,    /*Total Nb of int. due to covar. Fx with For/Dom*/
            found;          /*to test if Vol Exp date is on the time line   */

    char    ErrorMsg[MAXBUFF];  /*Error message to be sent to DR_Error      */


    /* Quick checks first */

    if (DomNbFac != 1 && DomNbFac != 2 && DomNbFac != 3)   
    {
        DR_Error("Invalid Number of domestic yield curve factors "
                 "in Hyb3_MultiFac_Spot_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if (ForNbFac != 1 && ForNbFac != 2 && ForNbFac != 3)   
    {
        DR_Error("Invalid Number of foreign yield curve factors "
                "in Hyb3_MultiFac_Spot_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }

    if (NbTP < 2)
    {
        DR_Error ("Invalid number of time points in Hyb3_MultiFac_Spot_FxVol: "
                    "should be at least 2\n");
        return (status);
    }

    if (NbVol < 1)
    {
        DR_Error ("Invalid number of Vol points in MultiFac_Spot_fxVol: "
                    "should be at least 1");
        return (status);
    }
    
    if (RhoDom == NULL && DomNbFac > 1)
    {
        DR_Error("Dom factor correlation array cannot be NULL"
                " when more than 1 factor\n");
        return(status);
    }
    if (RhoFor == NULL && ForNbFac > 1)
    {
        DR_Error("Foreign factor correlation array cannot be NULL"
                " when more than 1 factor\n");
        return(status);
    }

    if( DomSpotVol == NULL ||
        ForSpotVol == NULL ||
        DomAlpha   == NULL ||
        ForAlpha   == NULL ||
        DomFwdRate == NULL ||
        ForFwdRate == NULL ||
        RhoFxDom   == NULL ||
        RhoFxFor   == NULL ||
        RhoDomFor  == NULL ||
        VolDate    == NULL ||
        VolIntegrationDate == NULL ||
        FxVol      == NULL ||
        DomBeta    == NULL ||
        ForBeta    == NULL ||
        TPDate     == NULL)
    {   
        DR_Error("Invalid pointer inputs to MultiFac_SpotFxVol\n");
        return (status);
    }
    
    NbNoFxInt = (DomNbFac + ForNbFac) + 
                (DomNbFac + ForNbFac - 1) * (DomNbFac + ForNbFac) / 2;
    NbWithFxInt = (DomNbFac + ForNbFac);

    DomA =        (double*) DR_Array (DOUBLE,0,DomNbFac - 1);
    ForA =        (double*) DR_Array (DOUBLE,0,ForNbFac - 1);
    NoFxInt =     (double*) DR_Array (DOUBLE,0,NbNoFxInt - 1);
    WithFxInt =   (double*) DR_Array (DOUBLE,0,NbWithFxInt - 1);
    TPSpotFxVol = (double*) DR_Array (DOUBLE,0,NbTP - 2); 
    Length=       (double*) DR_Array (DOUBLE,0,NbTP - 2);
    
    if (DomA        == NULL ||
        ForA        == NULL ||
        NoFxInt     == NULL ||
        WithFxInt   == NULL ||
        TPSpotFxVol == NULL ||
        Length      == NULL)
    {
        goto RETURN;
    }

    
    for (i = 0 ; i < NbTP - 1;i++)
    {   
        if ( TPDate[i+1] <= TPDate[i])
        {
            sprintf(ErrorMsg,"TPDate[%d] and TPDate[%d] in input time line"
                      "are not in a strictly ascending order\n",i,i+1);
            DR_Error (ErrorMsg);
            return (status);
        }
    
        Length[i]=Daysact(TPDate[i],TPDate[i+1]) / 365.;
    }

    TPIntBefore = TPIntSofar = 0;
    TPMatSofar = 0;

    for (kVol = 0; kVol < NbVol; kVol++)
    {   
        /* check that IntegrationDate and VolDate are in the right order*/
        if (VolIntegrationDate[kVol] > VolDate[kVol])
        {
            sprintf(ErrorMsg, "Vol Integration date %ld is after"
                    "corresponding vol mat date (MultiFac_Spot_FxVol1)\n",
                    VolIntegrationDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
        if ( VolIntegrationDate[kVol] < TPDate[0])
        {
            DR_Error("VolIntegrationdate cannot be before TPDate[0].\n");
            goto RETURN;
        }
        /* make sure VolIntegrationDate is on the input time line*/

        found = 0;  /* initialise as "not on the time line"*/

        for ( ; TPIntSofar < NbTP; TPIntSofar++)
        {   
            found = (TPDate[TPIntSofar] == VolIntegrationDate[kVol]);

            if (found)
            {
                break;
            }
        }/*  for index */

        if (!found)
        {
            sprintf(ErrorMsg, "Vol Integration date %ld is not on the"
                    "input time line (or dates are not in ascending order)"
                    "(Hyb3_MultiFac_Spot_FxVol)\n",VolIntegrationDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
        /* check that Voldate is on input time line         */
        found = 0;  /* re-initialise as "not on the time line"*/

        for (TPMatSofar = TPIntSofar ; TPMatSofar < NbTP; TPMatSofar++)
        {   
            found = (TPDate[TPMatSofar] == VolDate[kVol]);

            if (found)
            {
                break;
            }
        }/*  for index */

        if (!found)
        {
            sprintf(ErrorMsg, "Vol date %ld is not on the"
                    "input time line (Hyb3_MultiFac_Spot_FxVol)\n",VolDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* the following should not normally happen as fx   */ 
        /* dates are input using New_Fx_Input_W,            */
        /* which checks for strictly ascending order and    */
        /* eliminates dates before FxBaseDate=TPDate[0]     */

        if (TPIntSofar == TPIntBefore)
        {
          continue;
        }
        
        /* first,some initialisations*/
        
        for (i = 0; i < DomNbFac;i++)
        {
            DomA[i] = 0.;

        } /* for i*/
        
        for (i = 0; i < ForNbFac;i++)
        {
            ForA[i] = 0.;

        } /* for i*/

        for (i = 0; i < NbNoFxInt; i++)
        {
            NoFxInt[i] = 0.;

        } /* for i*/
        
        for (i = 0; i < NbWithFxInt; i++)
        {
            WithFxInt[i] = 0.;

        } /* for i*/

        for (i = 0; i < 3; i++)
        {
            a[i] = 0.;

        } /* for i*/
        
        /* prepare A-factors    */

        for (i = TPMatSofar - 1; i>= TPIntSofar; i--)
        {
            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-DomBeta[l] * Length[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(DomBeta[l],Length[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0;l < ForNbFac;l++)
            {
                ForA[l] *= exp(-ForBeta[l] * Length[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(ForBeta[l],Length[i]);

            }/* for l */


        }/* for i */
        
        /* we can start integrating now */
        for (i = TPIntSofar - 1; i >= TPIntBefore; i--)
        {
            CurrentNbInt = -1;
            
            /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-DomBeta[l] * Length[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(DomBeta[l],Length[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0;l < ForNbFac;l++)
            {
                ForA[l] *= exp(-ForBeta[l] * Length[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(ForBeta[l],Length[i]);

            }/* for l */

            
            /* variance of domestic factors*/
            
            for (l = 0;l < DomNbFac;l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += DomAlpha[l] * DomAlpha[l] * 
                                         DomSpotVol[i] * DomSpotVol[i]
                                        * DomA[l] * DomA[l] * Length[i];
            }
            
            /* covariance between domestic factors*/

            if (DomNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* DomAlpha[0] * DomAlpha[1] * 
                                         DomA[0] * DomA[1]
                                        * DomSpotVol[i] * DomSpotVol[i] *
                                        RhoDom[0]*Length[i];

            }/* end of case: DomNbFac > 1 */
            
            if (DomNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * DomAlpha[0] * DomAlpha[2] *
                                         DomA[0] * DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] * 
                                        RhoDom[1] * Length[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * DomAlpha[1] * DomAlpha[2] * 
                                         DomA[1] * DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] * 
                                        RhoDom[2] * Length[i];

            }/* end of case DomNbFac > 2 */

            /* variance of foreign factors*/
            
            for (l = 0; l < ForNbFac;l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += ForAlpha[l] * ForAlpha[l] * 
                                         ForSpotVol[i] * ForSpotVol[i]
                                        * ForA[l] * ForA[l]*Length[i];

            } /* for l */

            /* covariance between foreign factors*/

            if (ForNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* ForAlpha[0] * ForAlpha[1] *
                                            ForA[0] * ForA[1]
                                            * ForSpotVol[i] * ForSpotVol[i] 
                                            * RhoFor[0]*Length[i];
            }
            
            if (ForNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* ForAlpha[0] * ForAlpha[2] * 
                                         ForA[0]*ForA[2]
                                        * ForSpotVol[i] * ForSpotVol[i] *
                                        RhoFor[1] * Length[i];
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * ForAlpha[1] * ForAlpha[2] * 
                                         ForA[1] * ForA[2]
                                         * ForSpotVol[i] * ForSpotVol[i] * 
                                         RhoFor[2] * Length[i];
            }
            
            /* covariance between domestic and foreign factors */

            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {
                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt] -= 
                                2.0 * DomAlpha[l] * DomA[l] * DomSpotVol[i] *
                                ForAlpha[m] * ForA[m] * ForSpotVol[i] * 
                                RhoDomFor[m + l * ForNbFac]
                                * Length[i];
                }/* for m */

            }/* for l */
            
            /*covariance between Fx and domestic*/

            for (l = 0; l < DomNbFac; l++)
            {
                WithFxInt[l] += 2.0 * DomSpotVol[i] * DomAlpha[l] * 
                                    DomA[l] * RhoFxDom[l]*Length[i];

            }/* for l*/


            /*covariance between fx and for*/

            for (l = 0; l < ForNbFac; l++)
            {
                WithFxInt[l+DomNbFac] -= 2.0 * ForSpotVol[i] * 
                                        ForAlpha[l] * ForA[l] * 
                                        RhoFxFor[l]*Length[i];

            }/* for l*/

        }   /* for i*/
        
        for (i = TPIntBefore - 1; i >= 0; i--)
        {               
            CurrentNbInt = -1;

          /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {
                DomA[l] *= exp(-DomBeta[l] * Length[i]);
                DomA[l] += DomFwdRate[i] * Hyb3_ExpDecay(DomBeta[l],Length[i]);

            }/* for l */
            
            /* foreign A's */

            for (l = 0; l < ForNbFac; l++)
            {
                ForA[l] *= exp(-ForBeta[l] * Length[i]);
                ForA[l] += ForFwdRate[i] * Hyb3_ExpDecay(ForBeta[l],Length[i]);

            }/* for l */  
                        
            /* variance of each domestic factor*/
            
            for (l = 0; l < DomNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += DomAlpha[l] * DomAlpha[l] * 
                                         DomSpotVol[i] * DomSpotVol[i]
                                        * DomA[l] * DomA[l] * Length[i];
            }
            
            /* covariance between domestic factors*/

            if (DomNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* DomAlpha[0] * DomAlpha[1] * 
                                            DomA[0] * DomA[1]
                                            * DomSpotVol[i] * DomSpotVol[i] * 
                                            RhoDom[0]*Length[i];
            }
            
            if (DomNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * DomAlpha[0] * DomAlpha[2] * 
                                         DomA[0]*DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] * 
                                        RhoDom[1] * Length[i];

                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * DomAlpha[1] * DomAlpha[2] * 
                                        DomA[1] * DomA[2]
                                        * DomSpotVol[i] * DomSpotVol[i] * 
                                        RhoDom[2] * Length[i];

            }/* for case d_Nbfcator > 2*/

            /* variance of foreign factors*/
            
            for (l = 0; l < ForNbFac; l++)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += ForAlpha[l] * ForAlpha[l] * 
                                            ForSpotVol[i]*ForSpotVol[i]
                                            * ForA[l] * ForA[l]*Length[i];

            } /* for l */
            
            /* covariance between foreign factors*/

            if (ForNbFac > 1)
            {
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.* ForAlpha[0] * ForAlpha[1] * 
                                         ForA[0] * ForA[1]
                                         * ForSpotVol[i] * ForSpotVol[i] *
                                         RhoFor[0]*Length[i];

            }
            
            if (ForNbFac > 2)
            {   
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2. * ForAlpha[0] * ForAlpha[2] *
                                         ForA[0]*ForA[2]
                                        * ForSpotVol[i] * ForSpotVol[i] * 
                                        RhoFor[1] * Length[i];
                CurrentNbInt++;
                NoFxInt[CurrentNbInt] += 2.0 * ForAlpha[1] * ForAlpha[2] * 
                                         ForA[1] * ForA[2]
                                         * ForSpotVol[i] * ForSpotVol[i] *
                                         RhoFor[2] * Length[i];
            }
           
            /* covariance between domestic and foreign factors */

            for (l = 0; l < DomNbFac; l++)
            {
                for (m = 0; m < ForNbFac; m++)
                {
                    CurrentNbInt++;
                    NoFxInt[CurrentNbInt] -= 2.0 * DomAlpha[l] * DomA[l] * 
                                       DomSpotVol[i] * ForSpotVol[i] *
                                       ForAlpha[m] * ForA[m] *
                                       RhoDomFor[m + l * ForNbFac] * Length[i];
                }/* for m */

            }/* for l */
                        
            /* constant coefficients due to covar fx/dom*/
            
            for (l = 0; l < DomNbFac; l++)
            {
                a[0] += 2.0 * DomSpotVol[i] * DomAlpha[l] * DomA[l] *
                                            TPSpotFxVol[i] 
                                           * RhoFxDom[l] * Length[i];
            }

            /*constant coeffs due to covar fx/for*/

            for (l = 0; l < ForNbFac; l++)
            {
                a[0] -= 2.0 * ForSpotVol[i] * ForAlpha[l] * 
                                            ForA[l] * TPSpotFxVol[i]
                                            * RhoFxFor[l] * Length[i];
            }

            /* constant coeffs due to fx variance*/

            a[0] += TPSpotFxVol[i] * TPSpotFxVol[i] * Length[i];

        }/* for i */

        for (l = 0; l < NbNoFxInt; l++)
        {
            a[0] += NoFxInt[l];
        }

        /* recall that we assumed TPDate[0] = FxBaseDate */
        
        
        a[0] -= FxVol[kVol] * FxVol[kVol] * 
                            Daysact(TPDate[0],VolIntegrationDate[kVol]) / 365.;

        for (l = 0; l < NbWithFxInt; l++)
        {
            a[1] += WithFxInt[l];
        }
        

        /* notice that TPDate[TPIntBefore] < VolIntegrationDate[kVol]   */
        /* i.e (a[2] non zero)                                          */
        /* because of earlier check on TPIntBefore==TPIntSofar          */

        a[2] = (double) Daysact(TPDate[TPIntBefore],VolIntegrationDate[kVol])
                        / 365.;

        discriminant = a[1] * a[1] - 4.* a[2] * a[0];

        VolTooLow = (discriminant < TINY);
        
        if (discriminant >= TINY)
        {            
            Sigma  = (-a[1] + sqrt(discriminant)) / 2. / a[2];

            if ((-a[1] - sqrt(discriminant)) >= TINY)
            {
                sprintf(ErrorMsg,"Unable to uniquely bootstrap spot fx vol "
                        "for expiry date = %ld \n",
                        VolDate[kVol]);
                DR_Error(ErrorMsg);
                goto RETURN;

            }
                            
            /* Evaluate whether spotvol satisfies appropriate criterion */

            if((FxCutOffLast) || (!FxCutOffFlag))
            {
                VolTooLow  = (FxVol[kVol] > FxCutOffLevel * Sigma);
            }
            else
            {   
                DR_Error("Incorrect FxCutOff mode \n");
                goto RETURN;
                /*
                VolTooLow = (Sigma < FxCutOffLevel);
                */
            }
        }
        /* now implement cutoff actions as the case may be */

        if (VolTooLow)
        {
            if (FxCutOffFlag)
            {
                if (FxCutOffLast)
                {
                    if(kVol == 0)
                    {
                        DR_Error("Unable to bootstrap at least one FX\n"
                                "vol point and therefore unable to\n"
                                "cut off at last spot vol level.\n");

                        goto RETURN;
                    }
                    else
                    {   
                        for (n = kVol; n < NbVol; n++)
                        {
                            SpotFxVol[n] = SpotFxVol[kVol-1];
                        }
                    }
                }/* end of FxCutOffLast=true */

                else
                {   
                    DR_Error("Incorrect Fx CutOff mode \n");
                    goto RETURN;
                    /*
                    for(n = kVol; n < NbVol; n++)
                    {                            
                        SpotFxVol[n] = FxCutOffLevel;
                    }
                    */
                }

                break;

            }/* end of FxCutOffFlag */

            else
            {
                sprintf(ErrorMsg,"Problem in calculating fx vol at %ld:less"
                    " than %4.2lf%% \n(level determined by ratio to fwd vol)"
                    ".\n(MultiFac_FX_Spot_Vol)\n",
                    VolDate[kVol],100.*FxVol[kVol]/FxCutOffLevel);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        } /* end of VolTooLow */

        SpotFxVol[kVol] = Sigma;

        for (i = TPIntBefore; i < TPIntSofar; i++)
        {
            TPSpotFxVol[i] = Sigma;

        }/* for i*/

        TPIntBefore  = TPIntSofar;

    
    }/*for kVol*/
    
    status = SUCCESS;

RETURN:

    Free_DR_Array (DomA,        DOUBLE,0,DomNbFac - 1);
    Free_DR_Array (ForA,        DOUBLE,0,ForNbFac - 1);
    Free_DR_Array (NoFxInt,     DOUBLE,0,NbNoFxInt-1);
    Free_DR_Array (WithFxInt,   DOUBLE,0,NbWithFxInt - 1);
    Free_DR_Array (TPSpotFxVol, DOUBLE,0,NbTP - 2);
    Free_DR_Array (Length,      DOUBLE,0,NbTP - 2);

    return(status);

}/* MultiFac_SpotFxVol*/

