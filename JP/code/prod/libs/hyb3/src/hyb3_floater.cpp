/****************************************************************************/
/*      Calculation of floater price in the lattice.                        */
/****************************************************************************/
/*      FLOATER.c                                                           */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"




/*****  Hyb3_Floater_Price  *******************************************************/
/*
*       Floater price with arbitrary index.
*/
int   	Hyb3_Floater_Price(
            TSLICE         FloaterPtr,    /* (I/O) Floater                   */
            int            IndexDim,      /* (I) Dimension of the Index slice*/
            TSLICE         IndexPtr,      /* (I) Floating index              */
            TSLICE         ZeroPtr,       /* (I) Zero maturing on nxt pmt    */
            long           FixingFlag,    /* (I) TRUE if index refixes       */
            double         DayCntFtn,     /* (I) Day count fraction          */
            double         FloatSpd,      /* (I) Spread                      */
            char	       CoS,           /* (I) 'C'ompounded, 'S'imple rate */
            char           Arrears,       /* (I) 'Y' if reset in arrears     */
            double	       Outstanding,   /* (I) Outstanding notional        */
            double         Principal,     /* (I) Principal payment           */
            long           PrincipalFlag, /* (I) TRUE if principal paid      */
            int            t,             /* (I) Current time period         */
            int            T,             /* (I) Total number of time points */
            int            DCurve,        /* (I) Curve for discounting(0,1,2)*/
            int            DMode,         /* (I) Dim of slices and of DEV    */
            HYB3_DEV_DATA      *dev_data,      /* (I) Data kept for DEV'ing       */
            HYB3_TREE_DATA     *tree_data)     /* (I) Structure of tree data      */
{


        double   IndexAtNode;

        double	*Floater = NULL;
        double	*Zero    = NULL;
	double  *IndexAtNodeSlice;

        int      offset;

        int
                Top1,      Bottom1,   /* Limits of the tree (1rst dimension) */
               *Top2,     *Bottom2,
              **Top3,    **Bottom3,

                i, j, k,              /* Node indices                        */
                status = FAILURE;     /* Error status = FAILURE initially    */

        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];	                


    /* Discounted expected value function */        
    if (Hyb3_Dev(FloaterPtr,         
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) != SUCCESS)
    {
        goto RETURN;
    }



    if (DMode == DISC_1D_NOCUPS)
    {
 
        if (IndexDim > 1)
        {
            DR_Error("Index dimension exceeds DEV dimension! "
                     "(Hyb3_Floater_Price)\n");
            goto RETURN;
        }

        if (PrincipalFlag)
        {
            Floater = (double *)FloaterPtr + Hyb3_Node_Offset(1,0,0,t,tree_data);
            for (i = Bottom1; i <= Top1; i ++)
            {
                Floater[i] += Principal;

            }  /* for i */	
        }  /* if */


        if (FixingFlag)	         /* Add fixing  */
        {
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            Floater = (double *)FloaterPtr + offset;
            Zero    = (double *)ZeroPtr    + offset;

            if (CoS == 'S')	     /* Simple rate */
            {
                if (Arrears == 'Y')
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                     i,0,0,t,tree_data);
                        Floater[i] += Outstanding * DayCntFtn * 
                                       (IndexAtNode + FloatSpd);
                    
                    }  /* for i */	
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                     i,0,0,t,tree_data);
                        Floater[i] += Outstanding * Zero[i] * DayCntFtn 
                                       * (IndexAtNode + FloatSpd);
                    
                    }  /* for i */	
                }  /* if then else */
            }
            else                    /* Compounded floating rate */
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                     i,0,0,t,tree_data);
                        Floater[i] += Outstanding * 
                         (pow (1. + (IndexAtNode + FloatSpd), DayCntFtn) - 1.);
                    
                    }  /* for i */	
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                     i,0,0,t,tree_data);
                        Floater[i] += Outstanding * Zero[i] * 
                         (pow (1. + (IndexAtNode + FloatSpd), DayCntFtn) - 1.);
                    
                    }  /* for i */	
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (DMode == DISC_2D_CUPS   ||
             DMode == DISC_2D_NOCUPS || 
             DMode == DISC_2D_1IR2F_NOCUPS )
    {
        
        

        if (IndexDim > 2)
        {
            DR_Error("Index dimension exceeds DEV dimension! "
                     "(Hyb3_Floater_Price)\n");
            goto RETURN;
        }


        if (PrincipalFlag)
        {
            
            for (i = Bottom1; i <= Top1; i ++)
            {
                Floater =(double *)FloaterPtr + Hyb3_Node_Offset(2,i,0,t,tree_data);
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Floater[j] += Principal;

                }  /* for j */
            }	
        }  /* if */


        if (FixingFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                        Floater = (double *)FloaterPtr + offset;
                        
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                         i,j,0,t,tree_data);
                            Floater[j] += Outstanding * DayCntFtn 
                                            * (IndexAtNode + FloatSpd);
                    
                        }  /* for j */	
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                        Floater = (double *)FloaterPtr + offset;
                        Zero    = (double *)ZeroPtr + offset;
                        
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode=Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                       i,j,0,t,tree_data);
                            Floater[j] += Outstanding * Zero[j] 
                                        * DayCntFtn * (IndexAtNode + FloatSpd);
                    
                        }  /* for j */	
                    }
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                        Floater = (double *)FloaterPtr + offset;
                        
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                         i,j,0,t,tree_data);
                            Floater[j] += Outstanding * (pow 
                              (1. + (IndexAtNode + FloatSpd), DayCntFtn) - 1.);

                        }  /* for j */
                    }	
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                        Floater = (double *)FloaterPtr + offset;
                        Zero    = (double *)ZeroPtr + offset;
                       
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                                         i,j,0,t,tree_data);
                            Floater[j] += Outstanding * Zero[j] * (pow 
                              (1. + (IndexAtNode + FloatSpd), DayCntFtn) - 1.);
                        
                        }  /* for j */
                    }	
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }


    else if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {

       

        if (IndexDim > 3)
        {
            DR_Error("Index dimension exceeds DEV dimension! "
                     "(Hyb3_Floater_Price)\n");
            goto RETURN;
        }


        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Floater = (double *)FloaterPtr 
                              + Hyb3_Node_Offset(3,i,j,t,tree_data);
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        Floater[k] += Principal;

                    }  /* for k */	
                }
            }
        }  /* if */


        if (FixingFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)	
                        {           
			    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                            Floater = (double *)FloaterPtr + offset;
			    IndexAtNodeSlice = (double *)IndexPtr + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                Floater[k] += Outstanding * DayCntFtn 
                                                   * (IndexAtNode + FloatSpd);
                            }  /* for k */	
                        }
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
			    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                            Floater = (double *)FloaterPtr + offset;
                            Zero    = (double *)ZeroPtr + offset;
			    IndexAtNodeSlice = (double *)IndexPtr + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                Floater[k] += Outstanding * Zero[k] 
                                        * DayCntFtn * (IndexAtNode + FloatSpd);
                            }  /* for k */
                        }
                    }	
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
			    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                            Floater = (double *)FloaterPtr + offset;
			    IndexAtNodeSlice = (double *)IndexPtr + offset;
                           
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                Floater[k] += Outstanding * (pow 
                                (1. + (IndexAtNode + FloatSpd),DayCntFtn) -1.);
                    
                            }  /* for k */	
                        }
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
			    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                            Floater = (double *)FloaterPtr + offset;
                            Zero =    (double *)ZeroPtr    + offset;
			    IndexAtNodeSlice = (double *)IndexPtr + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                Floater[k] += Outstanding * Zero[k]
                                * (pow (1. + (IndexAtNode + FloatSpd), 
                                              DayCntFtn) - 1.);
                        
                            }  /* for k */
                        }
                    }	
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }  /* if then else */

    status = SUCCESS;

  RETURN:

    return (status);



}  /* Hyb3_Floater_Price */








/*****  Hyb3_FwdFloater_t  ********************************************************/
/*
*       Forward floater price. Floating leg is arbitrary.
*/
int   Hyb3_FwdFloater_t(
          TSLICE      FwdFloater,       /* (O) Forward floater               */
          TSLICE      Floater,          /* (I) Underlying floater            */
          int         IndexDim,         /* (I) Dimension of FltIndex         */
          TSLICE      FltIndex,         /* (I) Last/current index rate       */
          int         ZeroDim,          /* (I) Dimension of FltZeroToPmt     */
          TSLICE      FltZeroToPmt,     /* (I) Zero w/mat on next flt coupon */
          long        FltResetFlag,     /* (I) Flt coupon flag               */
          double      Outstanding,      /* (I) Flt coupon outstanding        */
          double      FltCpnDcf,        /* (I) Flt coupon dcf                */
          long        FltCpnAccSt,      /* (I) Flt coupon accrued start      */
          char        FltDCConv,        /* (I) Flt coupon day count conv     */
          double      FloatSpd,         /* (I) Floating leg spread           */
          char        CoS,              /* (I) Flt coupon compound flag      */
          long        ExerFlag,         /* (I) Exercise flag                 */
          char        OptStub,          /* (I) Option stub rule              */
          char        ArreasReset,      /* (I) 'Y' if set-in-arreas          */
          long        CurrentDate,      /* (I) Current date                  */
          int         t,                /* (I) Current time period           */
          int         T,                /* (I) Total number of period        */
          int         DCurve,           /* (I) Discount curve                */
          int         DMode,            /* (I) Dimension of slices and  DEV  */
          HYB3_DEV_DATA   *dev_data,         /* (I) Hyb3_Dev data structure            */
          HYB3_TREE_DATA  *tree_data)        /* (I) Tree data structure           */
{



    double  *FwdFloaterL;           /* Local slice pointers               */
    double  *FloaterL;
    double  *FltZeroToPmtL;
    double  *IndexAtNodeSlice, *ZeroAtNodeSlice;


    double   IndexAtNode;
    double   ZeroAtNode;
    double   FltAccrued = 0.0;  /* Flt coupon accrued at current date */
    
    
    int              
              Top1,   Bottom1,  /* Tree limits (1rst dim)             */
             *Top2,  *Bottom2,  /* Tree limits (2nd dim)              */
            **Top3, **Bottom3,  /* Tree limits (3rd dim)              */
            i, j, k,            /* Node indices                       */
            offset,
            status = FAILURE;   /* Error status                       */

    
    

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /*  If this is not an exercise date discount and return. */
    if (!ExerFlag)                                              
    {
        if (Hyb3_Dev (FwdFloater,/* Disc expd value of fwd flaoter starting on */
                 t,         /* last exercise date                         */
                 T,
                 DCurve,
                 DMode,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
            
        }  /* if */
    
        status = SUCCESS;

        return (status);
        
    }  /* if */

   
    /* Must check if DMode is compatible with underlying tree */
    if (DMode == DISC_3D_CUPS|| DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if (tree_data->TreeType == TTYPE_2IR ||
            tree_data->TreeType == TTYPE_EQ1IR)
        {
            DR_Error("Slice dimension chosen is incompatible with tree type! "
                    "(Hyb3_FwdFloater_t)\n");
            goto RETURN;
        }
    }

    /* Check that index dimension is compatible with other variables */
    if ( IndexDim == 3    && 
        (tree_data->TreeType == TTYPE_2IR     ||
         tree_data->TreeType == TTYPE_EQ1IR ))
    {
        DR_Error("Index dimension chosen is incompatible with tree type! "
                    "(Hyb3_FwdFloater_t)\n");
        goto RETURN;
    }

    /* Check that zero dimension is compatible with other variables */
    if (  ZeroDim == 3 && 
        (tree_data->TreeType == TTYPE_2IR     ||
         tree_data->TreeType == TTYPE_EQ1IR ))
    {
        DR_Error("Zero dimension chosen is incompatible with tree type! "
                    "(Hyb3_FwdFloater_t)\n");
        goto RETURN;
    }

    /* NB: TTYPE_2IR and TTYPE_EQ1IR are the only */
    /*     non-3F types available in HYB3         */



    /*   Start from "basic" floater  */
    if (DMode == DISC_1D_NOCUPS)
    {
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        FwdFloaterL = (double *) FwdFloater + offset;
        FloaterL    = (double *) Floater + offset;
       

        for (i = Bottom1; i <= Top1; i ++)
        {
            FwdFloaterL[i] = FloaterL[i];
        }  

    }
    else if (DMode == DISC_2D_CUPS   || 
             DMode == DISC_2D_NOCUPS || 
             DMode == DISC_2D_1IR2F_NOCUPS )
    {

        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);  
            FwdFloaterL  = (double *)FwdFloater + offset; 
            FloaterL     = (double *)Floater + offset;    

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                FwdFloaterL[j] = FloaterL[j];
            }  
        }
    }
    else if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
     
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                                                            
                FwdFloaterL = (double *)FwdFloater + offset;          
                FloaterL    = (double *)Floater    + offset;          

                
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    FwdFloaterL[k] = FloaterL[k];
                }
            }
        }
        

    }



    /*  Calculate flt accrued interests if we are within accrual period. */
    if (FltCpnAccSt <= CurrentDate)
    {
        if (DrDayCountFraction (FltCpnAccSt,
                                CurrentDate,                                 
                                FltDCConv,
                                &FltAccrued) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for flt accrued (Hyb3_FwdFloater_t)");
            goto RETURN;
        }        
    }  /* if */
 

    if ((FltCpnAccSt < CurrentDate) && (DMode == DISC_1D_NOCUPS))
    {
        
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);  
                                                      
        FwdFloaterL   = (double *)FwdFloater   + offset;        
        FloaterL      = (double *)Floater      + offset;                
        FltZeroToPmtL = (double *)FltZeroToPmt + offset;  
        
              
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)   /* Subtract flt cpn just made  */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                     i,0,0,t,tree_data);
                        FwdFloaterL[i] -= (IndexAtNode + FloatSpd) * FltCpnDcf
                                        * Outstanding;

                    }          
                }      /* Swap,Bond stub : subtract PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                     i,0,0,t,tree_data);
                        ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                     i,0,0,t,tree_data);
                        FwdFloaterL[i] -= (IndexAtNode + FloatSpd * ZeroAtNode) 
                                        * FltAccrued * Outstanding ;
                        
                    }  
                }  /* if then else */   /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)       /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub: subtract accrued     */
                {                         /*            add paid PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                    i,0,0,t,tree_data);
                        ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                    i,0,0,t,tree_data);
                        FwdFloaterL[i] += (IndexAtNode + FloatSpd) *Outstanding 
                                        * (FltCpnDcf * ZeroAtNode - FltAccrued);
                        
                    }  
                }
                else if (OptStub == 'S') /* Swap stub: add PVed(cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                     i,0,0,t,tree_data);
                        ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                     i,0,0,t,tree_data);
                        FwdFloaterL[i] += (IndexAtNode + FloatSpd) * Outstanding 
                                        * (FltCpnDcf - FltAccrued) * ZeroAtNode;
                        
                    }  
                }
                else                  /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                     i,0,0,t,tree_data);
                        ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                     i,0,0,t,tree_data);
                        FwdFloaterL[i] += (IndexAtNode+FloatSpd) * Outstanding 
                                        * FltCpnDcf * ZeroAtNode;
                        
                    }  
                }  /* if then else */ 
            }
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)         /* Subtract flt cpn just paid */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                     i,0,0,t,tree_data);
                        FwdFloaterL[i] -= (pow(1. + (IndexAtNode + FloatSpd), 
                                                       FltCpnDcf) -1.)
                                        * Outstanding;
                        
                    }  
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("Accrued not supported for comp flt payment "
                             "(Hyb3_FwdFloater_t)");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)   /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Accrued not supported for comp flt payment"
                             " (Hyb3_FwdFloater_t)");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }


    else if ((FltCpnAccSt < CurrentDate) && 
             ((DMode == DISC_2D_CUPS)   ||
              (DMode == DISC_2D_NOCUPS) ||
              (DMode == DISC_2D_1IR2F_NOCUPS )))
    {
        

        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Subtract flt cpn just made   */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                       
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                        FwdFloaterL = (double *)FwdFloater + offset;          
                        
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                        i,j,0,t,tree_data);
                            FwdFloaterL[j] -= (IndexAtNode + FloatSpd) 
                                               * FltCpnDcf * Outstanding;

                        }          
                    }
                }                  /* Swap,Bond stub : subtract PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                        FwdFloaterL = (double *)FwdFloater + offset;          
                        
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                         i,j,0,t,tree_data);
                            ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                         i,j,0,t,tree_data);
                            FwdFloaterL[j] -= (IndexAtNode + FloatSpd 
                                 * ZeroAtNode) * FltAccrued * Outstanding ;
                        
                        } 
                    } 
                }  /* if then else */  /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)      /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub:subtract accrued    */
                {                         /*           add paid PVed coupon*/
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                        FwdFloaterL = (double *)FwdFloater + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                         i,j,0,t,tree_data);
                            ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                         i,j,0,t,tree_data);
                            FwdFloaterL[j] += (IndexAtNode+FloatSpd) *Outstanding 
                                       * (FltCpnDcf * ZeroAtNode - FltAccrued);
                            
                        }  
                    }
                }
                else if (OptStub == 'S')  /* Swap stub:add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                        FwdFloaterL = (double *)FwdFloater + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                         i,j,0,t,tree_data);
                            ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                         i,j,0,t,tree_data);
                            FwdFloaterL[j] += (IndexAtNode + FloatSpd) 
                            * Outstanding * (FltCpnDcf - FltAccrued) * ZeroAtNode;
                            
                        } 
                    } 
                }
                else               /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                        FwdFloaterL = (double *)FwdFloater + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                         i,j,0,t,tree_data);
                            ZeroAtNode  = Hyb3_GetValueAtNode(ZeroDim,FltZeroToPmt,
                                                         i,j,0,t,tree_data);
                            FwdFloaterL[j] += (IndexAtNode + FloatSpd) 
                                               * Outstanding 
                                               * FltCpnDcf * ZeroAtNode;
                            
                        } 
                    } 
                }  /* if then else */ 
            }
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)  /* Subtract flt cpn just made  */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {

                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                        FwdFloaterL = (double *)FwdFloater + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,FltIndex,
                                                         i,j,0,t,tree_data);
                            FwdFloaterL[j] -= (pow(1. + (IndexAtNode + FloatSpd), 
                                                  FltCpnDcf) -1.) * Outstanding;
                        
                        } 
                    } 
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("Accrued not supported for "
                             "comp flt payment (Hyb3_FwdFloater_t)");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)   /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Accrued not supported for "
                             "comp flt payment (Hyb3_FwdFloater_t)");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }


    else if ((FltCpnAccSt < CurrentDate) && 
             (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS))
    {
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)     /* Subtract flt cpn just made  */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data); 
                            FwdFloaterL = (double *)FwdFloater + offset;                
			    IndexAtNodeSlice = (double *)FltIndex + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                FwdFloaterL[k] -= (IndexAtNode + FloatSpd) 
                                                 * FltCpnDcf * Outstanding;

                            }  
                        }
                    }        
                }                  /* Swap,Bond stub : subtract PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);     
                            FwdFloaterL = (double *)FwdFloater + offset;               
			    IndexAtNodeSlice = (double *)FltIndex + offset;
			    ZeroAtNodeSlice = (double *)FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                ZeroAtNode  = ZeroAtNodeSlice[k];
                                FwdFloaterL[k] -= (IndexAtNode + FloatSpd * ZeroAtNode) 
                                                      * FltAccrued * Outstanding ;
                        
                            }
                        }
                    }  
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub : subtract accrued       */
                {                         /*             add paid PVed coupon   */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {

                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);     
                            FwdFloaterL = (double *)FwdFloater + offset;               
			    IndexAtNodeSlice = (double *)FltIndex + offset;
			    ZeroAtNodeSlice = (double *)FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                ZeroAtNode  = ZeroAtNodeSlice[k];
                                FwdFloaterL[k] += (IndexAtNode + FloatSpd) * Outstanding 
                                                * (FltCpnDcf * ZeroAtNode - FltAccrued);
                            } 
                        }
                    } 
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                            FwdFloaterL = (double *)FwdFloater + offset;          
			    IndexAtNodeSlice = (double *)FltIndex + offset;
			    ZeroAtNodeSlice = (double *)FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                ZeroAtNode  = ZeroAtNodeSlice[k];
                                FwdFloaterL[k] += (IndexAtNode + FloatSpd) * Outstanding 
                                                * (FltCpnDcf - FltAccrued) * ZeroAtNode;
                            
                            }  
                        }
                    }
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);     
                            FwdFloaterL = (double *)FwdFloater + offset;               
			    IndexAtNodeSlice = (double *)FltIndex + offset;
			    ZeroAtNodeSlice = (double *)FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                ZeroAtNode  = ZeroAtNodeSlice[k];
                                FwdFloaterL[k] += (IndexAtNode + FloatSpd) * Outstanding 
                                                 * FltCpnDcf * ZeroAtNode;
                            
                            } 
                        }
                    } 
                }  /* if then else */ 
            }
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)         /* Subtract flt cpn just paid  */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);     
                            FwdFloaterL = (double *)FwdFloater + offset;               
			    IndexAtNodeSlice = (double *)FltIndex + offset;
			    ZeroAtNodeSlice = (double *)FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexAtNode = IndexAtNodeSlice[k];
                                ZeroAtNode  = ZeroAtNodeSlice[k];
                                FwdFloaterL[k] -= (pow(1.+(IndexAtNode + FloatSpd), 
                                                     FltCpnDcf) -1.) * Outstanding;
                        
                            } 
                        }
                    } 
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("Accrued not supported for compounding flt "
                            "payment (Hyb3_FwdFloater_t)");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)   /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else           /* Stub : Not supported */
                {
                    DR_Error("Accrued not supported for comp flt payment "
                             "(Hyb3_FwdFloater_t)");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);



}  /* Hyb3_FwdFloater_t */







/*****  Hyb3_Collaret_t  ***********************************************************/
/*
*       Single collaret.
*/
int     Hyb3_Collaret_t (TSLICE      Collaret,     /* (I/O) Collaret              */
                    TSLICE      Index,        /* (I) Floating index          */
                    TSLICE      Zero,         /* (I) Zero to next payment    */
                    int         CollaretFlag, /* (I) Collare reset flag      */
                    double      CapRate,      /* (I) Cap rate                */
                    double      FloorRate,    /* (I) Floor rate              */
                    double      DayCntFtn,    /* (I) Day count fraction      */
                    char        CoS,          /* (I) 'C'ompound or 'S'imple  */
                    double      FloatSpd,     /* (I) Floating leg spread     */
                    char        Arrears,      /* (I) 'Y' if reset in arrears */
                    double      Notional,     /* (I) Caplet notional         */
                    int         t,            /* (I) Current time point      */
                    int         T,            /* (I) Last time point         */
                    int         DCurve,       /* (I) Discount curve          */
                    int         DMode,
                    HYB3_DEV_DATA    *dev_data,    /* (I) Hyb3_Dev data structure      */
                    HYB3_TREE_DATA   *tree_data)   /* (I) Tree data structure     */
{

    double  *CollaretL;                 /* Local slice pointers */
    double  *IndexL;
    double  *ZeroL;
        
    double  CompCapRate;
    double  CompFloorRate;

    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */

    int     i, j, k;                    /* Node indices           */
    int     offset;                     /* Node offset            */
    int     status = FAILURE;           /* Error status           */


    /* if not a reset, simply DEV the existing collaret slice and return */
    if (!CollaretFlag)
    {
        if (Hyb3_Dev(Collaret,
                t,
                T,
                DCurve,
                DMode,
                dev_data,
                tree_data) == FAILURE)
        {
            goto RETURN;
        }

        return SUCCESS;
    }

    /* calculate collaret value for the current slice */

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    CompCapRate   = pow (1. + CapRate, DayCntFtn) - 1.; 
    CompFloorRate = pow (1. + FloorRate, DayCntFtn) - 1.; 


    if (DMode == DISC_1D_NOCUPS)
    {
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        CollaretL = Collaret + offset;
        IndexL    = Index    + offset;
        ZeroL     = Zero     + offset;
    
        /* If there is a fixing we add a collaret */
        if (CoS == 'S')
        {
            if (Arrears == 'Y')                                         
            {                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CollaretL[i] = Notional * DayCntFtn 
                                 * MINMAX (IndexL[i] + FloatSpd, FloorRate, CapRate);
                }
            } 
            else
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CollaretL[i] = Notional * ZeroL[i] * DayCntFtn 
                                 * MINMAX (IndexL[i] + FloatSpd, FloorRate, CapRate);
                }
            }  /* if then else */
        }
        else
        {
            if (Arrears == 'Y')                                         
            {                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CollaretL[i] = Notional 
                                 * MINMAX ( pow (1. + IndexL[i] + FloatSpd, DayCntFtn) - 1., CompFloorRate, CompCapRate);
                }
            } 
            else
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CollaretL[i] = Notional * ZeroL[i] 
                                 * MINMAX ( pow (1. + IndexL[i] + FloatSpd, DayCntFtn) - 1., CompFloorRate, CompCapRate);
                }
            }  /* if then else */
        }  /* if then else */
    }
    else if (DMode == DISC_2D_CUPS   ||
             DMode == DISC_2D_NOCUPS || 
             DMode == DISC_2D_1IR2F_NOCUPS )
    {
        if (CoS == 'S')
        {
            if (Arrears == 'Y')                                         
            {                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    CollaretL = Collaret + offset;
                    IndexL    = Index    + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CollaretL[j] = Notional * DayCntFtn 
                                     * MINMAX (IndexL[j] + FloatSpd, FloorRate, CapRate);
                    }
                }  /* for i */
            } 
            else
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CollaretL = Collaret + offset;
                    IndexL    = Index    + offset;
                    ZeroL     = Zero     + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CollaretL[j] = Notional * ZeroL[j] * DayCntFtn 
                                     * MINMAX (IndexL[j] + FloatSpd, FloorRate, CapRate);
                    }
                }  /* for i */
            }  /* if then else */
        }
        else
        {
            if (Arrears == 'Y')                                         
            {                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    CollaretL = Collaret + offset;
                    IndexL    = Index    + offset;
 
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CollaretL[j] = Notional 
                                     * MINMAX ( pow (1. + IndexL[j] + FloatSpd, DayCntFtn) - 1., CompFloorRate, CompCapRate);
                    }
                }  /* for i */
            } 
            else
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CollaretL = Collaret + offset;
                    IndexL    = Index    + offset;
                    ZeroL     = Zero     + offset;
    
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CollaretL[j] = Notional * ZeroL[j] 
                                     * MINMAX ( pow (1. + IndexL[j] + FloatSpd, DayCntFtn) - 1., CompFloorRate, CompCapRate);
                    }
                }  /* for i */
            }  /* if then else */
        }  /* if then else */
    }
    else if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if (CoS == 'S')
        {
            if (Arrears == 'Y')                                         
            {                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
    
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CollaretL[k] = Notional * DayCntFtn 
                                         * MINMAX (IndexL[k] + FloatSpd, FloorRate, CapRate);
                        }
                    }  /* for j */
                }
            }
            else
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
                        ZeroL     = Zero     + offset;
    
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CollaretL[k] = Notional * ZeroL[k] * DayCntFtn 
                                         * MINMAX (IndexL[k] + FloatSpd, FloorRate, CapRate);
                        }
                    }  /* for j */
                } /* for j */
            }  /* if then else */
        }
        else
        {
            if (Arrears == 'Y')                                         
            {                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
    
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CollaretL[k] = Notional 
                                         * MINMAX ( pow (1. + IndexL[k] + FloatSpd, DayCntFtn) - 1., CompFloorRate, CompCapRate);
                        }
                    }  /* for j */
                }  /* for i */
            } 
            else
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
                        ZeroL     = Zero     + offset;
    
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CollaretL[k] = Notional * ZeroL[k] 
                                         * MINMAX ( pow (1. + IndexL[k] + FloatSpd, DayCntFtn) - 1., CompFloorRate, CompCapRate);
                        }
                    }  /* for j */
                }
            }  /* if then else */
        }  /* if then else */
    }  /* if then else */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Collaret_t */






/*****  Hyb3_Floater_Pmt *******************************************************/
/*
*       Price a single floater pmt.
*
*       Similar to Hyb3_Floater_Price, but no dev is done here. Also, anything 
*       already stored in FloaterPtr Slice will be overwritten.
*/
int   	Hyb3_Floater_Pmt(
            TSLICE         FloaterPtr,    /* (I/O) Floater                   */
            int            IndexDim,      /* (I) Dimension of the Index slice*/
            TSLICE         IndexPtr,      /* (I) Floating index              */
            TSLICE         ZeroPtr,       /* (I) Zero maturing on nxt pmt    */
            double         DayCntFtn,     /* (I) Day count fraction          */
            double         FloatSpd,      /* (I) Spread                      */
            double	       Outstanding,   /* (I) Outstanding notional        */
            double         Principal,     /* (I) Principal payment           */
            int            t,             /* (I) Current time period         */
            int            DMode,         /* (I) Dim of slices and of DEV    */
            HYB3_TREE_DATA     *tree_data)     /* (I) Structure of tree data      */
{


        double   IndexAtNode;
        double  *IndexAtNodeSlice;

        double	*Floater = NULL;
        double	*Zero    = NULL;

        int      offset;

        int
                Top1,      Bottom1,   /* Limits of the tree (1rst dimension) */
               *Top2,     *Bottom2,
              **Top3,    **Bottom3,

                i, j, k,              /* Node indices                        */
                status = FAILURE;     /* Error status = FAILURE initially    */

        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];	                



    if (DMode == DISC_1D_NOCUPS)
    {
        if (IndexDim > 1)
        {
            DR_Error("Index dimension exceeds DEV dimension! "
                     "(Hyb3_Floater_Price)\n");
            goto RETURN;
        }

        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        Floater = (double *)FloaterPtr + offset;
        Zero    = (double *)ZeroPtr    + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            Floater[i] = Principal;

            IndexAtNode = Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                         i,0,0,t,tree_data);
            Floater[i] += Outstanding * Zero[i] * DayCntFtn 
                           * (IndexAtNode + FloatSpd);
        
        }  /* for i */	
    }
    else if (DMode == DISC_2D_CUPS   ||
             DMode == DISC_2D_NOCUPS ||
             DMode == DISC_2D_1IR2F_NOCUPS )
    {
        
        

        if (IndexDim > 2)
        {
            DR_Error("Index dimension exceeds DEV dimension! "
                     "(Hyb3_Floater_Price)\n");
            goto RETURN;
        }


            
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
            Floater =(double *)FloaterPtr + offset;
            Zero    = (double *)ZeroPtr + offset;
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                Floater[j] = Principal;

                IndexAtNode=Hyb3_GetValueAtNode(IndexDim,IndexPtr,
                                           i,j,0,t,tree_data);
                Floater[j] += Outstanding * Zero[j] 
                            * DayCntFtn * (IndexAtNode + FloatSpd);
            }  /* for j */	
        }  /* if */
    }


    else if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {

        if (IndexDim > 3)
        {
            DR_Error("Index dimension exceeds DEV dimension! "
                     "(Hyb3_Floater_Price)\n");
            goto RETURN;
        }


        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                Floater = (double *)FloaterPtr + offset;
                Zero    = (double *)ZeroPtr + offset;
		IndexAtNodeSlice = (double *)IndexPtr + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {
                    Floater[k] = Principal;
		    
		    IndexAtNode = IndexAtNodeSlice[k];

                    Floater[k] += Outstanding * Zero[k] 
                            * DayCntFtn * (IndexAtNode + FloatSpd);
        
                }  /* for k */
            }
        }	
    }  /* if then else */

    status = SUCCESS;

  RETURN:

    return (status);



}  /* Hyb3_Floater_Pmt */


