/**************************************************************************/
/*      Ioucap                                                            */
/**************************************************************************/
/*      IOUCAP.c                                                          */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"



/*****  Ioucap_t  *********************************************************/
/*
*       Ioucap cap price.
*/
int     Ioucap_t 
           (double      **IouCap,     /* (I/O) Array of ioucap caps       */
            double      *Caplet,      /* (I) Current caplet/floorlet      */
            double      *Zero,        /* (I) Zero maturuing at next pmt   */
            long        CapletFlag,   /* (I) Reset flag                   */
            int         NbStates,     /* (I) Maximum number of caplets    */
            double      SMin,         /* (I) Minimum state variable       */
            double      SMax,         /* (I) Maximum state variable       */
            char        AoB,          /* (I) Above or Below ioucap        */
            char        Arrears,      /* (I) 'Y' if reset in arrears      */
            int         t,            /* (I) Current time point           */
            int         T,            /* (I) Last time point              */
            int         DCurve,       /* (I) Discount curve               */
            DEV_DATA    *dev_data,    /* (I) Dev data structure           */
            TREE_DATA   *tree_data)   /* (I) Tree data structure          */
{

    double  *IouCapL[MAXNBSTATES];    /* Local slice pointers */
    double  *CapletL;
    double  *ZeroL;
            
    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs     */
    double  E[MAXNBSTATES];           /* Precomputed linear coeffs        */
    double  State[MAXNBSTATES];       /* Levels of state variable (states)*/
    double  IC[MAXNBSTATES];          /* Intermediate values of ioucap    */
    double  deltaS;                   /* Increament between cons states   */
    double  NState;                   /* Next state                       */
    double  ICmin, ICmax;             /* Min and max of IC around sj      */
    
    int     s;                        /* State variable index             */
    int     sj;                       /* Interp state variable index      */
    int     q, l;                     /* Indices for quad, linear intrp   */
    
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Ioucap_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    if (SMax < SMin)
    {
        DR_Error("Ioucap_t: spent amount is bigger than level amount!");
        goto RETURN;
    }

    /* Discount all the state variables */
    
    for (s = 0; s < NbStates; s++)
    {                                  
        if (Dev (IouCap[s],
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }  /* for s */


    /* Ioucap state variable propagation */
    
    if (CapletFlag)
    {
        /* Prepare interpolation quadratic and linear polynomials */
        
        deltaS = (SMax - SMin)/(NbStates - 1);
        for (s = 0; s < NbStates; s++)
        {
            State[s] = SMin + deltaS * s;
        }
        for (s = 0; s < NbStates-2; s++)
        {
            D[s][0] = 1. / ((State[s]-State[s+1])*(State[s]-State[s+2]));
            D[s][1] = 1. / ((State[s+1]-State[s])*(State[s+1]-State[s+2]));
            D[s][2] = 1. / ((State[s+2]-State[s])*(State[s+2]-State[s+1]));
        }
        for (s = 0; s < NbStates-1; s++)
        {
            E[s] = 1. / (State[s]-State[s+1]);
        }            

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                IouCapL[s] = IouCap[s] + offset;
            }

            CapletL = Caplet + offset;
            ZeroL   = Zero   + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {                                    
                for (s = 0; s < NbStates ; s++)
                {
                    /* Determine new value of the state variable as a */
                    /* function of index s and caplet payment         */
                    /* Determine indices for interpolation sj         */

                    /* Remember to remove discount factor from caplet */
                    /* if reset-in-advance                            */
                    
                    if (Arrears == 'Y')
                    { 
                        NState = State[s] + CapletL[i];
                    }
                    else
                    {
                        NState = State[s] + CapletL[i] / ZeroL[i];
                    }
                    sj = MAX((int)((NState - SMin) / deltaS), 0);
                    
                    if (sj >= (NbStates - 1))
                    {
                        /* New value outside of [0,Level].      */      
                        /* Use cap with amount state[NbState-1] */

                        IC[s] = IouCapL[NbStates-1][i];
                    }
                    else
                    {
                        /* New value inside [0,Level].          */
                        /* Interpolate.                         */
                        
                        q = MIN(sj, NbStates - 3);
                    
                        IC[s] = (NState - State[q+1]) * (NState - State[q+2])
                                * D[q][0] * IouCapL[q][i]    
                              + (NState - State[q]) * (NState - State[q+2])
                                * D[q][1] * IouCapL[q+1][i]    
                              + (NState - State[q]) * (NState - State[q+1])
                                * D[q][2] * IouCapL[q+2][i];    

                        /* Check if quadratic interp inside [min,max] */

                        ICmin = MIN(    IouCapL[q][i], 
                                    MIN(IouCapL[q+1][i],
                                        IouCapL[q+2][i]));
                        ICmax = MAX(    IouCapL[q][i], 
                                    MAX(IouCapL[q+1][i],
                                        IouCapL[q+2][i]));

                        if ((IC[s] < ICmin) || (ICmax < IC[s]))
                        {
                            /* Quadratic interp failed. Use linear */

                            l = MIN(sj, NbStates - 2);
                            IC[s] = (NState - State[l+1]) 
                                      * E[l] * IouCapL[l][i]    
                                  + (State[l] - NState) 
                                      * E[l] * IouCapL[l+1][i];
                        }
                    }

                    if (Arrears == 'Y')
                    {
                        if (AoB == 'A')
                        {
                            IC[s] += MAX(CapletL[i] + State[s] - SMax, 0.);
                        }
                        else
                        {
                            IC[s] += MIN(CapletL[i], SMax - State[s]);
                        }
                    }
                    else
                    {
                        if (AoB == 'A')
                        {
                            IC[s] += MAX(ZeroL[i] * (State[s] - SMax)
                                         + CapletL[i], 0.);
                        }
                        else
                        {
                            IC[s] += MIN(ZeroL[i] * (SMax - State[s]),
                                         CapletL[i]);
                        }
                    } /* if Arrears */
                }  /* for s */

                /* Copy from intermadiate variable to ioucap */
                for (s = 0; s < NbStates ; s++)  
                {
                    IouCapL[s][i] = IC[s];
                }
            }  /* for i */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    IouCapL[s] = IouCap[s]+offset;
                }

                CapletL = Caplet + offset;
                ZeroL   = Zero   + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)      
                    {
                        /* Determine new value of the state variable as a */
                        /* function of index s and caplet payment         */
                        /* Determine indices for interpolation sj         */

                        /* Remember to remove discount factor from caplet */
                        /* if reset-in-advance                            */
                            
                        if (Arrears == 'Y')
                        { 
                            NState = State[s] + CapletL[j];
                        }
                        else
                        {
                            NState = State[s] + CapletL[j] / ZeroL[j];
                        }
                        sj = MAX((int)((NState - SMin) / deltaS), 0);
                        
                        if (sj >= (NbStates - 1))
                        {
                            /* New value outside of [0,Level].      */      
                            /* Use cap with amount state[NbState-1] */

                            IC[s] = IouCapL[NbStates-1][j];
                        }
                        else
                        {
                            /* New value inside [0,Level].          */
                            /* Interpolate.                         */
                            
                            q = MIN(sj, NbStates - 3);
                    
                            IC[s] = (NState - State[q+1])*(NState - State[q+2])
                                      * D[q][0] * IouCapL[q][j]    
                                  + (NState - State[q]) * (NState - State[q+2])
                                      * D[q][1] * IouCapL[q+1][j]    
                                  + (NState - State[q]) * (NState - State[q+1])
                                      * D[q][2] * IouCapL[q+2][j];    
                            
                            /* Check if quadratic interp inside [min,max] */

                            ICmin = MIN(    IouCapL[q][j], 
                                        MIN(IouCapL[q+1][j],
                                            IouCapL[q+2][j]));
                            ICmax = MAX(    IouCapL[q][j], 
                                        MAX(IouCapL[q+1][j],
                                            IouCapL[q+2][j]));

                            if ((IC[s] < ICmin) || (ICmax < IC[s]))
                            {
                                /* Quadratic interp failed. Use linear */
                                
                                l = MIN(sj, NbStates - 2);
                                IC[s] = (NState - State[l+1]) 
                                          * E[l] * IouCapL[l][j]    
                                      + (State[l] - NState) 
                                          * E[l] * IouCapL[l+1][j];
                            }
                        }

                        if (Arrears == 'Y')
                        {
                            if (AoB == 'A')
                            {
                                IC[s] += MAX(CapletL[j] + State[s] - SMax, 0.);
                            }
                            else
                            {
                                IC[s] += MIN(CapletL[j], SMax - State[s]);
                            }
                        }
                        else
                        {
                            if (AoB == 'A')
                            {
                                IC[s] += MAX(ZeroL[j] * (State[s] - SMax)
                                             + CapletL[j], 0.);
                            }
                            else
                            {
                                IC[s] += MIN(ZeroL[j] * (SMax - State[s]),
                                             CapletL[j]);
                            }
                        } /* if Arrears */
                    }  /* for s */
                    
                    /* Copy from intermadiate variable to ioucap */
                    for (s = 0; s < NbStates ; s++)      
                    {
                        IouCapL[s][j] = IC[s];
                    }
                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        IouCapL[s] = IouCap[s]+offset;
                    }

                    CapletL = Caplet + offset;
                    ZeroL   = Zero   + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s = 0; s < NbStates ; s++) 
                        {
                            /* Determine new value of the state variable as */
                            /* a function of index s and caplet payment     */
                            /* Determine indices for interpolation sj       */
                            
                            /* Remember to remove discount factor from */
                            /* caplet if reset-in-advance              */
                            
                            if (Arrears == 'Y')
                            { 
                                NState = State[s] + CapletL[k];
                            }
                            else
                            {
                                NState = State[s] 
                                       + CapletL[k] / ZeroL[k];
                            }
                            sj = MAX((int)((NState - SMin) / deltaS), 0);
                                                        
                            if (sj >= (NbStates - 1))
                            {
                                /* New value outside of [0,Level].      */      
                                /* Use cap with amount state[NbState-1] */

                                IC[s] = IouCapL[NbStates-1][k];
                            }
                            else
                            {
                                /* New value inside [0,Level].          */
                                /* Interpolate.                         */
                            
                                q = MIN(sj, NbStates - 3);
                                
                                IC[s] = (NState-State[q+1])*(NState-State[q+2])
                                          * D[q][0] * IouCapL[q][k]    
                                      + (NState - State[q])*(NState-State[q+2])
                                          * D[q][1] * IouCapL[q+1][k]   
                                      + (NState - State[q])*(NState-State[q+1])
                                          * D[q][2] * IouCapL[q+2][k];  
                            
                                /* Check if quadratic interp ins [min,max] */
                                
                                ICmin = MIN(    IouCapL[q][k], 
                                            MIN(IouCapL[q+1][k],
                                                IouCapL[q+2][k]));
                                ICmax = MAX(    IouCapL[q][k], 
                                            MAX(IouCapL[q+1][k],
                                                IouCapL[q+2][k]));

                                if ((IC[s] < ICmin) || (ICmax < IC[s]))
                                {
                                    /* Quadratic interp failed. Use linear */
                                
                                    l = MIN(sj, NbStates - 2);
                                    IC[s] = (NState - State[l+1]) 
                                              * E[l] * IouCapL[l][k]    
                                          + (State[l] - NState) 
                                              * E[l] * IouCapL[l+1][k];
                                }
                            }

                            if (Arrears == 'Y')
                            {
                                if (AoB == 'A')
                                {
                                    IC[s] += MAX(CapletL[k] 
                                                 + State[s] - SMax, 0.);
                                }
                                else
                                {
                                    IC[s] += MIN(CapletL[k],
                                                 SMax - State[s]);
                                }
                            }
                            else
                            {
                                if (AoB == 'A')
                                {
                                    IC[s] += MAX(ZeroL[k] * (State[s] - SMax)
                                                 + CapletL[k], 0.);
                                }
                                else
                                {
                                    IC[s] += MIN(ZeroL[k] * (SMax - State[s]),
                                                 CapletL[k]);
                                }
                            } /* if Arrears */

                        }  /* for s */
                        
                        /* Copy from intermadiate variable to ioucap */
                        for (s = 0; s < NbStates ; s++) 
                        {
                            IouCapL[s][k] = IC[s];
                        }
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
    }  /* if CapletFlag*/

    status = SUCCESS;
    
    RETURN:

    return (status);

}  /* Ioucap_t */



/*****  Iouswap_t  *********************************************************/
/*
*       Iouswap price.
*/
int     Iouswap_t 
           (double      **IouSwap,    /* (I/O) Array of iou swaps         */
            double      *Swaplet,     /* (I) Current swaplet/floorlet     */
            double      *Zero,        /* (I) Zero maturuing at next pmt   */
            long        SwapletFlag,  /* (I) Reset flag                   */
            int         NbSlices,     /* (I) Maximum number of slices     */
            double      SMin,         /* (I) Minimum state variable       */
            double      SMax,         /* (I) Maximum state variable       */
            double      LevelAmt,     /* (I) Level amount                 */
            double      **PrevState,  /* (I) Ptr to array of prev state   */ 
            char        AoB,          /* (I) Above or Below iouswap       */
            char        KoC,          /* (I) Knockout or contingent       */
            char        Arrears,      /* (I) 'Y' if reset in arrears      */
            int         t,            /* (I) Current time point           */
            int         T,            /* (I) Last time point              */
            int         DCurve,       /* (I) Discount curve               */
            DEV_DATA    *dev_data,    /* (I) Dev data structure           */
            TREE_DATA   *tree_data)   /* (I) Tree data structure          */
{

    double  *IouSwapL[MAXNBSTATES+1]; /* Local slice pointers */
    double  *SwapletL;
    double  *ZeroL;
            
    double  D[MAXNBSTATES][3];        /* Precomputed quadratic coeffs   */
    double  *State = NULL;            /* Current levels of state var    */
    double  *PState = NULL;           /* Previous levels of state var   */
    double  IC[MAXNBSTATES];          /* Intermediate values of iouswap */
    double  deltaS;                   /* Increament between cons states */
    double  NState;                   /* Next state                     */
    double  ICSpec=0.;
    double  IouPmt, Pmt;              /* Convenience variables          */
    
    int     NbStates;                 /* Number of iouswaps             */
    int     s;                        /* State variable index           */
    int     sj;                       /* Interp state variable index    */
    int     q;                        /* Index for quad, linear intrp   */
    
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* If KoC = K then NbSlices-1 is a regular swap */
    if (KoC == 'K')
    {
        NbStates = NbSlices - 1;
    }
    else
    {
        NbStates = NbSlices;
    }

    if (NbStates > MAXNBSTATES)
    {
        DR_Error("Iouswap_t: nb of state variables exceeds limit (200)!");
        goto RETURN;
    }

    /* Discount all the state variables */
    
    for (s = 0; s < NbSlices; s++)
    {                                  
        if (Dev (IouSwap[s],
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }  /* for l */


    /* Iouswap state variable propagation */
    
    if (SwapletFlag)
    {
        if (SMax < SMin)
        {
            DR_Error("Iouswap_t: max state smaller than min state!");
            goto RETURN;
        }
        
        /* Prepare interpolation quadratic and linear polynomials */
        
        deltaS = (SMax - SMin)/(NbStates - 1);

        State = (double *) DR_Array (DOUBLE, 0, NbStates - 1);
        if (State == NULL)
        {
            DR_Error("Iouswap_t: could not allocate memory for State!");
            goto RETURN;
        }

        for (s = 0; s < NbStates; s++)
        {
            State[s] = SMin + deltaS * s;
        }

        /* Find previous state var array */
        PState = *PrevState;

        if (PState != NULL)
        {
            /* PState exists, it is necessary to interpolate */
            for (s = 0; s < NbStates-2; s++)
            {
                D[s][0] = 1. / ((PState[s]-PState[s+1])
                               *(PState[s]-PState[s+2]));
                D[s][1] = 1. / ((PState[s+1]-PState[s])
                               *(PState[s+1]-PState[s+2]));
                D[s][2] = 1. / ((PState[s+2]-PState[s])
                               *(PState[s+2]-PState[s+1]));
            }
        }
        
        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbSlices ; s++)
            {
                IouSwapL[s] = IouSwap[s] + offset;
            }

            SwapletL = Swaplet + offset;
            ZeroL    = Zero    + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {                                    
                /* Treat special slice for knock iou swap */
                if (KoC == 'K')
                {
                    if (   ((LevelAmt > 0) && (AoB == 'A'))
                        || ((LevelAmt < 0) && (AoB == 'B'))  )
                    {
                        ICSpec = IouSwapL[NbStates][i] + SwapletL[i];
                    }
                    else if (   ((LevelAmt > 0) && (AoB == 'B'))
                             || ((LevelAmt < 0) && (AoB == 'A'))  )
                    {
                        ICSpec = 0.;
                    }
                }

                for (s = 0; s < NbStates; s++)
                {
                    /* Remember to remove discount factor from swaplet*/
                    /* if reset-in-advance                            */
                    
                    if (Arrears == 'Y')
                    { 
                        Pmt = SwapletL[i];
                    }
                    else
                    {
                        Pmt = SwapletL[i] / ZeroL[i];
                    }

                    /* Treat out-of-domain states for knock out iouswap first. If   */
                    /* contingent or knock out not applying, find next state values */
                    if (   ((KoC == 'K') && (LevelAmt > 0) && (State[s] >= LevelAmt))
                        || ((KoC == 'K') && (LevelAmt <= 0) && (State[s] <= LevelAmt))  )
                    {
                        IC[s] = ICSpec;
                    }
                    else 
                    {
                        
                        /* Decide if there are iouswap to interpolate off */
                        if (PState != NULL)
                        {
                            /* Interpolate */
                            
                            /* Determine new value of the state variable as a */
                            /* function of index s and swaplet payment        */
                            /* Determine indices for interpolation sj         */
                            
                            NState = State[s] + Pmt;
                        
                            /* Check if the next state does not knock out.    */
                            /* If not then interpolate                        */

                            if (   ((KoC == 'K') && (LevelAmt > 0) && (NState >= LevelAmt))
                                || ((KoC == 'K') && (LevelAmt <= 0) && (NState <= LevelAmt))  )
                            {
                                IC[s] = IouSwapL[NbStates][i];
                            }
                            else
                            {
                                sj = (int) floor((NState - PState[0])/(PState[1] - PState[0]));
                                
                                if (sj >= (NbStates - 1))
                                {
                                    /* New value outside of [0,Level].      */      
                                    /* Use swap with amount state[NbState-1] */
                                    
                                    IC[s] = IouSwapL[NbStates-1][i];
                                }
                                else if (sj < 0)
                                {
                                    /* New value outside of [0,Level].      */      
                                    /* Use swap with amount state[0]        */
                                    
                                    IC[s] = IouSwapL[0][i];
                                }
                                else
                                {
                                    /* New value inside [0,Level].          */
                                    /* Interpolate.                         */
                                    
                                    q = MIN(sj, NbStates - 3);
                    
                                    sqinterp(PState[q], PState[q+1], PState[q+2],
                                             IouSwapL[q][i],IouSwapL[q+1][i],IouSwapL[q+2][i],
                                             D[q][0], D[q][1], D[q][2],
                                             NState,
                                             &(IC[s]));
                                       
                                } /* if then else (sj) */
                            } /* if then else (KoC) */
                        }   
                        else
                        {
                            /* Nothing to interpolate off */
                            IC[s] = 0.;

                        } /* if then else (PState) */
                        
                        /* Find iou pmt */

                        if (AoB == 'A')
                        {
                            if (LevelAmt > State[s])
                            {
                                IouPmt = MAX(Pmt + State[s] - LevelAmt, 0.);
                            }
                            else
                            {
                                IouPmt = MAX(Pmt, LevelAmt - State[s]);
                            }
                        }
                        else
                        {
                            if (LevelAmt > State[s])
                            {
                                IouPmt = MIN(Pmt, LevelAmt - State[s]);
                            }
                            else
                            {
                                IouPmt = MIN(Pmt + State[s] - LevelAmt, 0.);
                            }
                        }
                        
                        if (Arrears == 'Y')
                        {
                            IC[s] += IouPmt;
                        }
                        else
                        {
                            IC[s] += ZeroL[i] * IouPmt;
                        }
                    } /* if then else (KoC) */
                }  /* for s */

                /* Treatment of limiting value of state variable */
                if (KoC == 'K')
                {
                    IouSwapL[NbStates][i] = ICSpec;
                } 

                /* Copy from intermadiate variable to iouswap.   */
                /* Use only NbStates IC values                   */
                for (s = 0; s < NbStates ; s++)  
                {
                    IouSwapL[s][i] = IC[s];
                }
            }  /* for i */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbSlices ; s++)
                {
                    IouSwapL[s] = IouSwap[s] + offset;
                }

                SwapletL = Swaplet + offset;
                ZeroL    = Zero    + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    /* Treat special slice for knock iou swap */
                    if (KoC == 'K')
                    {
                        if (   ((LevelAmt > 0) && (AoB == 'A'))
                            || ((LevelAmt < 0) && (AoB == 'B'))  )
                        {
                            ICSpec = IouSwapL[NbStates][j] + SwapletL[j];
                        }
                        else if (   ((LevelAmt > 0) && (AoB == 'B'))
                                 || ((LevelAmt < 0) && (AoB == 'A'))  )
                        {
                            ICSpec = 0.;
                        }
                    }
                    
                    for (s = 0; s < NbStates; s++)
                    {
                        /* Remember to remove discount factor from swaplet*/
                        /* if reset-in-advance                            */
                        
                        if (Arrears == 'Y')
                        { 
                            Pmt = SwapletL[j];
                        }
                        else
                        {
                            Pmt = SwapletL[j] / ZeroL[j];
                        }
                        
                        /* Treat out-of-domain states for knock out iouswap first. If   */
                        /* contingent or knock out not applying, find next state values */
                        if (   ((KoC == 'K') && (LevelAmt > 0) && (State[s] >= LevelAmt))
                            || ((KoC == 'K') && (LevelAmt <= 0) && (State[s] <= LevelAmt))  )
                        {
                            IC[s] = ICSpec;
                        }
                        else 
                        {
                            
                            /* Decide if there are iouswap to interpolate off */
                            if (PState != NULL)
                            {
                                /* Interpolate */
                                
                                /* Determine new value of the state variable as a */
                                /* function of index s and swaplet payment        */
                                /* Determine indices for interpolation sj         */
                                
                                NState = State[s] + Pmt;
                                
                                /* Check if the next state does not knock out.    */
                                /* If not then interpolate                        */
                                
                                if (   ((KoC == 'K') && (LevelAmt > 0) && (NState >= LevelAmt))
                                    || ((KoC == 'K') && (LevelAmt <= 0) && (NState <= LevelAmt))  )
                                {
                                    IC[s] = IouSwapL[NbStates][j];
                                }
                                else
                                {
                                    sj = (int) floor((NState - PState[0])
                                                     /(PState[1] - PState[0]));
                                    
                                    if (sj >= (NbStates - 1))
                                    {
                                        /* New value outside of [0,Level].      */      
                                        /* Use swap with amount state[NbState-1] */
                                        
                                        IC[s] = IouSwapL[NbStates-1][j];
                                    }
                                    else if (sj < 0)
                                    {
                                        /* New value outside of [0,Level].      */      
                                        /* Use swap with amount state[0]        */
                                    
                                        IC[s] = IouSwapL[0][j];
                                    }
                                    else
                                    {
                                        /* New value inside [0,Level].          */
                                        /* Interpolate.                         */
                                        
                                        q = MIN(sj, NbStates - 3);
                                        
                                        sqinterp(PState[q], PState[q+1], PState[q+2],
                                                 IouSwapL[q][j],IouSwapL[q+1][j],IouSwapL[q+2][j],
                                                 D[q][0], D[q][1], D[q][2],
                                                 NState,
                                                 &(IC[s]));

                                    } /* if then else (sj) */
                                } /* if then else (KoC) */
                            }   
                            else
                            {
                                /* Nothing to interpolate off */
                                IC[s] = 0.;
                                
                            } /* if then else (PState) */
                            
                            /* Find iou pmt */

                            if (AoB == 'A')
                            {
                                if (LevelAmt > State[s])
                                {
                                    IouPmt = MAX(Pmt + State[s] - LevelAmt, 0.);
                                }
                                else
                                {
                                    IouPmt = MAX(Pmt, LevelAmt - State[s]);
                                }
                            }
                            else
                            {
                                if (LevelAmt > State[s])
                                {
                                    IouPmt = MIN(Pmt, LevelAmt - State[s]);
                                }
                                else
                                {
                                    IouPmt = MIN(Pmt + State[s] - LevelAmt, 0.);
                                }
                            }
                            
                            if (Arrears == 'Y')
                            {
                                IC[s] += IouPmt;
                            }
                            else
                            {
                                IC[s] += ZeroL[j] * IouPmt;
                            }
                        } /* if then else (KoC) */
                    }  /* for s */
                    
                    /* Treatment of limiting value of state variable */
                    if (KoC == 'K')
                    {
                        IouSwapL[NbStates][j] = ICSpec;
                    } 
                    
                    /* Copy from intermadiate variable to iouswap.   */
                    /* Use only NbStates IC values                   */
                    for (s = 0; s < NbStates ; s++)  
                    {
                        IouSwapL[s][j] = IC[s];
                    }
                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbSlices ; s++)
                    {
                        IouSwapL[s] = IouSwap[s] + offset;
                    }

                    SwapletL = Swaplet + offset;
                    ZeroL    = Zero    + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        /* Treat special slice for knock iou swap */
                        if (KoC == 'K')
                        {
                            if (   ((LevelAmt > 0) && (AoB == 'A'))
                                || ((LevelAmt < 0) && (AoB == 'B'))  )
                            {
                                ICSpec = IouSwapL[NbStates][k] + SwapletL[k];
                            }
                            else if (   ((LevelAmt > 0) && (AoB == 'B'))
                                     || ((LevelAmt < 0) && (AoB == 'A'))  )
                            {
                                ICSpec = 0.;
                            }
                        }
                        
                        for (s = 0; s < NbStates; s++)
                        {
                            /* Remember to remove discount factor from swaplet*/
                            /* if reset-in-advance                            */
                            
                            if (Arrears == 'Y')
                            { 
                                Pmt = SwapletL[k];
                            }
                            else
                            {
                                Pmt = SwapletL[k] / ZeroL[k];
                            }
                            
                            /* Treat out-of-domain states for knock out iouswap first. If   */
                            /* contingent or knock out not applying, find next state values */
                            if (   ((KoC == 'K') && (LevelAmt > 0) && (State[s] >= LevelAmt))
                                || ((KoC == 'K') && (LevelAmt <= 0) && (State[s] <= LevelAmt))  )
                            {
                                IC[s] = ICSpec;
                            }
                            else 
                            {
                                
                                /* Decide if there are iouswap to interpolate off */
                                if (PState != NULL)
                                {
                                    /* Interpolate */
                                    
                                    /* Determine new value of the state variable as a */
                                    /* function of index s and swaplet payment        */
                                    /* Determine indices for interpolation sj         */
                                    
                                    NState = State[s] + Pmt;
                                    
                                    /* Check if the next state does not knock out.    */
                                    /* If not then interpolate                        */
                                    
                                    if (   ((KoC == 'K') && (LevelAmt > 0) && (NState >= LevelAmt))
                                        || ((KoC == 'K') && (LevelAmt <= 0) && (NState <= LevelAmt))  )
                                    {
                                        IC[s] = IouSwapL[NbStates][k];
                                    }
                                    else
                                    {
                                        sj = (int) floor((NState - PState[0])
                                                         /(PState[1] - PState[0]));

                                        if (sj >= (NbStates - 1))
                                        {
                                            /* New value outside of [0,Level].      */      
                                            /* Use swap with amount state[NbState-1] */
                                            
                                            IC[s] = IouSwapL[NbStates-1][k];
                                        }
                                        else if (sj < 0)
                                        {
                                            /* New value outside of [0,Level].      */      
                                            /* Use swap with amount state[0]        */
                                            
                                            IC[s] = IouSwapL[0][k];
                                        }
                                        else
                                        {
                                            /* New value inside [0,Level].          */
                                            /* Interpolate.                         */
                                            
                                            q = MIN(sj, NbStates - 3);
                                            
                                            sqinterp(PState[q], PState[q+1], PState[q+2],
                                                     IouSwapL[q][k],IouSwapL[q+1][k],IouSwapL[q+2][k],
                                                     D[q][0], D[q][1], D[q][2],
                                                     NState,
                                                     &(IC[s]));

                                        } /* if then else (sj) */
                                    } /* if then else (KoC) */
                                }   
                                else
                                {
                                    /* Nothing to interpolate off */
                                    IC[s] = 0.;
                                    
                                } /* if then else (PState) */
                                
                                /* Find iou pmt */

                                if (AoB == 'A')
                                {
                                    if (LevelAmt > State[s])
                                    {
                                        IouPmt = MAX(Pmt + State[s] - LevelAmt, 0.);
                                    }
                                    else
                                    {
                                        IouPmt = MAX(Pmt, LevelAmt - State[s]);
                                    }
                                }
                                else
                                {
                                    if (LevelAmt > State[s])
                                    {
                                        IouPmt = MIN(Pmt, LevelAmt - State[s]);
                                    }
                                    else
                                    {
                                        IouPmt = MIN(Pmt + State[s] - LevelAmt, 0.);
                                    }
                                }
                                
                                if (Arrears == 'Y')
                                {
                                    IC[s] += IouPmt;
                                }
                                else
                                {
                                    IC[s] += ZeroL[k] * IouPmt;
                                }
                            } /* if then else (KoC) */
                        }  /* for s */
                        
                        /* Treatment of limiting value of state variable */
                        if (KoC == 'K')
                        {
                            IouSwapL[NbStates][k] = ICSpec;
                        } 
                        
                        /* Copy from intermadiate variable to iouswap.   */
                        /* Use only NbStates IC values                   */
                        for (s = 0; s < NbStates ; s++)  
                        {
                            IouSwapL[s][k] = IC[s];
                        }
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        
        /* Free PState and assign State to PrevSate */
        if (Free_DR_Array (PState, DOUBLE, 0, NbStates - 1) == FAILURE)
        {
            goto RETURN;
        }
        *PrevState = State;

    }  /* if SwapletFlag*/

    status = SUCCESS;
    
    RETURN:

    if (status == FAILURE)
    {
        Free_DR_Array (State, DOUBLE, 0, NbStates-1);
    }

    return (status);

}  /* Iouswap_t */
