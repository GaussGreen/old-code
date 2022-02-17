/****************************************************************************/
/*      Standard cap and caplet.                                            */
/****************************************************************************/
/*      CAP.c                                                               */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"



/*****  Cap_t  **************************************************************/
/*
*       Cap price: add caplets one by one.
*/
int     Cap_t ( double     *CapFloor,   /* (I/O) Cap/floor                   */
                double     *Index,      /* (I) Index                         */
                double     *Zero,       /* (I) Zero to next payment date     */
                int         CoF,        /* (I) 1 for cap, -1 for floor       */
                long        CapletFlag, /* (I) Cap/Floor reset flag          */
                double      Strike,     /* (I) Strike                        */
                double      DayCntFtn,  /* (I) Day count fraction            */
                char        CoS,        /* (I) 'C'ompounded or 'S'imple rate */
                double      FloatSpd,   /* (I) Spread                        */
                char        Arrears,    /* (I) 'Y' if reset in arrears       */
                double      Notional,   /* (I) Caplet notional               */
                int         t,          /* (I) Current time point            */
                int         T,          /* (I) Last time point               */
                int         DCurve,     /* (I) Discount curve                */
                DEV_DATA    *dev_data,  /* (I) Dev data structure            */
                TREE_DATA   *tree_data) /* (I) Tree data structure           */
{

    double  *CapFloorL;                 /* Local slice pointers */
    double  *IndexL;
    double  *ZeroL;
        
    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */

    int     i, j, k;                    /* Node indices           */
    int     offset;                     /* Node offset            */
    int     status = FAILURE;           /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Dev (   CapFloor,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */

    
    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        CapFloorL = CapFloor + offset;
        IndexL    = Index    + offset;
        ZeroL     = Zero     + offset;
    
        if (CapletFlag)                                                         /* Add a caplet */
        {
            if (CoS == 'S')                                                     /* Simple rate */
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapFloorL[i] += Notional * DayCntFtn * MAX (CoF * (IndexL[i] - Strike), 0.);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapFloorL[i] += Notional * ZeroL[i] * DayCntFtn * MAX (CoF * (IndexL[i] - Strike), 0.);
                    }
                }  /* if then else */
            }
            else                                                                /* Compounded floating rate */
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapFloorL[i] += Notional * MAX (CoF * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) 
                                                             - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapFloorL[i] += Notional * ZeroL[i] * MAX (CoF * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) 
                                                                        - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                    }
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (CapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapFloorL = CapFloor + offset;
                        IndexL    = Index    + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapFloorL[j] += Notional * DayCntFtn * MAX (CoF * (IndexL[j] - Strike), 0.);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapFloorL = CapFloor + offset;
                        IndexL    = Index    + offset;
                        ZeroL     = Zero     + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapFloorL[j] += Notional * ZeroL[j] * DayCntFtn * MAX (CoF * (IndexL[j] - Strike), 0.);
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
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapFloorL = CapFloor + offset;
                        IndexL    = Index    + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapFloorL[j] += Notional * MAX (CoF * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) 
                                                                 - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapFloorL = CapFloor + offset;
                        IndexL    = Index    + offset;
                        ZeroL     = Zero     + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapFloorL[j] += Notional * ZeroL[j] * MAX (CoF * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) 
                                                                            - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                        }
                    }  /* for i */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (CapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapFloorL = CapFloor + offset;
                            IndexL    = Index    + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapFloorL[k] += Notional * DayCntFtn * MAX (CoF * (IndexL[k] - Strike), 0.);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapFloorL = CapFloor + offset;
                            IndexL    = Index    + offset;
                            ZeroL     = Zero     + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapFloorL[k] += Notional * ZeroL[k] * DayCntFtn * MAX (CoF * (IndexL[k] - Strike), 0.);
                            }
                        }  /* for j */
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapFloorL = CapFloor + offset;
                            IndexL    = Index    + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapFloorL[k] += Notional * MAX (CoF * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) 
                                                                     - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapFloorL = CapFloor + offset;
                            IndexL    = Index    + offset;
                            ZeroL     = Zero     + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapFloorL[k] += Notional * ZeroL[k] * MAX (CoF * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) 
                                                                                - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                            }
                        }  /* for j */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }  /* if then else */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Cap_t */



/*****  Caplet_t  ***********************************************************/
/*
*       Single caplet.
*/
int     Caplet_t (  double     *Caplet,     /* (I/O) Caplet/floorlet         */
                    double     *Index,      /* (I) Index                     */
                    double     *Zero,       /* (I) Zero to next payment date */
                    int         CoF,        /* (I) 1 for cap, -1 for floor   */
                    long        CapletFlag, /* (I) Cap/Floor reset flag      */
                    double      Strike,     /* (I) Strike                    */
                    double      DayCntFtn,  /* (I) Day count fraction        */
                    char        CoS,        /* (I) 'C'ompounded or 'S'imple  */
                    double      FloatSpd,   /* (I) Spread                    */
                    char        Arrears,    /* (I) 'Y' if reset in arrears   */
                    double      Notional,   /* (I) Caplet notional           */
                    int         t,          /* (I) Current time point        */
                    int         T,          /* (I) Last time point           */
                    int         DCurve,     /* (I) Discount curve            */
                    DEV_DATA    *dev_data,  /* (I) Dev data structure        */
                    TREE_DATA   *tree_data) /* (I) Tree data structure       */
{

    double  *CapletL;                   /* Local slice pointers */
    double  *IndexL;
    double  *ZeroL;
        
    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */

    int     i, j, k;                    /* Node indices           */
    int     offset;                     /* Node offset            */
    int     status = FAILURE;           /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Dev (   Caplet,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */

    
    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        CapletL = Caplet + offset;
        IndexL  = Index  + offset;
        ZeroL   = Zero   + offset;
    
        if (CapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapletL[i] = Notional * DayCntFtn * MAX (CoF * (IndexL[i] - Strike), 0.);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapletL[i] = Notional * ZeroL[i] * DayCntFtn * MAX (CoF * (IndexL[i] - Strike), 0.);
                    }
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapletL[i] = Notional * MAX (CoF * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) 
                                                          - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CapletL[i] = Notional * ZeroL[i] * MAX (CoF * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) 
                                                                     - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                    }
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (CapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapletL = Caplet + offset;
                        IndexL  = Index  + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapletL[j] = Notional * DayCntFtn * MAX (CoF * (IndexL[j] - Strike), 0.);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapletL = Caplet + offset;
                        IndexL  = Index  + offset;
                        ZeroL   = Zero   + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapletL[j] = Notional * ZeroL[j] * DayCntFtn * MAX (CoF * (IndexL[j] - Strike), 0.);
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
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapletL = Caplet + offset;
                        IndexL  = Index  + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapletL[j] = Notional * MAX (CoF * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) 
                                                              - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapletL = Caplet + offset;
                        IndexL  = Index  + offset;
                        ZeroL   = Zero   + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CapletL[j] = Notional * ZeroL[j] * MAX (CoF * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) 
                                                                         - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                        }
                    }  /* for i */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (CapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapletL = Caplet + offset;
                            IndexL  = Index  + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapletL[k] = Notional * DayCntFtn * MAX (CoF * (IndexL[k] - Strike), 0.);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapletL = Caplet + offset;
                            IndexL  = Index  + offset;
                            ZeroL   = Zero   + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapletL[k] = Notional * ZeroL[k] * DayCntFtn * MAX (CoF * (IndexL[k] - Strike), 0.);
                            }
                        }  /* for j */
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapletL = Caplet + offset;
                            IndexL  = Index  + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapletL[k] = Notional * MAX (CoF * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) 
                                                                  - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapletL = Caplet + offset;
                            IndexL  = Index  + offset;
                            ZeroL   = Zero   + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CapletL[k] = Notional * ZeroL[k] * MAX (CoF * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) 
                                                                             - pow (1. + (Strike + FloatSpd), DayCntFtn)), 0.);
                            }
                        }  /* for j */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }  /* if then else */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Caplet_t */


/*****  Collaret_t  ***********************************************************/
/*
*       Single collaret.
*/
int     Collaret_t (double     *Collaret,     /* (I/O) Collaret              */
                    double     *Index,        /* (I) Floating index          */
                    double     *Zero,         /* (I) Zero to next payment    */
                    long        CollaretFlag, /* (I) Collare reset flag      */
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
                    DEV_DATA    *dev_data,    /* (I) Dev data structure      */
                    TREE_DATA   *tree_data)   /* (I) Tree data structure     */
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


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Dev (   Collaret,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */


    CompCapRate   = pow (1. + CapRate, DayCntFtn) - 1.; 
    CompFloorRate = pow (1. + FloorRate, DayCntFtn) - 1.; 


    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        CollaretL = Collaret + offset;
        IndexL    = Index    + offset;
        ZeroL     = Zero     + offset;
    
        /* If there is a fixing we add a collaret */
        if (CollaretFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CollaretL[i] = Notional * DayCntFtn 
                                     * MIN( MAX (IndexL[i] + FloatSpd, FloorRate), CapRate);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CollaretL[i] = Notional * ZeroL[i] * DayCntFtn 
                                     * MIN( MAX (IndexL[i] + FloatSpd, FloorRate), CapRate);
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
                                     * MIN( MAX ( pow (1. + IndexL[i] + FloatSpd, DayCntFtn) - 1., CompFloorRate), CompCapRate);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CollaretL[i] = Notional * ZeroL[i] 
                                     * MIN( MAX ( pow (1. + IndexL[i] + FloatSpd, DayCntFtn) - 1., CompFloorRate), CompCapRate);
                    }
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (CollaretFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CollaretL[j] = Notional * DayCntFtn 
                                         * MIN( MAX (IndexL[j] + FloatSpd, FloorRate), CapRate);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
                        ZeroL     = Zero     + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CollaretL[j] = Notional * ZeroL[j] * DayCntFtn 
                                         * MIN( MAX (IndexL[j] + FloatSpd, FloorRate), CapRate);
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
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CollaretL[j] = Notional 
                                         * MIN( MAX ( pow (1. + IndexL[j] + FloatSpd, DayCntFtn) - 1., CompFloorRate), CompCapRate);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CollaretL = Collaret + offset;
                        IndexL    = Index    + offset;
                        ZeroL     = Zero     + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CollaretL[j] = Notional * ZeroL[j] 
                                         * MIN( MAX ( pow (1. + IndexL[j] + FloatSpd, DayCntFtn) - 1., CompFloorRate), CompCapRate);
                        }
                    }  /* for i */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (CollaretFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CollaretL = Collaret + offset;
                            IndexL    = Index    + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CollaretL[k] = Notional * DayCntFtn 
                                             * MIN( MAX (IndexL[k] + FloatSpd, FloorRate), CapRate);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CollaretL = Collaret + offset;
                            IndexL    = Index    + offset;
                            ZeroL     = Zero     + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CollaretL[k] = Notional * ZeroL[k] * DayCntFtn 
                                             * MIN( MAX (IndexL[k] + FloatSpd, FloorRate), CapRate);
                            }
                        }  /* for j */
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CollaretL = Collaret + offset;
                            IndexL    = Index    + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CollaretL[k] = Notional 
                                             * MIN( MAX ( pow (1. + IndexL[k] + FloatSpd, DayCntFtn) - 1., CompFloorRate), CompCapRate);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CollaretL = Collaret + offset;
                            IndexL    = Index    + offset;
                            ZeroL     = Zero     + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CollaretL[k] = Notional * ZeroL[k] 
                                             * MIN( MAX ( pow (1. + IndexL[k] + FloatSpd, DayCntFtn) - 1., CompFloorRate), CompCapRate);
                            }
                        }  /* for j */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }  /* if then else */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Collaret_t */
