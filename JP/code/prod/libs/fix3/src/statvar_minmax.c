/****************************************************************************/
/*      Extract the min/max values from a slice                             */
/****************************************************************************/
/*      minmax.c                                                            */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  Slice_MinMax **********************************************************/
/**
*       Extract the min/max values from a slice
*/
int     Slice_MinMax (  
                    double      *Slice,         /**< (I) Slice                */
                    double      *Min,           /**< (O) min value            */
                    double      *Max,           /**< (O) max value            */
                    int         t,              /**< (I) Current time point   */
                    FIX3_TREE_DATA   *tree_data) /* (I) Tree data structure   */
{

    double  *SliceL;                          /* Local slice pointer    */

    double  Vmin;
    double  Vmax;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */

    int     i, j, k;                        /* Node indices           */
    int     offset;                         /* Node offset            */
    int     status = FAILURE;               /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            SliceL   = Slice   + offset;
    
            Vmin = Vmax = SliceL[Top1];

            for (i = Bottom1; i < Top1; i ++)
            {
                if (SliceL[i] < Vmin)
                    Vmin = SliceL[i];
                else if (SliceL[i] > Vmax)
                    Vmax = SliceL[i]; 
            }
        }
        else if (tree_data->NbFactor == 2)
        {
            i = Top1; j = Top2[i];
 
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            SliceL   = Slice   + offset;
    
            Vmin = Vmax = SliceL[j];

            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SliceL   = Slice   + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    if (SliceL[j] < Vmin)
                        Vmin = SliceL[j];
                    else if (SliceL[j] > Vmax)
                        Vmax = SliceL[j]; 
                }
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            i = Top1; j = Top2[i]; k = Top3[i][j];

            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

            SliceL   = Slice   + offset;
    
            Vmin = Vmax = SliceL[k];

            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SliceL   = Slice   + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        if (SliceL[k] < Vmin)
                            Vmin = SliceL[k];
                        else if (SliceL[k] > Vmax)
                            Vmax = SliceL[k]; 
                    }
                }  /* for j */
        }  /* if then else */

    *Min = Vmin;
    *Max = Vmax;

    status = SUCCESS;

    return (status);

}  /* Slice_MinMax */

/*****  LadderRateRange ***************************************************/
/**
*       compute the ladder rate range given the previous range, sticky weight,
*       the precomputed H function (double Binary) and I function (leveraged
*       collar - Spread) of rates on a reset date
*/
static  int     LadderRateRange(
                    double  PrevCpnMin,
                    double  PrevCpnMax,
                    double  *CurrCpnMin,
                    double  *CurrCpnMax,
                    double  Binary,
                    double  Spread,
                    double  StickyWeight,
                    double  FlrLad,
                    double  CapLad,
                    char    AoM)
{
    double CurrCpn1;
    double CurrCpn2;

    if ( AoM == 'A')
    {
        CurrCpn1 = COLLAR(PrevCpnMin*StickyWeight + Spread + Binary,
                          CapLad, FlrLad);
        CurrCpn2 = COLLAR(PrevCpnMax*StickyWeight + Spread + Binary,
                          CapLad, FlrLad);
    }
    else if ( AoM == 'B')
    {
        CurrCpn1 = COLLAR(PrevCpnMin*StickyWeight + Spread * Binary,
                          CapLad, FlrLad);
        CurrCpn2 = COLLAR(PrevCpnMax*StickyWeight + Spread * Binary,
                          CapLad, FlrLad);
    }
    else
    {
        CurrCpn1 = COLLAR(PrevCpnMin*StickyWeight + Spread,
                          CapLad, FlrLad) * Binary;
        CurrCpn2 = COLLAR(PrevCpnMax*StickyWeight + Spread,
                          CapLad, FlrLad) * Binary;
    }

    if (CurrCpn1 < CurrCpn2)
    {
        *CurrCpnMin = CurrCpn1;
        *CurrCpnMax = CurrCpn2;
    }
    else
    {
        *CurrCpnMin = CurrCpn2;
        *CurrCpnMax = CurrCpn1;
    }

    return SUCCESS;
}

/*****  Ladder_Rate_MinMax ***************************************************/
/**
*       compute the min/max ladder rate on a reset date
*/
int     Ladder_Rate_MinMax (
                    double      *Binary,         /**< (I) Ladder Binary Slice */
                    double      *Spread,         /**< (I) Ladder Spread Slice */
                    double      StickyWeight,    /**< (I) Sticky Weight       */
                    double      FlrLad,          /**< (I) Ladder Rate Floor   */
                    double      CapLad,          /**< (I) Ladder Rate Cap     */
                    char        AoM,             /**< (I) Ladder calc flag    */
                    double      *Min,           /**< (I/O) min value          */
                    double      *Max,           /**< (I/O) max value          */
                    int         t,              /**< (I) Current time point   */
                    FIX3_TREE_DATA   *tree_data) /* (I) Tree data structure   */
{

    double  *BinaryL;                          /* Local slice pointer    */
    double  *SpreadL;                          /* Local slice pointer    */

    double  PrevCpnMin = *Min;
    double  PrevCpnMax = *Max;

    double  CurrCpnMin;
    double  CurrCpnMax;

    double  Vmin = CapLad;
    double  Vmax = FlrLad;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */

    int     i, j, k;                        /* Node indices           */
    int     offset;                         /* Node offset            */
    int     status = FAILURE;               /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


        if (tree_data->NbFactor == 1)
        {
            offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

            BinaryL   = Binary   + offset;
            SpreadL   = Spread   + offset;

            if (AoM == 'M') 
            {
                Vmin *= BinaryL[Top1];
                Vmax *= BinaryL[Top1];
            }

            for (i = Bottom1; i <= Top1; i ++)
            {
                LadderRateRange(PrevCpnMin,
                                PrevCpnMax,
                                &CurrCpnMin,
                                &CurrCpnMax,
                                BinaryL[i],
                                SpreadL[i],
                                StickyWeight,
                                FlrLad,
                                CapLad,
                                AoM);
                if (CurrCpnMin < Vmin) Vmin = CurrCpnMin;
                if (CurrCpnMax > Vmax) Vmax = CurrCpnMax;

            }
        }
        else if (tree_data->NbFactor == 2)
        {
            if (AoM == 'M') 
            {
                i = Top1; j = Top2[i];
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                BinaryL   = Binary   + offset;

                Vmin *= BinaryL[j];
                Vmax *= BinaryL[j];
            }

            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                BinaryL   = Binary   + offset;
                SpreadL   = Spread   + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    LadderRateRange(PrevCpnMin,
                                PrevCpnMax,
                                &CurrCpnMin,
                                &CurrCpnMax,
                                BinaryL[j],
                                SpreadL[j],
                                StickyWeight,
                                FlrLad,
                                CapLad,
                                AoM);
                    if (CurrCpnMin < Vmin) Vmin = CurrCpnMin;
                    if (CurrCpnMax > Vmax) Vmax = CurrCpnMax;

                }
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            if (AoM == 'M') 
            {
                i = Top1; j = Top2[i]; k = Top3[i][j];
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                BinaryL   = Binary   + offset;

                Vmin *= BinaryL[k];
                Vmax *= BinaryL[k];
            }

            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    BinaryL   = Binary   + offset;
                    SpreadL   = Spread   + offset;

    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        LadderRateRange(PrevCpnMin,
                                PrevCpnMax,
                                &CurrCpnMin,
                                &CurrCpnMax,
                                BinaryL[k],
                                SpreadL[k],
                                StickyWeight,
                                FlrLad,
                                CapLad,
                                AoM);
                        if (CurrCpnMin < Vmin) Vmin = CurrCpnMin;
                        if (CurrCpnMax > Vmax) Vmax = CurrCpnMax;

                    }
                }  /* for j */
        }  /* if then else */

    *Min = Vmin;
    *Max = Vmax;

    status = SUCCESS;

    return (status);

}  /* Slice_MinMax */
