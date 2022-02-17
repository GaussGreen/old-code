/****************************************************************************/
/*      Calculate discounted expected value in the tree.                    */
/****************************************************************************/
/*      DEV.C                                                               */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"


    



/*****  Hyb3_Dev   ****************************************************************/
/*
 *       Discounted expected value in variable dimension (i.e. 1, 2 or 3). 
 */
int     Hyb3_Dev (TSLICE          Price,         /* (I/O) Slice to be discounted   */
                  int             t,             /* (I) Current time period        */
                  int             T,             /* (I) Total number of periods    */
                  int             DCurve,        /* (I) Disc curve to use (0,1,2)  */
                  int             DMode,         /* (I) Nb of dims of dev required */
                  HYB3_DEV_DATA  *dev_data,      /* (I) Probas, shifts, limits, etc*/
                  HYB3_TREE_DATA *tree_data)     /* (I) Structure of tree data     */
{





        /* Locals for convenience in addressing */
        double       *DiscountL;     /* Local discount slice       */
        double       *Discount_1D;   /* used if discounting in 1-D */
        double       *Discount_2D;   /* used if discounting in 2-D */
        double       *Discount_3D;   /* used if discounting in 3-D */

        double       *PriceL;
        double       *NewPriceL;
        
        double       *Price0;  
        double       *Price1;
        double       *Price2;
        double       *Price00;
        double       *Price01;
        double       *Price02;
        double       *Price10; 
        double       *Price11;
        double       *Price12;
        double       *Price20; 
        double       *Price21; 
        double       *Price22;

        /*  1-D probability slices    */  
        TPROB_0  *p;

        /*  2-D probability slices    */
        TPROB_0  *q;

        double   *quu,  *qu0,  *qud;
        double   *q0u,  *q00,  *q0d;
        double   *qdu,  *qd0,  *qdd;

        /*  3-D probability slices    */
        TPROB_0  *r;  

        double   *tx;

        /* 1-D probabilities (values) */
        double     Pu,    P0,   Pd;

        /* 2-D probabilities (values) */
        double    Quu,   Qu0,  Qud;
        double    Q0u,   Q00,  Q0d;
        double    Qdu,   Qd0,  Qdd;

        double DiscountLi, DiscountLij;
      
        int

            Top1,      Bottom1,    /* Limits of the tree (1rst dimension)                */
           *Top2,     *Bottom2,    /* Limits of the tree (2nd dimension)                 */
          **Top3,    **Bottom3,    /* Limits of the tree (3rd dimension)                 */
            NTop1,     NBottom1,   /* Limits of the tree at next period (1rst dimension) */
           *NTop2,    *NBottom2,   /* Limits of the tree at next period (2nd dimension)  */
          **NTop3,   **NBottom3,   /* Limits of the tree at next period (3rd dimension)  */
            OutTop1,   OutBottom1, /* Limits of the tree at next period (1rst dimension) */
           *OutTop2,  *OutBottom2, /* Limits of the tree at next period (2nd dimension)  */
          **OutTop3, **OutBottom3, /* Limits of the tree at next period (3rd dimension)  */


           *Shift1,  *Shift2,  *Shift3,

            i, i0, i1, i2,         /* Node indices                                       */
            j, j0, j1, j2,
            k, k0, k1, k2;
           
        int
            offset,

            jLowerNext,
            jUpperNext,

            kLowerNext,
            kUpperNext;
            
        int
            status = FAILURE;
        


        if (t == T) /* Nothing to do at the back of the tree */
        {
            return(SUCCESS);
        }

        
        Top1    = tree_data->Top1[t];	        
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];	    
        Bottom1 = tree_data->Bottom1[t];            
        Bottom2 = tree_data->Bottom2[t];    
        Bottom3 = tree_data->Bottom3[t];
                
        NTop1    = tree_data->Top1[t+1];
        NTop2    = tree_data->Top2[t+1];
        NTop3    = tree_data->Top3[t+1];
        NBottom1 = tree_data->Bottom1[t+1];
        NBottom2 = tree_data->Bottom2[t+1];
        NBottom3 = tree_data->Bottom3[t+1];

        OutTop1    = tree_data->OutTop1[t+1];	        
        OutTop2    = tree_data->OutTop2[t+1];
        OutTop3    = tree_data->OutTop3[t+1];	    
        OutBottom1 = tree_data->OutBottom1[t+1];            
        OutBottom2 = tree_data->OutBottom2[t+1];    
        OutBottom3 = tree_data->OutBottom3[t+1];


        /* Validate DEV mode */
        if (   DMode != DISC_1D_NOCUPS       &&  DMode != DISC_2D_NOCUPS
            && DMode != DISC_2D_CUPS         &&  DMode != DISC_3D_CUPS
            && DMode != DISC_2D_1IR2F_NOCUPS &&  DMode != DISC_3D_2IR2F1D_CUPS)
        {
            DR_Error("Invalid DEV mode requested! (Hyb3_Dev)\n");
            goto RETURN;
        }


        /* Validate DEV dimension for given tree type    */
        /* Do we need to check all possible combinaison? */
        /* Currently, the following improper combinaison are allowed: */
        /* - EQ1IR & 2D_CUPS       */
        /* - 1IR (CET) & 2D_CUPS   */




        if (DMode == DISC_1D_NOCUPS)
        {
            if (tree_data->TreeType == TTYPE_1IR2F||
                tree_data->TreeType == TTYPE_2IR2F1D)
            {
                DR_Error("Hyb3_Dev mode requested (1D_NOCUPS) not valid\n"
                        "for two factors foreign IR tree type! (Hyb3_Dev)\n");
                goto RETURN;
            }
        }
        else if (DMode == DISC_2D_CUPS)
        {
            if (tree_data->TreeType == TTYPE_1IR2F||
                tree_data->TreeType == TTYPE_2IR2F1D)
            {
                DR_Error("Hyb3_Dev mode requested (2D_CUPS) not valid\n"
                        "for two factor foreign IR tree type! (Hyb3_Dev)\n");
                goto RETURN;
            }
        }
        else if (DMode == DISC_2D_NOCUPS)
        {
            if (tree_data->TreeType != TTYPE_EQ1IR)
            {
                DR_Error("Hyb3_Dev mode requested (2D_NOCUPS) only valid\n"
                        "for equity plus interest rate tree type! (Hyb3_Dev)\n");
                goto RETURN;
            }
        }
        else if (DMode == DISC_2D_1IR2F_NOCUPS)
        {
            if (tree_data->TreeType != TTYPE_1IR2F&&
                tree_data->TreeType != TTYPE_2IR2F1D)
            {
                DR_Error("Hyb3_Dev mode requested (2D_1IR2F_NOCUPS) only valid\n"
                        "for two factors foreign IR tree type! (Hyb3_Dev)\n");
                goto RETURN;
            }
        }
        else if (DMode == DISC_3D_2IR2F1D_CUPS)
        {
            if (tree_data->TreeType != TTYPE_2IR2F1D)
            {
                DR_Error("Hyb3_Dev mode requested (3D_2IR2F1D_CUPS) only valid\n"
                        "for TTYPE_2IR2F1D tree type! (Hyb3_Dev)\n");
                goto RETURN;
            }
        }
        else
        {
            /* Must be runing in 3-D, therefore TreeType must be 3-D */
            if ((tree_data->TreeType != TTYPE_FX2IR)   &&
                (tree_data->TreeType != TTYPE_EQF2IR)  &&
                (tree_data->TreeType != TTYPE_EQD2IR)  &&
                (tree_data->TreeType != TTYPE_EQC2IR))
            {
                DR_Error("DEV mode chosen would need 3-D tree type! (Hyb3_Dev)\n");
                goto RETURN;
            }

        }


        /* Choice of the discount curve  */
        if ((DCurve == 0) || (DCurve == 1) || (DCurve == 2))
        {
            Discount_1D = dev_data->Discount_1D[DCurve];
            Discount_2D = dev_data->Discount_2D[DCurve];
            Discount_3D = dev_data->Discount_3D[DCurve];

        }
        else
        {
            DR_Error ("Incorrect specification for discount curve! (Hyb3_Dev)");
            goto RETURN;
        }





        /*** One factor DEV ***/
        if (DMode == DISC_1D_NOCUPS)
        {

            PriceL    = (double *)Price + Hyb3_Node_Offset(1,0,0,t+1,tree_data);

            /* Values of variable are assumed to be flat in the region */
            /* between the inner and the outer ellipses (which are de- */
            /* limited by NTop, OutTop, NBottom and OutBottom)         */

            for (i=OutBottom1; i<NBottom1; i++)
            {
                PriceL[i] = PriceL[NBottom1];
            }
            for (i=NTop1+1; i<=OutTop1; i++)
            {
                PriceL[i] = PriceL[NTop1];
            }

            /* Discounted expected value using the non-CUPS probs and shift */
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            Shift1 = dev_data->Shift4 + offset;
            p      = dev_data->s + offset;

            DiscountL = Discount_1D + offset;
            NewPriceL = dev_data->Aux1D + offset;
            

            for (i=Bottom1; i<=Top1; i++)
            {
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                NewPriceL[i] = p[i].u * PriceL[i0] +
                               p[i].m * PriceL[i1] +
                               p[i].d * PriceL[i2];

                NewPriceL[i] *= DiscountL[i];
            }

            /* Put the prices back in the original slice */
            PriceL = (double *)Price + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                PriceL[i] = NewPriceL[i];
            }
        }



        /*** Two factor DEV for CUPS mode ***/
        else if (DMode == DISC_2D_CUPS)
        {
            
            double    PFlat; 
            
            
            /* Fill values for nodes between limiting ellipses */
            PriceL    = (double *)Price + Hyb3_Node_Offset(2,NBottom1,0,t+1,tree_data);
            PFlat = PriceL[NTop2[NBottom1]];
            for (i=OutBottom1; i<NBottom1; i++)
            {
                PriceL = (double *)Price + Hyb3_Node_Offset(2,i,0,t+1,tree_data);

                for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            for (i=NBottom1; i<=NTop1; i++)
            {
                PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
                PFlat = PriceL[NBottom2[i]];

                for (j=OutBottom2[i]; j<NBottom2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }

                PFlat = PriceL[NTop2[i]];
                for (j=NTop2[i]+1; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            PriceL = (double *)Price + Hyb3_Node_Offset(2,NTop1,0,t+1,tree_data);
            PFlat = PriceL[NTop2[NTop1]];
            for (i=NTop1+1; i<=OutTop1; i++)
            {

                PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
             
                for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            /* Discounted expected value */
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

            Shift1 = dev_data->Shift1 + offset;
            p = dev_data->p + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                q = dev_data->q + offset;

                Shift2    = dev_data->Shift2 + offset;
                DiscountL = Discount_2D + offset;
                NewPriceL = dev_data->Aux2D + offset;

                Price0 = (double *)Price + Hyb3_Node_Offset(2,i0,0,t+1,tree_data);
                Price1 = (double *)Price + Hyb3_Node_Offset(2,i1,0,t+1,tree_data);
                Price2 = (double *)Price + Hyb3_Node_Offset(2,i2,0,t+1,tree_data);

                jLowerNext =                 OutBottom2[i0] ;
                jLowerNext = MAX(jLowerNext, OutBottom2[i1]);
                jLowerNext = MAX(jLowerNext, OutBottom2[i2]);

                jUpperNext =                 OutTop2[i0] ;
                jUpperNext = MIN(jUpperNext, OutTop2[i1]);
                jUpperNext = MIN(jUpperNext, OutTop2[i2]);

                tx = dev_data->t2 - jLowerNext;

                Pu = p[i].u;
                P0 = p[i].m;
                Pd = p[i].d;

                for (j=jLowerNext; j<=jUpperNext; ++j)
                {
                    tx[j] = Pu * Price0[j] +
                            P0 * Price1[j] +
                            Pd * Price2[j];
                }

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;

                    NewPriceL[j] = q[j].u * tx[j0] +
                                   q[j].m * tx[j1] +
                                   q[j].d * tx[j2];

                    NewPriceL[j] *=DiscountL[j];
                }
            }

            /* Finally transfer values from auxiliary slice to variable slice*/
            for (i=Bottom1; i<=Top1; i++)
            {
                
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                PriceL    = (double *)Price           + offset;
                NewPriceL = dev_data->Aux2D + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    PriceL[j] = NewPriceL[j];
                }
            }
        }


        
        /*** Two dim DEV for 1IR 2 factor NOCUPS mode   ***/
        /*   2D tree, 2D Discount, NOCUPS prob            */
        /*   Only change with DISC_2D_CUPS is prob. used  */
        /***  can be easily merged                      ***/
        else if (DMode == DISC_2D_1IR2F_NOCUPS)
        {
            
            double    PFlat; 
            
            
            /* Fill values for nodes between limiting ellipses */
            PriceL    = (double *)Price + Hyb3_Node_Offset(2,NBottom1,0,t+1,tree_data);
            PFlat = PriceL[NTop2[NBottom1]];
            for (i=OutBottom1; i<NBottom1; i++)
            {
                PriceL = (double *)Price + Hyb3_Node_Offset(2,i,0,t+1,tree_data);

                for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            for (i=NBottom1; i<=NTop1; i++)
            {
                PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
                PFlat = PriceL[NBottom2[i]];

                for (j=OutBottom2[i]; j<NBottom2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }

                PFlat = PriceL[NTop2[i]];
                for (j=NTop2[i]+1; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            PriceL = (double *)Price + Hyb3_Node_Offset(2,NTop1,0,t+1,tree_data);
            PFlat = PriceL[NTop2[NTop1]];
            for (i=NTop1+1; i<=OutTop1; i++)
            {

                PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
             
                for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            /* Discounted expected value */
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

            Shift1 = dev_data->Shift4 + offset;     /* NOCUPS prob.*/
            p = dev_data->s + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                q = dev_data->t + offset;              /*NOCUPS prob.*/
                Shift2    = dev_data->Shift5 + offset;

                DiscountL = Discount_2D + offset;
                NewPriceL = dev_data->Aux2D + offset;

                Price0 = (double *)Price + Hyb3_Node_Offset(2,i0,0,t+1,tree_data);
                Price1 = (double *)Price + Hyb3_Node_Offset(2,i1,0,t+1,tree_data);
                Price2 = (double *)Price + Hyb3_Node_Offset(2,i2,0,t+1,tree_data);

                jLowerNext =                 OutBottom2[i0] ;
                jLowerNext = MAX(jLowerNext, OutBottom2[i1]);
                jLowerNext = MAX(jLowerNext, OutBottom2[i2]);

                jUpperNext =                 OutTop2[i0] ;
                jUpperNext = MIN(jUpperNext, OutTop2[i1]);
                jUpperNext = MIN(jUpperNext, OutTop2[i2]);

                tx = dev_data->t2 - jLowerNext;

                Pu = p[i].u;
                P0 = p[i].m;
                Pd = p[i].d;

                for (j=jLowerNext; j<=jUpperNext; ++j)
                {
                    tx[j] = Pu * Price0[j] +
                            P0 * Price1[j] +
                            Pd * Price2[j];
                }

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;

                    NewPriceL[j] = q[j].u * tx[j0] +
                                   q[j].m * tx[j1] +
                                   q[j].d * tx[j2];

                    NewPriceL[j] *=DiscountL[j];
                }
            }

            /* Finally transfer values from auxiliary slice to variable slice*/
            for (i=Bottom1; i<=Top1; i++)
            {
                
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                PriceL    = (double *)Price           + offset;
                NewPriceL = dev_data->Aux2D + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    PriceL[j] = NewPriceL[j];
                }
            }
        }
        /*** Two factor DEV for EQ1IR mode ***/
        else if (DMode == DISC_2D_NOCUPS)
        {
            
            double    PFlat; 
            
            
            /* Fill values for nodes between limiting ellipses */
            PriceL    = (double *)Price + Hyb3_Node_Offset(2,NBottom1,0,t+1,tree_data);
            PFlat = PriceL[NTop2[NBottom1]];
            for (i=OutBottom1; i<NBottom1; i++)
            {
                PriceL = (double *)Price + Hyb3_Node_Offset(2,i,0,t+1,tree_data);

                for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            for (i=NBottom1; i<=NTop1; i++)
            {
                PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
                PFlat = PriceL[NBottom2[i]];

                for (j=OutBottom2[i]; j<NBottom2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }

                PFlat = PriceL[NTop2[i]];
                for (j=NTop2[i]+1; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            PriceL = (double *)Price + Hyb3_Node_Offset(2,NTop1,0,t+1,tree_data);
            PFlat = PriceL[NTop2[NTop1]];
            for (i=NTop1+1; i<=OutTop1; i++)
            {

                PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
             
                for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
                {
                    PriceL[j] = PFlat;
                }
            }

            /* Discounted expected value */
            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

            DiscountL = Discount_1D + offset;
            Shift1 = dev_data->Shift1 + offset;
            p = dev_data->p + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                q = dev_data->q + offset;

                Shift2    = dev_data->Shift2 + offset;
                NewPriceL = dev_data->Aux2D + offset;

                Price0 = (double *)Price + Hyb3_Node_Offset(2,i0,0,t+1,tree_data);
                Price1 = (double *)Price + Hyb3_Node_Offset(2,i1,0,t+1,tree_data);
                Price2 = (double *)Price + Hyb3_Node_Offset(2,i2,0,t+1,tree_data);

                jLowerNext =                 OutBottom2[i0] ;
                jLowerNext = MAX(jLowerNext, OutBottom2[i1]);
                jLowerNext = MAX(jLowerNext, OutBottom2[i2]);

                jUpperNext =                 OutTop2[i0] ;
                jUpperNext = MIN(jUpperNext, OutTop2[i1]);
                jUpperNext = MIN(jUpperNext, OutTop2[i2]);

                tx = dev_data->t2 - jLowerNext;

                Pu = p[i].u;
                P0 = p[i].m;
                Pd = p[i].d;

                DiscountLi = DiscountL[i];

                for (j=jLowerNext; j<=jUpperNext; ++j)
                {
                    tx[j] = Pu * Price0[j] +
                            P0 * Price1[j] +
                            Pd * Price2[j];
                }

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;

                    NewPriceL[j] = q[j].u * tx[j0] +
                                   q[j].m * tx[j1] +
                                   q[j].d * tx[j2];

                    NewPriceL[j] *= DiscountLi;
                }
            }

            /* Finally transfer values from auxiliary slice to variable slice*/
            for (i=Bottom1; i<=Top1; i++)
            {
                
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                PriceL    = (double *)Price           + offset;
                NewPriceL = dev_data->Aux2D + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    PriceL[j] = NewPriceL[j];
                }
            }
        }
        /*** Three dim DEV                  ***/
        /*   3D discount factor, CUPS prob. ***/
        else if (DMode == DISC_3D_2IR2F1D_CUPS)
        {
            /* Discounted expected value */
            Shift1 = dev_data->Shift1 + Hyb3_Node_Offset(1,0,0,t,tree_data);
            for (i = Bottom1; i <= Top1; i++)
            {
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                quu = dev_data->quu + offset; 
                q0u = dev_data->q0u + offset;
                qdu = dev_data->qdu + offset;
                qu0 = dev_data->qu0 + offset;
                q00 = dev_data->q00 + offset;
                qd0 = dev_data->qd0 + offset;
                qud = dev_data->qud + offset;
                q0d = dev_data->q0d + offset;
                qdd = dev_data->qdd + offset;

                Shift2 = dev_data->Shift2 + offset;



                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;


                                    
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    DiscountL = Discount_3D + offset;

                    Quu = quu[j]; Qu0 = qu0[j]; Qud = qud[j];
                    Q0u = q0u[j]; Q00 = q00[j]; Q0d = q0d[j];
                    Qdu = qdu[j]; Qd0 = qd0[j]; Qdd = qdd[j];

                    r = dev_data->r + offset; 

                    Shift3    = dev_data->Shift3 + offset;
                    NewPriceL = dev_data->Aux3D  + offset;

                    Price00 = (double *)Price + Hyb3_Node_Offset(3,i0,j0,t+1,tree_data);
                    Price01 = (double *)Price + Hyb3_Node_Offset(3,i0,j1,t+1,tree_data);
                    Price02 = (double *)Price + Hyb3_Node_Offset(3,i0,j2,t+1,tree_data);
                    Price10 = (double *)Price + Hyb3_Node_Offset(3,i1,j0,t+1,tree_data);
                    Price11 = (double *)Price + Hyb3_Node_Offset(3,i1,j1,t+1,tree_data);
                    Price12 = (double *)Price + Hyb3_Node_Offset(3,i1,j2,t+1,tree_data);
                    Price20 = (double *)Price + Hyb3_Node_Offset(3,i2,j0,t+1,tree_data);
                    Price21 = (double *)Price + Hyb3_Node_Offset(3,i2,j1,t+1,tree_data);
                    Price22 = (double *)Price + Hyb3_Node_Offset(3,i2,j2,t+1,tree_data);

                    kLowerNext =                 OutBottom3[i0][j0] ;
                    kLowerNext = MAX(kLowerNext, OutBottom3[i0][j1]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i0][j2]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i1][j0]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i1][j1]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i1][j2]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i2][j0]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i2][j1]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i2][j2]);

                    kUpperNext =                 OutTop3[i0][j0] ;
                    kUpperNext = MIN(kUpperNext, OutTop3[i0][j1]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i0][j2]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i1][j0]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i1][j1]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i1][j2]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i2][j0]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i2][j1]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i2][j2]);                    

                    tx = dev_data->t3 - kLowerNext;

                    for (k=kLowerNext; k<=kUpperNext; ++k)
                    {
                        tx[k] = Quu * Price00[k] +
                                Qu0 * Price01[k] +
                                Qud * Price02[k] +
                                Q0u * Price10[k] +
                                Q00 * Price11[k] +
                                Q0d * Price12[k] +
                                Qdu * Price20[k] +
                                Qd0 * Price21[k] +
                                Qdd * Price22[k];
                    }
                        
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {        
                        k0 = (k1 = (k2 = k + Shift3[k] - 1) + 1) + 1;

                        NewPriceL[k] = r[k].u * tx[k0] +
                                       r[k].m * tx[k1] +
                                       r[k].d * tx[k2];

                        NewPriceL[k] *= DiscountL[k];
                    }  /* For k */
                } /* For j */
            } /* For i */

            for (i = Bottom1; i <= Top1; i++)
            {        
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);

                    PriceL    = (double *)Price + offset;
                    NewPriceL = dev_data->Aux3D + offset;
                        
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {        
                        PriceL[k] = NewPriceL[k];
                    } 
                }
            }

        }
        else if (DMode == DISC_3D_CUPS)
        {
            /* Discounted expected value */
            Shift1 = dev_data->Shift1 + Hyb3_Node_Offset(1,0,0,t,tree_data);
            for (i = Bottom1; i <= Top1; i++)
            {
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                quu = dev_data->quu + offset; 
                q0u = dev_data->q0u + offset;
                qdu = dev_data->qdu + offset;
                qu0 = dev_data->qu0 + offset;
                q00 = dev_data->q00 + offset;
                qd0 = dev_data->qd0 + offset;
                qud = dev_data->qud + offset;
                q0d = dev_data->q0d + offset;
                qdd = dev_data->qdd + offset;

                Shift2 = dev_data->Shift2 + offset;

                DiscountL = Discount_2D + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;

                    DiscountLij = DiscountL[j];
                                    
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);

                    Quu = quu[j]; Qu0 = qu0[j]; Qud = qud[j];
                    Q0u = q0u[j]; Q00 = q00[j]; Q0d = q0d[j];
                    Qdu = qdu[j]; Qd0 = qd0[j]; Qdd = qdd[j];

                    r = dev_data->r + offset; 

                    Shift3    = dev_data->Shift3 + offset;
                    NewPriceL = dev_data->Aux3D  + offset;

                    Price00 = (double *)Price + Hyb3_Node_Offset(3,i0,j0,t+1,tree_data);
                    Price01 = (double *)Price + Hyb3_Node_Offset(3,i0,j1,t+1,tree_data);
                    Price02 = (double *)Price + Hyb3_Node_Offset(3,i0,j2,t+1,tree_data);
                    Price10 = (double *)Price + Hyb3_Node_Offset(3,i1,j0,t+1,tree_data);
                    Price11 = (double *)Price + Hyb3_Node_Offset(3,i1,j1,t+1,tree_data);
                    Price12 = (double *)Price + Hyb3_Node_Offset(3,i1,j2,t+1,tree_data);
                    Price20 = (double *)Price + Hyb3_Node_Offset(3,i2,j0,t+1,tree_data);
                    Price21 = (double *)Price + Hyb3_Node_Offset(3,i2,j1,t+1,tree_data);
                    Price22 = (double *)Price + Hyb3_Node_Offset(3,i2,j2,t+1,tree_data);

                    kLowerNext =                 OutBottom3[i0][j0] ;
                    kLowerNext = MAX(kLowerNext, OutBottom3[i0][j1]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i0][j2]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i1][j0]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i1][j1]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i1][j2]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i2][j0]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i2][j1]);
                    kLowerNext = MAX(kLowerNext, OutBottom3[i2][j2]);

                    kUpperNext =                 OutTop3[i0][j0] ;
                    kUpperNext = MIN(kUpperNext, OutTop3[i0][j1]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i0][j2]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i1][j0]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i1][j1]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i1][j2]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i2][j0]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i2][j1]);
                    kUpperNext = MIN(kUpperNext, OutTop3[i2][j2]);                    

                    tx = dev_data->t3 - kLowerNext;

                    for (k=kLowerNext; k<=kUpperNext; ++k)
                    {
                        tx[k] = Quu * Price00[k] +
                                Qu0 * Price01[k] +
                                Qud * Price02[k] +
                                Q0u * Price10[k] +
                                Q00 * Price11[k] +
                                Q0d * Price12[k] +
                                Qdu * Price20[k] +
                                Qd0 * Price21[k] +
                                Qdd * Price22[k];
                    }
                        
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {        
                        k0 = (k1 = (k2 = k + Shift3[k] - 1) + 1) + 1;

                        NewPriceL[k] = r[k].u * tx[k0] +
                                       r[k].m * tx[k1] +
                                       r[k].d * tx[k2];

                        NewPriceL[k] *= DiscountLij;
                    }  /* For k */
                } /* For j */
            } /* For i */

            for (i = Bottom1; i <= Top1; i++)
            {        
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);

                    PriceL    = (double *)Price + offset;
                    NewPriceL = dev_data->Aux3D + offset;
                        
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {        
                        PriceL[k] = NewPriceL[k];
                    } 
                }
            }

        }       /* END of 1-2-3 Factor Hyb3_Dev if */

        status = SUCCESS;

      RETURN:

        return(status);

}  /* Hyb3_Dev */


/*****  Hyb3_Ev  *******************************************************************/
/*
 *       Expected value function.(no discounting)
 *       Used for stats on payoff
 *
 *       Negative probs p_i get modified as follows: 
 *       q_i = | p_i | / Sum_i | p_i |
 *       This procedure removes numerical instabilities 
 *       in the calculation of exercise probabilities and fugits 
 *
 */

int Hyb3_Ev(TSLICE           Price,       /* (I/O) Values to be EV'd       */
            int              t,           /* (I) Current time point        */
            int              T,           /* (I) Last time point           */
            int              DMode,       /* (I) Nb of dims of dev required*/
            HYB3_DEV_DATA   *dev_data,    /* (I) Hyb3_Dev data structure   */
            HYB3_TREE_DATA  *tree_data)   /* (I) Tree data structure       */
{

   /*  Locals for convenience in addressing */
    

    double       *PriceL;
    double       *NewPriceL;

    double       *Price0;  
    double       *Price1;
    double       *Price2;
    double       *Price00;
    double       *Price01;
    double       *Price02;
    double       *Price10; 
    double       *Price11;
    double       *Price12;
    double       *Price20; 
    double       *Price21; 
    double       *Price22;

    /*  1-D probability slices    */  
    TPROB_0  *p;

    /*  2-D probability slices    */
    double   *quu,  *qu0,  *qud;
    double   *q0u,  *q00,  *q0d;
    double   *qdu,  *qd0,  *qdd;

    TPROB_0  *q;  

    /*  3-D probability slices    */
    TPROB_0  *r;  

    double   *tx;

    /* 1-D probabilities (values) */
    double    Pu,    P0,    Pd;

    /* 2-D probabilities (values) */
    double    quuL,   qu0L,  qudL,
              q0uL,   q00L,  q0dL,
              qduL,   qd0L,  qddL;

    double    Quu,   Qu0,  Qud,
              Q0u,   Q00,  Q0d,
              Qdu,   Qd0,  Qdd;

    /* 3-D probabilities (values) */
    double    Ru,    R0,    Rd;

    /* normalisation factor used in case of negative probs */
	double    NormFac;
    
	
    int
        Top1,      Bottom1,  /* Limits of the tree (1rst dimension)          */
        *Top2,     *Bottom2, /* Limits of the tree (2nd dimension)           */
        **Top3,    **Bottom3,/* Limits of the tree (3rd dimension)           */
        NTop1,     NBottom1, /* Limits of the tree at next period (1rst dim) */
        *NTop2,    *NBottom2,/* Limits of the tree at next period (2nd dim)  */
        **NTop3,   **NBottom3,/* Limits of the tree at next period (3rd dim) */
        OutTop1,   OutBottom1,/* Limits of the tree at next period (1rst dim)*/
        *OutTop2,  *OutBottom2,/* Limits of the tree at next period (2nd dim)*/
        **OutTop3, **OutBottom3,/* Limits of the tree at next period (3rd dim)*/


        *Shift1,  *Shift2,  *Shift3,

        i, i0, i1, i2, /* Node indices                                 */
        j, j0, j1, j2,
        k, k0, k1, k2;

    int
            offset,
            
            kLowerNext,
            kUpperNext;
            
    int
            status = FAILURE;
        


    if (t == T) /* Nothing to do at the back of the tree */
    {
        return(SUCCESS);
    }

        
    Top1    = tree_data->Top1[t];	        
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];	    
    Bottom1 = tree_data->Bottom1[t];            
    Bottom2 = tree_data->Bottom2[t];    
    Bottom3 = tree_data->Bottom3[t];
                
    NTop1    = tree_data->Top1[t+1];
    NTop2    = tree_data->Top2[t+1];
    NTop3    = tree_data->Top3[t+1];
    NBottom1 = tree_data->Bottom1[t+1];
    NBottom2 = tree_data->Bottom2[t+1];
    NBottom3 = tree_data->Bottom3[t+1];

    OutTop1    = tree_data->OutTop1[t+1];	        
    OutTop2    = tree_data->OutTop2[t+1];
    OutTop3    = tree_data->OutTop3[t+1];	    
    OutBottom1 = tree_data->OutBottom1[t+1];            
    OutBottom2 = tree_data->OutBottom2[t+1];    
    OutBottom3 = tree_data->OutBottom3[t+1];


    /* Validate DEV mode */
    if (   DMode != DISC_1D_NOCUPS  &&  DMode != DISC_2D_NOCUPS
        && DMode != DISC_2D_CUPS    &&  DMode != DISC_3D_CUPS   )
    {

        DR_Error("Invalid DEV mode requested! (Hyb3_Dev)\n");
        goto RETURN;
    }


    /* Validate DEV dimension for given tree type */
    if (DMode == DISC_1D_NOCUPS)
    {
        if (tree_data->TreeType == TTYPE_1IR2F||
            tree_data->TreeType == TTYPE_2IR2F1D)
        {
            DR_Error("Hyb3_Dev mode requested (1D_NOCUPS) not valid\n"
                    "for two factors foreign IR tree type! (Hyb3_Ev)\n");
            goto RETURN;
        }
    }
    else if (DMode == DISC_2D_CUPS)
    {
        if (tree_data->TreeType == TTYPE_1IR2F||
            tree_data->TreeType == TTYPE_2IR2F1D)
        {
            DR_Error("Hyb3_Dev mode requested (2D_CUPS) not valid\n"
                     "for two factor foreign IR tree type! (Hyb3_Ev)\n");
            goto RETURN;
        }

    }
    else if (DMode == DISC_2D_NOCUPS)
    {
        if (tree_data->TreeType != TTYPE_EQ1IR)
        {
                DR_Error("Hyb3_Dev mode requested (2D_NOCUPS) only valid\n"
                        "for equity plus interest rate tree type! (Hyb3_Ev)\n");
                goto RETURN;
        }
    }
    else if (DMode == DISC_2D_1IR2F_NOCUPS)
    {
        if (tree_data->TreeType != TTYPE_1IR2F&&
            tree_data->TreeType != TTYPE_2IR2F1D)
        {
            DR_Error("Hyb3_Dev mode requested (2D_1IR2F_NOCUPS) only valid\n"
                    "for two factors foreign IR tree type! (Hyb3_Ev)\n");
                goto RETURN;
        }
    }
    else if (DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if (tree_data->TreeType != TTYPE_2IR2F1D)
        {
            DR_Error("Hyb3_Dev mode requested (3D_2IR2F1D_CUPS) only valid\n"
                     "for TTYPE_2IR2F1D tree type! (Hyb3_Ev)\n");
            goto RETURN;
        }
    }
    else
    {

        /* Must be runing in 3-D, therefore TreeType must be 3-D */
        if ((tree_data->TreeType != TTYPE_FX2IR)   &&
            (tree_data->TreeType != TTYPE_EQF2IR)  &&
            (tree_data->TreeType != TTYPE_EQD2IR)  &&
            (tree_data->TreeType != TTYPE_EQC2IR))
        {
                DR_Error("DEV mode chosen would need 3-D tree type! (Hyb3_Dev)\n");
                goto RETURN;
        }

    }

    /*** One factor EV ***/
    if (DMode == DISC_1D_NOCUPS)
    {


        PriceL    = (double *)Price + Hyb3_Node_Offset(1,0,0,t+1,tree_data);

        /* Values of variable are assumed to be flat in the region */
        /* between the inner and the outer ellipses (which are de- */
        /* limited by NTop, OutTop, NBottom and OutBottom)         */

        for (i=OutBottom1; i<NBottom1; i++)
        {
                PriceL[i] = PriceL[NBottom1];
        }
        for (i=NTop1+1; i<=OutTop1; i++)
        {

            PriceL[i] = PriceL[NTop1];
        }

            /* expected value using the non-CUPS probs and shift */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        Shift1 = dev_data->Shift4 + offset;
        p      = dev_data->p + offset;
        
        NewPriceL = dev_data->Aux1D + offset;
            
        for (i=Bottom1; i<=Top1; i++)
        {
            i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

			/* if any transition probabilities are negative, */
			/* change signs and normalise (sum = 1)          */ 			
		    if ( (p[i].u + TINY < 0.0) ||
	  	         (p[i].d + TINY < 0.0) ||
		    	 (p[i].m + TINY < 0.0) )
			{
                  Pu = fabs(p[i].u); Pd = fabs(p[i].d); P0 = fabs(p[i].m);
                  NormFac = Pu + Pd + P0;
 	 	          NormFac = 1.0/NormFac;
		          Pu *= NormFac; Pd *= NormFac; P0 *= NormFac;
			}
            else
			{
                 Pu = p[i].u; Pd = p[i].d; P0 = p[i].m;
			}

            NewPriceL[i] = Pu * PriceL[i0] +
                           P0 * PriceL[i1] +
                           Pd * PriceL[i2];
        }

        /* Put the prices back in the original slice */
        PriceL = (double *)Price + offset;

        for (i=Bottom1; i<=Top1; i++)
        {
            PriceL[i] = NewPriceL[i];
        }
    }

    /*** Two factor EV ***/
    else if ((DMode == DISC_2D_CUPS)  ||
             (DMode == DISC_2D_NOCUPS)||
             (DMode == DISC_2D_1IR2F_NOCUPS))
    {

            
        double    PFlat; 
            
            
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

        Shift1 = dev_data->Shift1 + offset;

        /* Choose between cups and no-cups probability */
        if (DMode == DISC_2D_1IR2F_NOCUPS)
        {
            p = dev_data->s + offset;
        }
        else /* Cups probability */
        {
            p = dev_data->p + offset;
        }

        /* Fill values for nodes between limiting ellipses */
        PriceL    = (double *)Price + Hyb3_Node_Offset(2,NBottom1,0,t+1,tree_data);
        PFlat = PriceL[NTop2[NBottom1]];
        for (i=OutBottom1; i<NBottom1; i++)
        {

            PriceL = (double *)Price + Hyb3_Node_Offset(2,i,0,t+1,tree_data);

            for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }

        for (i=NBottom1; i<=NTop1; i++)
        {
            PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
            PFlat = PriceL[NBottom2[i]];

            for (j=OutBottom2[i]; j<NBottom2[i]; j++)
            {
                PriceL[j] = PFlat;
            }

            PFlat = PriceL[NTop2[i]];
            for (j=NTop2[i]+1; j<=OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }

        PriceL = (double *)Price + Hyb3_Node_Offset(2,NTop1,0,t+1,tree_data);
        PFlat = PriceL[NTop2[NTop1]];
        for (i=NTop1+1; i<=OutTop1; i++)
        {
            PriceL = (double *)Price + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
             
            for (j=OutBottom2[i]; j<=OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }

        /* expected value */
        for (i=Bottom1; i<=Top1; i++)
        {
            i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

            /* Choose between cups and no-cups probability */
            if (DMode == DISC_2D_1IR2F_NOCUPS)
            {
                q = dev_data->t + offset;
            }
            else /* Cups prob. */
            {
                q = dev_data->q + offset;
            }

            Shift2    = dev_data->Shift2 + offset;            
            NewPriceL = dev_data->Aux2D + offset;

            Price0 = (double *)Price + Hyb3_Node_Offset(2,i0,0,t+1,tree_data);
            Price1 = (double *)Price + Hyb3_Node_Offset(2,i1,0,t+1,tree_data);
            Price2 = (double *)Price + Hyb3_Node_Offset(2,i2,0,t+1,tree_data);

            Pu = p[i].u;
            P0 = p[i].m;
            Pd = p[i].d;

            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;

                quuL = Pu * q[j].u;
                qu0L = Pu * q[j].m;
                qudL = Pu * q[j].d;
                q0uL = P0 * q[j].u;
                q00L = P0 * q[j].m;
                q0dL = P0 * q[j].d;
                qduL = Pd * q[j].u;
                qd0L = Pd * q[j].m;
                qddL = Pd * q[j].d;

				/* if any transition probabilities are negative, */
				/* change signs and normalise (sum = 1)          */ 
				if (  (quuL + TINY < 0.0) ||
                      (qu0L + TINY < 0.0) ||
					  (qudL + TINY < 0.0) ||
					  (q0uL + TINY < 0.0) ||
					  (q00L + TINY < 0.0) ||
					  (q0dL + TINY < 0.0) ||
					  (qduL + TINY < 0.0) ||
					  (qd0L + TINY < 0.0) ||
					  (qddL + TINY < 0.0) )
				{
                      Quu = fabs(quuL); Qu0 = fabs(qu0L); Qud = fabs(qudL);
				      Q0u = fabs(q0uL); Q00 = fabs(q00L); Q0d = fabs(q0dL);
                      Qdu = fabs(qduL); Qd0 = fabs(qd0L); Qdd = fabs(qddL);

			          NormFac = Quu + Qu0 + Qud +
                                Q0u + Q00 + Q0d +
                                Qdu + Qd0 + Qdd;

				      NormFac = 1.0/NormFac;
			          Quu *= NormFac; Qu0 *= NormFac; Qud *= NormFac;
			          Q0u *= NormFac; Q00 *= NormFac; Q0d *= NormFac;
			          Qdu *= NormFac; Qd0 *= NormFac; Qdd *= NormFac;
				}
				else
				{
                      Quu = quuL; Qu0 = qu0L; Qud = qudL;
				      Q0u = q0uL; Q00 = q00L; Q0d = q0dL;
                      Qdu = qduL; Qd0 = qd0L; Qdd = qddL;
				}

                NewPriceL[j] = Quu*Price0[j0] + Qu0*Price0[j1] + Qud*Price0[j2] +
                               Q0u*Price1[j0] + Q00*Price1[j1] + Q0d*Price1[j2] +
                               Qdu*Price2[j0] + Qd0*Price2[j1] + Qdd*Price2[j2];                
            }
        }

        /* Finally transfer values from auxiliary slice to variable slice*/
        for (i=Bottom1; i<=Top1; i++)
        {
            
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

            PriceL    = (double *)Price           + offset;
            NewPriceL = dev_data->Aux2D + offset;

            for (j=Bottom2[i]; j<=Top2[i]; j++)
            {
                PriceL[j] = NewPriceL[j];
            }
        }
    }
        
    /*** Three factor EV ***/
    else if (DMode == DISC_3D_CUPS||
             DMode == DISC_3D_2IR2F1D_CUPS)
    {


        /* expected value */
        Shift1 = dev_data->Shift1 + Hyb3_Node_Offset(1,0,0,t,tree_data);
        for (i = Bottom1; i <= Top1; i++)
        {
            i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

            quu = dev_data->quu + offset; 
            q0u = dev_data->q0u + offset;
            qdu = dev_data->qdu + offset;
            qu0 = dev_data->qu0 + offset;
            q00 = dev_data->q00 + offset;
            qd0 = dev_data->qd0 + offset;
            qud = dev_data->qud + offset;
            q0d = dev_data->q0d + offset;
            qdd = dev_data->qdd + offset;

            Shift2 = dev_data->Shift2 + offset;                
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {                    
                j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;
                                    
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);

				/* if any transition probabilities are negative, */
				/* change signs and normalise (sum = 1)          */ 
				if (  (quu[j] + TINY < 0.0) ||
                      (qu0[j] + TINY < 0.0) ||
					  (qud[j] + TINY < 0.0) ||
					  (q0u[j] + TINY < 0.0) ||
					  (q00[j] + TINY < 0.0) ||
					  (q0d[j] + TINY < 0.0) ||
					  (qdu[j] + TINY < 0.0) ||
					  (qd0[j] + TINY < 0.0) ||
					  (qdd[j] + TINY < 0.0) )
				{
                      Quu = fabs(quu[j]); Qu0 = fabs(qu0[j]); Qud = fabs(qud[j]);
				      Q0u = fabs(q0u[j]); Q00 = fabs(q00[j]); Q0d = fabs(q0d[j]);
                      Qdu = fabs(qdu[j]); Qd0 = fabs(qd0[j]); Qdd = fabs(qdd[j]);

			          NormFac = Quu + Qu0 + Qud +
                                Q0u + Q00 + Q0d +
                                Qdu + Qd0 + Qdd;

				      NormFac = 1.0/NormFac;
			          Quu *= NormFac; Qu0 *= NormFac; Qud *= NormFac;
			          Q0u *= NormFac; Q00 *= NormFac; Q0d *= NormFac;
			          Qdu *= NormFac; Qd0 *= NormFac; Qdd *= NormFac;
				}
				else
				{
                      Quu = quu[j]; Qu0 = qu0[j]; Qud = qud[j];
				      Q0u = q0u[j]; Q00 = q00[j]; Q0d = q0d[j];
                      Qdu = qdu[j]; Qd0 = qd0[j]; Qdd = qdd[j];
				}

                r = dev_data->r + offset; 

                Shift3    = dev_data->Shift3 + offset;
                NewPriceL = dev_data->Aux3D  + offset;
                    
                Price00 = (double *)Price + Hyb3_Node_Offset(3,i0,j0,t+1,tree_data);
                Price01 = (double *)Price + Hyb3_Node_Offset(3,i0,j1,t+1,tree_data);
                Price02 = (double *)Price + Hyb3_Node_Offset(3,i0,j2,t+1,tree_data);
                Price10 = (double *)Price + Hyb3_Node_Offset(3,i1,j0,t+1,tree_data);
                Price11 = (double *)Price + Hyb3_Node_Offset(3,i1,j1,t+1,tree_data);
                Price12 = (double *)Price + Hyb3_Node_Offset(3,i1,j2,t+1,tree_data);
                Price20 = (double *)Price + Hyb3_Node_Offset(3,i2,j0,t+1,tree_data);
                Price21 = (double *)Price + Hyb3_Node_Offset(3,i2,j1,t+1,tree_data);
                Price22 = (double *)Price + Hyb3_Node_Offset(3,i2,j2,t+1,tree_data);

                kLowerNext =                 OutBottom3[i0][j0] ;
                kLowerNext = MAX(kLowerNext, OutBottom3[i0][j1]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i0][j2]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i1][j0]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i1][j1]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i1][j2]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i2][j0]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i2][j1]);
                kLowerNext = MAX(kLowerNext, OutBottom3[i2][j2]);

                kUpperNext =                 OutTop3[i0][j0] ;
                kUpperNext = MIN(kUpperNext, OutTop3[i0][j1]);
                kUpperNext = MIN(kUpperNext, OutTop3[i0][j2]);
                kUpperNext = MIN(kUpperNext, OutTop3[i1][j0]);
                kUpperNext = MIN(kUpperNext, OutTop3[i1][j1]);
                kUpperNext = MIN(kUpperNext, OutTop3[i1][j2]);
                kUpperNext = MIN(kUpperNext, OutTop3[i2][j0]);
                kUpperNext = MIN(kUpperNext, OutTop3[i2][j1]);
                kUpperNext = MIN(kUpperNext, OutTop3[i2][j2]);                    

                tx = dev_data->t3 - kLowerNext;

                for (k=kLowerNext; k<=kUpperNext; ++k)
                {
                    tx[k] = Quu * Price00[k] +
                            Qu0 * Price01[k] +
                            Qud * Price02[k] +
                            Q0u * Price10[k] +
                            Q00 * Price11[k] +
                            Q0d * Price12[k] +
                            Qdu * Price20[k] +
                            Qd0 * Price21[k] +
                            Qdd * Price22[k];
                }
                    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {        
                    k0 = (k1 = (k2 = k + Shift3[k] - 1) + 1) + 1;

					if ( (r[k].u + TINY < 0.0) ||
						 (r[k].d + TINY < 0.0) ||
						 (r[k].m + TINY < 0.0) )
					{
                        Ru = fabs(r[k].u); Rd = fabs(r[k].d); R0 = fabs(r[k].m);
  			            NormFac  = Ru + Rd + R0;					
                        NormFac  = 1.0/NormFac;
			            Ru  *= NormFac; Rd  *= NormFac; R0  *= NormFac;
					}
					else 
					{
                        Ru = r[k].u; Rd = r[k].d; R0 = r[k].m;
					}

                    NewPriceL[k] = Ru * tx[k0] +
                                   R0 * tx[k1] +
                                   Rd * tx[k2];
                }  /* For k */
            } /* For j */
        } /* For i */

        for (i = Bottom1; i <= Top1; i++)
        {        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);

                PriceL    = (double *)Price + offset;
                NewPriceL = dev_data->Aux3D + offset;
                    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {        
                    PriceL[k] = NewPriceL[k];
                } 
            }
        }

    } /* END of 1-2-3 Factor Hyb3_Dev if */

        

        status = SUCCESS;

      RETURN:

        return(status);



    
}  /* NormEv */
