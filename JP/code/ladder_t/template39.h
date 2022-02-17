/****************************************************************************/
/*      Function templates for 3 factor model of callable ladder            */
/****************************************************************************/
/*      TEMPLATE39.h                                                        */
/****************************************************************************/

/*
$Header$
*/


#include "ladder.h"


#ifdef __cplusplus
int     Ladder_Manager 
            (T_CURVE         *t_curve,       /* (O) Zero curve data        */
             MKTVOL_DATA     *mktvol_data,   /* (O) Vol data               */
             LADDER_DATA     *ladder_data,   /* (O) Deal data              */
             FIX3_TREE_DATA  *tree_data,     /* (O) Tree data              */
             std::istream&   is);

#else
int     Ladder_Manager 
            (T_CURVE         *t_curve,       /* (O) Zero curve data        */
             MKTVOL_DATA     *mktvol_data,   /* (O) Vol data               */
             LADDER_DATA     *ladder_data,   /* (O) Deal data              */
             FIX3_TREE_DATA  *tree_data);    /* (O) Tree data              */
#endif



#ifdef __cplusplus
extern "C" {
#endif

/* calc2.c */
int     Calc_Ladder 
             (MKTVOL_DATA         *mktvol_data,    /* (I) Vol data          */
              LADDER_DATA         *ladder_data,    /* (I) Deal data         */
              FIX3_TREE_DATA      *tree_data,      /* (I) Tree data         */
              OPT_OUT_DATA        *opt_out_data);  /* (I) Output data       */

/* manager2.c */
int     Ladder_PreProcess 
            (T_CURVE         *t_curve,        /* (O) Zero curve data        */
             MKTVOL_DATA     *mktvol_data,    /* (O) Vol data               */
             long            ValueDate,       /* (I) Value date             */
             LADDER_DATA     *ladder_data,    /* (O) Deal data              */
             FIX3_TREE_DATA  *tree_data);     /* (O) Tree data              */

int     Ladder_Check 
           (LADDER_DATA      *ladder_data,    /* (I) Ladder data        */
            long             ValueDate,       /* (I) Value date             */
            FIX3_TREE_DATA   *tree_data);     /* (I) Structure of tree data */

int     Print_Ladder 
            (T_CURVE          *t_curve,           /* (O) Zero curve data    */
             MKTVOL_DATA      *mktvol_data,       /* (O) Vol  data          */
             LADDER_DATA      *ladder_data,       /* (O) Deal data          */
             FIX3_TREE_DATA   *tree_data);        /* (O) Tree data          */

int     Print_Flows_Deal 
            (MKTVOL_DATA      *mktvol_data,
             LADDER_DATA      *ladder_data,
             FIX3_TREE_DATA   *tree_data);

void    Ladder_Data_Init(
             LADDER_DATA      *ladder_data);      /* (I/O) Deal data        */

void Ladder_Data_Free(
             LADDER_DATA      *ladder_data);      /* (I) Deal data          */

/* time2.c */
int     Ladder_Schedule 
            (long              ValueDate,         /* (I) Value date         */
             T_CURVE          *t_curve,           /* (I) Zero curve data    */
             MKTVOL_DATA      *mktvol_data,       /* (I) volatility data    */
             LADDER_DATA      *ladder_data,       /* (I) Structure of deal  */
             FIX3_TREE_DATA   *tree_data);        /* (O) Tree data          */

#ifdef  __cplusplus
}
#endif


