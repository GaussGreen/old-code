/****************************************************************************/
/*      Function templates for 3 factor model.                              */
/****************************************************************************/
/*      TEMPLATE15.H                                                         */
/****************************************************************************/

/*
$Header$
*/


#include "cupslib.h"
#include "fxkoseries_t.h"



/* main15.c */
int     main (void);

/* manager15.c */
int     Fxkoseries_Manager(T_CURVE [2][3],MKTVOL_DATA *,FX_DATA *,FXKOSERIES_DATA *,HYB3_TREE_DATA *);
int     Fxkoseries_Check (FXKOSERIES_DATA *,FX_DATA *,HYB3_TREE_DATA *);
int     Print_Fxkoseries(T_CURVE  [2][3],HYB3_TREE_DATA  *);


/* time15.c */
int     Fxkoseries_Schedule(long,FXKOSERIES_DATA *,HYB3_TREE_DATA *);

/* calc15.c */
int     Calc_Fxkoseries(MKTVOL_DATA *,FXKOSERIES_DATA *,HYB3_TREE_DATA *,OPT_OUT_DATA *);

