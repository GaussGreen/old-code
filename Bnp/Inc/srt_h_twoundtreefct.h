#ifndef SRT_H_TWOUNDTREFCT_H
#define SRT_H_TWOUNDTREFCT_H

/* ------------------------------------------------------------------------------------------
                                     PART I
                           CONSTRUCTION OF THE TREE
   ------------------------------------------------------------------------------------------
 */
/* Set up the tree limits at each step */
Err twound_trelim(SrtStpPtr stp, SrtDiagTwoFacTreeInfo *maxtrinf,
                  SrtMdlType mdl_type1, SrtMdlType mdl_type2);

/* ----------------------------------------------------------------------------------
                                    PART II
                             POPULATING A NODE IN THE TREE
   ----------------------------------------------------------------------------------
 */
void populate_twound_tree_node(
    SrtMdlType mdl_type1, SrtMdlType mdl_type2, SrtStpPtr stp,
    SrtPentoTreeNodeInfo *node, /* output will be returned in node->forward */
    double dt,                  /* delta t between this step and next */
    SrtDiagTwoFacTreeInfo
        *nxttrinf); /* Info about tree structure at next time step */

/* ------------------------------------------------------------------------------------------
                                     PART III
                           DISCOUNTING IN THE TREE
   ------------------------------------------------------------------------------------------
 */
void twound_tree_expectation(
    SrtPentoTreeNodeInfo *node, /* Node info */
    double ***next_assets,      /* Next assets [x1][x2][i] */
    double *cur_assets,         /* Current assets at node */
    long nb_assets);            /* Number of assets */

#endif
