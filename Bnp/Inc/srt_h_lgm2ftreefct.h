#ifndef SRT_H_LGM2FTREEFCT_H
#define SRT_H_LGM2FTREEFCT_H

/* --------------------------------------------------------------------------------------
                                                                                PART I
                                                                CONSTRUCTION OF
   THE TREE
   --------------------------------------------------------------------------------------
 */
/*	Populates trinf and maxtrinf
    At each time step:
        - The dual basis is computed: eigenvectors of local cov matrix * sqrt
   (eigenvalues)
        - The grid is constructed in dual basis with spacing NUMSTDEVINSPACING
   (in srt_h_trestruct.h)
        - The grid is trimmed at CUTTREEATSTDEV (srt_h_trestruct.h) GLOBAL
   standard deviations in each direction
*/
Err lgm2f_trelim(SrtStpPtr stp, SrtDiagTwoFacTreeInfo *maxtrinf);

/* --------------------------------------------------------------------------------------
                                                                                PART II
                                                        POPULATING A NODE IN THE
   TREE
   --------------------------------------------------------------------------------------
 */
void populate_lgm2f_tree_node(
    SrtIRMTmInf
        *tminf, /* time info      , ie. sigma      , tau and rho at step */
    SrtPentoTreeNodeInfo *node, /* output will be returned in node->forward */
    double dt,                  /* delta t between this step and next */
    SrtDiagTwoFacTreeInfo
        *nxttrinf); /* Info about tree structure at next time step */

/* ------------------------------------------------------------------------------------------
                                     PART III
                           DISCOUNTING IN THE TREE
   ------------------------------------------------------------------------------------------
 */
void lgm2f_tree_expectation(SrtPentoTreeNodeInfo *node, /* Node info */
                            double ***next_assets, /* Next assets [x1][x2][i] */
                            double *cur_assets,    /* Current assets at node */
                            long nb_assets);       /* Number of assets */

#endif