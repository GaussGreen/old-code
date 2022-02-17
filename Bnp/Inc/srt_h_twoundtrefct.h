#ifndef SRT_H_TWOUNDTREFCT_H
#define SRT_H_TWOUNDTREFCT_H

/* Functions used in a nonomial tree */

/* Calculate delta_x in both x and y directions.  Assume that the same step
   sizes are used throughout the tree.  Results are stored in delta_x[0]
   (for x) and delta_x[1] (for y) */

void twound_calcdeltax(SrtStpPtr top, double *delta_x, SrtMdlType mdl_type1,
                       SrtMdlType mdl_type2);

/* Find the node which is closest to the expected value from the current
   node.  Results stored in node->mid_sons_index[0] and
   node->mid_sons_index[1] */

void twoundtre_xindex(SrtStpPtr stp,   /* current step */
                      double **prev_x, /* state vars at previous step */
                      SrtNonoTreeNodeInfo *node, /* current node in tree */
                      int closestmeth /* CLOSESTINX or CLOSESTINLNX */
);

/* Set up the tree limits at each step */

Err twound_trelim(SrtStpPtr stp, SrtTwoFacTreInf *maxtrinf,
                  SrtMdlType mdl_type1, SrtMdlType mdl_type2);

#endif
