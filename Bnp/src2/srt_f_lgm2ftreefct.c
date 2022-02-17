/* ----------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & VE

   DATE : JULY 98

   FILENAME:  srt_f_lgm2ftreefct.c

   PURPOSE:  utility functions for the implementation the 2 factor lgm tree in
   Grfn

   MODIFICATION:

   ----------------------------------------------------------------------------------
 */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgm2ftreefct.h"

#define TREE_MESH_SPACING                                                      \
  1.73205080756888 /* theoretically > 4  (1.414213562373) */
#define TREE_MESH_SPACING_SQUARE 3.0
#define ONE_QUARTER 0.25
#define ONE_HALF 0.5
#define TREE_LIM_IN_STDEV 5
#define NODE_UP 2
#define NODE_MID 1
#define NODE_DOWN 0

/* ----------------------------------------------------------------------------------*/

/* --------------------------------------------------------------------------------------
                                                                                PART I
                                                                CONSTRUCTION OF
   THE TREE
   --------------------------------------------------------------------------------------
 */
/*	Populates trinf and maxtrinf
    At each time step:
        -	The dual basis is computed: eigenvectors of local cov matrix *
   sqrt (eigenvalues) -	The grid is constructed in dual basis with spacing
   NUMSTDEVINSPACING (in srt_h_trestruct.h) -	The grid is trimmed at
   CUTTREEATSTDEV (srt_h_trestruct.h) GLOBAL standard deviations in each
   direction
*/
Err lgm2f_trelim(SrtStpPtr stp, SrtDiagTwoFacTreeInfo *maxtrinf) {
  Err err;
  SrtStpPtr cur;                /* Current time step */
  SrtDiagTwoFacTreeInfo *trinf; /* trinf at current step */
  SrtIRMTmInf *prvtminf,        /* tminf at previous step */
      *tminf;                   /* tminf at current step */
  double global_var[2], /* Cov matrix seen from today for current step */
      global_cov,
      global_var_dir[2], /* Global variance in direction of eigenvectors */
      local_var[2], /* Cov matrix seen from previopus step for current step */
      local_cov, local_correl, det, /* Det of dual basis used for inversion */
      std, lim, temp1, temp2, temp3, temp4, temp5, temp6,
      temp7, /* Temps used for diagonalisation */
      V1, V2, C;
  int i, max_trunc_index, /* Trimming */
      max_max_index[2],   /* Largest indexation in the tree */
      max_num_nodes;

  /* Go to first step */
  cur = gototop(stp);

  /* Allocates tree info pointer ( 1 stands for  trinf (0=tminf)  , 0: index of
   * the underlying ) */
  if (err = srtstpalloc(cur, sizeof(SrtDiagTwoFacTreeInfo), 1, 0))
    return err;

  /* First step */
  cur = gototop(stp);
  trinf = (SrtDiagTwoFacTreeInfo *)cur->trinf;

  /* No change of basis on first step: Transfer matrices are identioty */
  trinf->dual_basis[0][0] = trinf->dual_basis[1][1] = 1.0;
  trinf->dual_basis[0][0] = trinf->dual_basis[1][1] = 0.0;
  trinf->inverse_basis[1][1] = trinf->inverse_basis[0][0] = 1.0;
  trinf->inverse_basis[0][1] = trinf->inverse_basis[1][0] = 0.0;

  /* Set Spacing and max_index to 0 */
  trinf->spacing[0] = trinf->spacing[1] = 0.0;
  trinf->max_index[0] = trinf->max_index[1] = 0;

  /* Init maximum indexation */
  max_max_index[0] = max_max_index[1] = 0;

  /* Next step  , Loop on steps */
  cur = cur->next;
  while (cur) {

    /* Get current trinf and tminf and previous tminf */
    trinf = (SrtDiagTwoFacTreeInfo *)cur->trinf;
    tminf = (SrtIRMTmInf *)cur->tminf[0];
    prvtminf = (SrtIRMTmInf *)cur->prev->tminf[0];

    /* Compute Global and Local covariances in original basis */
    global_var[0] = sam_get(tminf->fwd_sam, 0, PHI1);
    global_var[1] = sam_get(tminf->fwd_sam, 0, PHI2);
    global_cov = sam_get(tminf->fwd_sam, 0, CROSSPHI);

    /* Local */
    local_var[0] = prvtminf->ev.twof[0].sig2 *
                   ((fabs(prvtminf->ev.twof[0].lambda) > EPS)
                        ? (1.0 - exp(-2 * prvtminf->ev.twof[0].lambda *
                                     cur->prev->delta_t)) /
                              (2 * prvtminf->ev.twof[0].lambda)
                        : cur->prev->delta_t);
    local_var[1] = prvtminf->ev.twof[1].sig2 *
                   ((fabs(prvtminf->ev.twof[1].lambda) > EPS)
                        ? (1.0 - exp(-2 * prvtminf->ev.twof[1].lambda *
                                     cur->prev->delta_t)) /
                              (2 * prvtminf->ev.twof[1].lambda)
                        : cur->prev->delta_t);
    local_cov =
        prvtminf->correl_x * prvtminf->ev.twof[0].sig *
        prvtminf->ev.twof[1].sig *
        ((fabs(prvtminf->ev.twof[0].lambda + prvtminf->ev.twof[1].lambda) > EPS)
             ? (1.0 - exp(-(prvtminf->ev.twof[0].lambda +
                            prvtminf->ev.twof[1].lambda) *
                          cur->prev->delta_t)) /
                   (prvtminf->ev.twof[0].lambda + prvtminf->ev.twof[1].lambda)
             : cur->prev->delta_t);

    /* Check for numerical errors */
    local_correl = local_cov / sqrt(local_var[0] * local_var[1]);
    if (local_correl > 1.0 - EPS) {
      smessage("Warning: local correl > 1 in Nono Tree");
      local_correl = 1.0 - EPS;
      local_cov = local_correl * sqrt(local_var[0] * local_var[1]);
    } else if (local_correl < -1.0 + EPS) {
      smessage("Warning: local correl < -1 in Nono Tree");
      local_correl = -1.0 + EPS;
      local_cov = local_correl * sqrt(local_var[0] * local_var[1]);
    }

    /* Diagonalise Local var matrix  , compute sqrt(cov) */
    V1 = local_var[0];
    V2 = local_var[1];
    C = local_cov;
    /* Covariance is significative: need to diagonalise and multiply by
     * eigenvalues */
    if (fabs(C) > EPS) {
      temp1 = sqrt(4 * C * C + (V2 - V1) * (V2 - V1));
      temp2 = sqrt(1.0 + pow(V2 - V1 + temp1, 2) / (4 * C * C));
      temp3 = 2 * sqrt(2) * C * temp2;
      temp4 = sqrt(2) * temp2;
      temp5 = sqrt(1.0 + pow(V2 - V1 - temp1, 2) / (4 * C * C));
      temp6 = 2 * sqrt(2) * C * temp5;
      temp7 = sqrt(2) * temp5;
      trinf->dual_basis[0][0] =
          sqrt(V1 + V2 - temp1) * (V1 - V2 - temp1) / temp3;
      trinf->dual_basis[1][0] = sqrt(V1 + V2 - temp1) / temp4;
      trinf->dual_basis[0][1] =
          sqrt(V1 + V2 + temp1) * (V1 - V2 + temp1) / temp6;
      trinf->dual_basis[1][1] = sqrt(V1 + V2 + temp1) / temp7;
    } else
    /*	covariance is already diagonal */
    {
      trinf->dual_basis[0][0] = sqrt(V1);
      trinf->dual_basis[0][1] = 0.0;
      trinf->dual_basis[1][0] = 0.0;
      trinf->dual_basis[1][1] = sqrt(V2);
    }

    /* compute Inverse basis = coordinates of x1 and x2 in terms of eigenvectors
     */
    det = trinf->dual_basis[0][0] * trinf->dual_basis[1][1] -
          trinf->dual_basis[0][1] * trinf->dual_basis[1][0];
    trinf->inverse_basis[0][0] = trinf->dual_basis[1][1] / det;
    trinf->inverse_basis[1][0] = -trinf->dual_basis[1][0] / det;
    trinf->inverse_basis[0][1] = -trinf->dual_basis[0][1] / det;
    trinf->inverse_basis[1][1] = trinf->dual_basis[0][0] / det;

    /* Compute Spacing */
    trinf->spacing[0] = trinf->spacing[1] = TREE_MESH_SPACING;

    /* Trimming */
    for (i = 0; i < 2; i++) {
      /* Compute global variance in the direction of local eigenvectors */
      global_var_dir[i] = global_var[0] * trinf->inverse_basis[i][0] *
                              trinf->inverse_basis[i][0] +
                          global_var[1] * trinf->inverse_basis[i][1] *
                              trinf->inverse_basis[i][1] +
                          2.0 * global_cov * trinf->inverse_basis[i][0] *
                              trinf->inverse_basis[i][1];

      /* Compute standard deviations and limits (TREE_LIM_IN_STDEV * std) */
      std = sqrt(global_var_dir[i]);
      lim = TREE_LIM_IN_STDEV * std;
      /* Compute max index */
      max_trunc_index = (int)DTOL(lim / trinf->spacing[i]) + 1;

      /* Store the information at the time step tree info */
      trinf->max_index[i] = max_trunc_index;

      /* Update maximum indexation ever met */
      max_max_index[i] = IMAX(max_max_index[i], max_trunc_index);
    }

    cur = cur->next;
  }

  /* Make maxtrinf */
  *maxtrinf = *trinf;
  maxtrinf->max_index[0] = max_max_index[0];
  maxtrinf->max_index[1] = max_max_index[1];

  /* Return error if maxnode is exceeded */
  max_num_nodes = (2 * max_max_index[0] + 1) * (2 * max_max_index[1] + 1);
  if (max_num_nodes > MAXNODE)
    return serror("--- MAXNODE exceeded in nono tree ---");

  return NULL;

} /* END Err lgm2f_trelim(...) */

/* --------------------------------------------------------------------------------------
                                                                                PART II
                                                        POPULATING A NODE IN THE
   TREE
   --------------------------------------------------------------------------------------
 */
void populate_lgm2f_tree_node(
    SrtIRMTmInf *tminf, /* time info  , ie. sigma  , tau and rho at step */
    SrtPentoTreeNodeInfo *node, /* output will be returned in node->forward */
    double dt,                  /* delta t between this step and next */
    SrtDiagTwoFacTreeInfo
        *nxttrinf) /* Info about tree structure at next time step */
{
  int i;
  double x, y, xx, yy;

  /* Compute forwards of original statevars */
  xx = sam_get(node->cur_sam, 0, X1) * exp(-tminf->ev.twof[0].lambda * dt);
  yy = sam_get(node->cur_sam, 0, X2) * exp(-tminf->ev.twof[1].lambda * dt);

  /* Transfer in the dual basis: compute forwards of dual statevars */
  node->forward[0] =
      nxttrinf->inverse_basis[0][0] * xx + nxttrinf->inverse_basis[0][1] * yy;
  node->forward[1] =
      nxttrinf->inverse_basis[1][0] * xx + nxttrinf->inverse_basis[1][1] * yy;

  for (i = 0; i <= 1; i++) {
    /* Connects mid to the closest node to forward and rebuilds the statevar
     * accordingly */
    xx = (node->forward[i] < 0.0) ? -0.5 : 0.5;
    node->son_index[i][NODE_MID] =
        (int)(node->forward[i] / nxttrinf->spacing[i] + xx);
    node->son_statevar[i][NODE_MID] =
        node->son_index[i][NODE_MID] * nxttrinf->spacing[i];
    node->son_index[i][NODE_MID] =
        DMAX(min(node->son_index[i][NODE_MID], nxttrinf->max_index[i]),
             -nxttrinf->max_index[i]);

    node->son_index[i][NODE_UP] =
        (node->son_index[i][NODE_MID] < nxttrinf->max_index[i])
            ? node->son_index[i][NODE_MID] + 1
            : nxttrinf->max_index[i];
    node->son_statevar[i][NODE_UP] =
        node->son_statevar[i][NODE_MID] + nxttrinf->spacing[i];

    node->son_index[i][NODE_DOWN] =
        (node->son_index[i][NODE_MID] > -nxttrinf->max_index[i])
            ? node->son_index[i][NODE_MID] - 1
            : -nxttrinf->max_index[i];
    node->son_statevar[i][NODE_DOWN] =
        node->son_statevar[i][NODE_MID] - nxttrinf->spacing[i];
  }

  /* Computes the 5 connection probabilities */
  x = node->forward[0] - node->son_statevar[0][NODE_MID];
  y = node->forward[1] - node->son_statevar[1][NODE_MID];
  xx = 1.0 + x * x;
  yy = 1.0 + y * y;

  /* Points are like this : 	0: center ; 1 : right (x1 Up/x2 Mid) ; 2 : up ;
   * 3 : left ; 4 : down */
  /*	node->p[1] = ONE_QUARTER*(xx + x*TREE_MESH_SPACING);
          node->p[2] = ONE_QUARTER*(yy + y*TREE_MESH_SPACING);
          node->p[3] = ONE_QUARTER*(xx - x*TREE_MESH_SPACING);
          node->p[4] = ONE_QUARTER*(yy - y*TREE_MESH_SPACING); */

  node->p[1] =
      ONE_HALF * (xx + x * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;
  node->p[3] =
      ONE_HALF * (xx - x * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;
  node->p[2] =
      ONE_HALF * (yy + y * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;
  node->p[4] =
      ONE_HALF * (yy - y * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;

  node->p[0] = 1.0 - node->p[1] - node->p[2] - node->p[3] - node->p[4];

} /* END void populate_lgm2f_tree_node(...) */

/* ----------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------------
                                     PART III
                           DISCOUNTING IN THE TREE
   ------------------------------------------------------------------------------------------
 */
void lgm2f_tree_expectation(SrtPentoTreeNodeInfo *node, /* Node info */
                            double ***next_assets, /* Next assets [x1][x2][i] */
                            double *cur_assets,    /* Current assets at node */
                            long nb_assets)        /* Number of assets */
{
  int i;

  for (i = 0; i < nb_assets; i++)
    cur_assets[i] =
        node->df * (node->p[0] * next_assets[node->son_index[0][NODE_MID]]
                                            [node->son_index[1][NODE_MID]][i] +
                    node->p[1] * next_assets[node->son_index[0][NODE_UP]]
                                            [node->son_index[1][NODE_MID]][i] +
                    node->p[2] * next_assets[node->son_index[0][NODE_MID]]
                                            [node->son_index[1][NODE_UP]][i] +
                    node->p[3] * next_assets[node->son_index[0][NODE_DOWN]]
                                            [node->son_index[1][NODE_MID]][i] +
                    node->p[4] * next_assets[node->son_index[0][NODE_MID]]
                                            [node->son_index[1][NODE_DOWN]][i]);

} /* END void lgm2f_tree_expectation(...) */

#undef TREE_LIM_IN_STDEV
#undef TREE_MESH_SPACING_SQUARE
#undef TREE_MESH_SPACING
#undef ONE_HALF
#undef ONE_QUARTER
#undef NODE_UP
#undef NODE_MID
#undef NODE_DOWN

/* ========= END OF FILE =================================================== */