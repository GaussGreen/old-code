/* ----------------------------------------------------------------------------------------------
   FILENAME:   srt_f_twoundtreefct.cxx


   PURPOSE:  general functions	for lgm/bs trees
   ----------------------------------------------------------------------------------------------
 */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_twoundtreefct.h"

#define TREE_MESH_SPACING                                                      \
  1.73205080756888 /* theoretically > 4  1.73205080756888 */
#define TREE_MESH_SPACING_SQUARE 3.0
#define ONE_QUARTER 0.25
#define ONE_HALF 0.5
#define TREE_LIM_IN_STDEV 5
#define NODE_UP 2
#define NODE_MID 1
#define NODE_DOWN 0

static long srt_minimum(long a, long b) {
  if (a > b)
    return b;
  else
    return a;
}
static long srt_maximum(long a, long b) {
  if (a > b)
    return a;
  else
    return b;
}

/* ----------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------
                                     PART I
                           CONSTRUCTION OF THE TREE
   ------------------------------------------------------------------------------------------
 */
/*	Populates trinf and maxtrinf */
/*  At each time step:
        -	The dual basis is computed: eigenvectors of local cov matrix *
   sqrt (eigenvalues) GLOBAL standard deviations in each direction
*/
Err twound_trelim(SrtStpPtr stp, SrtDiagTwoFacTreeInfo *maxtrinf,
                  SrtMdlType mdl_type1, SrtMdlType mdl_type2) {
  Err err;
  SrtStpPtr cur;                /* Current time step */
  SrtDiagTwoFacTreeInfo *trinf; /* trinf at current step */
  double global_var[2],
      global_cov,        /* Cov matrix seen from today for current step */
      global_var_dir[2], /* Global variance in direction of eigenvectors */
      local_var[2], local_cov,
      local_correl, /* Cov matrix seen from previous step for current step */
      det,          /* Det of dual basis used for inversion */
      std, lim, temp1, temp2, temp3, temp4, temp5, temp6,
      temp7, /* Temps used for diagonalisation */
      V1, V2, C, rho, dt;
  int i, max_trunc_index, /* Trimming */
      max_max_index[2],   /* Largest indexation in the tree */
      max_num_nodes;
  double GG;
  SrtIRMTmInf *irmtminf1;
  SrtIRMTmInf *irmprvtminf1;
  SrtIRMTmInf *irmtminf2;
  SrtIRMTmInf *irmprvtminf2;
  SrtLogTmInf *logtminf2;
  SrtLogTmInf *logprvtminf2;

  /* Go to first step */
  cur = gototop(stp);

  /* Allocates tree info pointer ( 1 stands for  trinf (0=tminf)        , 0:
   * index of the underlying ) */
  if (err = srtstpalloc(cur, sizeof(SrtDiagTwoFacTreeInfo), 1, 0))
    return err;

  /* First step */
  cur = gototop(stp);
  trinf = (SrtDiagTwoFacTreeInfo *)cur->trinf;

  /* Initialisation of cumulative variance variables */
  global_var[0] = global_var[1] = global_cov = GG = 0.0;

  /* No change of basis on first step: Transfer matrices are identioty */
  trinf->dual_basis[0][0] = trinf->dual_basis[1][1] = 1.0;
  trinf->dual_basis[0][1] = trinf->dual_basis[1][0] = 0.0;
  trinf->inverse_basis[0][0] = trinf->inverse_basis[1][1] = 1.0;
  trinf->inverse_basis[0][1] = trinf->inverse_basis[1][0] = 0.0;

  /* Set Spacing and max_index to 0 */
  trinf->spacing[0] = trinf->spacing[1] = 0.0;
  trinf->max_index[0] = trinf->max_index[1] = 0;

  /* Init maximum indexation */
  max_max_index[0] = max_max_index[1] = 0;

  /* Next steps: Loop on steps */
  cur = cur->next;
  while (cur) {
    dt = cur->prev->delta_t;

    /* Get current trinf and tminf and previous tminf for each underlying */
    trinf = (SrtDiagTwoFacTreeInfo *)cur->trinf;

    /* The correlation at this time step */
    rho = cur->prev->correl[0][1];

    /* Compute Global and Local covariances in original basis */
    if (mdl_type1 == LGM) {
      irmtminf1 = (SrtIRMTmInf *)(cur->tminf[0]);
      irmprvtminf1 = (SrtIRMTmInf *)(cur->prev->tminf[0]);

      /* LGM variance from the previous step to this one */
      local_var[0] = irmprvtminf1->ev.twof[0].sig2 * dt;

      /* For LGM        , the global covariance of the short rate if Phi ... */
      global_var[0] = sam_get(irmtminf1->fwd_sam, 0, PHI);

      if (mdl_type2 == LGM) {
        irmtminf2 = (SrtIRMTmInf *)(cur->tminf[1]);
        irmprvtminf2 = (SrtIRMTmInf *)(cur->prev->tminf[1]);

        /* LGM variance from the previous step to this one */
        local_var[1] = irmprvtminf2->ev.twof[0].sig2 * dt;

        /* For LGM        , the global covariance of the short rate if Phi ...
         */
        global_var[1] = sam_get(irmtminf2->fwd_sam, 1, PHI);

        /* The local Covariance LGM/LGM from the previous step to this one */
        local_cov = rho * irmprvtminf1->ev.twof[0].sig *
                    irmprvtminf2->ev.twof[0].sig * dt;

        /* For the Computation of the Global Covariance */
        GG += rho *
              (irmprvtminf1->ev.twof[0].sig * irmprvtminf2->ev.twof[0].sig) /
              (irmprvtminf1->ev.onef.F * irmprvtminf2->ev.onef.F) * dt;

        /* The Global Covariance : the previous sum multiplied by the product of
         * F's */
        global_cov = irmtminf1->ev.onef.F * irmtminf2->ev.onef.F * GG;
      } else if ((mdl_type2 == BLACK_SCHOLES) ||
                 (mdl_type2 == EQ_STOCH_RATES)) {
        logtminf2 = (SrtLogTmInf *)(cur->tminf[1]);
        logprvtminf2 = (SrtLogTmInf *)(cur->prev->tminf[1]);

        /* For a BS like model : the Local Variance from the previous step to
         * this one */
        local_var[1] = logprvtminf2->int_sig2_dt;

        /* The global covariance is the simple integral of sigma^2 */
        global_var[1] = logtminf2->vol2_cum;

        /* The local covariance LGM/BS from the previous step to this one */
        local_cov =
            rho * irmprvtminf1->ev.twof[0].sig * logprvtminf2->int_sig_dt;

        /* For the computation of the Global covariance */
        GG += rho * (irmprvtminf1->ev.twof[0].sig * logprvtminf2->int_sig_dt) /
              irmprvtminf1->ev.onef.F;

        /* The Global Covariance : the previous sum multiplied by the LGM F */
        global_cov = irmtminf1->ev.onef.F * GG;
      }

    } else if ((mdl_type1 == BLACK_SCHOLES) || (mdl_type1 == EQ_STOCH_RATES)) {
      SrtLogTmInf *logtminf1 = (SrtLogTmInf *)(cur->tminf[0]);
      SrtLogTmInf *logprvtminf1 = (SrtLogTmInf *)(cur->prev->tminf[0]);

      /* For a BS like model : the Local Variance from this step to the next one
       */
      local_var[0] = logprvtminf1->int_sig2_dt;

      /* The global covariance is the simple integral of sigma^2 */
      global_var[0] = logtminf1->vol2_cum;

      if (mdl_type2 == LGM) {
        SrtIRMTmInf *tminf2 = (SrtIRMTmInf *)(cur->tminf[1]);
        SrtIRMTmInf *prvtminf2 = (SrtIRMTmInf *)(cur->prev->tminf[1]);

        /* LGM variance from the previous step to this one */
        local_var[1] = prvtminf2->ev.twof[0].sig2 * dt;

        /* For LGM        , the global covariance of the short rate if Phi ...
         */
        global_var[0] = sam_get(irmtminf2->fwd_sam, 1, PHI);

        /* The local Covariance BS/LGM from the previous step to this one */
        local_cov = rho * logprvtminf1->int_sig_dt * prvtminf2->ev.twof[0].sig;

        /* For the computation of the Global covariance */
        GG += rho * (prvtminf2->ev.twof[0].sig * logprvtminf1->int_sig_dt) /
              prvtminf2->ev.onef.F;

        /* The Global Covariance : the previous sum multiplied by the LGM F */
        global_cov = tminf2->ev.onef.F * GG;
      } else if ((mdl_type2 == BLACK_SCHOLES) ||
                 (mdl_type1 == EQ_STOCH_RATES)) {
        SrtLogTmInf *tminf2 = (SrtLogTmInf *)(cur->tminf[1]);
        SrtLogTmInf *prvtminf2 = (SrtLogTmInf *)(cur->prev->tminf[1]);

        /* For a BS like model : the Local Variance from this step to the next
         * one */
        local_var[1] = prvtminf2->int_sig2_dt;

        /* The global covariance is the simple integral of sigma^2 */
        global_var[1] += local_var[1];

        /* The local covariance BS/BS from the previous step to this one */
        local_cov = rho * logprvtminf1->int_sig_dt * prvtminf2->int_sig_dt / dt;

        /* The global variance is just the integral of the local variance (no
         * mean reversion) */
        global_cov += local_cov;
      }

    } /* END if ((mdl_type1 == BLACK_SCHOLES)  || (mdl_type1 == EQ_STOCH_RATES))
       */

    /* Check for numerical errors */
    local_correl = local_cov / sqrt(local_var[0] * local_var[1]);
    if (local_correl > 1.0 - EPS) {
      smessage("Warning: local correl > 1 in 2d Tree");
      local_correl = 1.0 - EPS;
      local_cov = local_correl * sqrt(local_var[0] * local_var[1]);
    } else if (local_correl < -1.0 + EPS) {
      smessage("Warning: local correl < -1 in 2d Tree");
      local_correl = -1.0 + EPS;
      local_cov = local_correl * sqrt(local_var[0] * local_var[1]);
    }

    /* Diagonalise Local var matrix        , compute sqrt(cov) */
    V1 = local_var[0];
    V2 = local_var[1];
    C = local_cov;

    /* Covariance is significative: diagonalise and multiply by eigenvalues */
    if (fabs(C) > EPS) {
      temp1 = sqrt(4.0 * C * C + (V2 - V1) * (V2 - V1));
      temp2 = sqrt(1.0 + pow(V2 - V1 + temp1, 2.0) / (4.0 * C * C));
      temp3 = 2.0 * sqrt(2.0) * C * temp2;
      temp4 = sqrt(2.0) * temp2;
      temp5 = sqrt(1.0 + pow(V2 - V1 - temp1, 2.0) / (4.0 * C * C));
      temp6 = 2.0 * sqrt(2.0) * C * temp5;
      temp7 = sqrt(2.0) * temp5;
      trinf->dual_basis[0][0] =
          sqrt(V1 + V2 - temp1) * (V1 - V2 - temp1) / temp3;
      trinf->dual_basis[1][0] = sqrt(V1 + V2 - temp1) / temp4;
      trinf->dual_basis[0][1] =
          sqrt(V1 + V2 + temp1) * (V1 - V2 + temp1) / temp6;
      trinf->dual_basis[1][1] = sqrt(V1 + V2 + temp1) / temp7;
    } else
    /*	Covariance is already diagonal */
    {
      trinf->dual_basis[0][0] = sqrt(V1);
      trinf->dual_basis[0][1] = 0.0;
      trinf->dual_basis[1][0] = 0.0;
      trinf->dual_basis[1][1] = sqrt(V2);
    }

    /* Compute Inverse basis = coordinates of x1 and x2 in terms of eigenvectors
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

      /* Compute standard deviations and TREE_LIM_IN_STDEVs (TREE_LIM_IN_STDEV *
       * std) */
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
  maxtrinf->max_index[0] = max_max_index[0];
  maxtrinf->max_index[1] = max_max_index[1];

  /* Return error if maxnode is exceeded */
  max_num_nodes = (2 * max_max_index[0] + 1) * (2 * max_max_index[1] + 1);
  if (max_num_nodes > MAXNODE)
    return serror("--- MAXNODE exceeded in 2d tree ---");

  return NULL;

} /* END Err twound_trelim(...) */

/* ----------------------------------------------------------------------------------*/

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
        *nxttrinf) /* Info about tree structure at next time step */
{
  int i;
  double x, y, xx, yy;

  /* Compute forwards of original statevars */
  if (mdl_type1 == LGM)
    xx = sam_get(node->cur_sam, 0, STATEVAR) *
         exp(-((SrtIRMTmInf *)stp->tminf[0])->ev.onef.lambda * dt);
  else {
    if (mdl_type2 == LGM) {
      SrtIRMTmInf *tminf_irm = (SrtIRMTmInf *)stp->tminf[1];
      xx = sam_get(node->cur_sam, 0, STATEVAR) +
           (sam_get(node->cur_sam, 1, SHORT_RATE) -
            sam_get(tminf_irm->fwd_sam, 1, F_0_t)) *
               dt;
    } else
      xx = sam_get(node->cur_sam, 0, STATEVAR);
  }

  if (mdl_type2 == LGM) {
    yy = (sam_get(node->cur_sam, 1, STATEVAR) + node->quanto_fwd) *
         exp(-((SrtIRMTmInf *)stp->tminf[1])->ev.onef.lambda * dt);
  } else {
    if (mdl_type1 == LGM) {
      SrtIRMTmInf *tminf_irm = (SrtIRMTmInf *)stp->tminf[0];
      yy = sam_get(node->cur_sam, 1, STATEVAR) +
           (sam_get(node->cur_sam, 0, SHORT_RATE) -
            sam_get(tminf_irm->fwd_sam, 0, F_0_t)) *
               dt;
    } else
      yy = sam_get(node->cur_sam, 1, STATEVAR);
    yy += node->quanto_fwd;
  }

  /* Transfer in the dual basis: compute forwards of dual statevars */
  node->forward[0] =
      nxttrinf->inverse_basis[0][0] * xx + nxttrinf->inverse_basis[0][1] * yy;
  node->forward[1] =
      nxttrinf->inverse_basis[1][0] * xx + nxttrinf->inverse_basis[1][1] * yy;

  /* Loops on the two dimensions */
  for (i = 0; i <= 1; i++) {
    /* Connects mid to the closest node to forward and rebuilds the statevar
     * accordingly */
    xx = (node->forward[i] < 0.0) ? -0.5 : 0.5;

    node->son_index[i][NODE_MID] =
        (int)(node->forward[i] / nxttrinf->spacing[i] + xx);
    node->son_statevar[i][NODE_MID] =
        node->son_index[i][1] * nxttrinf->spacing[i];
    node->son_index[i][NODE_MID] = srt_maximum(
        srt_minimum(node->son_index[i][NODE_MID], nxttrinf->max_index[i]),
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
  node->p[1] =
      ONE_HALF * (xx + x * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;
  node->p[3] =
      ONE_HALF * (xx - x * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;
  node->p[2] =
      ONE_HALF * (yy + y * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;
  node->p[4] =
      ONE_HALF * (yy - y * TREE_MESH_SPACING) / TREE_MESH_SPACING_SQUARE;

  node->p[0] = 1.0 - node->p[1] - node->p[2] - node->p[3] - node->p[4];

} /* END void populate_twound_tree_node(...) */

/* ----------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------------
                                     PART III
                           DISCOUNTING IN THE TREE
   ------------------------------------------------------------------------------------------
 */

void twound_tree_expectation(
    SrtPentoTreeNodeInfo *node, /* Node info */
    double ***next_assets,      /* Next assets [x1][x2][i] */
    double *cur_assets,         /* Current assets at node */
    long nb_assets)             /* Number of assets */
{
  int i;

  for (i = 0; i < nb_assets; i++) {
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
  }
}

/* ----------------------------------------------------------------------------------*/

#undef TREE_LIM_IN_STDEV
#undef TREE_MESH_SPACING_SQUARE
#undef TREE_MESH_SPACING
#undef ONE_HALF
#undef ONE_QUARTER
#undef NODE_UP
#undef NODE_MID
#undef NODE_DOWN

/* ========= END OF FILE =================================================== */
