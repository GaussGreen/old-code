/* ==============================================================================

   FILE NAME:      srt_f_cheybetatree.cxx

   MAIN FUNCTION:  srt_f_CheyBetaTree

   OBJECT:         main function for local vol Cheyette tree

  ===============================================================================
*/

/* ------------------------------------------------------------------------------
 */

/* ------------------------- Include files
 * -------------------------------------- */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheybetadynamics.h"
#include "srt_h_cheybetatree.h"
#include "srt_h_cheybetatreefnc.h"
#include "srt_h_irmfct.h"

/* ------------------------------------------------------------------------------
 */

/* ------------------------ Macro definitions
 * ----------------------------------- */

/* ------------------------------------------------------------------------------
 */
#define PSWAP(X, Y)                                                            \
  {                                                                            \
    void *_VPTR = X;                                                           \
    X = Y;                                                                     \
    Y = _VPTR;                                                                 \
  }

#define FREE_CHEY_BETA_TREE_MEMORY                                             \
  {                                                                            \
    free_f3tensor(cur_assets, maxtrinf.min_state_var_index,                    \
                  maxtrinf.max_state_var_index, maxtrinf.min_phi_index,        \
                  maxtrinf.max_phi_index, 0, grfnparam->node_dim - 1);         \
    free_f3tensor(next_assets, maxtrinf.min_state_var_index,                   \
                  maxtrinf.max_state_var_index, maxtrinf.min_phi_index,        \
                  maxtrinf.max_phi_index, 0, grfnparam->node_dim - 1);         \
    free_dvector(state_var_grid, maxtrinf.min_state_var_index,                 \
                 maxtrinf.max_state_var_index);                                \
    free_dvector(cur_phi_grid, maxtrinf.min_phi_index,                         \
                 maxtrinf.max_phi_index);                                      \
    free_dvector(next_phi_grid, maxtrinf.min_phi_index,                        \
                 maxtrinf.max_phi_index);                                      \
  }

/* ------------------------------------------------------------------------------
 */

/* ------------------------- Main function
 * -------------------------------------- */

Err srt_f_CheyBetaTree(SrtUndPtr und, /* Market pointer: the underlying */
                       SrtGrfnParam *grfnparam, /* Model Parameters */
                       SrtStpPtr stp, /* Step pointer (time discretisation)  */
                       GrfnDeal *gd,  /* GRFN Deal description */
                       EvalEventFct evalcf,  /* Cashflow evalution function*/
                       SrtIOStruct *iolist,  /* List of requests */
                       SrtUndInfo *und_info) /* Underlying info */

{

  /* ------------------------- 1.- Declarations
   * ----------------------------------- */

  Err err = NULL;
  int i;
  int j;
  double zc_yield;
  double cashflow;
  /* Time steps: the first one and the last one */
  SrtStpPtr top;
  SrtStpPtr bot;

  /* Information for one node in the tree (for local volatility) */
  SrtTrinTreNdInf node;

  /* Information for the tree at one particular time step (spacing        ,...)
   */
  SrtLocTreInf maxtrinf;
  SrtLocTreInf *trinf;

  /* Model information for a particular time step (vol        ,...) */
  SrtIRMTmInf *tminf;

  /* Assets computed through the tree at cur and next time step
          3 dimensions:	1- R
                                          2- Phi
                                          3- Asset */
  double ***next_assets;
  double ***cur_assets;

  /* State variables at current time step and next time step */
  double *state_var_grid;
  double *cur_phi_grid;
  double *next_phi_grid;
  double next_phi;

  TermStruct *ts;
  SrtMdlType mdl_type;

  /* ------------------------ 2.- Initialisations
   * --------------------------------- */

  /* Get the volatility Term Struct from the underlying */
  err = get_underlying_ts(und, &ts);
  if (err) {
    return err;
  }

  /* Get the mdl type */
  err = get_underlying_mdltype(und, &mdl_type);
  if (err) {
    return err;
  }

  /* Set top (1st step)        , bot (last step) and stp (current step        ,
   * init to top)
   */
  top = stp = gototop(stp);
  bot = gotobot(stp);

  /* Allocate space for one tminf per step (only one underlying) */
  if (err = srtstptminfalloc(top, 1)) {
    return err;
  }

  /* Initialises within stp (all steps in tminf) using info in mkt */
  err = srt_f_irministp(stp, und, 0, und, und_info);
  if (err) {
    return err;
  }

  /* Initialises the tree information: geometry & limits (in maxtrinf) */
  if (err = srt_f_CheyBetaTree_init_trinf(stp, grfnparam, und, &maxtrinf)) {
    return err;
  }

  /* Memory Allocation */
  next_assets = f3tensor(maxtrinf.min_state_var_index,
                         maxtrinf.max_state_var_index, maxtrinf.min_phi_index,
                         maxtrinf.max_phi_index, 0, grfnparam->node_dim - 1);
  cur_assets = f3tensor(maxtrinf.min_state_var_index,
                        maxtrinf.max_state_var_index, maxtrinf.min_phi_index,
                        maxtrinf.max_phi_index, 0, grfnparam->node_dim - 1);
  state_var_grid =
      dvector(maxtrinf.min_state_var_index, maxtrinf.max_state_var_index);

  cur_phi_grid = dvector(maxtrinf.min_phi_index, maxtrinf.max_phi_index);

  next_phi_grid = dvector(maxtrinf.min_phi_index, maxtrinf.max_phi_index);

  if (!next_assets || !cur_assets || !state_var_grid || !cur_phi_grid ||
      !next_phi_grid) {
    FREE_CHEY_BETA_TREE_MEMORY;
    return serror("Memory allocation failure in srt_f_locvoltree");
  }

  /* Construct the grid in the short rate (r) */
  err = srt_f_CheyBetaTree_make_state_var_grid(maxtrinf, state_var_grid);
  if (err) {
    FREE_CHEY_BETA_TREE_MEMORY;
    return err;
  }

  /* -------------------- 3.- TREE: Backward Induction
   * --------------------------- */

  /* For each time stp        , starting from the botom one */
  for (stp = bot; stp != NULL; stp = stp->prev) {
    /* Get trinf and tminf */
    trinf = (SrtLocTreInf *)stp->trinf;
    tminf = (SrtIRMTmInf *)stp->tminf[0];

    /* Create current phi grid */
    err = srt_f_CheyBetaTree_make_phi_grid(trinf, cur_phi_grid);
    if (err) {
      FREE_CHEY_BETA_TREE_MEMORY;
      return err;
    }

    /* Loop on r */
    for (i = trinf->min_state_var_index; i <= trinf->max_state_var_index; i++) {
      /* Fill node info */
      sam_get(node.cxxur_sam, 0, STATEVAR) = state_var_grid[i];
      sam_get(node.cxxur_sam, 0, SHORT_RATE) =
          state_var_grid[i] + sam_get(tminf->fwd_sam, 0, F_0_t);

      /* Loop on phi */
      for (j = trinf->min_phi_index; j <= trinf->max_phi_index; j++) {
        /* Fill node info */
        sam_get(node.cxxur_sam, 0, PHI) = cur_phi_grid[j];

        /* If not at the bottom: discount */
        if (stp->next) {
          /* Compute DF for that stp (from t to t + Dt) */
          Y_T_at_t_compute(1, &node.cxxur_sam, &tminf->yp, &zc_yield, 0,
                           ONE_FAC, mdl_type);
          node.df = exp(-zc_yield);

          /* Compute expected values for t+delta_t at this node */
          err = srt_f_CheyBeta_drift_at_sam(mdl_type, stp, &node.cxxur_sam,
                                            &node.drift_sam, 0);
          if (err)
            return err;
          node.state_var_expectation = sam_get(node.drift_sam, 0, STATEVAR);

          /* Compute the local variance at this node */
          err = srt_f_CheyBeta_var_at_sam(mdl_type, stp, &node.cxxur_sam,
                                          &node.state_var_variance, 0);
          if (err)
            return err;

          /* Recombine in r : find the closest node to the fwd + up & do */
          err = srt_f_CheyBetaTree_recombine_in_r(&node, state_var_grid,
                                                  stp->next->trinf);
          if (err) {
            FREE_CHEY_BETA_TREE_MEMORY;
            return err;
          }

          /* Compute probabilities */
          err = srt_f_CheyBetaTree_calc_prob(&node);
          if (err) {
            FREE_CHEY_BETA_TREE_MEMORY;
            return err;
          }

          /* Compute discounted sum of next cashflows **/
          next_phi = sam_get(node.drift_sam, 0, PHI);
          err = srt_f_CheyBetaTree_discount_assets(
              cur_assets[i][j], next_assets, &node, grfnparam->node_dim,
              stp->next->trinf, state_var_grid, next_phi, next_phi_grid);
          if (err) {
            FREE_CHEY_BETA_TREE_MEMORY;
            return err;
          }
        } /* End if (stp->next) (i.e. if not last step) */

        /* Compute current cash flow */
        err = evalcf((GrfnEvent *)stp->e, &node.cxxur_sam, gd,
                     (double *)cur_assets[i][j],
                     (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info,
                     &cashflow);
        if (err) {
          return err;
        }

      } /* END of loop on phi */

    } /* END of loop on r */

    /* Swaps the assets and phi grids (the new grid becomes the old one) */
    PSWAP(cur_assets, next_assets);
    PSWAP(cur_phi_grid, next_phi_grid);

  } /* END of loop on time steps */

  /* Stores the premium in the Input/Output list */
  err = srt_f_IOstructsetpremium(
      iolist, SRT_NO, next_assets[0][0][grfnparam->node_dim - 1], "");
  if (err) {
    FREE_CHEY_BETA_TREE_MEMORY;
    return (err);
  }

  /* Stores all the columns PV in the I/O list */
  err = srt_f_IOstructsetcolpvs(iolist, SRT_NO, (double *)next_assets[0][0],
                                grfnparam->node_dim, "");
  if (err) {
    FREE_CHEY_BETA_TREE_MEMORY;
    return err;
  }

  /* Return a success message  */
  FREE_CHEY_BETA_TREE_MEMORY;
  return NULL;

} /* END of Err srt_f_locvoltree (...) */

/* ==============================================================================
 */
