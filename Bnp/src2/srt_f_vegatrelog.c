/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT  , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_LOGTREE.C                                       */
/*                                                                            */
/*      PURPOSE:        Functions to compute bs tree                          */
/*                                                                            */
/*      AUTHORS:        Whoever (KLC  ,OVE...)			              */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "srt_h_all.h"

#define SWAP(X, Y)                                                             \
  {                                                                            \
    void *_VPTR;                                                               \
    _VPTR = X;                                                                 \
    X = Y;                                                                     \
    Y = _VPTR;                                                                 \
  }

/* -------------------------------------------------------------------------- */

static Err
srt_f_logtrelim(SrtStpPtr stp, SrtTreStpInf *maxtrinf,
                double *old_u /* previous value of state var spacing u*/
) {
  Err err;
  SrtStpPtr cur;
  SrtTreStpInf *trinf;
  SrtLogTmInf *tminf;
  double delta_x;

  cur = gototop(stp);

  /* allocate space */

  if (err = srtstpalloc(cur, sizeof(SrtTreStpInf), 1, 0))
    return err;

  if (*old_u ==
      -1) { /* assume here that there is only one underlying : und_ind = 0 */
    srt_f_trintredeltax(cur, &delta_x, 0);
    *old_u = delta_x;
  } else {
    delta_x = *old_u;
  }

  while (cur->next) {
    trinf = cur->trinf;
    tminf = cur->tminf[0];

    trinf->min_x_index = -cur->index;
    trinf->max_x_index = cur->index;
    trinf->u = delta_x;
    trinf->xmin = -cur->index * delta_x;
    /*sam_get(tminf->fwd_sam  ,0  ,STATEVAR)  */
    cur = cur->next;
  }

  trinf = cur->trinf;
  tminf = cur->tminf[0];

  trinf->min_x_index = -cur->index;
  trinf->max_x_index = cur->index;
  trinf->u = delta_x;
  trinf->xmin = -cur->index * delta_x;
  /*sam_get(tminf->fwd_sam  ,0  ,STATEVAR) */
  *maxtrinf = *trinf;
  return NULL;
}

/* -------------------------------------------------------------------------- */

/* This function is inspired by srt_f_vegashifttrelgm1d */
Err srt_f_logtree(SrtUndPtr und, SrtGrfnParam *grfnparam, SrtStpPtr stp,
                  GrfnDeal *gd, EvalEventFct evalcf, SrtIOStruct *iolist,
                  SrtUndInfo *und_info) {

  /* ----- declare stuff ---- */
  Err err = NULL;
  int i;
  double cashflow;
  SrtStpPtr top, bot;

  /* information about the  particular  node within the tree we are at */
  SrtTrinTreNdInf node;

  /* information about the dimensions of this tree */
  SrtTreStpInf maxtrinf, *trinf;

  /* information needed for the model at a particular time + next time*/
  SrtLogTmInf *tminf;

  /* val of assets being computed through the tree at cur and prev time step */
  /* 2D arrays: x dimension  , grfn_cols dimension*/
  double **prev_assets, **cur_assets;

  /* state variable values at current time step and previous time step */
  double *cur_st_var, *prev_st_var;

  /* amt of space for assets */
  long sz;

  /* current and constant value for rate spacing in the tree  , for the prices
  and all the schifted prices */
  double old_u = -1;

  /* number of extra calculation = 1 (for the price) +
                                   number of shifted values requested.
                                  The price can always be computed  */

  /* in order to know which shift corresponds to which calc.value */

  SrtIOVal *io_request_val_p;
  SrtListAtom *io_request;

  /* --------------- end of declarations ---------------------------------------
   */

  top = stp = gototop(stp);
  bot = gotobot(stp);

  /* determine how much space we need to allocate  ,etc;
     attach tree geometry info to steps */

  if (err = srtstptminfalloc(stp, 1))
    return err;

  if (err = srt_f_loginistp(stp, und, 0, und, und_info))
    return err;

  if (err = srt_f_logtrelim(stp, &maxtrinf, &old_u))
    return err;

  /* Memory Allocation -------------------------------------*/

  sz = (maxtrinf.max_x_index - maxtrinf.min_x_index + 1) * grfnparam->node_dim;
  prev_assets = dmatrix(maxtrinf.min_x_index, maxtrinf.max_x_index, 0,
                        grfnparam->node_dim - 1);
  cur_assets = dmatrix(maxtrinf.min_x_index, maxtrinf.max_x_index, 0,
                       grfnparam->node_dim - 1);
  prev_st_var = dvector(maxtrinf.min_x_index, maxtrinf.max_x_index);
  cur_st_var = dvector(maxtrinf.min_x_index, maxtrinf.max_x_index);

  if (!prev_assets || !cur_assets || !prev_st_var || !cur_st_var)
    return serror("allocation failure srt_f_logtree.c");

  memset(&cur_assets[maxtrinf.min_x_index][0], 0, sz * sizeof(double));
  /* Set to zero all the cell of the array */

  /* Tree (going backwards in time) ---------------------------*/

  for (stp = bot; stp != NULL; stp = stp->prev) {
    trinf = (SrtTreStpInf *)(stp->trinf);
    tminf = (SrtLogTmInf *)(stp->tminf[0]);

    /** create vector of state variables **/

    cur_st_var[trinf->min_x_index] = trinf->xmin;

    for (i = trinf->min_x_index + 1; i <= trinf->max_x_index; i++)
      cur_st_var[i] = cur_st_var[i - 1] + trinf->u;

    if (stp->next) {
      node.df = ((SrtLogTmInf *)(stp->next->tminf[0]))->df / tminf->df;
    }
    /** For each value of the state variable == at each node of the tree step
     * **/
    for (i = trinf->min_x_index; i <= trinf->max_x_index; i++) {

      io_request = iolist->head;
      tminf = (SrtLogTmInf *)(stp->tminf[0]);

      sam_get(node.cur_sam, 0, STATEVAR) = cur_st_var[i];
      sam_get(node.cur_sam, 0, SPOT) = exp(cur_st_var[i]) * tminf->init_fwd_val;

      /** if not at the end: **/
      if (stp->next) {

        /** calculate drift at this node **/

        sam_get(node.drift_sam, 0, STATEVAR) =
            sam_get(node.cur_sam, 0, STATEVAR) - 0.5 * tminf->int_sig2_dt;

        /** calculate x index to center on in previous level **/

        srt_f_trintrexindex(stp, prev_st_var, &node, grfnparam, CLOSESTINX);

        /** calculate variance at this node **/

        node.var_at_sam = tminf->int_sig2_dt;

        /** calculate probabilities **/
        if (i == trinf->min_x_index) {
          srt_f_trintreprob(prev_st_var, &node, grfnparam);
        }
        /** calculate discounted prob weighted sum of previous cashflows **/

        /*CHECK*/ srt_f_trintredcntvec(stp, cur_assets[i], &node, prev_assets,
                                       grfnparam->node_dim);
      }

      /** call cashflow function **/

      err = evalcf(
          (GrfnEvent *)stp->e, &node.cur_sam, gd, (double *)cur_assets[i],
          (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info, &cashflow);
      if (err)
        return err;

    } /*end of loop for every node in step*/

    SWAP(cur_assets, prev_assets);
    SWAP(cur_st_var, prev_st_var);
  }

  /*Sets the premium in the IO list */
  for (io_request = iolist->head; io_request != NULL;
       io_request = io_request->next) {
    io_request_val_p = (SrtIOVal *)(io_request->element->val.pval);
    strncpy(io_request_val_p->computation_origin, "Tree_Log",
            strlen("Tree_Log"));
    if (io_request_val_p->type == IO_PREMIUM) {
      io_request_val_p->dval = prev_assets[0][grfnparam->node_dim - 1];
    }
  }

  /* Stores all the columns PV in the I/O list */
  err = srt_f_IOstructsetcolpvs((iolist), SRT_NO, (double *)prev_assets[0],
                                grfnparam->node_dim, "");

  free_dmatrix(prev_assets, maxtrinf.min_x_index, maxtrinf.max_x_index, 0,
               grfnparam->node_dim - 1);
  free_dmatrix(cur_assets, maxtrinf.min_x_index, maxtrinf.max_x_index, 0,
               grfnparam->node_dim - 1);

  free_dvector(prev_st_var, maxtrinf.min_x_index, maxtrinf.max_x_index);
  free_dvector(cur_st_var, maxtrinf.min_x_index, maxtrinf.max_x_index);

  /* return a success message */
  return err;
}

/* Problemo:
probs are computed each time although sometimes it is useless.
*/

/* -------------------------------------------------------------------------- */
