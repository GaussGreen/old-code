/* ---------------------------------------------------------------------

   FILENAME : srt_f_lsmcsv.c

   PURPOSE : implementation of LSM algo for the Cheyette beta stoch vol

   AUTHOR : Eric Fournie

   DATE: SEPTEMBRE 1999

 ----------------------------------------------------------------------- */
#include "grf_h_all.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lsmcsv.h"
#include "srt_h_lsmcsvfct.h"

/* --------------------------------------------------------------------- */

#define SWAP(X, Y)                                                             \
  {                                                                            \
    void *_VPTR;                                                               \
    _VPTR = X;                                                                 \
    X = Y;                                                                     \
    Y = _VPTR;                                                                 \
  }

#define FREE_LSMCSV_MEMORY                                                     \
  {                                                                            \
    for (nt = 0; nt < num_time_pts; nt++)                                      \
      free(grid[nt]);                                                          \
    free(grid);                                                                \
    free_f3tensor(bwd_assets, 0, num_time_pts - 1, 0,                          \
                  grfnparams->num_MCarlo_paths - 1, 0,                         \
                  grfnparams->node_dim - 1);                                   \
    free_f3tensor(fwd_assets, 0, num_time_pts - 1, 0,                          \
                  grfnparams->num_MCarlo_paths - 1, 0,                         \
                  grfnparams->node_dim - 1);                                   \
    free_f3tensor(cur_assets, 0, num_time_pts - 1, 0,                          \
                  grfnparams->num_MCarlo_paths - 1, 0,                         \
                  grfnparams->node_dim - 1);                                   \
    free_dvector(cvalue, 0, grfnparams->node_dim - 1);                         \
    free_imatrix(datex, 0, grfnparams->num_MCarlo_paths - 1, 0,                \
                 grfnparams->node_dim - 1);                                    \
    free(cur_val);                                                             \
  }

/* -------------------------------------------------------------------------------
 */

Err srt_f_csv_lsm(SrtGrfnParam *grfnparams, /*grfn parameters*/
                  SrtStpPtr stp,            /* step pointer*/
                  GrfnDeal *gd,             /* GRFN deal description*/
                  EvalEventFct evalcf,      /* cash-flow evaluation function*/
                  SrtIOStruct *iolist,      /* list of requests*/
                  SrtUndInfo *und_info      /* underlying info*/
) {
  Err err = NULL;
  SrtUndPtr und1;
  SrtStpPtr top, bot;
  SrtSample sam;                  /* sample for evalcf */
  SrtIRMTmInf *tminf = NULL;      /* info for a particular time step*/
  SrtListAtom *io_request = NULL; /* For the storage of the premium ... */
  SrtIOVal *io_request_val_p = NULL;

  long nt, i, k, num_time_pts;
  double cash_flow;

  /* Arrays used to store the Grfn columns for all paths in forward-backward */
  double ***fwd_assets = NULL, ***cur_assets = NULL, ***bwd_assets = NULL;
  double *cvalue;
  int **datex;

  /* A massive grid to store all paths */
  MCTreePointCheyBetaStochVol **grid = NULL;

  SrtMdlType mdl_type1;
  int nreg = 13; /* basis size */
  double *cur_val = (double *)srt_calloc(grfnparams->node_dim, sizeof(double));

  /* -------------------- STEP 1: INITIALISATION ---------------------------- */

  if (und_info->no_of_underlyings != 1)
    return serror("Wrong number of und for americam mc");

  /* Set top (1st step)  , bot (last step) and stp (set to top) and create index
   */
  top = stp = gototop(stp);
  bot = gotobot(stp);
  num_time_pts = create_index(top) + 1;

  /* Allocate space for the time info to each step */
  err = srtstptminfalloc(top, und_info->no_of_underlyings);
  if (err)
    return err;

  /* Allocates and populates time_info at each step.
     Gets the models from the FIRST (== DOMESTIC) underlying
     and initialise the time info in the steps accordingly */
  und1 = lookup_und(und_info->und_data[0].und_name);
  err = get_underlying_mdltype(und1, &mdl_type1);
  err = srt_f_irministp(top, und1, 0, und1, und_info);
  if (err)
    return err;

  /* -------------- STEP 2 : PATHS GENERATION AND STORAGE ------------------- */

  /* Memory allocation for the asset grids ( npaths * Grfn dim ) */
  fwd_assets =
      f3tensor(0, num_time_pts - 1, 0, grfnparams->num_MCarlo_paths - 1, 0,
               grfnparams->node_dim - 1);
  cur_assets =
      f3tensor(0, num_time_pts - 1, 0, grfnparams->num_MCarlo_paths - 1, 0,
               grfnparams->node_dim - 1);
  bwd_assets =
      f3tensor(0, num_time_pts - 1, 0, grfnparams->num_MCarlo_paths - 1, 0,
               grfnparams->node_dim - 1);
  cvalue = dvector(0, grfnparams->node_dim - 1);
  datex =
      imatrix(0, grfnparams->num_MCarlo_paths - 1, 0, grfnparams->node_dim - 1);

  /* Memory allocation for the grid */
  grid = (MCTreePointCheyBetaStochVol **)srt_calloc(
      num_time_pts, sizeof(MCTreePointCheyBetaStochVol *));
  for (nt = 0; nt < num_time_pts; nt++)
    grid[nt] = (MCTreePointCheyBetaStochVol *)calloc(
        grfnparams->num_MCarlo_paths, sizeof(MCTreePointCheyBetaStochVol));

  if (!bwd_assets || !fwd_assets || !cur_assets || !grid)
    return serror("Big memory allocation failure in srt_f_lsmcsv");

  /* Generate Monte Carlo paths  , and store everything */
  err =
      cheybetastochvol_simul(grfnparams, top, und_info, grid, num_time_pts - 1);
  if (err) {
    srt_free(grid);
    return err;
  }

  /* -------------------- STEP 3 : FORWARD  --------------------------- */

  if ((grfnparams->lsm == SRT_FOR) || (grfnparams->lsm == SRT_FORBACK)) {
    for (i = 0; i < grfnparams->num_MCarlo_paths; i++) {

      for (nt = 0, stp = top; nt < num_time_pts; nt++, stp = stp->next) {
        /* Evaluate cash flows in Grfn tableau using passed function evalcf */
        sam_get(sam, 0, STATEVAR) = grid[nt][i].x;
        sam_get(sam, 0, PHI) = grid[nt][i].phi;
        sam_get(sam, 0, SHORT_RATE) = grid[nt][i].short_rate;
        sam.numeraire = grid[nt][i].numeraire;

        err = grfn_eval_event_ammc(
            0, (GrfnEvent *)stp->e, &sam, gd, (double *)cur_assets[nt][i],
            (double *)fwd_assets[nt][i], (double *)NULL,
            (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info, &cash_flow);
        /*
        err = grfn_eval_event((GrfnEvent *) stp->e  ,
                        &sam  , gd  , (double *) cur_assets[nt][i]  ,
                        (EvalEventDfsFct) srt_f_calc_grfn_event_dfs  ,
                        und_info  , &cash_flow);
        */

        if (err) {
          FREE_LSMCSV_MEMORY;
          return err;
        }

        for (k = 0; k < grfnparams->node_dim; k++)
          fwd_assets[nt][i][k] =
              (nt == 0) ? cur_assets[nt][i][k]
                        : fwd_assets[nt - 1][i][k] / grid[nt][i].df +
                              cur_assets[nt][i][k];

      } /* END of loop on time steps */

      for (k = 0; k < grfnparams->node_dim; k++) {
        bwd_assets[num_time_pts - 1][i][k] = fwd_assets[num_time_pts - 1][i][k];
        cur_val[k] += bwd_assets[num_time_pts - 1][i][k] / sam.numeraire;
      }

    } /* END of loop on paths */

    for (k = 0; k < grfnparams->node_dim; k++)
      cur_val[k] /= grfnparams->num_MCarlo_paths;
  }

  /* -------------------- STEP 4 : BACKWARD   --------------------------- */

  if ((grfnparams->lsm == SRT_BACK) || (grfnparams->lsm == SRT_FORBACK)) {
    srt_f_alloc_lsmcsv_resources(nreg, grfnparams->node_dim);

    /* LOOP ON TIME STEPS */
    for (stp = bot, nt = num_time_pts - 1; stp->prev != NULL;
         stp = stp->prev, nt--) {
      /* LSM ALORITHM */
      if (stp->next != NULL) {
        double temp1, temp2, temp3, temp4, temp5, temp6, temp7, x1, x2,
            m1 = 0.0, m2 = 0.0, v11 = 0.0, v12 = 0.0, v22 = 0.0, d11, d12, d21,
            d22, s11, s12, s21, s22, det;

        /* Compute mean-variance */
        for (i = 0; i < grfnparams->num_MCarlo_paths; i++) {
          x1 = grid[nt][i].x;
          x2 = grid[nt][i].sigma;
          m1 += x1;
          m2 += x2;
          v11 += x1 * x1;
          v12 += x1 * x2;
          v22 += x2 * x2;
        }
        m1 /= grfnparams->num_MCarlo_paths;
        m2 /= grfnparams->num_MCarlo_paths;
        v11 /= grfnparams->num_MCarlo_paths;
        v12 /= grfnparams->num_MCarlo_paths;
        v22 /= grfnparams->num_MCarlo_paths;
        v11 = (v11 - m1 * m1);
        v12 = (v12 - m1 * m2);
        v22 = (v22 - m2 * m2);

        temp1 = sqrt(4 * v12 * v12 + (v22 - v11) * (v22 - v11));
        temp2 = sqrt(1.0 + pow(v22 - v11 + temp1, 2) / (4 * v12 * v12));
        temp3 = 2 * sqrt(2) * v12 * temp2;
        temp4 = sqrt(2) * temp2;
        temp5 = sqrt(1.0 + pow(v22 - v11 - temp1, 2) / (4 * v12 * v12));
        temp6 = 2 * sqrt(2) * v12 * temp5;
        temp7 = sqrt(2) * temp5;
        d11 = sqrt(v11 + v22 - temp1) * (v11 - v22 - temp1) / temp3;
        d21 = sqrt(v11 + v22 - temp1) / temp4;
        d12 = sqrt(v11 + v22 + temp1) * (v11 - v22 + temp1) / temp6;
        d22 = sqrt(v11 + v22 + temp1) / temp7;

        det = d11 * d22 - d12 * d21;
        s11 = d22 / det;
        s21 = -d21 / det;
        s12 = -d12 / det;
        s22 = d11 / det;

        /* Compute discounted asset prices at current node */
        if (grfnparams->lsm == SRT_FORBACK) {
          err = srt_f_lsmcsv_regr(nreg, grfnparams->node_dim,
                                  grfnparams->num_MCarlo_paths, nt, grid,
                                  bwd_assets, m1, m2, s11, s12, s21, s22);
        } else {
          err = srt_f_lsmcsv_regr(nreg, grfnparams->node_dim,
                                  grfnparams->num_MCarlo_paths, nt, grid,
                                  cur_assets, m1, m2, s11, s12, s21, s22);
        }
        if (err) {
          FREE_LSMCSV_MEMORY
          srt_f_free_lsmcsv_resources(nreg, grfnparams->node_dim);
          return err;
        }
      }

      /* LOOP FOR THE CONSTRAIN PART */
      for (i = 0; i < grfnparams->num_MCarlo_paths; i++) {
        /* Evaluate cash flows in Grfn tableau */
        sam_get(sam, 0, STATEVAR) = grid[nt][i].x;
        sam_get(sam, 0, PHI) = grid[nt][i].phi;
        sam_get(sam, 0, SHORT_RATE) = grid[nt][i].short_rate;
        sam.numeraire = grid[nt][i].numeraire;

        if (grfnparams->lsm == SRT_FORBACK) {
          if (stp != bot)
            for (k = 0; k < grfnparams->node_dim; k++)
              cvalue[k] = bwd_assets[nt][i][k];

          err = grfn_eval_event_ammc(
              1, (GrfnEvent *)stp->e, &sam, gd, (double *)bwd_assets[nt][i],
              (double *)fwd_assets[nt][i], (double *)cur_assets[nt][i],
              (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info, &cash_flow);

          if (stp != bot) {
            for (k = 0; k < grfnparams->node_dim; k++)
              if (bwd_assets[nt][i][k] != cvalue[k])
                datex[i][k] = nt;
          } else
            for (k = 0; k < grfnparams->node_dim; k++)
              datex[i][k] = nt;
        } else {
          err = grfn_eval_event(
              (GrfnEvent *)stp->e, &sam, gd, (double *)cur_assets[nt][i],
              (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info, &cash_flow);
        }

        if (err) {
          FREE_LSMCSV_MEMORY
          srt_f_free_lsmcsv_resources(nreg, grfnparams->node_dim);
          return err;
        }
      }

    } /* END of backward induction */

    for (k = 0; k < grfnparams->node_dim; k++) {
      cur_val[k] = 0.0;
      if (grfnparams->lsm == SRT_FORBACK) {
        for (i = 0; i < grfnparams->num_MCarlo_paths; i++)
          cur_val[k] +=
              bwd_assets[datex[i][k]][i][k] / grid[datex[i][k]][i].numeraire;
      } else {
        for (i = 0; i < grfnparams->num_MCarlo_paths; i++)
          cur_val[k] += (grid[0][i].df * cur_assets[1][i][k]);
      }
      cur_val[k] /= grfnparams->num_MCarlo_paths;
    }

    srt_f_free_lsmcsv_resources(nreg, grfnparams->node_dim);

  } /* End if backward */

  /* -------------------- STEP 4 : RETURN RESULTS  ---------------------------
   */

  /* Stores the premium in the Input/Output list */
  err = srt_f_IOstructsetpremium(iolist, SRT_NO,
                                 cur_val[grfnparams->node_dim - 1], "");
  if (err) {
    FREE_LSMCSV_MEMORY
    return (err);
  }

  /* Stores all the columns PV in the I/O list */
  err = srt_f_IOstructsetcolpvs(iolist, SRT_NO, (double *)cur_val,
                                grfnparams->node_dim, "");
  if (err) {
    FREE_LSMCSV_MEMORY
    return err;
  }

  /* Free the memory allocated for the arrays of assets  , the mc grid */
  FREE_LSMCSV_MEMORY
  return err;

} /* END of	*/

#undef FREE_LSMCSV_MEMORY
#undef SWAP

/*
for ( i = 0; i < grfnparams->num_MCarlo_paths; i++)
        cur_val[k] +=  (grid[0][i].df * bwd_assets[1][i][k]);
*/

/*--------------------------------- End of File
 * -------------------------------------*/