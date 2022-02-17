/*******************************************************************************
**
**              Copyright (c) 1994 PARIBAS Capital Markets Group
**
********************************************************************************
**
**      MODULE NAME:    SRT_F_MC_CORE.C
**
**      SYSTEM:         SRT     SORT        , Fixed Income 2020 Addins
**      SUB_SYSTEM:     GRF     GRFN model in swap_tools
**
**      PURPOSE:        Monte-Carlo for GRFN interest rate model.
**
**      AUTHORS:        Rishin ROY        ,
**                      Adam LITKE
**						Eric AULD
**                      Jasbir MALHI
**
**      DATE:           14th October        , 1994
**
*******************************************************************************/

/* ==========================================================================
   include files
   ========================================================================== */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_grfn_undinfo.h"
#include "srt_h_mc_evolve.h"
#include "srt_h_powermodel.h"
#include "srt_h_stpcorata.h"

#define shift 0.01

#define FREE_MC_CORE                                                           \
  {                                                                            \
    if (grfnparam->calib == SRT_NO)                                            \
      SrtMCRanSrc_end(&mcrs);                                                  \
    if (dfs)                                                                   \
      free_dvector(dfs, 0, num_time_pts);                                      \
    if (sum_disc_cf)                                                           \
      free_dvector(sum_disc_cf, 0, num_time_pts);                              \
    if (Z_sum)                                                                 \
      free_dvector(Z_sum, 0, num_time_pts);                                    \
    if (Pt_by_Zt)                                                              \
      free_dvector(Pt_by_Zt, 0, num_time_pts);                                 \
    if (sam)                                                                   \
      srt_free(sam);                                                           \
    if (col_pvs)                                                               \
      free_dvector(col_pvs, 0, grfnparam->node_dim - 1);                       \
    if (path_cash_flows)                                                       \
      free_dmatrix(path_cash_flows, 0, num_time_pts - 1, 0,                    \
                   grfnparam->node_dim - 1);                                   \
  }

#define FREE_MC_CORE_EXFR                                                      \
  {                                                                            \
    if (grfnparam->exfrontierlast_time == SRT_YES)                             \
      SrtMCRanSrc_end(&mcrs);                                                  \
    if (dfs)                                                                   \
      free_dvector(dfs, 0, num_time_pts);                                      \
    if (sum_disc_cf)                                                           \
      free_dvector(sum_disc_cf, 0, num_time_pts);                              \
    if (Z_sum)                                                                 \
      free_dvector(Z_sum, 0, num_time_pts);                                    \
    if (sam)                                                                   \
      srt_free(sam);                                                           \
    if (col_pvs)                                                               \
      free_dvector(col_pvs, 0, grfnparam->node_dim - 1);                       \
    if (path_cash_flows)                                                       \
      free_dmatrix(path_cash_flows, 0, num_time_pts - 1, 0,                    \
                   grfnparam->node_dim - 1);                                   \
    if (Pt_by_Zt)                                                              \
      free_dvector(Pt_by_Zt, 0, num_time_pts);                                 \
  }

/* ----------------------------------------------------------------
  FUNCNAME        : McCore
  AUTHOR          : J.Malhi        , R.Roy        , E.Auld        , K.Chau
  DESCRIPTION     :computes prices of assets/contingent claims
   over the (n x num_time_pts)-dimensional
   VECTOR of sample structures        , and is called each time for iterating
  through each one of the (num_mcarlo_paths). MODIFIES        : gd        ,price
  ,mcrs CALL            :

  ---------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

static SrtMCRanSrc mcrs;

Err free_MCrandSrc() {
  if (&mcrs)
    SrtMCRanSrc_end(&mcrs);
  return NULL;
}

Err MCCore(SrtGrfnParam *grfnparam, SrtStpPtr step, GrfnDeal *gd,
           EvalEventFct evalcf, void *iolist, SrtUndInfo *und_info) {
  int i, j,             /* index for underlying */
      st_index = 0;     /* index for stochastic underlying */
  int dom_ind, for_ind; /* index for domestic and foreign underlying when quanto
                     adjustment */
  long t,               /* index for time */
      n,                /* index for path */
      num_time_pts;
  double cash_flow, sum_price;
  double *dfs = NULL, *sum_disc_cf = NULL, *Z_sum = NULL, *Pt_by_Zt = NULL;
  SrtSample *sam = NULL;
  SrtIRMTmInf *tminf = NULL;
  Err err = NULL;
  SrtUndPtr und = NULL;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;

  /* For FX_STOCH_RATES_jumping */

  String dom_und_name;
  String for_und_name;
  String numeraire_ccy;
  SrtUndPtr dom_und, for_und, numeraire_und;
  SrtMdlType dom_mdl_type, for_mdl_type;

  SrtStpPtr top;
  long seed = 0;
  long last_index;

  SrtIOStruct *iorequest = iolist;
  double exp_disc_pay = 0.0, exp_disc_pay_sqr = 0.0;
  double stdev = 0.0, temp = 0.0;
  double cash_flow_pv;
  double **deal_cash_flows_pv = NULL;
  double **path_cash_flows = NULL;
  double *col_pvs = NULL;

  /* A few initialisation */
  top = gototop(step);
  last_index = create_index(top);

  seed = grfnparam->rand_seed;
  num_time_pts = create_index(top) + 1;

  /* Initialises mcrs according to inputs        , makes memory allocation or
     even generates fully the non correlated random numbers */
  if ((grfnparam->calib == SRT_NO) || (grfnparam->first_time == SRT_YES)) {
    err = SrtMCRanSrc_start(&mcrs, 0, grfnparam->num_MCarlo_paths - 1, 1,
                            last_index, und_info->no_of_brownians,
                            grfnparam->sample_type, seed, top);
    if (err) {
      FREE_MC_CORE;
      return err;
    }

    if (grfnparam->calib == SRT_YES) {
      err = SrtMCRanSrc_init(&mcrs);
      if (err) {
        FREE_MC_CORE;
        return err;
      }
    }
    grfnparam->first_time = SRT_NO;
  } else
    mcrs.cxxur_path = mcrs.pl - 1;

  /* Initialise the sample structure */
  sam = srt_calloc(num_time_pts, sizeof(SrtSample));
  for (t = 0, step = top; t < num_time_pts; t++, step = step->next) {
    sam[t].numeraire_index = und_info->numeraire_index;
  }

  /* Need to create tminf structure in SrtStpPtr first */
  err = srtstptminfalloc(top, und_info->no_of_underlyings);
  if (err) {
    FREE_MC_CORE;
    return err;
  }

  /* Keep track of the Main Currency (for Quanto adjustments)  (the One of the
   * Numeraire) */
  numeraire_und = lookup_und(und_info->und_data[0].und_name);
  numeraire_ccy = get_underlying_ccy(numeraire_und);

  /* Initialises the steps: stacks all relevant information */
  for (i = 0; i < und_info->no_of_underlyings; i++) {
    /* Gets the Relevant Underlying (in the order of the Grfn Tableau) */
    und = lookup_und(und_info->und_data[i].und_name);

    if (und_info->und_data[i].stochastic == SRT_YES) {
      /* Gets the Model type and dimension */
      err = get_underlying_mdltype(und, &mdl_type);
      if (err) {
        FREE_MC_CORE;
        return err;
      }

      err = get_underlying_mdldim(und, &mdl_dim);
      if (err) {
        FREE_MC_CORE;
        return err;
      }

      /* Initialisation for Black Sholes like diffusions */
      if ((mdl_type == BLACK_SCHOLES) || (mdl_type == NORMAL_BS) ||
          (mdl_type == FX_STOCH_RATES) || (mdl_type == EQ_STOCH_RATES) ||
          (mdl_type == EQ_STOCH_RATES_SRVGS)) {

        if (mdl_type == FX_STOCH_RATES) {

          dom_und_name = get_domname_from_fxund(und);
          dom_und = lookup_und(dom_und_name);
          if (!dom_und)
            return serror("Could not find %s underlying", dom_und_name);
          err = get_underlying_mdltype(dom_und, &dom_mdl_type);
          if (err) {
            FREE_MC_CORE;
            return err;
          }

          for_und_name = get_forname_from_fxund(und);
          for_und = lookup_und(for_und_name);
          if (!for_und)
            return serror("Could not find %s underlying", for_und_name);
          err = get_underlying_mdltype(for_und, &for_mdl_type);
          if (err) {
            FREE_MC_CORE;
            return err;
          }

          /* Force the Jumping if possible  and set in srt_f_grfn_und*/
          if ((dom_mdl_type == LGM) && (for_mdl_type == LGM) &&
              (und_info->jumping == SRT_YES)) {
            err = srt_f_fx_initstp(top, und_info, und, dom_und, for_und, i,
                                   numeraire_und);
            if (err) {
              FREE_MC_CORE;
              return err;
            }
          } else {
            err = srt_f_loginistp(top, und, i, numeraire_und, und_info);
            if (err) {
              FREE_MC_CORE;
              return err;
            }
          } /* END if not jumping with FX STOCH RATES */

        } /* END if(mdl_type == FX_STOCH_RATES) */

        else {
          err = srt_f_loginistp(top, und, i, numeraire_und, und_info);
          if (err) {
            FREE_MC_CORE;
            return err;
          }

          if (mdl_type == EQ_STOCH_RATES_SRVGS)
            mcrs.nbr = und_info->no_of_underlyings;
        }

      } /* END of initilisation for Black Scholes like diffusions */
      else {
        if (err = srt_f_irministp(top, und, i, numeraire_und, und_info)) {
          FREE_MC_CORE;
          return err;
        }
      }

      /* This is required for Y_T_at_t_compute to work properly */

      if ((mdl_type == LGM) || (mdl_type == LGM_STOCH_VOL)) {
        for (t = 0, step = top; t < num_time_pts; t++, step = step->next) {
          tminf = (SrtIRMTmInf *)step->tminf[i];
          switch (mdl_dim) {
          case ONE_FAC:
            sam_get(sam[t], i, PHI) = sam_get(tminf->fwd_sam, i, PHI);
            break;
          case TWO_FAC:
            sam_get(sam[t], i, PHI1) = sam_get(tminf->fwd_sam, i, PHI1);
            sam_get(sam[t], i, PHI2) = sam_get(tminf->fwd_sam, i, PHI2);
            sam_get(sam[t], i, CROSSPHI) = sam_get(tminf->fwd_sam, i, CROSSPHI);
            break;
          default:
            sam_get(sam[t], i, PHI) = 0.0;
          } /* END switch declaration  */

        } /* END of for t loop   - for LGM only */

      } /* END if mdl type = LGM */

    } /* END if stochastic == SRT_YES */

    else /* if stochastic == SRT_NO */
    {
      if (err = srt_f_basicinistp(top, und, i)) {
        FREE_MC_CORE;
        return err;
      }

      /* we only need to set up sam[t].numeriare once if it is deterministic */

      monte_DETERMINISTIC_evolve(top, sam);

    } /* END if stochastic == SRT_NO */

  } /* END for i loop on underlying number*/

  /* Attaches the right correlation coefficients to the right steps        ,
   * when needed  */
  if ((und_info->no_of_underlyings > 1) &&
      (und_info->two_factor_model == SRT_NO) &&
      (und_info->use_stochastic_vol == SRT_NO)) {
    mcrs.need_to_correl = SRT_YES;

    if (err = srt_f_attach_correl_to_stp(top, und_info->corr_ts)) {
      FREE_MC_CORE;
      return err;
    }

  } else {
    mcrs.need_to_correl = SRT_NO;
  }

  dfs = dvector(0, num_time_pts);
  sum_disc_cf = dvector(0, num_time_pts);
  Z_sum = dvector(0, num_time_pts);
  Pt_by_Zt = dvector(0, num_time_pts);
  if ((dfs == NULL) || (sum_disc_cf == NULL) || (Z_sum == NULL) ||
      (Pt_by_Zt == NULL)) {
    FREE_MC_CORE;
    return serror("Memory allocation error inMC_Core");
  }

  /* Pricing the security: sum_price=PRICE */

  /* Initialization */

  cash_flow = 0.0;
  memset(dfs, 0, num_time_pts * sizeof(double));
  memset(sum_disc_cf, 0, num_time_pts * sizeof(double));
  memset(Z_sum, 0, num_time_pts * sizeof(double));
  memset(Pt_by_Zt, 0, num_time_pts * sizeof(double));

  /* If we just want the intrinsic value i.e. no brownian is defined        , we
    just need one path to evaluate the answer.  Therefore we modified the no. of
    Monte Carlo paths below if this is the case.  */
  if (und_info->no_of_brownians == 0)
    grfnparam->num_MCarlo_paths = 1;

  /* Compute discount factors from 0 to the step date for discounting using dom
   * und */
  t = 0;
  step = top;
  und = lookup_und(und_info->und_data[0].und_name);
  err = get_underlying_mdltype(und, &mdl_type);
  if (err) {
    FREE_MC_CORE;
    return err;
  }

  while (step != NULL) {
    switch (mdl_type) {
      /* Treats FX stoch rates separately */
    case FX_STOCH_RATES:
      if (und_info->jumping == SRT_YES) {
        dfs[t] = ((SrtFXTmInf *)step->tminf[0])->df;
        break;
      }
    case EQ_STOCH_RATES:
    case BLACK_SCHOLES:
    case NORMAL_BS:
    case EQ_STOCH_RATES_SRVGS:
      dfs[t] = ((SrtLogTmInf *)step->tminf[0])->df;
      break;
    case NONE:
      dfs[t] = ((SrtBasicTmInf *)step->tminf[0])->df;
      break;
    default:
      dfs[t] = ((SrtIRMTmInf *)step->tminf[0])->df;
      break;
    }
    step = step->next;
    t = t + 1;
  }

  /* Memory allocation for a vector where all the Columns PV's will be stored */
  col_pvs = dvector(0, grfnparam->node_dim - 1);
  deal_cash_flows_pv = dmatrix(0, num_time_pts - 1, 0, grfnparam->node_dim - 1);
  path_cash_flows = dmatrix(0, num_time_pts - 1, 0, grfnparam->node_dim - 1);
  memset(&deal_cash_flows_pv[0][0], 0,
         (int)(num_time_pts * grfnparam->node_dim * sizeof(double)));
  memset(&path_cash_flows[0][0], 0,
         (int)(num_time_pts * grfnparam->node_dim * sizeof(double)));

  /* Starts loop on number of paths */
  for (n = 0; n < grfnparam->num_MCarlo_paths; n++) {

    /* Generates or gets the next MC sample path and puts it in mcrs.r[cur_path]
     */
    if (grfnparam->calib == SRT_NO) {
      if (err = SrtMCRanSrc_next(&mcrs)) {
        FREE_MC_CORE;
        return err;
      }
    } else
      mcrs.cxxur_path++;

    /* Correlates (if necessary) the n(numb of undelyings)
           random (independent) Brownian paths mcrs.r[cur_path][][]        ,
           for all the steps with time dependent correlation matrixes */
    if (grfnparam->calib == SRT_NO) {
      if (err = SrtMCRanSrc_correl(&mcrs, top)) {
        FREE_MC_CORE;
        return err;
      }
    }

    /* Evolves each underlying        , one by one
      We have to evolve first the equity/interest rate underlyings
     then the forex ones (this is due to the stochatic rates
    needed for quanto adjustment. */
    for (j = 0; j < 2; j++) {
      for (i = 0; i < und_info->no_of_underlyings; i++) {
        und = lookup_und(und_info->und_data[i].und_name);
        err = get_underlying_mdltype(und, &mdl_type);
        if (err) {
          FREE_MC_CORE;
          return err;
        }

        if (((j == 0) &&
             (!((mdl_type == EQ_STOCH_RATES) ||
                (mdl_type == EQ_STOCH_RATES_SRVGS))) &&
             (!ISUNDTYPE(und, FOREX_UND))) ||
            ((j == 1) && (ISUNDTYPE(und, FOREX_UND))) ||
            ((j == 1) && (mdl_type == EQ_STOCH_RATES)) ||
            ((j == 1) && (mdl_type == EQ_STOCH_RATES_SRVGS))) {

          /* if underlying is stochastic        , set up sam[t]; otherwise
          sam[t]. pv_money_mkt should already be set up above because it is
          determinstic */

          if (und_info->und_data[i].stochastic == SRT_YES) {
            err = get_underlying_mdltype(und, &mdl_type);
            if (err) {
              FREE_MC_CORE;
              return err;
            }

            err = get_underlying_mdldim(und, &mdl_dim);
            if (err) {
              FREE_MC_CORE;
              return err;
            }

            switch (mdl_type) {
            case NORMAL_BS:
              monte_NORMALBS_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam, i);
              break;
            case BLACK_SCHOLES:
              monte_BLACKSCHOLES_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                        i);
              break;
            case EQ_STOCH_RATES:
              err = get_index_from_und_info(
                  und_info, get_discname_from_underlying(und), &dom_ind);
              if (err) {
                FREE_MC_CORE;
                return err;
              }
              monte_EQ_STOCH_RATES_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                          i, dom_ind);
              break;
            case EQ_STOCH_RATES_SRVGS:
              err = get_index_from_und_info(
                  und_info, get_discname_from_underlying(und), &dom_ind);
              if (err) {
                FREE_MC_CORE;
                return err;
              }

              monte_EQ_STOCH_RATES_SRVGS(
                  mcrs.r[mcrs.cxxur_path][i],
                  mcrs.r[mcrs.cxxur_path][und_info->no_of_brownians - 1], top,
                  sam, i, dom_ind);
              break;
            case FX_STOCH_RATES: /* DONE ONLY FOR FOREX UNDERLYING */
              err = get_index_from_und_info(
                  und_info, get_domname_from_fxund(und), &dom_ind);
              if (err) {
                FREE_MC_CORE;
                return err;
              }
              err = get_index_from_und_info(
                  und_info, get_forname_from_fxund(und), &for_ind);
              if (err) {
                FREE_MC_CORE;
                return err;
              }

              if ((mdl_type == FX_STOCH_RATES) && (dom_mdl_type == LGM) &&
                  (for_mdl_type == LGM) && (und_info->jumping == SRT_YES)) {

                monte_FX_STOCH_RATES_Jumping_evolve(
                    mcrs.r[mcrs.cxxur_path][i], top, sam, i, dom_ind, for_ind);
              } else

                monte_FX_STOCH_RATES_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i, dom_ind, for_ind);

              break;
            case ETABETA:
              if (mdl_dim == TWO_FAC) {
                FREE_MC_CORE;
                return serror("Two Factor ETABETA not implemented yet...");
              } else {
                monte_BETAETA_1f_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                        i);
              }
              break;
            case LGM:
              if (mdl_dim == ONE_FAC) {
                if (und_info->jumping == SRT_YES)
                  monte_LGM_1f_Jumping_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                              sam, i);
                else if (und_info->jumping == SRT_NO)
                  monte_LGM_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i);
              } else if (mdl_dim == TWO_FAC) {
                monte_LGM_2f_evolve(mcrs.r[mcrs.cxxur_path], top, sam, i);
              }
              break;
            case NEWLGM:
              monte_NEWLGM_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                           i);
              break;
            case VASICEK:
              monte_vasicek_1f_euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i);
              break;
            case NEWCHEYBETA:
              monte_NEWCHEYBETA_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                                sam, i);
              break;
            case LGM_STOCH_VOL:
              monte_LGM_1f_STOCHVOL_evolve(mcrs.r[mcrs.cxxur_path], top, sam,
                                           Pt_by_Zt, i);
              break;
            case CHEY_STOCH_VOL:
              monte_CHE_1f_STOCHVOL_evolve(mcrs.r[mcrs.cxxur_path], top, sam,
                                           Pt_by_Zt, i);
              break;
            case CHEY_BETA_STOCH_VOL:
              monte_CHEYBETA_1f_STOCHVOL_evolve(mcrs.r[mcrs.cxxur_path], top,
                                                sam, Pt_by_Zt, i);
              break;
            case CHEY_BETA:
              if (mdl_dim == ONE_FAC) {
                monte_CHEYBETA_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                               sam, i);
              } else if (mdl_dim == TWO_FAC) {
                monte_CHEYBETA_2f_Euler_evolve(mcrs.r[mcrs.cxxur_path], top,
                                               sam, i);
              }
              break;
            case MIXED_BETA:
              if (mdl_dim == ONE_FAC) {
                FREE_MC_CORE;
                return serror("One Factor MIXEDBETA is nonsense");
              } else if (mdl_dim == TWO_FAC) {
                monte_MIXEDBETA_2f_Euler_evolve(mcrs.r[mcrs.cxxur_path], top,
                                                sam, i);
              }
              break;
            case CHEY:
            default:
              if (mdl_dim == ONE_FAC) {
                if (grfnparam->difference_scheme == MILSHTEIN) {
                  monte_CHE_1f_Milshtein_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                                sam, i);
                } else /* EULER */
                {
                  monte_CHE_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i);
                }
              } else if (mdl_dim == TWO_FAC) {
                monte_CHE_2f_Euler_evolve(mcrs.r[mcrs.cxxur_path], top, sam, i);
              }
              break;

            } /* END switch declaration */

          } /* END if stochastic == SRT_YES */

        } /* FOREX or NOT FOREX underlying */

      } /* END for i loop on underlying number*/

    } /* END of the loop on underlying types*/

    /* Loops on all the time steps of this path to compute the cash flow */
    temp = 0.0;
    for (t = 0, step = top; t < num_time_pts; step = step->next, t++) {

      sam[t].pathnum = n;

      /* Call to the main Grfn COMLL evaluation function */

      err = evalcf(
          (GrfnEvent *)step->e, &sam[t], gd, (double *)(path_cash_flows[t]),
          (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info, &cash_flow);
      if (err) {
        FREE_MC_CORE;
        return err;
      }

      /* Here we divide by the numeraire */
      if (mdl_type == VASICEK) {
        cash_flow_pv = cash_flow / sam[t].numeraire;
        sum_disc_cf[t] += cash_flow_pv;
        temp += cash_flow_pv;
      } else {
        /* Here we multiply by dfs[t] that has been removed from the numeraire
           in HJM type model        ,
           (it has been computed with exp(r - f(0        ,t))  or exp(A(t)+h(t
           ,x))
         */
        cash_flow_pv = cash_flow / (sam[t].numeraire / dfs[t]);
        sum_disc_cf[t] += cash_flow_pv;
        temp += cash_flow_pv;
      }

      /* Store the information for all the columns        , in a global deal pv
       * calculation */
      for (i = 0; i < grfnparam->node_dim; i++) {
        if (mdl_type == VASICEK) {
          deal_cash_flows_pv[t][i] += path_cash_flows[t][i] / sam[t].numeraire;
        } else {
          deal_cash_flows_pv[t][i] +=
              path_cash_flows[t][i] / (sam[t].numeraire / dfs[t]);
        }
        /* Resets the path_cash_flow to 0.0 */
        path_cash_flows[t][i] = 0.0;
      }

      /* The zero-coupon price computed for this time step (sum on all paths)
       * (without its forward) */
      Z_sum[t] += 1.0 / sam[t].numeraire;

    } /* END for t loop on all time steps of current path */

    /* Increments the current sums of all cashflows and cash flow squares for
     * all paths */
    exp_disc_pay += temp;
    exp_disc_pay_sqr += temp * temp;

  } /* END for n loop on path number*/

  /* Frees the memory allocated for random numbers generation */
  if (grfnparam->calib == SRT_NO)
    SrtMCRanSrc_end(&mcrs);

  /* Computes the renormalisation factor for the price */
  if (grfnparam->mc_renormalize) {
    /* The correction factor is (at each step) the zero-coupon bond price
     * computed by full MC simulation  */
    for (t = 0; t < num_time_pts; t++)
      Pt_by_Zt[t] = 1.0 / Z_sum[t];

  } else {
    /* The corrector is just the number paths */
    for (t = 0; t < num_time_pts; t++)
      Pt_by_Zt[t] = 1.0 / (double)grfnparam->num_MCarlo_paths;
  }

  /* The Real price is the price obtained by Monte Carlo        , renormalised
   * step by step */
  sum_price = 0.0;
  memset(&col_pvs[0], 0, sizeof(double) * grfnparam->node_dim);
  for (t = 0; t < num_time_pts; t++) {
    sum_price += (Pt_by_Zt[t] * sum_disc_cf[t]);
    for (i = 0; i < grfnparam->node_dim; i++)
      col_pvs[i] += deal_cash_flows_pv[t][i] * Pt_by_Zt[t];
  }

  /* Stores the premium in the Input/Output list */
  err =
      srt_f_IOstructsetpremium((SrtIOStruct *)iorequest, SRT_NO, sum_price, "");
  if (err) {
    FREE_MC_CORE;
    return (err);
  }

  /* Computes and stores the standard deviation as well */
  exp_disc_pay_sqr *=
      1.0 + DBL_EPSILON *
                grfnparam->num_MCarlo_paths; /* To prevent rounding errors */
  stdev = (sqrt(exp_disc_pay_sqr / grfnparam->num_MCarlo_paths -
                exp_disc_pay * exp_disc_pay / grfnparam->num_MCarlo_paths /
                    grfnparam->num_MCarlo_paths)) /
          sqrt(grfnparam->num_MCarlo_paths);
  err = srt_f_IOstructsetstdev((SrtIOStruct *)iorequest, SRT_NO, stdev, "");
  if (err) {
    FREE_MC_CORE;
    return (err);
  }

  /* Stores all the columns PV in the I/O list */
  err = srt_f_IOstructsetcolpvs((SrtIOStruct *)iorequest, SRT_NO, col_pvs,
                                grfnparam->node_dim, "");
  if (err) {
    FREE_MC_CORE;
    return err;
  }

  /* Free all allocated memory */
  FREE_MC_CORE;

  /* Return a success message */
  return NULL;

} /* END srt_f_mccore */

Err MCOptimizeExFrontier(SrtGrfnParam *grfnparam, SrtStpPtr step, GrfnDeal *gd,
                         EvalEventFct evalcf, void *iolist,
                         SrtUndInfo *und_info)

{
  int numex; /*number of exercize dates        , equal to the length of the
                appropriate aux*/
  int exfr;  /*the index of the auxiliary containing the exercise frontier*/
  int n;
  int j;
  double *storeexfr;
  long num_paths;
  SrtIOStruct *iorequest = iolist;
  Err err = NULL;

  /*initialization*/

  exfr = grfnparam->exfrontier;
  numex = gd->auxlen[exfr];
  num_paths = grfnparam->num_MCarlo_paths;
  grfnparam->num_MCarlo_paths = 100;

  /* Memory allocation*/
  storeexfr = dvector(0, numex - 1);

  /* storage and no previous exercise condition*/
  for (j = 0; j < numex; j++) {
    storeexfr[j] = gd->aux[exfr][j];
    gd->aux[exfr][j] = 999;
  }

  /* optimization of the exfrontier using MCCoreExFrontier*/
  for (n = numex - 1; n >= 0; n--) {
    gd->aux[exfr][n] = storeexfr[n];

    if (err = MCCoreExFrontier(grfnparam, step, gd, evalcf, iolist, und_info)) {
      smessage("Error in MCCoreExFrontier");
      free_dvector(storeexfr, 0, numex - 1);
      return err;
    }
  }

  /*Price computation*/
  grfnparam->num_MCarlo_paths = num_paths;

  if (err = MCCore(grfnparam, step, gd, evalcf, iolist, und_info)) {
    smessage("Error in MCCore");
    free_dvector(storeexfr, 0, numex - 1);
    return err;
  }

  /*free*/
  free_dvector(storeexfr, 0, numex - 1);
  return err;
}

Err MCCoreExFrontier(SrtGrfnParam *grfnparam, SrtStpPtr step, GrfnDeal *gd,
                     EvalEventFct evalcf, void *iolist, SrtUndInfo *und_info) {
  int i, j,             /* index for underlying */
      st_index = 0;     /* index for stochastic underlying */
  int dom_ind, for_ind; /* index for domestic and foreign underlying when quanto
                     adjustment */
  long t,               /* index for time */
      n,                /* index for path */
      num_time_pts;
  double cash_flow, sum_price;
  double *dfs = NULL, *sum_disc_cf = NULL, *Z_sum = NULL, *Pt_by_Zt = NULL;
  SrtSample *sam = NULL;
  SrtIRMTmInf *tminf = NULL;
  Err err = NULL;
  SrtUndPtr und = NULL;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;

  /* For FX_STOCH_RATES_jumping */

  String dom_und_name;
  String for_und_name;
  String numeraire_ccy;
  SrtUndPtr dom_und, for_und, numeraire_und;
  SrtMdlType dom_mdl_type, for_mdl_type;

  SrtStpPtr top;
  long seed = 0;
  long last_index;

  SrtIOStruct *iorequest = iolist;
  double exp_disc_pay = 0.0, exp_disc_pay_sqr = 0.0;
  double stdev = 0.0, temp = 0.0;
  double cash_flow_pv;
  double **deal_cash_flows_pv = NULL;
  double **path_cash_flows = NULL;
  double *col_pvs = NULL;
  long nstart, nend;
  /* A few initialisation */
  top = gototop(step);
  last_index = create_index(top);

  seed = grfnparam->rand_seed;
  num_time_pts = create_index(top) + 1;

  /* Initialises mcrs according to inputs        , makes memory allocation or
     even generates fully the non correlated random numbers */
  if (grfnparam->exfrontierfirst_time == SRT_YES) {
    err = SrtMCRanSrc_start(&mcrs, 0, grfnparam->num_MCarlo_paths - 1, 1,
                            last_index, und_info->no_of_brownians,
                            grfnparam->sample_type, seed, top);
    if (err) {
      FREE_MC_CORE_EXFR;
      return err;
    }

    if (grfnparam->calib == SRT_YES) {
      err = SrtMCRanSrc_init(&mcrs);
      if (err) {
        FREE_MC_CORE_EXFR;
        return err;
      }
    }
  } else
    mcrs.cxxur_path = mcrs.pl - 1;

  /* Initialise the sample structure */
  sam = srt_calloc(num_time_pts, sizeof(SrtSample));
  for (t = 0, step = top; t < num_time_pts; t++, step = step->next) {
    sam[t].numeraire_index = und_info->numeraire_index;
  }

  /* Need to create tminf structure in SrtStpPtr first */
  err = srtstptminfalloc(top, und_info->no_of_underlyings);
  if (err) {
    FREE_MC_CORE_EXFR;
    return err;
  }

  /* Keep track of the Main Currency (for Quanto adjustments)  (the One of the
   * Numeraire) */
  numeraire_und = lookup_und(und_info->und_data[0].und_name);
  numeraire_ccy = get_underlying_ccy(numeraire_und);

  /* Initialises the steps: stacks all relevant information */
  for (i = 0; i < und_info->no_of_underlyings; i++) {
    /* Gets the Relevant Underlying (in the order of the Grfn Tableau) */
    und = lookup_und(und_info->und_data[i].und_name);

    if (und_info->und_data[i].stochastic == SRT_YES) {
      /* Gets the Model type and dimension */
      err = get_underlying_mdltype(und, &mdl_type);
      if (err) {
        FREE_MC_CORE_EXFR;
        return err;
      }

      err = get_underlying_mdldim(und, &mdl_dim);
      if (err) {
        FREE_MC_CORE_EXFR;
        return err;
      }

      /* Initialisation for Black Sholes like diffusions */
      if ((mdl_type == BLACK_SCHOLES) || (mdl_type == NORMAL_BS) ||
          (mdl_type == FX_STOCH_RATES) || (mdl_type == EQ_STOCH_RATES)) {

        if (mdl_type == FX_STOCH_RATES) {

          dom_und_name = get_domname_from_fxund(und);
          dom_und = lookup_und(dom_und_name);
          if (!dom_und)
            return serror("Could not find %s underlying", dom_und_name);
          err = get_underlying_mdltype(dom_und, &dom_mdl_type);
          if (err) {
            FREE_MC_CORE_EXFR;
            return err;
          }

          for_und_name = get_forname_from_fxund(und);
          for_und = lookup_und(for_und_name);
          if (!for_und)
            return serror("Could not find %s underlying", for_und_name);
          err = get_underlying_mdltype(for_und, &for_mdl_type);
          if (err) {
            FREE_MC_CORE_EXFR;
            return err;
          }

          /* Force the Jumping if possible  and set in srt_f_grfn_und*/
          if ((dom_mdl_type == LGM) && (for_mdl_type == LGM) &&
              (und_info->jumping == SRT_YES)) {
            err = srt_f_fx_initstp(top, und_info, und, dom_und, for_und, i,
                                   numeraire_und);
            if (err) {
              FREE_MC_CORE_EXFR;
              return err;
            }
          } else {
            err = srt_f_loginistp(top, und, i, numeraire_und, und_info);
            if (err) {
              FREE_MC_CORE_EXFR;
              return err;
            }
          } /* END if not jumping with FX STOCH RATES */

        } /* END if(mdl_type == FX_STOCH_RATES) */

        else {
          err = srt_f_loginistp(top, und, i, numeraire_und, und_info);
          if (err) {
            FREE_MC_CORE_EXFR;
            return err;
          }
        }

      } /* END of initilisation for Black Scholes like diffusions */
      else {
        if (err = srt_f_irministp(top, und, i, numeraire_und, und_info)) {
          FREE_MC_CORE_EXFR;
          return err;
        }
      }

      /* This is required for Y_T_at_t_compute to work properly */

      if ((mdl_type == LGM) || (mdl_type == LGM_STOCH_VOL)) {
        for (t = 0, step = top; t < num_time_pts; t++, step = step->next) {
          tminf = (SrtIRMTmInf *)step->tminf[i];
          switch (mdl_dim) {
          case ONE_FAC:
            sam_get(sam[t], i, PHI) = sam_get(tminf->fwd_sam, i, PHI);
            break;
          case TWO_FAC:
            sam_get(sam[t], i, PHI1) = sam_get(tminf->fwd_sam, i, PHI1);
            sam_get(sam[t], i, PHI2) = sam_get(tminf->fwd_sam, i, PHI2);
            sam_get(sam[t], i, CROSSPHI) = sam_get(tminf->fwd_sam, i, CROSSPHI);
            break;
          default:
            sam_get(sam[t], i, PHI) = 0.0;
          } /* END switch declaration  */

        } /* END of for t loop   - for LGM only */

      } /* END if mdl type = LGM */

    } /* END if stochastic == SRT_YES */

    else /* if stochastic == SRT_NO */
    {
      if (err = srt_f_basicinistp(top, und, i)) {
        FREE_MC_CORE_EXFR;
        return err;
      }

      /* we only need to set up sam[t].numeriare once if it is deterministic */

      monte_DETERMINISTIC_evolve(top, sam);

    } /* END if stochastic == SRT_NO */

  } /* END for i loop on underlying number*/

  /* Attaches the right correlation coefficients to the right steps        ,
   * when needed  */
  if ((und_info->no_of_underlyings > 1) &&
      (und_info->two_factor_model == SRT_NO) &&
      (und_info->use_stochastic_vol == SRT_NO)) {
    mcrs.need_to_correl = SRT_YES;

    if (err = srt_f_attach_correl_to_stp(top, und_info->corr_ts)) {
      FREE_MC_CORE_EXFR;
      return err;
    }

  } else {
    mcrs.need_to_correl = SRT_NO;
  }

  dfs = dvector(0, num_time_pts);
  sum_disc_cf = dvector(0, num_time_pts);
  Z_sum = dvector(0, num_time_pts);
  Pt_by_Zt = dvector(0, num_time_pts);
  if ((dfs == NULL) || (sum_disc_cf == NULL) || (Z_sum == NULL) ||
      (Pt_by_Zt == NULL)) {
    FREE_MC_CORE_EXFR;
    return serror("Memory allocation error inMC_Core");
  }

  /* Pricing the security: sum_price=PRICE */

  /* Initialization */

  cash_flow = 0.0;
  memset(dfs, 0, num_time_pts * sizeof(double));
  memset(sum_disc_cf, 0, num_time_pts * sizeof(double));
  memset(Z_sum, 0, num_time_pts * sizeof(double));
  memset(Pt_by_Zt, 0, num_time_pts * sizeof(double));

  /* If we just want the intrinsic value i.e. no brownian is defined        , we
    just need one path to evaluate the answer.  Therefore we modified the no. of
    Monte Carlo paths below if this is the case.  */
  if (und_info->no_of_brownians == 0)
    grfnparam->num_MCarlo_paths = 1;

  /* Compute discount factors from 0 to the step date for discounting using dom
   * und */
  t = 0;
  step = top;
  und = lookup_und(und_info->und_data[0].und_name);
  err = get_underlying_mdltype(und, &mdl_type);
  if (err) {
    FREE_MC_CORE_EXFR;
    return err;
  }

  while (step != NULL) {
    switch (mdl_type) {
      /* Treats FX stoch rates separately */
    case FX_STOCH_RATES:
      if (und_info->jumping == SRT_YES) {
        dfs[t] = ((SrtFXTmInf *)step->tminf[0])->df;
        break;
      }
    case EQ_STOCH_RATES:
    case BLACK_SCHOLES:
    case NORMAL_BS:
      dfs[t] = ((SrtLogTmInf *)step->tminf[0])->df;
      break;
    case NONE:
      dfs[t] = ((SrtBasicTmInf *)step->tminf[0])->df;
      break;
    default:
      dfs[t] = ((SrtIRMTmInf *)step->tminf[0])->df;
      break;
    }
    step = step->next;
    t = t + 1;
  }

  /* Memory allocation for a vector where all the Columns PV's will be stored */
  col_pvs = dvector(0, grfnparam->node_dim - 1);
  deal_cash_flows_pv = dmatrix(0, num_time_pts - 1, 0, grfnparam->node_dim - 1);
  path_cash_flows = dmatrix(0, num_time_pts - 1, 0, grfnparam->node_dim - 1);
  memset(&deal_cash_flows_pv[0][0], 0,
         (int)(num_time_pts * grfnparam->node_dim * sizeof(double)));
  memset(&path_cash_flows[0][0], 0,
         (int)(num_time_pts * grfnparam->node_dim * sizeof(double)));

  if (grfnparam->exfrontierlast_time == SRT_YES) {
    nstart = 0;
    nend = grfnparam->num_MCarlo_paths;
  }

  else {

    nstart = (grfnparam->exfrontiercounter) * 200;
    nend = (grfnparam->exfrontiercounter) * 200 + 200;
  }

  /* Starts loop on number of paths */
  for (n = nstart; n < nend; n++) {

    /* Generates or gets the next MC sample path and puts it in mcrs.r[cur_path]
     */
    /*if (grfnparam->calib == SRT_NO)
    {
            if ( err = SrtMCRanSrc_next(&mcrs) )
            {
                    FREE_MC_CORE_EXFR;
                    return err ;
            }
    }
    else */
    mcrs.cxxur_path = n;

    /* Correlates (if necessary) the n(numb of undelyings)
           random (independent) Brownian paths mcrs.r[cur_path][][]        ,
           for all the steps with time dependent correlation matrixes */
    if (grfnparam->calib == SRT_NO) {
      if (err = SrtMCRanSrc_correl(&mcrs, top)) {
        FREE_MC_CORE_EXFR;
        return err;
      }
    }

    /* Evolves each underlying        , one by one
      We have to evolve first the equity/interest rate underlyings
     then the forex ones (this is due to the stochatic rates
    needed for quanto adjustment. */
    for (j = 0; j < 2; j++) {
      for (i = 0; i < und_info->no_of_underlyings; i++) {
        und = lookup_und(und_info->und_data[i].und_name);
        err = get_underlying_mdltype(und, &mdl_type);
        if (err) {
          FREE_MC_CORE_EXFR;
          return err;
        }

        if (((j == 0) && (!(mdl_type == EQ_STOCH_RATES)) &&
             (!ISUNDTYPE(und, FOREX_UND))) ||
            ((j == 1) && (ISUNDTYPE(und, FOREX_UND))) ||
            ((j == 1) && (mdl_type == EQ_STOCH_RATES))) {

          /* if underlying is stochastic        , set up sam[t]; otherwise
          sam[t]. pv_money_mkt should already be set up above because it is
          determinstic */

          if (und_info->und_data[i].stochastic == SRT_YES) {
            err = get_underlying_mdltype(und, &mdl_type);
            if (err) {
              FREE_MC_CORE_EXFR;
              return err;
            }

            err = get_underlying_mdldim(und, &mdl_dim);
            if (err) {
              FREE_MC_CORE_EXFR;
              return err;
            }

            switch (mdl_type) {
            case NORMAL_BS:
              monte_NORMALBS_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam, i);
              break;
            case BLACK_SCHOLES:
              monte_BLACKSCHOLES_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                        i);
              break;
            case EQ_STOCH_RATES:
              err = get_index_from_und_info(
                  und_info, get_discname_from_underlying(und), &dom_ind);
              if (err) {
                FREE_MC_CORE_EXFR;
                return err;
              }
              monte_EQ_STOCH_RATES_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                          i, dom_ind);
              break;
            case FX_STOCH_RATES: /* DONE ONLY FOR FOREX UNDERLYING */
              err = get_index_from_und_info(
                  und_info, get_domname_from_fxund(und), &dom_ind);
              if (err) {
                FREE_MC_CORE_EXFR;
                return err;
              }
              err = get_index_from_und_info(
                  und_info, get_forname_from_fxund(und), &for_ind);
              if (err) {
                FREE_MC_CORE_EXFR;
                return err;
              }

              if ((mdl_type == FX_STOCH_RATES) && (dom_mdl_type == LGM) &&
                  (for_mdl_type == LGM) && (und_info->jumping == SRT_YES)) {

                monte_FX_STOCH_RATES_Jumping_evolve(
                    mcrs.r[mcrs.cxxur_path][i], top, sam, i, dom_ind, for_ind);
              } else

                monte_FX_STOCH_RATES_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i, dom_ind, for_ind);

              break;
            case ETABETA:
              if (mdl_dim == TWO_FAC) {
                FREE_MC_CORE_EXFR;
                return serror("Two Factor ETABETA not implemented yet...");
              } else {
                monte_BETAETA_1f_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                        i);
              }
              break;
            case LGM:
              if (mdl_dim == ONE_FAC) {
                if (und_info->jumping == SRT_YES)
                  monte_LGM_1f_Jumping_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                              sam, i);
                else if (und_info->jumping == SRT_NO)
                  monte_LGM_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i);
              } else if (mdl_dim == TWO_FAC) {
                monte_LGM_2f_evolve(mcrs.r[mcrs.cxxur_path], top, sam, i);
              }
              break;
            case NEWLGM:
              monte_NEWLGM_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top, sam,
                                           i);
              break;
            case NEWCHEYBETA:
              monte_NEWCHEYBETA_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                                sam, i);
              break;
            case LGM_STOCH_VOL:
              monte_LGM_1f_STOCHVOL_evolve(mcrs.r[mcrs.cxxur_path], top, sam,
                                           Pt_by_Zt, i);
              break;
            case CHEY_STOCH_VOL:
              monte_CHE_1f_STOCHVOL_evolve(mcrs.r[mcrs.cxxur_path], top, sam,
                                           Pt_by_Zt, i);
              break;
            case CHEY_BETA_STOCH_VOL:
              monte_CHEYBETA_1f_STOCHVOL_evolve(mcrs.r[mcrs.cxxur_path], top,
                                                sam, Pt_by_Zt, i);
              break;
            case CHEY_BETA:
              if (mdl_dim == ONE_FAC) {
                monte_CHEYBETA_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                               sam, i);
              } else if (mdl_dim == TWO_FAC) {
                monte_CHEYBETA_2f_Euler_evolve(mcrs.r[mcrs.cxxur_path], top,
                                               sam, i);
              }
              break;
            case MIXED_BETA:
              if (mdl_dim == ONE_FAC) {
                FREE_MC_CORE_EXFR;
                return serror("One Factor MIXEDBETA is nonsense");
              } else if (mdl_dim == TWO_FAC) {
                monte_MIXEDBETA_2f_Euler_evolve(mcrs.r[mcrs.cxxur_path], top,
                                                sam, i);
              }
              break;
            case CHEY:
            default:
              if (mdl_dim == ONE_FAC) {
                if (grfnparam->difference_scheme == MILSHTEIN) {
                  monte_CHE_1f_Milshtein_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                                sam, i);
                } else /* EULER */
                {
                  monte_CHE_1f_Euler_evolve(mcrs.r[mcrs.cxxur_path][i], top,
                                            sam, i);
                }
              } else if (mdl_dim == TWO_FAC) {
                monte_CHE_2f_Euler_evolve(mcrs.r[mcrs.cxxur_path], top, sam, i);
              }
              break;

            } /* END switch declaration */

          } /* END if stochastic == SRT_YES */

        } /* FOREX or NOT FOREX underlying */

      } /* END for i loop on underlying number*/

    } /* END of the loop on underlying types*/

    /* Loops on all the time steps of this path to compute the cash flow */
    temp = 0.0;
    for (t = 0, step = top; t < num_time_pts; step = step->next, t++) {

      sam[t].pathnum = n;

      /* Call to the main Grfn COMLL evaluation function */

      err = evalcf(
          (GrfnEvent *)step->e, &sam[t], gd, (double *)(path_cash_flows[t]),
          (EvalEventDfsFct)srt_f_calc_grfn_event_dfs, und_info, &cash_flow);
      if (err) {
        FREE_MC_CORE_EXFR;
        return err;
      }

      /* Here we multiply by dfs[t] that has been removed from the numeraire ,
         (it has been computed with exp(r - f(0        ,t))  or exp(A(t)+h(t
         ,x))  */
      cash_flow_pv = cash_flow / (sam[t].numeraire / dfs[t]);
      sum_disc_cf[t] += cash_flow_pv;
      temp += cash_flow_pv;

      /* Store the information for all the columns        , in a global deal pv
       * calculation */
      for (i = 0; i < grfnparam->node_dim; i++) {
        deal_cash_flows_pv[t][i] +=
            path_cash_flows[t][i] / (sam[t].numeraire / dfs[t]);
        /* Resets the path_cash_flow to 0.0 */
        path_cash_flows[t][i] = 0.0;
      }

      /* The zero-coupon price computed for this time step (sum on all paths)
       * (without its forward) */
      Z_sum[t] += 1.0 / sam[t].numeraire;

    } /* END for t loop on all time steps of current path */

    /* Increments the current sums of all cashflows and cash flow squares for
     * all paths */
    exp_disc_pay += temp;
    exp_disc_pay_sqr += temp * temp;

  } /* END for n loop on path number*/

  /* Frees the memory allocated for random numbers generation */
  if (grfnparam->exfrontierlast_time == SRT_YES)
    SrtMCRanSrc_end(&mcrs);

  /* Computes the renormalisation factor for the price */
  if (grfnparam->mc_renormalize) {
    /* The correction factor is (at each step) the zero-coupon bond price
     * computed by full MC simulation  */
    for (t = 0; t < num_time_pts; t++)
      Pt_by_Zt[t] = 1.0 / Z_sum[t];

  } else {
    /* The corrector is just the number paths */
    for (t = 0; t < num_time_pts; t++)
      Pt_by_Zt[t] = 1.0 / (double)grfnparam->num_MCarlo_paths;
  }

  /* The Real price is the price obtained by Monte Carlo        , renormalised
   * step by step */
  sum_price = 0.0;
  memset(&col_pvs[0], 0, sizeof(double) * grfnparam->node_dim);
  for (t = 0; t < num_time_pts; t++) {
    sum_price += (Pt_by_Zt[t] * sum_disc_cf[t]);
    for (i = 0; i < grfnparam->node_dim; i++)
      col_pvs[i] += deal_cash_flows_pv[t][i] * Pt_by_Zt[t];
  }

  /* Stores the premium in the Input/Output list */
  err =
      srt_f_IOstructsetpremium((SrtIOStruct *)iorequest, SRT_NO, sum_price, "");
  if (err) {
    FREE_MC_CORE_EXFR;
    return (err);
  }

  /* Computes and stores the standard deviation as well */
  exp_disc_pay_sqr *=
      1.0 + DBL_EPSILON *
                grfnparam->num_MCarlo_paths; /* To prevent rounding errors */
  stdev = (sqrt(exp_disc_pay_sqr / grfnparam->num_MCarlo_paths -
                exp_disc_pay * exp_disc_pay / grfnparam->num_MCarlo_paths /
                    grfnparam->num_MCarlo_paths)) /
          sqrt(grfnparam->num_MCarlo_paths);
  err = srt_f_IOstructsetstdev((SrtIOStruct *)iorequest, SRT_NO, stdev, "");
  if (err) {
    FREE_MC_CORE_EXFR;
    return (err);
  }

  /* Stores all the columns PV in the I/O list */
  err = srt_f_IOstructsetcolpvs((SrtIOStruct *)iorequest, SRT_NO, col_pvs,
                                grfnparam->node_dim, "");
  if (err) {
    FREE_MC_CORE_EXFR;
    return err;
  }

  /* Free all allocated memory */
  FREE_MC_CORE_EXFR;

  /* Return a success message */
  return NULL;

} /* END McCoreExFrontier */
