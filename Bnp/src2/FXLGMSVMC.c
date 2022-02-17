/* ==========================================================================
   FILE_NAME:	FXLGMSVMC.c

   PURPOSE:		Monte Carlo FX Black-Scholes / IR Dom LGMSV2F / IR for
   LGMSV2F

   DATE:		01/10/03
   ========================================================================== */

#include "MCEBOptimisation.h"
#include "RandomGen.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

Err fxlgmsv_mc_balsam(/*	Time Information  */
                      int iNbTime, int iNbEvent, double *dTime, double *dDate,

                      int iNumPaths,

                      // Domestic
                      /*	Model data Information	*/
                      double ddomLambdaX1, double ddomLambdaX2,

                      double *ddomSigma, double *ddomAlphaLGM,
                      double *ddomRhoLGM,

                      double *ddomAlpha, double *ddomLambdaEps,
                      double *ddomLvlEps, double *ddomRho, double *ddomRho2,

                      double *domzcvol1_star, double *domzcvol2_star,

                      /* Parameters for DF(t  ,T*) reconstruction */
                      double *ddomff_star, double *domgam1_star,
                      double *domgam2_star, double *domgam1_2_star,
                      double *domgam2_2_star, double *domgam12_star,

                      // Foreign
                      /*	Model data Information	*/
                      double dforLambdaX1, double dforLambdaX2,

                      double *dforSigma, double *dforAlphaLGM,
                      double *dforRhoLGM,

                      double *dforAlpha, double *dforLambdaEps,
                      double *dforLvlEps, double *dforRho, double *dforRho2,

                      double *forzcvol1_star, double *forzcvol2_star,

                      /* Parameters for DF(t  ,T*) reconstruction */
                      double *dforff_star, double *forgam1_star,
                      double *forgam2_star, double *forgam1_2_star,
                      double *forgam2_2_star, double *forgam12_star,

                      // FX
                      double fwd_fx_TStar, double *fx_vol,

                      //	Correlation
                      double ***CorrMatrix,
                      /*	0 : Dom1
                              1 : Dom2
                              2 : DomSV
                              3 : For1
                              4 : For2
                              5 : ForSV
                              6 : FX	*/

                      /*	Product data */
                      void **func_parm_tab, int *EvalEvent,

                      /* for Optimisation of exercise boundary */
                      int do_optimisation, int *optimise, MCEBPARAMS params,

                      /*	Initialisation function to be called at the
                         beggining of each path or NULL if none */
                      void (*init_func)(),

                      /*	Payoff function */
                      Err (*payoff_func)(long path_index, double evt_date,
                                         double evt_time, void *func_parm,

                                         // Domestic
                                         /* Model data	*/
                                         double domft1, double domft2,
                                         double domphi1, double domphi2,
                                         double domphi12, double domv,

                                         // Foreign
                                         /* Model data	*/
                                         double forft1, double forft2,
                                         double forphi1, double forphi2,
                                         double forphi12, double forv,

                                         // FX
                                         double fx_spot,

                                         /* Vector of results to be updated */
                                         int nprod,
                                         /* Result	*/
                                         double *prod_val, int *stop_path),
                      /*	Result */
                      int iNbProduct, double **res) {
  Err err = NULL;
  clock_t time1, time2;
  long i, j, k, l, m, n;

  double *sum_payoff = NULL, *sum_payoff2 = NULL, *path_payoff = NULL,
         *res_evt = NULL, **matrix = NULL,

         *dom_sigf1 = NULL, *dom_sigf2 = NULL, *dom_sigpsi1 = NULL,
         *dom_sigpsi2 = NULL, *dom_sigpsi12 = NULL, *dom_muv1 = NULL,
         *dom_muv2 = NULL, *dom_sigv = NULL,

         *for_sigf1 = NULL, *for_sigf2 = NULL, *for_sigpsi1 = NULL,
         *for_sigpsi2 = NULL, *for_sigpsi12 = NULL, *for_muv1 = NULL,
         *for_muv2 = NULL, *for_sigv = NULL,

         *fx_sig = NULL,

         **CorrMatrix_loc = NULL,

         ***cholmat = NULL,

         ***save_values = NULL,

         *gauss = NULL;

  int stop_path;
  int evtindex;
  int NbFactor, NbFactorReduction;

  double dt, sqdt;
  double *dom_ft1 = NULL, *dom_ft2 = NULL, *dom_v = NULL, dom_sqv,
         *dom_psi1 = NULL, *dom_psi2 = NULL, *dom_psi12 = NULL;
  double *for_ft1 = NULL, *for_ft2 = NULL, *for_v = NULL, for_sqv,
         *for_psi1 = NULL, *for_psi2 = NULL, *for_psi12 = NULL;
  double *fwdfx = NULL;
  double *brow = NULL;
  double dom_df, for_df;
  double domzcvol1, domzcvol2;
  double forzcvol1, forzcvol2;

  long seed = -123456789;
  int NbInfoSave = 13;

  SRandomGen rg;

  /* --------------------------------------------------------------------------------------------------
                                                                                  Transformation into 5 Factor
     --------------------------------------------------------------------------------------------------
   */
  if (ddomAlphaLGM[1] < 1.0e-6 && dforAlphaLGM[1] < 1.0e-6) {
    NbFactor = 5;
    NbFactorReduction = 1;
  } else {
    NbFactor = 7;
    NbFactorReduction = 0;
  }

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Initialisation
     --------------------------------------------------------------------------------------------------
   */
  memset(&rg, 0, sizeof(SRandomGen));

  /* For computational time calculation				 */
  time1 = clock();

  /* odd number of paths */
  iNumPaths = (int)(iNumPaths / 2) * 2 + 1;

  brow = dvector(0, 6);
  // matrix = dmatrix (0  , iNumPaths - 1  , 0  , NbFactor * (iNbTime - 1) - 1);
  gauss = dvector(0, NbFactor);
  matrix = dmatrix(0, iNumPaths - 1, 0, NbInfoSave * iNbEvent - 1);

  res_evt = dvector(0, iNbProduct - 1);
  path_payoff = dvector(0, iNbProduct - 1);
  sum_payoff = dvector(0, iNbProduct - 1);
  sum_payoff2 = dvector(0, iNbProduct - 1);

  /* for precalculations */
  dom_sigf1 = dvector(0, iNbTime - 1);
  dom_sigf2 = dvector(0, iNbTime - 1);
  dom_sigpsi1 = dvector(0, iNbTime - 1);
  dom_sigpsi2 = dvector(0, iNbTime - 1);
  dom_sigpsi12 = dvector(0, iNbTime - 1);
  dom_muv1 = dvector(0, iNbTime - 1);
  dom_muv2 = dvector(0, iNbTime - 1);
  dom_sigv = dvector(0, iNbTime - 1);

  for_sigf1 = dvector(0, iNbTime - 1);
  for_sigf2 = dvector(0, iNbTime - 1);
  for_sigpsi1 = dvector(0, iNbTime - 1);
  for_sigpsi2 = dvector(0, iNbTime - 1);
  for_sigpsi12 = dvector(0, iNbTime - 1);
  for_muv1 = dvector(0, iNbTime - 1);
  for_muv2 = dvector(0, iNbTime - 1);
  for_sigv = dvector(0, iNbTime - 1);

  dom_ft1 = dvector(0, iNumPaths - 1);
  dom_ft2 = dvector(0, iNumPaths - 1);
  dom_v = dvector(0, iNumPaths - 1);
  dom_psi1 = dvector(0, iNumPaths - 1);
  dom_psi2 = dvector(0, iNumPaths - 1);
  dom_psi12 = dvector(0, iNumPaths - 1);

  for_ft1 = dvector(0, iNumPaths - 1);
  for_ft2 = dvector(0, iNumPaths - 1);
  for_v = dvector(0, iNumPaths - 1);
  for_psi1 = dvector(0, iNumPaths - 1);
  for_psi2 = dvector(0, iNumPaths - 1);
  for_psi12 = dvector(0, iNumPaths - 1);

  fwdfx = dvector(0, iNumPaths - 1);

  fx_sig = dvector(0, iNbTime - 1);

  CorrMatrix_loc = dmatrix(0, 6, 0, 6);

  cholmat = f3tensor(0, iNbTime - 1, 0, 6, 0, 6);

  if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt ||
      !dom_sigf1 || !dom_sigf2 || !dom_sigpsi1 || !dom_sigpsi2 ||
      !dom_sigpsi12 || !dom_muv1 || !dom_muv2 || !dom_sigv || !for_sigf1 ||
      !for_sigf2 || !for_sigpsi1 || !for_sigpsi2 || !for_sigpsi12 ||
      !for_muv1 || !for_muv2 || !for_sigv || !fx_sig || !dom_ft1 || !dom_ft2 ||
      !dom_v || !dom_psi1 || !dom_psi2 || !dom_psi12 || !fwdfx || !for_ft1 ||
      !for_ft2 || !for_v || !for_psi1 || !for_psi2 || !for_psi12 || !cholmat) {
    err = "Memory allocation failure in fxlgmsv_mc_balsam";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(iNumPaths, iNbEvent, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  memset(sum_payoff, 0, iNbProduct * sizeof(double));
  memset(sum_payoff2, 0, iNbProduct * sizeof(double));

  /* All the needed precalculations */
  for (j = 1; j < iNbTime; j++) {
    dt = dTime[j] - dTime[j - 1];
    sqdt = sqrt(dt);

    dom_sigf1[j] = ddomSigma[j] * sqdt;
    dom_sigf2[j] = dom_sigf1[j] * ddomAlphaLGM[j];
    dom_sigpsi1[j] = dom_sigf1[j] * dom_sigf1[j];
    dom_sigpsi2[j] = dom_sigf2[j] * dom_sigf2[j];
    dom_sigpsi12[j] = ddomRhoLGM[j] * dom_sigf1[j] * dom_sigf2[j];

    dom_sigv[j] = ddomAlpha[j] * sqdt;
    dom_muv1[j] = ddomLambdaEps[j] * dt;
    dom_muv2[j] = ddomLvlEps[j] * dt;

    for_sigf1[j] = dforSigma[j] * sqdt;
    for_sigf2[j] = for_sigf1[j] * dforAlphaLGM[j];
    for_sigpsi1[j] = for_sigf1[j] * for_sigf1[j];
    for_sigpsi2[j] = for_sigf2[j] * for_sigf2[j];
    for_sigpsi12[j] = dforRhoLGM[j] * for_sigf1[j] * for_sigf2[j];

    for_sigv[j] = dforAlpha[j] * sqdt;
    for_muv1[j] = dforLambdaEps[j] * dt;
    for_muv2[j] = dforLvlEps[j] * dt;

    fx_sig[j] = fx_vol[j] * sqdt;

    for (i = 0; i < 7; ++i) {
      for (k = 0; k < 7; ++k) {
        CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
      }
    }

    err = choldc(7, CorrMatrix_loc, cholmat[j]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -BalSam generation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* fill the Brownian matrix */
  /*
  err = balsam_generation (iNumPaths  , NbFactor * (iNbTime - 1)  , matrix);
  if (err)
  {
          goto FREE_RETURN;
  }
  */

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Convolution
     --------------------------------------------------------------------------------------------------
   */
  err = ABS_Init(&rg, seed, iNumPaths, NbFactor, 0);
  if (err)
    goto FREE_RETURN;

  /* initialisation */
  evtindex = 0;

  for (i = 0; i < iNumPaths; i++) {
    dom_v[i] = 1.0;
    for_v[i] = 1.0;
    fwdfx[i] = log(fwd_fx_TStar);
  }

  if (func_parm_tab[0] && EvalEvent[0]) {
    for (i = 0; i < iNumPaths; i++) {
      matrix[i][evtindex * NbInfoSave + 0] = dom_ft1[i];
      matrix[i][evtindex * NbInfoSave + 1] = dom_ft2[i];
      matrix[i][evtindex * NbInfoSave + 2] = dom_psi1[i];
      matrix[i][evtindex * NbInfoSave + 3] = dom_psi2[i];
      matrix[i][evtindex * NbInfoSave + 4] = dom_psi12[i];
      matrix[i][evtindex * NbInfoSave + 5] = dom_v[i];
      matrix[i][evtindex * NbInfoSave + 6] = for_ft1[i];
      matrix[i][evtindex * NbInfoSave + 7] = for_ft2[i];
      matrix[i][evtindex * NbInfoSave + 8] = for_psi1[i];
      matrix[i][evtindex * NbInfoSave + 9] = for_psi2[i];
      matrix[i][evtindex * NbInfoSave + 10] = for_psi12[i];
      matrix[i][evtindex * NbInfoSave + 11] = for_v[i];
      matrix[i][evtindex * NbInfoSave + 12] = fwdfx[i];
    }

    evtindex++;
  }

  // Generate all paths simultaneously:
  for (j = 1; j < iNbTime; j++) {
    for (i = 0; i < iNumPaths; i++) {
      for (k = 0; k < NbFactor; k++) {
        err = rg.Gauss(&rg, &gauss[k]);
        if (err)
          goto FREE_RETURN;
      }

      for (k = 0; k < 7; k++) {
        brow[k] = 0;
        if (!NbFactorReduction) {
          for (l = 0; l < 7; l++) {
            brow[k] += cholmat[j][k][l] * gauss[l];
          }
        } else {
          n = 0;
          for (l = 0; l < 7; l++) {
            if (l != 1 && l != 4) {
              brow[k] += cholmat[j][k][l] * gauss[n];
              n++;
            }
          }
        }
      }

      dom_sqv = sqrt(dom_v[i]);
      for_sqv = sqrt(for_v[i]);

      domzcvol1 = domzcvol1_star[j] * dom_sigf1[j] * dom_sqv;
      domzcvol2 = domzcvol2_star[j] * dom_sigf2[j] * dom_sqv;
      forzcvol1 = forzcvol1_star[j] * for_sigf1[j] * for_sqv;
      forzcvol2 = forzcvol2_star[j] * for_sigf2[j] * for_sqv;

      dom_ft1[i] += dom_sigf1[j] * dom_sqv * brow[0];
      dom_ft2[i] += dom_sigf2[j] * dom_sqv * brow[1];
      dom_psi1[i] += dom_sigpsi1[j] * dom_v[i];
      dom_psi2[i] += dom_sigpsi2[j] * dom_v[i];
      dom_psi12[i] += dom_sigpsi12[j] * dom_v[i];
      dom_v[i] = max(dom_v[i] + dom_muv2[j] - dom_muv1[j] * dom_v[i] +
                         dom_sigv[j] * dom_sqv * brow[2],
                     0.0);

      for_ft1[i] += for_sigf1[j] * for_sqv * brow[3];
      for_ft1[i] += -for_sigf1[j] * for_sqv * CorrMatrix[j][3][6] * fx_sig[j];
      for_ft1[i] += -for_sigf1[j] * for_sqv *
                    (forzcvol1 + CorrMatrix[j][3][4] * forzcvol2);
      for_ft1[i] +=
          +for_sigf1[j] * for_sqv *
          (CorrMatrix[j][3][0] * domzcvol1 + CorrMatrix[j][3][1] * domzcvol2);

      for_ft2[i] += for_sigf2[j] * for_sqv * brow[4];
      for_ft2[i] += -for_sigf2[j] * for_sqv * CorrMatrix[j][4][6] * fx_sig[j];
      for_ft2[i] += -for_sigf2[j] * for_sqv *
                    (CorrMatrix[j][4][3] * forzcvol1 + forzcvol2);
      for_ft2[i] +=
          +for_sigf2[j] * for_sqv *
          (CorrMatrix[j][4][0] * domzcvol1 + CorrMatrix[j][4][1] * domzcvol2);

      for_psi1[i] += for_sigpsi1[j] * for_v[i];
      for_psi2[i] += for_sigpsi2[j] * for_v[i];
      for_psi12[i] += for_sigpsi12[j] * for_v[i];
      for_v[i] += for_muv2[j] - for_muv1[j] * for_v[i] +
                  for_sigv[j] * for_sqv * brow[5];
      for_v[i] += -for_sigv[j] * for_sqv * CorrMatrix[j][5][6] * fx_sig[j];
      for_v[i] +=
          -for_sigv[j] * for_sqv *
          (CorrMatrix[j][5][3] * forzcvol1 + CorrMatrix[j][5][4] * forzcvol2);
      for_v[i] +=
          +for_sigv[j] * for_sqv *
          (CorrMatrix[j][5][0] * domzcvol1 + CorrMatrix[j][5][1] * domzcvol2);
      for_v[i] = max(0.0, for_v[i]);

      fwdfx[i] +=
          (fx_sig[j] * brow[6] - domzcvol1 * brow[0] - domzcvol2 * brow[1] +
           forzcvol1 * brow[3] + forzcvol2 * brow[4]);

      fwdfx[i] += -0.5 * (fx_sig[j] * fx_sig[j] + domzcvol1 * domzcvol1 +
                          domzcvol2 * domzcvol2 + forzcvol1 * forzcvol1 +
                          forzcvol2 * forzcvol2);

      fwdfx[i] += -(-fx_sig[j] * domzcvol1 * CorrMatrix[j][6][0] -
                    fx_sig[j] * domzcvol2 * CorrMatrix[j][6][1] +
                    fx_sig[j] * forzcvol1 * CorrMatrix[j][6][3] +
                    fx_sig[j] * forzcvol2 * CorrMatrix[j][6][4]

                    + domzcvol1 * domzcvol2 * CorrMatrix[j][0][1] -
                    domzcvol1 * forzcvol1 * CorrMatrix[j][0][3] -
                    domzcvol1 * forzcvol2 * CorrMatrix[j][0][4]

                    - domzcvol2 * forzcvol1 * CorrMatrix[j][1][3] -
                    domzcvol2 * forzcvol2 * CorrMatrix[j][1][4]

                    + forzcvol1 * forzcvol2 * CorrMatrix[j][3][4]);

      /*	Case of evaluation events */
      if (EvalEvent[j]) {

        matrix[i][evtindex * NbInfoSave + 0] = dom_ft1[i];
        matrix[i][evtindex * NbInfoSave + 1] = dom_ft2[i];
        matrix[i][evtindex * NbInfoSave + 2] = dom_psi1[i];
        matrix[i][evtindex * NbInfoSave + 3] = dom_psi2[i];
        matrix[i][evtindex * NbInfoSave + 4] = dom_psi12[i];
        matrix[i][evtindex * NbInfoSave + 5] = dom_v[i];
        matrix[i][evtindex * NbInfoSave + 6] = for_ft1[i];
        matrix[i][evtindex * NbInfoSave + 7] = for_ft2[i];
        matrix[i][evtindex * NbInfoSave + 8] = for_psi1[i];
        matrix[i][evtindex * NbInfoSave + 9] = for_psi2[i];
        matrix[i][evtindex * NbInfoSave + 10] = for_psi12[i];
        matrix[i][evtindex * NbInfoSave + 11] = for_v[i];
        matrix[i][evtindex * NbInfoSave + 12] = fwdfx[i];

        if (i == iNumPaths - 1) {
          evtindex++;
        }
      }
    }
  }

  for (i = 0; i < iNumPaths; i++) {
    /* initialisation */
    evtindex = 0;
    stop_path = 0;

    memset(path_payoff, 0, iNbProduct * sizeof(double));

    if (func_parm_tab[0] && EvalEvent[0]) {
      dom_df = exp(ddomff_star[0]);
      for_df = exp(dforff_star[0]);

      /* Event at time 0.0 */
      err = payoff_func(0, dDate[0], dTime[0], func_parm_tab[0],

                        matrix[i][evtindex * NbInfoSave + 0],
                        matrix[i][evtindex * NbInfoSave + 1],
                        matrix[i][evtindex * NbInfoSave + 2],
                        matrix[i][evtindex * NbInfoSave + 3],
                        matrix[i][evtindex * NbInfoSave + 4],
                        matrix[i][evtindex * NbInfoSave + 5],

                        matrix[i][evtindex * NbInfoSave + 6],
                        matrix[i][evtindex * NbInfoSave + 7],
                        matrix[i][evtindex * NbInfoSave + 8],
                        matrix[i][evtindex * NbInfoSave + 9],
                        matrix[i][evtindex * NbInfoSave + 10],
                        matrix[i][evtindex * NbInfoSave + 11],

                        exp(matrix[i][evtindex * NbInfoSave + 12]) * dom_df /
                            for_df,

                        iNbProduct, res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < iNbProduct; k++) {
        path_payoff[k] += res_evt[k] / dom_df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[0], res_evt, i, dom_df,
                                       params);
      }

      evtindex++;
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < iNbTime; j++) {
      /*	Case of evaluation events */
      if (EvalEvent[j]) {
        dom_df = exp(
            ddomff_star[evtindex] +
            domgam1_star[evtindex] * matrix[i][evtindex * NbInfoSave + 0] +
            domgam2_star[evtindex] * matrix[i][evtindex * NbInfoSave + 1] +
            domgam1_2_star[evtindex] * matrix[i][evtindex * NbInfoSave + 2] +
            domgam2_2_star[evtindex] * matrix[i][evtindex * NbInfoSave + 3] +
            domgam12_star[evtindex] * matrix[i][evtindex * NbInfoSave + 4]);

        for_df = exp(
            dforff_star[evtindex] +
            forgam1_star[evtindex] * matrix[i][evtindex * NbInfoSave + 6] +
            forgam2_star[evtindex] * matrix[i][evtindex * NbInfoSave + 7] +
            forgam1_2_star[evtindex] * matrix[i][evtindex * NbInfoSave + 8] +
            forgam2_2_star[evtindex] * matrix[i][evtindex * NbInfoSave + 9] +
            forgam12_star[evtindex] * matrix[i][evtindex * NbInfoSave + 10]);

        /* Modification of the Payoff at t */
        err = payoff_func(i, dDate[j], dTime[j], func_parm_tab[j],

                          matrix[i][evtindex * NbInfoSave + 0],
                          matrix[i][evtindex * NbInfoSave + 1],
                          matrix[i][evtindex * NbInfoSave + 2],
                          matrix[i][evtindex * NbInfoSave + 3],
                          matrix[i][evtindex * NbInfoSave + 4],
                          matrix[i][evtindex * NbInfoSave + 5],

                          matrix[i][evtindex * NbInfoSave + 6],
                          matrix[i][evtindex * NbInfoSave + 7],
                          matrix[i][evtindex * NbInfoSave + 8],
                          matrix[i][evtindex * NbInfoSave + 9],
                          matrix[i][evtindex * NbInfoSave + 10],
                          matrix[i][evtindex * NbInfoSave + 11],

                          exp(matrix[i][evtindex * NbInfoSave + 12]) * dom_df /
                              for_df,

                          iNbProduct, res_evt, &stop_path);

        if (err)
          goto FREE_RETURN;

        for (k = 0; k < iNbProduct; k++) {
          path_payoff[k] += res_evt[k] / dom_df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                         dom_df, params);
        }

        evtindex++;
      }
    }

    for (k = 0; k < iNbProduct; k++) {
      sum_payoff[k] += path_payoff[k] / iNumPaths;
      sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
    }

    if (do_optimisation && params->iKnockInCol) {
      /* we recopy in the col pay the pv of the column */
      for (j = 0; j < iNbEvent; j++) {
        if (optimise[j]) {
          save_values[j][params->iNbIndex][i] =
              path_payoff[(int)(save_values[j][params->iNbIndex][i] + 0.5)];
        }
      }
    }
  }

  for (k = 0; k < iNbProduct; k++) {
    res[k][0] = sum_payoff[k];
    res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

    if (res[k][1] > 0.0) {
      res[k][1] = sqrt(res[k][1]);
    } else {
      res[k][1] = 0.0;
    }
  }

  /* Convolution time display */
  time1 = clock();
  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(time1 - time2) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    /* Free the big matrix of memory first */
    if (matrix)
      free_dmatrix(matrix, 0, iNumPaths - 1, 0, NbInfoSave * iNbEvent - 1);
    matrix = NULL;

    time1 = clock();

    err = find_and_optimise_boundary(save_values, iNbEvent, iNumPaths, optimise,
                                     params, &(res[iNbProduct][0]),
                                     &(res[iNbProduct][1]));

    if (err)
      goto FREE_RETURN;

    time2 = clock();
    smessage("Phase 3 -optimisation  , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (matrix)
    free_dmatrix(matrix, 0, iNumPaths - 1, 0, NbInfoSave * iNbEvent - 1);
  if (res_evt)
    free_dvector(res_evt, 0, iNbProduct - 1);
  if (path_payoff)
    free_dvector(path_payoff, 0, iNbProduct - 1);
  if (sum_payoff)
    free_dvector(sum_payoff, 0, iNbProduct - 1);
  if (sum_payoff2)
    free_dvector(sum_payoff2, 0, iNbProduct - 1);

  if (dom_sigf1)
    free_dvector(dom_sigf1, 0, iNbTime - 1);
  if (dom_sigf2)
    free_dvector(dom_sigf2, 0, iNbTime - 1);
  if (dom_sigpsi1)
    free_dvector(dom_sigpsi1, 0, iNbTime - 1);
  if (dom_sigpsi2)
    free_dvector(dom_sigpsi2, 0, iNbTime - 1);
  if (dom_sigpsi12)
    free_dvector(dom_sigpsi12, 0, iNbTime - 1);
  if (dom_muv1)
    free_dvector(dom_muv1, 0, iNbTime - 1);
  if (dom_muv2)
    free_dvector(dom_muv2, 0, iNbTime - 1);
  if (dom_sigv)
    free_dvector(dom_sigv, 0, iNbTime - 1);

  if (for_sigf1)
    free_dvector(for_sigf1, 0, iNbTime - 1);
  if (for_sigf2)
    free_dvector(for_sigf2, 0, iNbTime - 1);
  if (for_sigpsi1)
    free_dvector(for_sigpsi1, 0, iNbTime - 1);
  if (for_sigpsi2)
    free_dvector(for_sigpsi2, 0, iNbTime - 1);
  if (for_sigpsi12)
    free_dvector(for_sigpsi12, 0, iNbTime - 1);
  if (for_muv1)
    free_dvector(for_muv1, 0, iNbTime - 1);
  if (for_muv2)
    free_dvector(for_muv2, 0, iNbTime - 1);
  if (for_sigv)
    free_dvector(for_sigv, 0, iNbTime - 1);

  if (fx_sig)
    free_dvector(fx_sig, 0, iNbTime - 1);

  if (brow)
    free_dvector(brow, 0, 6);

  if (cholmat)
    free_f3tensor(cholmat, 0, iNbTime - 1, 0, 6, 0, 6);

  if (CorrMatrix_loc)
    free_dmatrix(CorrMatrix_loc, 0, 6, 0, 6);

  mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

  if (dom_ft1)
    free_dvector(dom_ft1, 0, iNumPaths - 1);
  if (dom_ft2)
    free_dvector(dom_ft2, 0, iNumPaths - 1);
  if (dom_v)
    free_dvector(dom_v, 0, iNumPaths - 1);
  if (dom_psi1)
    free_dvector(dom_psi1, 0, iNumPaths - 1);
  if (dom_psi2)
    free_dvector(dom_psi2, 0, iNumPaths - 1);
  if (dom_psi12)
    free_dvector(dom_psi12, 0, iNumPaths - 1);

  if (for_ft1)
    free_dvector(for_ft1, 0, iNumPaths - 1);
  if (for_ft2)
    free_dvector(for_ft2, 0, iNumPaths - 1);
  if (for_v)
    free_dvector(for_v, 0, iNumPaths - 1);
  if (for_psi1)
    free_dvector(for_psi1, 0, iNumPaths - 1);
  if (for_psi2)
    free_dvector(for_psi2, 0, iNumPaths - 1);
  if (for_psi12)
    free_dvector(for_psi12, 0, iNumPaths - 1);

  if (fwdfx)
    free_dvector(fwdfx, 0, iNumPaths - 1);

  if (gauss)
    free_dvector(gauss, 0, NbFactor);

  ABS_Free(&rg);

  /* Return the error message */
  return err;
}

Err fxlgmsv_mc_balsam2(
    /*	Time Information  */
    int iNbTime, int iNbEvent, double *dTime, double *dDate,

    int iNumPaths,

    // Domestic
    /*	Model data Information	*/
    double ddomLambdaX1, double ddomLambdaX2,

    double *ddomSigma, double *ddomAlphaLGM, double *ddomRhoLGM,

    double *ddomAlpha, double *ddomLambdaEps, double *ddomLvlEps,
    double *ddomRho, double *ddomRho2,

    double *domzcvol1_star, double *domzcvol2_star,

    /* Parameters for DF(t  ,T*) reconstruction */
    double *ddomff_star, double *domgam1_star, double *domgam2_star,
    double *domgam1_2_star, double *domgam2_2_star, double *domgam12_star,

    // Foreign
    /*	Model data Information	*/
    double dforLambdaX1, double dforLambdaX2,

    double *dforSigma, double *dforAlphaLGM, double *dforRhoLGM,

    double *dforAlpha, double *dforLambdaEps, double *dforLvlEps,
    double *dforRho, double *dforRho2,

    double *forzcvol1_star, double *forzcvol2_star,

    /* Parameters for DF(t  ,T*) reconstruction */
    double *dforff_star, double *forgam1_star, double *forgam2_star,
    double *forgam1_2_star, double *forgam2_2_star, double *forgam12_star,

    // FX
    double fwd_fx_TStar, double *fx_vol,

    //	Correlation
    double ***CorrMatrix,
    /*	0 : Dom1
            1 : Dom2
            2 : DomSV
            3 : For1
            4 : For2
            5 : ForSV
            6 : FX	*/

    /*	Product data */
    void **func_parm_tab, int *EvalEvent,

    /* for Optimisation of exercise boundary */
    int do_optimisation, int *optimise, MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or
       NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(long path_index, double evt_date, double evt_time,
                       void *func_parm,

                       // Domestic
                       /* Model data	*/
                       double domft1, double domft2, double domphi1,
                       double domphi2, double domphi12, double domv,

                       // Foreign
                       /* Model data	*/
                       double forft1, double forft2, double forphi1,
                       double forphi2, double forphi12, double forv,

                       // FX
                       double fx_spot,

                       /* Vector of results to be updated */
                       int nprod,
                       /* Result	*/
                       double *prod_val, int *stop_path),
    /*	Result */
    int iNbProduct, double **res) {
  Err err = NULL;
  clock_t time1, time2;
  long i, j, k, l, m, n;

  double *sum_payoff = NULL, *sum_payoff2 = NULL, *path_payoff = NULL,
         *res_evt = NULL, **matrix = NULL, *matrixi = NULL,

         *dom_sigf1 = NULL, *dom_sigf2 = NULL, *dom_sigpsi1 = NULL,
         *dom_sigpsi2 = NULL, *dom_sigpsi12 = NULL, *dom_muv1 = NULL,
         *dom_muv2 = NULL, *dom_sigv = NULL,

         *for_sigf1 = NULL, *for_sigf2 = NULL, *for_sigpsi1 = NULL,
         *for_sigpsi2 = NULL, *for_sigpsi12 = NULL, *for_muv1 = NULL,
         *for_muv2 = NULL, *for_sigv = NULL,

         *fx_sig = NULL,

         **CorrMatrix_loc = NULL,

         ***cholmat = NULL,

         ***save_values = NULL;

  int stop_path;
  int evtindex;
  int NbFactor, NbFactorReduction;

  double dt, sqdt;
  double dom_ft1, dom_ft2, dom_v, dom_sqv, dom_psi1, dom_psi2, dom_psi12;
  double for_ft1, for_ft2, for_v, for_sqv, for_psi1, for_psi2, for_psi12;
  double fwdfx;
  double *brow;
  double dom_df, for_df;
  double domzcvol1, domzcvol2;
  double forzcvol1, forzcvol2;

  long seed = -123456789;
  SRandomGen rg;
  double *gauss = NULL;
  /* --------------------------------------------------------------------------------------------------
                                                                                  Transformation into 5 Factor
     --------------------------------------------------------------------------------------------------
   */
  if (ddomAlphaLGM[1] < 1.0e-6 && dforAlphaLGM[1] < 1.0e-6) {
    NbFactor = 5;
    NbFactorReduction = 1;
  } else {
    NbFactor = 7;
    NbFactorReduction = 0;
  }

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Initialisation
     --------------------------------------------------------------------------------------------------
   */

  /* For computational time calculation				 */
  time1 = clock();

  /* odd number of paths */
  iNumPaths = (int)(iNumPaths / 2) * 2 + 1;

  brow = dvector(0, 6);
  gauss = dvector(0, NbFactor);
  matrix = dmatrix(0, iNumPaths - 1, 0, NbFactor * (iNbTime - 1) - 1);
  res_evt = dvector(0, iNbProduct - 1);
  path_payoff = dvector(0, iNbProduct - 1);
  sum_payoff = dvector(0, iNbProduct - 1);
  sum_payoff2 = dvector(0, iNbProduct - 1);

  /* for precalculations */
  dom_sigf1 = dvector(0, iNbTime - 1);
  dom_sigf2 = dvector(0, iNbTime - 1);
  dom_sigpsi1 = dvector(0, iNbTime - 1);
  dom_sigpsi2 = dvector(0, iNbTime - 1);
  dom_sigpsi12 = dvector(0, iNbTime - 1);
  dom_muv1 = dvector(0, iNbTime - 1);
  dom_muv2 = dvector(0, iNbTime - 1);
  dom_sigv = dvector(0, iNbTime - 1);

  for_sigf1 = dvector(0, iNbTime - 1);
  for_sigf2 = dvector(0, iNbTime - 1);
  for_sigpsi1 = dvector(0, iNbTime - 1);
  for_sigpsi2 = dvector(0, iNbTime - 1);
  for_sigpsi12 = dvector(0, iNbTime - 1);
  for_muv1 = dvector(0, iNbTime - 1);
  for_muv2 = dvector(0, iNbTime - 1);
  for_sigv = dvector(0, iNbTime - 1);

  fx_sig = dvector(0, iNbTime - 1);

  CorrMatrix_loc = dmatrix(0, 6, 0, 6);

  cholmat = f3tensor(0, iNbTime - 1, 0, 6, 0, 6);

  if (!gauss || !matrix || !path_payoff || !sum_payoff || !sum_payoff2 ||
      !res_evt || !dom_sigf1 || !dom_sigf2 || !dom_sigpsi1 || !dom_sigpsi2 ||
      !dom_sigpsi12 || !dom_muv1 || !dom_muv2 || !dom_sigv || !for_sigf1 ||
      !for_sigf2 || !for_sigpsi1 || !for_sigpsi2 || !for_sigpsi12 ||
      !for_muv1 || !for_muv2 || !for_sigv || !fx_sig || !cholmat) {
    err = "Memory allocation failure in fxlgmsv_mc_balsam";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(iNumPaths, iNbEvent, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  memset(sum_payoff, 0, iNbProduct * sizeof(double));
  memset(sum_payoff2, 0, iNbProduct * sizeof(double));

  /* All the needed precalculations */
  for (j = 1; j < iNbTime; j++) {
    dt = dTime[j] - dTime[j - 1];
    sqdt = sqrt(dt);

    dom_sigf1[j] = ddomSigma[j] * sqdt;
    dom_sigf2[j] = dom_sigf1[j] * ddomAlphaLGM[j];
    dom_sigpsi1[j] = dom_sigf1[j] * dom_sigf1[j];
    dom_sigpsi2[j] = dom_sigf2[j] * dom_sigf2[j];
    dom_sigpsi12[j] = ddomRhoLGM[j] * dom_sigf1[j] * dom_sigf2[j];

    dom_sigv[j] = ddomAlpha[j] * sqdt;
    dom_muv1[j] = ddomLambdaEps[j] * dt;
    dom_muv2[j] = ddomLvlEps[j] * dt;

    for_sigf1[j] = dforSigma[j] * sqdt;
    for_sigf2[j] = for_sigf1[j] * dforAlphaLGM[j];
    for_sigpsi1[j] = for_sigf1[j] * for_sigf1[j];
    for_sigpsi2[j] = for_sigf2[j] * for_sigf2[j];
    for_sigpsi12[j] = dforRhoLGM[j] * for_sigf1[j] * for_sigf2[j];

    for_sigv[j] = dforAlpha[j] * sqdt;
    for_muv1[j] = dforLambdaEps[j] * dt;
    for_muv2[j] = dforLvlEps[j] * dt;

    fx_sig[j] = fx_vol[j] * sqdt;

    for (i = 0; i < 7; ++i) {
      for (k = 0; k < 7; ++k) {
        CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
      }
    }

    err = choldc(7, CorrMatrix_loc, cholmat[j]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -BalSam generation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* fill the Brownian matrix */
  err = balsam_generation(iNumPaths, NbFactor * (iNbTime - 1), matrix);
  if (err) {
    goto FREE_RETURN;
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Convolution
     --------------------------------------------------------------------------------------------------
   */
  err = ABS_Init(&rg, seed, iNumPaths, NbFactor, 0);
  if (err)
    goto FREE_RETURN;

  for (i = 0; i < iNumPaths; i++) {
    /* initialisation */
    evtindex = 0;
    stop_path = 0;

    dom_ft1 = 0;
    dom_ft2 = 0;
    dom_v = 1;
    dom_psi1 = 0;
    dom_psi2 = 0;
    dom_psi12 = 0;

    for_ft1 = 0;
    for_ft2 = 0;
    for_v = 1;
    for_psi1 = 0;
    for_psi2 = 0;
    for_psi12 = 0;

    fwdfx = log(fwd_fx_TStar);

    matrixi = matrix[i];
    memset(path_payoff, 0, iNbProduct * sizeof(double));

    if (func_parm_tab[0] && EvalEvent[0]) {
      dom_df = exp(ddomff_star[0]);
      for_df = exp(dforff_star[0]);

      /* Event at time 0.0 */
      err = payoff_func(0, dDate[0], dTime[0], func_parm_tab[0],

                        dom_ft1, dom_ft2, dom_psi1, dom_psi2, dom_psi12, dom_v,

                        for_ft1, for_ft2, for_psi1, for_psi2, for_psi12, for_v,

                        exp(fwdfx) * dom_df / for_df,

                        iNbProduct, res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < iNbProduct; k++) {
        path_payoff[k] += res_evt[k] / dom_df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[0], res_evt, i, dom_df,
                                       params);
      }

      evtindex++;
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < iNbTime; j++) {

      for (k = 0; k < 7; k++) {
        brow[k] = 0;
        err = rg.Gauss(&rg, &gauss[k]);

        if (!NbFactorReduction) {
          for (l = 0; l < 7; l++) {
            brow[k] += cholmat[j][k][l] * matrixi[m + l];
          }
        } else {
          n = 0;
          for (l = 0; l < 7; l++) {
            if (l != 1 && l != 4) {
              brow[k] += cholmat[j][k][l] * matrixi[m + n];
              n++;
            }
          }
        }
      }

      dom_sqv = sqrt(dom_v);
      for_sqv = sqrt(for_v);

      domzcvol1 = domzcvol1_star[j] * dom_sigf1[j] * dom_sqv;
      domzcvol2 = domzcvol2_star[j] * dom_sigf2[j] * dom_sqv;
      forzcvol1 = forzcvol1_star[j] * for_sigf1[j] * for_sqv;
      forzcvol2 = forzcvol2_star[j] * for_sigf2[j] * for_sqv;

      dom_ft1 += dom_sigf1[j] * dom_sqv * brow[0];
      dom_ft2 += dom_sigf2[j] * dom_sqv * brow[1];
      dom_psi1 += dom_sigpsi1[j] * dom_v;
      dom_psi2 += dom_sigpsi2[j] * dom_v;
      dom_psi12 += dom_sigpsi12[j] * dom_v;
      dom_v = max(dom_v + dom_muv2[j] - dom_muv1[j] * dom_v +
                      dom_sigv[j] * dom_sqv * brow[2],
                  0.0);

      for_ft1 += for_sigf1[j] * for_sqv * brow[3];
      for_ft1 += -for_sigf1[j] * for_sqv * CorrMatrix[j][3][6] * fx_sig[j];
      for_ft1 += -for_sigf1[j] * for_sqv *
                 (forzcvol1 + CorrMatrix[j][3][4] * forzcvol2);
      for_ft1 +=
          +for_sigf1[j] * for_sqv *
          (CorrMatrix[j][3][0] * domzcvol1 + CorrMatrix[j][3][1] * domzcvol2);

      for_ft2 += for_sigf2[j] * for_sqv * brow[4];
      for_ft2 += -for_sigf2[j] * for_sqv * CorrMatrix[j][4][6] * fx_sig[j];
      for_ft2 += -for_sigf2[j] * for_sqv *
                 (CorrMatrix[j][4][3] * forzcvol1 + forzcvol2);
      for_ft2 +=
          +for_sigf2[j] * for_sqv *
          (CorrMatrix[j][4][0] * domzcvol1 + CorrMatrix[j][4][1] * domzcvol2);

      for_psi1 += for_sigpsi1[j] * for_v;
      for_psi2 += for_sigpsi2[j] * for_v;
      for_psi12 += for_sigpsi12[j] * for_v;
      for_v +=
          for_muv2[j] - for_muv1[j] * for_v + for_sigv[j] * for_sqv * brow[5];
      for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][5][6] * fx_sig[j];
      for_v +=
          -for_sigv[j] * for_sqv *
          (CorrMatrix[j][5][3] * forzcvol1 + CorrMatrix[j][5][4] * forzcvol2);
      for_v +=
          +for_sigv[j] * for_sqv *
          (CorrMatrix[j][5][0] * domzcvol1 + CorrMatrix[j][5][1] * domzcvol2);
      for_v = max(0.0, for_v);

      fwdfx +=
          (fx_sig[j] * brow[6] - domzcvol1 * brow[0] - domzcvol2 * brow[1] +
           forzcvol1 * brow[3] + forzcvol2 * brow[4]);

      fwdfx += -0.5 * (fx_sig[j] * fx_sig[j] + domzcvol1 * domzcvol1 +
                       domzcvol2 * domzcvol2 + forzcvol1 * forzcvol1 +
                       forzcvol2 * forzcvol2);

      fwdfx += -(-fx_sig[j] * domzcvol1 * CorrMatrix[j][6][0] -
                 fx_sig[j] * domzcvol2 * CorrMatrix[j][6][1] +
                 fx_sig[j] * forzcvol1 * CorrMatrix[j][6][3] +
                 fx_sig[j] * forzcvol2 * CorrMatrix[j][6][4]

                 + domzcvol1 * domzcvol2 * CorrMatrix[j][0][1] -
                 domzcvol1 * forzcvol1 * CorrMatrix[j][0][3] -
                 domzcvol1 * forzcvol2 * CorrMatrix[j][0][4]

                 - domzcvol2 * forzcvol1 * CorrMatrix[j][1][3] -
                 domzcvol2 * forzcvol2 * CorrMatrix[j][1][4]

                 + forzcvol1 * forzcvol2 * CorrMatrix[j][3][4]);

      /*	Case of evaluation events */
      if (EvalEvent[j]) {
        dom_df = exp(ddomff_star[evtindex] + domgam1_star[evtindex] * dom_ft1 +
                     domgam2_star[evtindex] * dom_ft2 +
                     domgam1_2_star[evtindex] * dom_psi1 +
                     domgam2_2_star[evtindex] * dom_psi2 +
                     domgam12_star[evtindex] * dom_psi12);

        for_df = exp(dforff_star[evtindex] + forgam1_star[evtindex] * for_ft1 +
                     forgam2_star[evtindex] * for_ft2 +
                     forgam1_2_star[evtindex] * for_psi1 +
                     forgam2_2_star[evtindex] * for_psi2 +
                     forgam12_star[evtindex] * for_psi12);

        /* Modification of the Payoff at t */
        err =
            payoff_func(i, dDate[j], dTime[j], func_parm_tab[j], dom_ft1,
                        dom_ft2, dom_psi1, dom_psi2, dom_psi12, dom_v,

                        for_ft1, for_ft2, for_psi1, for_psi2, for_psi12, for_v,

                        exp(fwdfx) * dom_df / for_df,

                        iNbProduct, res_evt, &stop_path);

        if (err)
          goto FREE_RETURN;

        for (k = 0; k < iNbProduct; k++) {
          path_payoff[k] += res_evt[k] / dom_df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                         dom_df, params);
        }

        evtindex++;
      }

      m += NbFactor;
    }

    for (k = 0; k < iNbProduct; k++) {
      sum_payoff[k] += path_payoff[k] / iNumPaths;
      sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
    }

    if (do_optimisation && params->iKnockInCol) {
      /* we recopy in the col pay the pv of the column */
      for (j = 0; j < iNbEvent; j++) {
        if (optimise[j]) {
          save_values[j][params->iNbIndex][i] =
              path_payoff[(int)(save_values[j][params->iNbIndex][i] + 0.5)];
        }
      }
    }
  }

  for (k = 0; k < iNbProduct; k++) {
    res[k][0] = sum_payoff[k];
    res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

    if (res[k][1] > 0.0) {
      res[k][1] = sqrt(res[k][1]);
    } else {
      res[k][1] = 0.0;
    }
  }

  /* Convolution time display */
  time1 = clock();
  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(time1 - time2) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    /* Free the big matrix of memory first */
    if (matrix)
      free_dmatrix(matrix, 0, iNumPaths - 1, 0, NbFactor * (iNbTime - 1) - 1);
    matrix = NULL;

    time1 = clock();

    err = find_and_optimise_boundary(save_values, iNbEvent, iNumPaths, optimise,
                                     params, &(res[iNbProduct][0]),
                                     &(res[iNbProduct][1]));

    if (err)
      goto FREE_RETURN;

    time2 = clock();
    smessage("Phase 3 -optimisation  , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (matrix)
    free_dmatrix(matrix, 0, iNumPaths - 1, 0, NbFactor * (iNbTime - 1) - 1);
  if (res_evt)
    free_dvector(res_evt, 0, iNbProduct - 1);
  if (path_payoff)
    free_dvector(path_payoff, 0, iNbProduct - 1);
  if (sum_payoff)
    free_dvector(sum_payoff, 0, iNbProduct - 1);
  if (sum_payoff2)
    free_dvector(sum_payoff2, 0, iNbProduct - 1);

  if (dom_sigf1)
    free_dvector(dom_sigf1, 0, iNbTime - 1);
  if (dom_sigf2)
    free_dvector(dom_sigf2, 0, iNbTime - 1);
  if (dom_sigpsi1)
    free_dvector(dom_sigpsi1, 0, iNbTime - 1);
  if (dom_sigpsi2)
    free_dvector(dom_sigpsi2, 0, iNbTime - 1);
  if (dom_sigpsi12)
    free_dvector(dom_sigpsi12, 0, iNbTime - 1);
  if (dom_muv1)
    free_dvector(dom_muv1, 0, iNbTime - 1);
  if (dom_muv2)
    free_dvector(dom_muv2, 0, iNbTime - 1);
  if (dom_sigv)
    free_dvector(dom_sigv, 0, iNbTime - 1);

  if (for_sigf1)
    free_dvector(for_sigf1, 0, iNbTime - 1);
  if (for_sigf2)
    free_dvector(for_sigf2, 0, iNbTime - 1);
  if (for_sigpsi1)
    free_dvector(for_sigpsi1, 0, iNbTime - 1);
  if (for_sigpsi2)
    free_dvector(for_sigpsi2, 0, iNbTime - 1);
  if (for_sigpsi12)
    free_dvector(for_sigpsi12, 0, iNbTime - 1);
  if (for_muv1)
    free_dvector(for_muv1, 0, iNbTime - 1);
  if (for_muv2)
    free_dvector(for_muv2, 0, iNbTime - 1);
  if (for_sigv)
    free_dvector(for_sigv, 0, iNbTime - 1);

  if (fx_sig)
    free_dvector(fx_sig, 0, iNbTime - 1);

  if (brow)
    free_dvector(brow, 0, 6);

  if (cholmat)
    free_f3tensor(cholmat, 0, iNbTime - 1, 0, 6, 0, 6);

  if (CorrMatrix_loc)
    free_dmatrix(CorrMatrix_loc, 0, 6, 0, 6);

  mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

  ABS_Free(&rg);

  /* Return the error message */
  return err;
}

Err fxlgmsv_mc_balsam3(
    /*	Time Information  */
    int iNbTime, int iNbEvent, double *dTime, double *dDate,

    int iNumPaths,

    // Domestic
    /*	Model data Information	*/
    double ddomLambdaX1, double ddomLambdaX2,

    double *ddomSigma, double *ddomAlphaLGM, double *ddomRhoLGM,

    double *ddomAlpha, double *ddomLambdaEps, double *ddomLvlEps,
    double *ddomRho, double *ddomRho2,

    double *domzcvol1_star, double *domzcvol2_star,

    /* Parameters for DF(t  ,T*) reconstruction */
    double *ddomff_star, double *domgam1_star, double *domgam2_star,
    double *domgam1_2_star, double *domgam2_2_star, double *domgam12_star,

    // Foreign
    /*	Model data Information	*/
    double dforLambdaX1, double dforLambdaX2,

    double *dforSigma, double *dforAlphaLGM, double *dforRhoLGM,

    double *dforAlpha, double *dforLambdaEps, double *dforLvlEps,
    double *dforRho, double *dforRho2,

    double *forzcvol1_star, double *forzcvol2_star,

    /* Parameters for DF(t  ,T*) reconstruction */
    double *dforff_star, double *forgam1_star, double *forgam2_star,
    double *forgam1_2_star, double *forgam2_2_star, double *forgam12_star,

    // FX
    double fwd_fx_TStar, double *fx_vol,

    //	Correlation
    double ***CorrMatrix,
    /*	0 : Dom1
            1 : Dom2
            2 : DomSV
            3 : For1
            4 : For2
            5 : ForSV
            6 : FX	*/

    /*	Product data */
    void **func_parm_tab, int *EvalEvent,

    /* for Optimisation of exercise boundary */
    int do_optimisation, int *optimise, MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or
       NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(long path_index, double evt_date, double evt_time,
                       void *func_parm,

                       // Domestic
                       /* Model data	*/
                       double domft1, double domft2, double domphi1,
                       double domphi2, double domphi12, double domv,

                       // Foreign
                       /* Model data	*/
                       double forft1, double forft2, double forphi1,
                       double forphi2, double forphi12, double forv,

                       // FX
                       double fx_spot,

                       /* Vector of results to be updated */
                       int nprod,
                       /* Result	*/
                       double *prod_val, int *stop_path),
    /*	Result */
    int iNbProduct, double **res) {
  Err err = NULL;
  clock_t time1, time2;
  long i, j, k, l, m, n;

  double *sum_payoff = NULL, *sum_payoff2 = NULL, *path_payoff = NULL,
         *res_evt = NULL, **matrix = NULL, *matrixi = NULL,

         *dom_sigf1 = NULL, *dom_sigf2 = NULL, *dom_sigpsi1 = NULL,
         *dom_sigpsi2 = NULL, *dom_sigpsi12 = NULL, *dom_muv1 = NULL,
         *dom_muv2 = NULL, *dom_sigv = NULL,

         *for_sigf1 = NULL, *for_sigf2 = NULL, *for_sigpsi1 = NULL,
         *for_sigpsi2 = NULL, *for_sigpsi12 = NULL, *for_muv1 = NULL,
         *for_muv2 = NULL, *for_sigv = NULL,

         *fx_sig = NULL,

         **CorrMatrix_loc = NULL,

         ***cholmat = NULL,

         ***save_values = NULL;

  int stop_path;
  int evtindex;
  int NbFactor, NbFactorReduction;

  double dt, sqdt;
  double dom_ft1, dom_ft2, dom_v, dom_sqv, dom_psi1, dom_psi2, dom_psi12;
  double for_ft1, for_ft2, for_v, for_sqv, for_psi1, for_psi2, for_psi12;
  double fwdfx;
  double *brow;
  double dom_df, for_df;
  double domzcvol1, domzcvol2;
  double forzcvol1, forzcvol2;

  long seed = -123456789;

  /* --------------------------------------------------------------------------------------------------
                                                                                  Transformation into 5 Factor
     --------------------------------------------------------------------------------------------------
   */
  if (ddomAlphaLGM[1] < 1.0e-6 && dforAlphaLGM[1] < 1.0e-6) {
    NbFactor = 5;
    NbFactorReduction = 1;
  } else {
    NbFactor = 7;
    NbFactorReduction = 0;
  }

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Initialisation
     --------------------------------------------------------------------------------------------------
   */

  /* For computational time calculation				 */
  time1 = clock();

  /* odd number of paths */
  iNumPaths = (int)(iNumPaths / 2) * 2 + 1;

  brow = dvector(0, 6);
  matrix = dmatrix(0, iNumPaths - 1, 0, NbFactor * (iNbTime - 1) - 1);
  res_evt = dvector(0, iNbProduct - 1);
  path_payoff = dvector(0, iNbProduct - 1);
  sum_payoff = dvector(0, iNbProduct - 1);
  sum_payoff2 = dvector(0, iNbProduct - 1);

  /* for precalculations */
  dom_sigf1 = dvector(0, iNbTime - 1);
  dom_sigf2 = dvector(0, iNbTime - 1);
  dom_sigpsi1 = dvector(0, iNbTime - 1);
  dom_sigpsi2 = dvector(0, iNbTime - 1);
  dom_sigpsi12 = dvector(0, iNbTime - 1);
  dom_muv1 = dvector(0, iNbTime - 1);
  dom_muv2 = dvector(0, iNbTime - 1);
  dom_sigv = dvector(0, iNbTime - 1);

  for_sigf1 = dvector(0, iNbTime - 1);
  for_sigf2 = dvector(0, iNbTime - 1);
  for_sigpsi1 = dvector(0, iNbTime - 1);
  for_sigpsi2 = dvector(0, iNbTime - 1);
  for_sigpsi12 = dvector(0, iNbTime - 1);
  for_muv1 = dvector(0, iNbTime - 1);
  for_muv2 = dvector(0, iNbTime - 1);
  for_sigv = dvector(0, iNbTime - 1);

  fx_sig = dvector(0, iNbTime - 1);

  CorrMatrix_loc = dmatrix(0, 6, 0, 6);

  cholmat = f3tensor(0, iNbTime - 1, 0, 6, 0, 6);

  if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt ||
      !dom_sigf1 || !dom_sigf2 || !dom_sigpsi1 || !dom_sigpsi2 ||
      !dom_sigpsi12 || !dom_muv1 || !dom_muv2 || !dom_sigv || !for_sigf1 ||
      !for_sigf2 || !for_sigpsi1 || !for_sigpsi2 || !for_sigpsi12 ||
      !for_muv1 || !for_muv2 || !for_sigv || !fx_sig || !cholmat) {
    err = "Memory allocation failure in fxlgmsv_mc_balsam";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(iNumPaths, iNbEvent, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  memset(sum_payoff, 0, iNbProduct * sizeof(double));
  memset(sum_payoff2, 0, iNbProduct * sizeof(double));

  /* All the needed precalculations */
  for (j = 1; j < iNbTime; j++) {
    dt = dTime[j] - dTime[j - 1];
    sqdt = sqrt(dt);

    dom_sigf1[j] = ddomSigma[j] * sqdt;
    dom_sigf2[j] = dom_sigf1[j] * ddomAlphaLGM[j];
    dom_sigpsi1[j] = dom_sigf1[j] * dom_sigf1[j];
    dom_sigpsi2[j] = dom_sigf2[j] * dom_sigf2[j];
    dom_sigpsi12[j] = ddomRhoLGM[j] * dom_sigf1[j] * dom_sigf2[j];

    dom_sigv[j] = ddomAlpha[j] * sqdt;
    dom_muv1[j] = ddomLambdaEps[j] * dt;
    dom_muv2[j] = ddomLvlEps[j] * dt;

    for_sigf1[j] = dforSigma[j] * sqdt;
    for_sigf2[j] = for_sigf1[j] * dforAlphaLGM[j];
    for_sigpsi1[j] = for_sigf1[j] * for_sigf1[j];
    for_sigpsi2[j] = for_sigf2[j] * for_sigf2[j];
    for_sigpsi12[j] = dforRhoLGM[j] * for_sigf1[j] * for_sigf2[j];

    for_sigv[j] = dforAlpha[j] * sqdt;
    for_muv1[j] = dforLambdaEps[j] * dt;
    for_muv2[j] = dforLvlEps[j] * dt;

    fx_sig[j] = fx_vol[j] * sqdt;

    for (i = 0; i < 7; ++i) {
      for (k = 0; k < 7; ++k) {
        CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
      }
    }

    err = choldc(7, CorrMatrix_loc, cholmat[j]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -BalSam generation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* fill the Brownian matrix */
  err = balsam_generation(iNumPaths, NbFactor * (iNbTime - 1), matrix);
  if (err) {
    goto FREE_RETURN;
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Convolution
     --------------------------------------------------------------------------------------------------
   */

  for (i = 0; i < iNumPaths; i++) {
    /* initialisation */
    evtindex = 0;
    stop_path = 0;

    dom_ft1 = 0;
    dom_ft2 = 0;
    dom_v = 1;
    dom_psi1 = 0;
    dom_psi2 = 0;
    dom_psi12 = 0;

    for_ft1 = 0;
    for_ft2 = 0;
    for_v = 1;
    for_psi1 = 0;
    for_psi2 = 0;
    for_psi12 = 0;

    fwdfx = log(fwd_fx_TStar);

    matrixi = matrix[i];
    memset(path_payoff, 0, iNbProduct * sizeof(double));

    if (func_parm_tab[0] && EvalEvent[0]) {
      dom_df = exp(ddomff_star[0]);
      for_df = exp(dforff_star[0]);

      /* Event at time 0.0 */
      err = payoff_func(0, dDate[0], dTime[0], func_parm_tab[0],

                        dom_ft1, dom_ft2, dom_psi1, dom_psi2, dom_psi12, dom_v,

                        for_ft1, for_ft2, for_psi1, for_psi2, for_psi12, for_v,

                        exp(fwdfx) * dom_df / for_df,

                        iNbProduct, res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < iNbProduct; k++) {
        path_payoff[k] += res_evt[k] / dom_df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                       dom_df, params);
      }

      evtindex++;
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < iNbTime; j++) {

      for (k = 0; k < 7; k++) {
        brow[k] = 0;
        if (!NbFactorReduction) {
          for (l = 0; l < 7; l++) {
            brow[k] += cholmat[j][k][l] * matrixi[m + l];
          }
        } else {
          n = 0;
          for (l = 0; l < 7; l++) {
            if (l != 1 && l != 4) {
              brow[k] += cholmat[j][k][l] * matrixi[m + n];
              n++;
            }
          }
        }
      }

      dom_sqv = sqrt(dom_v);
      for_sqv = sqrt(for_v);

      domzcvol1 = domzcvol1_star[j] * dom_sigf1[j] * dom_sqv;
      domzcvol2 = domzcvol2_star[j] * dom_sigf2[j] * dom_sqv;
      forzcvol1 = forzcvol1_star[j] * for_sigf1[j] * for_sqv;
      forzcvol2 = forzcvol2_star[j] * for_sigf2[j] * for_sqv;

      dom_ft1 += dom_sigf1[j] * dom_sqv * brow[0];
      dom_ft2 += dom_sigf2[j] * dom_sqv * brow[1];
      dom_psi1 += dom_sigpsi1[j] * dom_v;
      dom_psi2 += dom_sigpsi2[j] * dom_v;
      dom_psi12 += dom_sigpsi12[j] * dom_v;
      dom_v = max(dom_v + dom_muv2[j] - dom_muv1[j] * dom_v +
                      dom_sigv[j] * dom_sqv * brow[2],
                  0.0);

      for_ft1 += for_sigf1[j] * for_sqv * brow[3];
      for_ft1 += -for_sigf1[j] * for_sqv * CorrMatrix[j][3][6] * fx_sig[j];
      for_ft1 += -for_sigf1[j] * for_sqv *
                 (forzcvol1 + CorrMatrix[j][3][4] * forzcvol2);
      for_ft1 +=
          +for_sigf1[j] * for_sqv *
          (CorrMatrix[j][3][0] * domzcvol1 + CorrMatrix[j][3][1] * domzcvol2);

      for_ft2 += for_sigf2[j] * for_sqv * brow[4];
      for_ft2 += -for_sigf2[j] * for_sqv * CorrMatrix[j][4][6] * fx_sig[j];
      for_ft2 += -for_sigf2[j] * for_sqv *
                 (CorrMatrix[j][4][3] * forzcvol1 + forzcvol2);
      for_ft2 +=
          +for_sigf2[j] * for_sqv *
          (CorrMatrix[j][4][0] * domzcvol1 + CorrMatrix[j][4][1] * domzcvol2);

      for_psi1 += for_sigpsi1[j] * for_v;
      for_psi2 += for_sigpsi2[j] * for_v;
      for_psi12 += for_sigpsi12[j] * for_v;
      for_v +=
          for_muv2[j] - for_muv1[j] * for_v + for_sigv[j] * for_sqv * brow[5];
      for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][5][6] * fx_sig[j];
      for_v +=
          -for_sigv[j] * for_sqv *
          (CorrMatrix[j][5][3] * forzcvol1 + CorrMatrix[j][5][4] * forzcvol2);
      for_v +=
          +for_sigv[j] * for_sqv *
          (CorrMatrix[j][5][0] * domzcvol1 + CorrMatrix[j][5][1] * domzcvol2);
      for_v = max(0.0, for_v);

      fwdfx +=
          (fx_sig[j] * brow[6] - domzcvol1 * brow[0] - domzcvol2 * brow[1] +
           forzcvol1 * brow[3] + forzcvol2 * brow[4]);

      fwdfx += -0.5 * (fx_sig[j] * fx_sig[j] + domzcvol1 * domzcvol1 +
                       domzcvol2 * domzcvol2 + forzcvol1 * forzcvol1 +
                       forzcvol2 * forzcvol2);

      fwdfx += -(-fx_sig[j] * domzcvol1 * CorrMatrix[j][6][0] -
                 fx_sig[j] * domzcvol2 * CorrMatrix[j][6][1] +
                 fx_sig[j] * forzcvol1 * CorrMatrix[j][6][3] +
                 fx_sig[j] * forzcvol2 * CorrMatrix[j][6][4]

                 + domzcvol1 * domzcvol2 * CorrMatrix[j][0][1] -
                 domzcvol1 * forzcvol1 * CorrMatrix[j][0][3] -
                 domzcvol1 * forzcvol2 * CorrMatrix[j][0][4]

                 - domzcvol2 * forzcvol1 * CorrMatrix[j][1][3] -
                 domzcvol2 * forzcvol2 * CorrMatrix[j][1][4]

                 + forzcvol1 * forzcvol2 * CorrMatrix[j][3][4]);

      /*	Case of evaluation events */
      if (EvalEvent[j]) {
        dom_df = exp(ddomff_star[evtindex] + domgam1_star[evtindex] * dom_ft1 +
                     domgam2_star[evtindex] * dom_ft2 +
                     domgam1_2_star[evtindex] * dom_psi1 +
                     domgam2_2_star[evtindex] * dom_psi2 +
                     domgam12_star[evtindex] * dom_psi12);

        for_df = exp(dforff_star[evtindex] + forgam1_star[evtindex] * for_ft1 +
                     forgam2_star[evtindex] * for_ft2 +
                     forgam1_2_star[evtindex] * for_psi1 +
                     forgam2_2_star[evtindex] * for_psi2 +
                     forgam12_star[evtindex] * for_psi12);

        /* Modification of the Payoff at t */
        err =
            payoff_func(i, dDate[j], dTime[j], func_parm_tab[j], dom_ft1,
                        dom_ft2, dom_psi1, dom_psi2, dom_psi12, dom_v,

                        for_ft1, for_ft2, for_psi1, for_psi2, for_psi12, for_v,

                        exp(fwdfx) * dom_df / for_df,

                        iNbProduct, res_evt, &stop_path);

        if (err)
          goto FREE_RETURN;

        for (k = 0; k < iNbProduct; k++) {
          path_payoff[k] += res_evt[k] / dom_df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                         dom_df, params);
        }

        evtindex++;
      }

      m += NbFactor;
    }

    for (k = 0; k < iNbProduct; k++) {
      sum_payoff[k] += path_payoff[k] / iNumPaths;
      sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
    }

    if (do_optimisation && params->iKnockInCol) {
      /* we recopy in the col pay the pv of the column */
      for (j = 0; j < iNbEvent; j++) {
        if (optimise[j]) {
          save_values[j][params->iNbIndex][i] =
              path_payoff[(int)(save_values[j][params->iNbIndex][i] + 0.5)];
        }
      }
    }
  }

  for (k = 0; k < iNbProduct; k++) {
    res[k][0] = sum_payoff[k];
    res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

    if (res[k][1] > 0.0) {
      res[k][1] = sqrt(res[k][1]);
    } else {
      res[k][1] = 0.0;
    }
  }

  /* Convolution time display */
  time1 = clock();
  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(time1 - time2) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    /* Free the big matrix of memory first */
    if (matrix)
      free_dmatrix(matrix, 0, iNumPaths - 1, 0, NbFactor * (iNbTime - 1) - 1);
    matrix = NULL;

    time1 = clock();

    err = find_and_optimise_boundary(save_values, iNbEvent, iNumPaths, optimise,
                                     params, &(res[iNbProduct][0]),
                                     &(res[iNbProduct][1]));

    if (err)
      goto FREE_RETURN;

    time2 = clock();
    smessage("Phase 3 -optimisation  , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (matrix)
    free_dmatrix(matrix, 0, iNumPaths - 1, 0, NbFactor * (iNbTime - 1) - 1);
  if (res_evt)
    free_dvector(res_evt, 0, iNbProduct - 1);
  if (path_payoff)
    free_dvector(path_payoff, 0, iNbProduct - 1);
  if (sum_payoff)
    free_dvector(sum_payoff, 0, iNbProduct - 1);
  if (sum_payoff2)
    free_dvector(sum_payoff2, 0, iNbProduct - 1);

  if (dom_sigf1)
    free_dvector(dom_sigf1, 0, iNbTime - 1);
  if (dom_sigf2)
    free_dvector(dom_sigf2, 0, iNbTime - 1);
  if (dom_sigpsi1)
    free_dvector(dom_sigpsi1, 0, iNbTime - 1);
  if (dom_sigpsi2)
    free_dvector(dom_sigpsi2, 0, iNbTime - 1);
  if (dom_sigpsi12)
    free_dvector(dom_sigpsi12, 0, iNbTime - 1);
  if (dom_muv1)
    free_dvector(dom_muv1, 0, iNbTime - 1);
  if (dom_muv2)
    free_dvector(dom_muv2, 0, iNbTime - 1);
  if (dom_sigv)
    free_dvector(dom_sigv, 0, iNbTime - 1);

  if (for_sigf1)
    free_dvector(for_sigf1, 0, iNbTime - 1);
  if (for_sigf2)
    free_dvector(for_sigf2, 0, iNbTime - 1);
  if (for_sigpsi1)
    free_dvector(for_sigpsi1, 0, iNbTime - 1);
  if (for_sigpsi2)
    free_dvector(for_sigpsi2, 0, iNbTime - 1);
  if (for_sigpsi12)
    free_dvector(for_sigpsi12, 0, iNbTime - 1);
  if (for_muv1)
    free_dvector(for_muv1, 0, iNbTime - 1);
  if (for_muv2)
    free_dvector(for_muv2, 0, iNbTime - 1);
  if (for_sigv)
    free_dvector(for_sigv, 0, iNbTime - 1);

  if (fx_sig)
    free_dvector(fx_sig, 0, iNbTime - 1);

  if (brow)
    free_dvector(brow, 0, 6);

  if (cholmat)
    free_f3tensor(cholmat, 0, iNbTime - 1, 0, 6, 0, 6);

  if (CorrMatrix_loc)
    free_dmatrix(CorrMatrix_loc, 0, 6, 0, 6);

  mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

  /* Return the error message */
  return err;
}

Err qtolgmsv2f_mc_balsam(
    //	Time Information
    int iNbTime, int iNbEvent, double *dTime, double *dDate,

    int iNumPaths,

    // Domestic
    //	Model data Information
    double ddomLambdaX,

    double *ddomSigma,

    double *domzcvol_star,

    // Parameters for DF(t  ,T*) reconstruction
    double *ddomff_star, double *domgam1_star, double *domgam1_2_star,

    // Foreign
    //	Model data Information
    double dforLambdaX1, double dforLambdaX2,

    double *dforSigma, double *dforAlphaLGM, double *dforRhoLGM,

    double *dforAlpha, double *dforLambdaEps, double *dforLvlEps,
    double *dforRho, double *dforRho2,

    double *forzcvol1_star, double *forzcvol2_star,

    // Parameters for DF(t  ,T*) reconstruction
    double *dforff_star, double *forgam1_star, double *forgam2_star,
    double *forgam1_2_star, double *forgam2_2_star, double *forgam12_star,

    // FX Vol
    double *fx_vol,

    //	Correlation
    double ***CorrMatrix,
    //	0 : Dom
    //	1 : For1
    //	2 : For2
    //	3 : ForSV
    //	4 : FX

    //	Product data
    void **func_parm_tab, int *EvalEvent,

    // for Optimisation of exercise boundary
    int do_optimisation, int *optimise, MCEBPARAMS params,

    //	Initialisation function to be called at the beggining of each path or
    //NULL if none
    void (*init_func)(),

    //	Payoff function
    Err (*payoff_func)(long path_index, double evt_date, double evt_time,
                       void *func_parm,

                       // Domestic
                       // Model data
                       double domft, double domphi,

                       // Foreign
                       // Model data
                       double forft1, double forft2, double forphi1,
                       double forphi2, double forphi12, double forv,

                       // Vector of results to be updated
                       int nprod,
                       // Result
                       double *prod_val, int *stop_path),
    // Result
    int iNbProduct, double **res) {
  Err err = NULL;
  clock_t time1, time2;
  long i, j, k, l, m;

  double *sum_payoff = NULL, *sum_payoff2 = NULL, *path_payoff = NULL,
         *res_evt = NULL, **matrix = NULL, *matrixi = NULL,

         *dom_sigf = NULL, *dom_sigpsi = NULL,

         *for_sigf1 = NULL, *for_sigf2 = NULL, *for_sigpsi1 = NULL,
         *for_sigpsi2 = NULL, *for_sigpsi12 = NULL, *for_muv1 = NULL,
         *for_muv2 = NULL, *for_sigv = NULL,

         *fx_sig = NULL,

         **CorrMatrix_loc = NULL,

         ***cholmat = NULL,

         ***save_values = NULL;

  int stop_path;
  int evtindex;

  double dt, sqdt;
  double dom_ft, dom_psi;
  double for_ft1, for_ft2, for_v, for_sqv, for_psi1, for_psi2, for_psi12;
  double *brow;
  double dom_df, for_df;
  double domzcvol;
  double forzcvol1, forzcvol2;

  long seed = -123456789;

  //-------------------------------------------------------------------------------------------------
  //										Initialisation
  //-------------------------------------------------------------------------------------------------

  // For computational time calculation
  time1 = clock();

  // odd number of paths
  iNumPaths = (int)(iNumPaths / 2) * 2 + 1;

  brow = dvector(0, 3);
  matrix = dmatrix(0, iNumPaths - 1, 0, 4 * (iNbTime - 1) - 1);
  res_evt = dvector(0, iNbProduct - 1);
  path_payoff = dvector(0, iNbProduct - 1);
  sum_payoff = dvector(0, iNbProduct - 1);
  sum_payoff2 = dvector(0, iNbProduct - 1);

  // for precalculations
  dom_sigf = dvector(0, iNbTime - 1);
  dom_sigpsi = dvector(0, iNbTime - 1);

  for_sigf1 = dvector(0, iNbTime - 1);
  for_sigf2 = dvector(0, iNbTime - 1);
  for_sigpsi1 = dvector(0, iNbTime - 1);
  for_sigpsi2 = dvector(0, iNbTime - 1);
  for_sigpsi12 = dvector(0, iNbTime - 1);
  for_muv1 = dvector(0, iNbTime - 1);
  for_muv2 = dvector(0, iNbTime - 1);
  for_sigv = dvector(0, iNbTime - 1);

  fx_sig = dvector(0, iNbTime - 1);

  CorrMatrix_loc = dmatrix(0, 4, 0, 4);

  cholmat = f3tensor(0, iNbTime - 1, 0, 3, 0, 3);

  if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt ||
      !dom_sigf || !dom_sigpsi || !for_sigf1 || !for_sigf2 || !for_sigpsi1 ||
      !for_sigpsi2 || !for_sigpsi12 || !for_muv1 || !for_muv2 || !for_sigv ||
      !fx_sig || !cholmat) {
    err = "Memory allocation failure in qtolgmsv2f_mc_balsam";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(iNumPaths, iNbEvent, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  memset(sum_payoff, 0, iNbProduct * sizeof(double));
  memset(sum_payoff2, 0, iNbProduct * sizeof(double));

  // All the needed precalculations
  for (j = 1; j < iNbTime; j++) {
    dt = dTime[j] - dTime[j - 1];
    sqdt = sqrt(dt);

    dom_sigf[j] = ddomSigma[j] * sqdt;
    dom_sigpsi[j] = dom_sigf[j] * dom_sigf[j];

    for_sigf1[j] = dforSigma[j] * sqdt;
    for_sigf2[j] = for_sigf1[j] * dforAlphaLGM[j];
    for_sigpsi1[j] = for_sigf1[j] * for_sigf1[j];
    for_sigpsi2[j] = for_sigf2[j] * for_sigf2[j];
    for_sigpsi12[j] = dforRhoLGM[j] * for_sigf1[j] * for_sigf2[j];

    for_sigv[j] = dforAlpha[j] * sqdt;
    for_muv1[j] = dforLambdaEps[j] * dt;
    for_muv2[j] = dforLvlEps[j] * dt;

    fx_sig[j] = fx_vol[j] * sqdt;

    for (i = 0; i < 5; ++i) {
      for (k = 0; k < 5; ++k) {
        CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
      }
    }

    err = choldc(4, CorrMatrix_loc, cholmat[j]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  // Initialisation time display
  time2 = clock();
  smessage("Phase 1 -BalSam generation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  // fill the Brownian matrix
  err = balsam_generation(iNumPaths, 4 * (iNbTime - 1), matrix);

  if (err) {
    goto FREE_RETURN;
  }

  // Initialisation time display
  time2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  // --------------------------------------------------------------------------------------------------
  //												Convolution
  // --------------------------------------------------------------------------------------------------

  for (i = 0; i < iNumPaths; i++) {
    // initialisation
    evtindex = 0;
    stop_path = 0;

    dom_ft = 0;
    dom_psi = 0;

    for_ft1 = 0;
    for_ft2 = 0;
    for_v = 1;
    for_psi1 = 0;
    for_psi2 = 0;
    for_psi12 = 0;

    matrixi = matrix[i];
    memset(path_payoff, 0, iNbProduct * sizeof(double));

    if (func_parm_tab[0] && EvalEvent[0]) {
      dom_df = exp(ddomff_star[0]);
      for_df = exp(dforff_star[0]);

      // Event at time 0.0
      err = payoff_func(0, dDate[0], dTime[0], func_parm_tab[0],

                        dom_ft, dom_psi,

                        for_ft1, for_ft2, for_psi1, for_psi2, for_psi12, for_v,

                        iNbProduct, res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < iNbProduct; k++) {
        path_payoff[k] += res_evt[k] / dom_df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                       dom_df, params);
      }

      evtindex++;
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < iNbTime; j++) {

      for (k = 0; k < 4; k++) {
        brow[k] = 0;
        for (l = 0; l < 4; l++) {
          brow[k] += cholmat[j][k][l] * matrixi[m + l];
        }
      }

      dom_ft += dom_sigf[j] * brow[0];
      dom_psi += dom_sigpsi[j];

      for_sqv = sqrt(for_v);

      domzcvol = domzcvol_star[j] * dom_sigf[j];
      forzcvol1 = forzcvol1_star[j] * for_sigf1[j] * for_sqv;
      forzcvol2 = forzcvol2_star[j] * for_sigf2[j] * for_sqv;

      for_ft1 += for_sigf1[j] * for_sqv * brow[1];
      for_ft1 += -for_sigf1[j] * for_sqv * CorrMatrix[j][1][4] * fx_sig[j];
      for_ft1 += -for_sigf1[j] * for_sqv *
                 (forzcvol1 + CorrMatrix[j][1][2] * forzcvol2);
      for_ft1 += +for_sigf1[j] * for_sqv * (CorrMatrix[j][1][0] * domzcvol);

      for_ft2 += for_sigf2[j] * for_sqv * brow[2];
      for_ft2 += -for_sigf2[j] * for_sqv * CorrMatrix[j][2][4] * fx_sig[j];

      // Previous code:
      // for_ft2 += -for_sigf2[j] * for_sqv * ( CorrMatrix[j][2][3] * forzcvol1
      // + forzcvol2 );
      for_ft2 += -for_sigf2[j] * for_sqv *
                 (CorrMatrix[j][2][1] * forzcvol1 + forzcvol2);
      for_ft2 += +for_sigf2[j] * for_sqv * (CorrMatrix[j][2][0] * domzcvol);

      for_psi1 += for_sigpsi1[j] * for_v;
      for_psi2 += for_sigpsi2[j] * for_v;
      for_psi12 += for_sigpsi12[j] * for_v;
      for_v +=
          for_muv2[j] - for_muv1[j] * for_v + for_sigv[j] * for_sqv * brow[3];
      for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][3][4] * fx_sig[j];
      for_v +=
          -for_sigv[j] * for_sqv *
          (CorrMatrix[j][3][1] * forzcvol1 + CorrMatrix[j][3][2] * forzcvol2);
      for_v += +for_sigv[j] * for_sqv * (CorrMatrix[j][3][0] * domzcvol);
      for_v = max(0.0, for_v);

      //	Case of evaluation events
      if (EvalEvent[j]) {
        dom_df = exp(ddomff_star[evtindex] + domgam1_star[evtindex] * dom_ft +
                     domgam1_2_star[evtindex] * dom_psi);

        for_df = exp(dforff_star[evtindex] + forgam1_star[evtindex] * for_ft1 +
                     forgam2_star[evtindex] * for_ft2 +
                     forgam1_2_star[evtindex] * for_psi1 +
                     forgam2_2_star[evtindex] * for_psi2 +
                     forgam12_star[evtindex] * for_psi12);

        // Modification of the Payoff at t
        err =
            payoff_func(i, dDate[j], dTime[j], func_parm_tab[j],

                        dom_ft, dom_psi,

                        for_ft1, for_ft2, for_psi1, for_psi2, for_psi12, for_v,

                        iNbProduct, res_evt, &stop_path);

        if (err)
          goto FREE_RETURN;

        for (k = 0; k < iNbProduct; k++) {
          path_payoff[k] += res_evt[k] / dom_df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                         dom_df, params);
        }

        evtindex++;
      }

      m += 4;
    }

    for (k = 0; k < iNbProduct; k++) {
      sum_payoff[k] += path_payoff[k] / iNumPaths;
      sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
    }

    if (do_optimisation && params->iKnockInCol) {
      // we recopy in the col pay the pv of the column
      for (j = 0; j < iNbEvent; j++) {
        if (optimise[j]) {
          save_values[j][params->iNbIndex][i] =
              path_payoff[(int)(save_values[j][params->iNbIndex][i] + 0.5)];
        }
      }
    }
  }

  for (k = 0; k < iNbProduct; k++) {
    res[k][0] = sum_payoff[k];
    res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

    if (res[k][1] > 0.0) {
      res[k][1] = sqrt(res[k][1]);
    } else {
      res[k][1] = 0.0;
    }
  }

  // Convolution time display
  time1 = clock();
  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(time1 - time2) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    // Free the big matrix of memory first
    if (matrix)
      free_dmatrix(matrix, 0, iNumPaths - 1, 0, 4 * (iNbTime - 1) - 1);
    matrix = NULL;

    time1 = clock();

    err = find_and_optimise_boundary(save_values, iNbEvent, iNumPaths, optimise,
                                     params, &(res[iNbProduct][0]),
                                     &(res[iNbProduct][1]));

    if (err)
      goto FREE_RETURN;

    time2 = clock();
    smessage("Phase 3 -optimisation  , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (matrix)
    free_dmatrix(matrix, 0, iNumPaths - 1, 0, 4 * (iNbTime - 1) - 1);
  if (res_evt)
    free_dvector(res_evt, 0, iNbProduct - 1);
  if (path_payoff)
    free_dvector(path_payoff, 0, iNbProduct - 1);
  if (sum_payoff)
    free_dvector(sum_payoff, 0, iNbProduct - 1);
  if (sum_payoff2)
    free_dvector(sum_payoff2, 0, iNbProduct - 1);

  if (dom_sigf)
    free_dvector(dom_sigf, 0, iNbTime - 1);
  if (dom_sigpsi)
    free_dvector(dom_sigpsi, 0, iNbTime - 1);

  if (for_sigf1)
    free_dvector(for_sigf1, 0, iNbTime - 1);
  if (for_sigf2)
    free_dvector(for_sigf2, 0, iNbTime - 1);
  if (for_sigpsi1)
    free_dvector(for_sigpsi1, 0, iNbTime - 1);
  if (for_sigpsi2)
    free_dvector(for_sigpsi2, 0, iNbTime - 1);
  if (for_sigpsi12)
    free_dvector(for_sigpsi12, 0, iNbTime - 1);
  if (for_muv1)
    free_dvector(for_muv1, 0, iNbTime - 1);
  if (for_muv2)
    free_dvector(for_muv2, 0, iNbTime - 1);
  if (for_sigv)
    free_dvector(for_sigv, 0, iNbTime - 1);

  if (fx_sig)
    free_dvector(fx_sig, 0, iNbTime - 1);

  if (brow)
    free_dvector(brow, 0, 3);

  if (cholmat)
    free_f3tensor(cholmat, 0, iNbTime - 1, 0, 3, 0, 3);

  if (CorrMatrix_loc)
    free_dmatrix(CorrMatrix_loc, 0, 4, 0, 4);

  mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

  // Return the error message
  return err;
}

Err qtolgmsv1f_mc_balsam(
    /*	Time Information  */
    int iNbTime, int iNbEvent, double *dTime, double *dDate,

    int iNumPaths,

    // Domestic
    /*	Model data Information	*/
    double ddomLambdaX,

    double *ddomSigma,

    double *domzcvol_star,

    /* Parameters for DF(t  ,T*) reconstruction */
    double *ddomff_star, double *domgam1_star, double *domgam1_2_star,

    // Foreign
    /*	Model data Information	*/
    double dforLambdaX1,

    double *dforSigma,

    double *dforAlpha, double *dforLambdaEps, double *dforLvlEps,
    double *dforRho,

    double *forzcvol_star,

    /* Parameters for DF(t  ,T*) reconstruction */
    double *dforff_star, double *forgam1_star, double *forgam1_2_star,

    // FX Vol
    double *fx_vol,

    //	Correlation
    double ***CorrMatrix,
    /*	0 : Dom
            1 : For1
            2 : ForSV
            3 : Fx	*/

    /*	Product data */
    void **func_parm_tab, int *EvalEvent,

    /* for Optimisation of exercise boundary */
    int do_optimisation, int *optimise, MCEBPARAMS params,

    /*	Initialisation function to be called at the beggining of each path or
       NULL if none */
    void (*init_func)(),

    /*	Payoff function */
    Err (*payoff_func)(long path_index, double evt_date, double evt_time,
                       void *func_parm,

                       // Domestic
                       /* Model data	*/
                       double domft, double domphi,

                       // Foreign
                       /* Model data	*/
                       double forft, double forphi, double forv,

                       /* Vector of results to be updated */
                       int nprod,
                       /* Result	*/
                       double *prod_val, int *stop_path),
    /*	Result */
    int iNbProduct, double **res) {
  Err err = NULL;
  clock_t time1, time2;
  long i, j, k, l, m;

  double *sum_payoff = NULL, *sum_payoff2 = NULL, *path_payoff = NULL,
         *res_evt = NULL, **matrix = NULL, *matrixi = NULL,

         *dom_sigf = NULL, *dom_sigpsi = NULL,

         *for_sigf = NULL, *for_sigpsi = NULL, *for_muv1 = NULL,
         *for_muv2 = NULL, *for_sigv = NULL,

         *fx_sig = NULL,

         **CorrMatrix_loc = NULL,

         ***cholmat = NULL,

         ***save_values = NULL;

  int stop_path;
  int evtindex;

  double dt, sqdt;
  double dom_ft, dom_psi;
  double for_ft, for_v, for_sqv, for_psi;
  double *brow;
  double dom_df, for_df;
  double domzcvol;
  double forzcvol;

  long seed = -123456789;

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Initialisation
     --------------------------------------------------------------------------------------------------
   */

  /* For computational time calculation				 */
  time1 = clock();

  /* odd number of paths */
  iNumPaths = (int)(iNumPaths / 2) * 2 + 1;

  brow = dvector(0, 2);
  matrix = dmatrix(0, iNumPaths - 1, 0, 3 * (iNbTime - 1) - 1);
  res_evt = dvector(0, iNbProduct - 1);
  path_payoff = dvector(0, iNbProduct - 1);
  sum_payoff = dvector(0, iNbProduct - 1);
  sum_payoff2 = dvector(0, iNbProduct - 1);

  /* for precalculations */
  dom_sigf = dvector(0, iNbTime - 1);
  dom_sigpsi = dvector(0, iNbTime - 1);

  for_sigf = dvector(0, iNbTime - 1);
  for_sigpsi = dvector(0, iNbTime - 1);
  for_muv1 = dvector(0, iNbTime - 1);
  for_muv2 = dvector(0, iNbTime - 1);
  for_sigv = dvector(0, iNbTime - 1);

  fx_sig = dvector(0, iNbTime - 1);

  CorrMatrix_loc = dmatrix(0, 3, 0, 3);

  cholmat = f3tensor(0, iNbTime - 1, 0, 2, 0, 2);

  if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt ||
      !dom_sigf || !dom_sigpsi || !for_sigf || !for_sigpsi || !for_muv1 ||
      !for_muv2 || !for_sigv || !fx_sig || !cholmat) {
    err = "Memory allocation failure in qtolgmsv1f_mc_balsam";
    goto FREE_RETURN;
  }

  if (do_optimisation) {
    err = mceb_allocate_savevalues_for_GRFN(iNumPaths, iNbEvent, params,
                                            &save_values);

    if (err)
      goto FREE_RETURN;
  }

  memset(sum_payoff, 0, iNbProduct * sizeof(double));
  memset(sum_payoff2, 0, iNbProduct * sizeof(double));

  /* All the needed precalculations */
  for (j = 1; j < iNbTime; j++) {
    dt = dTime[j] - dTime[j - 1];
    sqdt = sqrt(dt);

    dom_sigf[j] = ddomSigma[j] * sqdt;
    dom_sigpsi[j] = dom_sigf[j] * dom_sigf[j];

    for_sigf[j] = dforSigma[j] * sqdt;
    for_sigpsi[j] = for_sigf[j] * for_sigf[j];

    for_sigv[j] = dforAlpha[j] * sqdt;
    for_muv1[j] = dforLambdaEps[j] * dt;
    for_muv2[j] = dforLvlEps[j] * dt;

    fx_sig[j] = fx_vol[j] * sqdt;

    for (i = 0; i < 4; ++i) {
      for (k = 0; k < 4; ++k) {
        CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
      }
    }

    err = choldc(3, CorrMatrix_loc, cholmat[j]);
    if (err) {
      goto FREE_RETURN;
    }
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -BalSam generation  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* fill the Brownian matrix */
  err = balsam_generation(iNumPaths, 3 * (iNbTime - 1), matrix);
  if (err) {
    goto FREE_RETURN;
  }

  /* Initialisation time display */
  time2 = clock();
  smessage("Phase 1 -preprocessing  , time in sec: %.2f",
           (double)(time2 - time1) / CLOCKS_PER_SEC);

  /* --------------------------------------------------------------------------------------------------
                                                                                                  Convolution
     --------------------------------------------------------------------------------------------------
   */

  for (i = 0; i < iNumPaths; i++) {
    /* initialisation */
    evtindex = 0;
    stop_path = 0;

    dom_ft = 0;
    dom_psi = 0;

    for_ft = 0;
    for_v = 1;
    for_psi = 0;

    matrixi = matrix[i];
    memset(path_payoff, 0, iNbProduct * sizeof(double));

    if (func_parm_tab[0] && EvalEvent[0]) {
      dom_df = exp(ddomff_star[0]);
      for_df = exp(dforff_star[0]);

      /* Event at time 0.0 */
      err = payoff_func(0, dDate[0], dTime[0], func_parm_tab[0],

                        dom_ft, dom_psi,

                        for_ft, for_psi, for_v,

                        iNbProduct, res_evt, &stop_path);

      if (err)
        goto FREE_RETURN;

      for (k = 0; k < iNbProduct; k++) {
        path_payoff[k] += res_evt[k] / dom_df;
      }

      if (do_optimisation) {
        mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                       dom_df, params);
      }

      evtindex++;
    }

    m = 0;

    for (j = 1; stop_path == 0 && j < iNbTime; j++) {

      for (k = 0; k < 3; k++) {
        brow[k] = 0;
        for (l = 0; l < 3; l++) {
          brow[k] += cholmat[j][k][l] * matrixi[m + l];
        }
      }

      dom_ft += dom_sigf[j] * brow[0];
      dom_psi += dom_sigpsi[j];

      for_sqv = sqrt(for_v);

      domzcvol = domzcvol_star[j] * dom_sigf[j];
      forzcvol = forzcvol_star[j] * for_sigf[j] * for_sqv;

      for_ft += for_sigf[j] * for_sqv * brow[1];
      for_ft += -for_sigf[j] * for_sqv * CorrMatrix[j][1][3] * fx_sig[j];
      for_ft += -for_sigf[j] * for_sqv * forzcvol;
      for_ft += +for_sigf[j] * for_sqv * CorrMatrix[j][1][0] * domzcvol;

      for_psi += for_sigpsi[j] * for_v;
      for_v +=
          for_muv2[j] - for_muv1[j] * for_v + for_sigv[j] * for_sqv * brow[2];
      for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][2][3] * fx_sig[j];
      for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][2][1] * forzcvol;
      for_v += +for_sigv[j] * for_sqv * CorrMatrix[j][2][0] * domzcvol;
      for_v = max(0.0, for_v);

      /*	Case of evaluation events */
      if (EvalEvent[j]) {
        dom_df = exp(ddomff_star[evtindex] + domgam1_star[evtindex] * dom_ft +
                     domgam1_2_star[evtindex] * dom_psi);

        for_df = exp(dforff_star[evtindex] + forgam1_star[evtindex] * for_ft +
                     forgam1_2_star[evtindex] * for_psi);

        /* Modification of the Payoff at t */
        err = payoff_func(i, dDate[j], dTime[j], func_parm_tab[j],

                          dom_ft, dom_psi,

                          for_ft, for_psi, for_v,

                          iNbProduct, res_evt, &stop_path);

        if (err)
          goto FREE_RETURN;

        for (k = 0; k < iNbProduct; k++) {
          path_payoff[k] += res_evt[k] / dom_df;
        }

        if (do_optimisation) {
          mceb_fill_savevalues_from_GRFN(save_values[evtindex], res_evt, i,
                                         dom_df, params);
        }

        evtindex++;
      }

      m += 3;
    }

    for (k = 0; k < iNbProduct; k++) {
      sum_payoff[k] += path_payoff[k] / iNumPaths;
      sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
    }

    if (do_optimisation && params->iKnockInCol) {
      /* we recopy in the col pay the pv of the column */
      for (j = 0; j < iNbEvent; j++) {
        if (optimise[j]) {
          save_values[j][params->iNbIndex][i] =
              path_payoff[(int)(save_values[j][params->iNbIndex][i] + 0.5)];
        }
      }
    }
  }

  for (k = 0; k < iNbProduct; k++) {
    res[k][0] = sum_payoff[k];
    res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

    if (res[k][1] > 0.0) {
      res[k][1] = sqrt(res[k][1]);
    } else {
      res[k][1] = 0.0;
    }
  }

  /* Convolution time display */
  time1 = clock();
  smessage("Phase 2 -convolution  , time in sec: %.2f",
           (double)(time1 - time2) / CLOCKS_PER_SEC);

  if (do_optimisation) {
    /* Free the big matrix of memory first */
    if (matrix)
      free_dmatrix(matrix, 0, iNumPaths - 1, 0, 3 * (iNbTime - 1) - 1);
    matrix = NULL;

    time1 = clock();

    err = find_and_optimise_boundary(save_values, iNbEvent, iNumPaths, optimise,
                                     params, &(res[iNbProduct][0]),
                                     &(res[iNbProduct][1]));

    if (err)
      goto FREE_RETURN;

    time2 = clock();
    smessage("Phase 3 -optimisation  , time in sec: %.2f",
             (double)(time2 - time1) / CLOCKS_PER_SEC);
  }

FREE_RETURN:

  if (matrix)
    free_dmatrix(matrix, 0, iNumPaths - 1, 0, 3 * (iNbTime - 1) - 1);
  if (res_evt)
    free_dvector(res_evt, 0, iNbProduct - 1);
  if (path_payoff)
    free_dvector(path_payoff, 0, iNbProduct - 1);
  if (sum_payoff)
    free_dvector(sum_payoff, 0, iNbProduct - 1);
  if (sum_payoff2)
    free_dvector(sum_payoff2, 0, iNbProduct - 1);

  if (dom_sigf)
    free_dvector(dom_sigf, 0, iNbTime - 1);
  if (dom_sigpsi)
    free_dvector(dom_sigpsi, 0, iNbTime - 1);

  if (for_sigf)
    free_dvector(for_sigf, 0, iNbTime - 1);
  if (for_sigpsi)
    free_dvector(for_sigpsi, 0, iNbTime - 1);
  if (for_muv1)
    free_dvector(for_muv1, 0, iNbTime - 1);
  if (for_muv2)
    free_dvector(for_muv2, 0, iNbTime - 1);
  if (for_sigv)
    free_dvector(for_sigv, 0, iNbTime - 1);

  if (fx_sig)
    free_dvector(fx_sig, 0, iNbTime - 1);

  if (brow)
    free_dvector(brow, 0, 2);

  if (cholmat)
    free_f3tensor(cholmat, 0, iNbTime - 1, 0, 2, 0, 2);

  if (CorrMatrix_loc)
    free_dmatrix(CorrMatrix_loc, 0, 3, 0, 3);

  mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

  /* Return the error message */
  return err;
}

// REWRITING THE QTOLGMSV2F_MC_BALSAM
/*
Err	 qtolgmsv2f_mc_balsam(

                                        int			iNbTime  ,
                                        int			iNbEvent  ,
                                        double		*dTime  ,
                                        double		*dDate  ,

                                        int			iNumPaths  ,

                                        //Domestic

                                        double		ddomLambdaX  ,

                                        double		*ddomSigma  ,

                                        double		*domzcvol_star  ,


                                        double		*ddomff_star  ,
                                        double		*domgam1_star  ,
                                        double		*domgam1_2_star  ,

                                        //Foreign

                                        double		dforLambdaX1  ,
                                        double		dforLambdaX2  ,

                                        double		*dforSigma  ,
                                        double		*dforAlphaLGM  ,
                                        double		*dforRhoLGM  ,

                                        double		*dforAlpha  ,
                                        double		*dforLambdaEps  ,
                                        double		*dforLvlEps  ,
                                        double		*dforRho  ,
                                        double		*dforRho2  ,

                                        double		*forzcvol1_star  ,
                                        double		*forzcvol2_star  ,


                                        double		*dforff_star  ,
                                        double		*forgam1_star  ,
                                        double		*forgam2_star  ,
                                        double		*forgam1_2_star  ,
                                        double		*forgam2_2_star  ,
                                        double		*forgam12_star  ,

                                        //FX Vol
                                        double		*fx_vol  ,

                                        //	Correlation
                                        double		***CorrMatrix  ,
                                        //	0 : Dom
                                        //	1 : For1
                                        //	2 : For2
                                        //	3 : ForSV
                                        //	4 : FX


                                        void		**func_parm_tab  ,
                                        int			*EvalEvent  ,


                                        int do_optimisation  ,
                                        int			*optimise  ,
                                        MCEBPARAMS	params  ,


                                        void		(*init_func)()  ,

                                        // Payoff function of QUANTOLGMSV2F
                                        Err (*payoff_func)(	long
path_index  , double	evt_date  , double	evt_time  ,
                                                                                void	*func_parm  ,

                                                                                //Domestic

                                                                                double	domft  ,
                                                                                double	domphi  ,

                                                                                //Foreign

                                                                                double	forft1  ,
                                                                                double	forft2  ,
                                                                                double	forphi1  ,
                                                                                double	forphi2  ,
                                                                                double	forphi12  ,
                                                                                double	forv  ,


                                                                                int		nprod  ,

                                                                                double	*prod_val  ,
                                                                                int		*stop_path)  ,

                                        int			iNbProduct  ,
                                        double		**res)
{
Err		err = NULL;
clock_t	time1  , time2;
long	i  , j  , k  , l  , m;

double	*sum_payoff			= NULL  ,
                *sum_payoff2		= NULL  ,
                *path_payoff		= NULL  ,
                *res_evt			= NULL  ,
                **matrix			= NULL  ,
                *matrixi			= NULL  ,

                *dom_sigf1			= NULL  ,
                *dom_sigpsi1		= NULL  ,

                *for_sigf1			= NULL  ,
                *for_sigf2			= NULL  ,
                *for_sigpsi1		= NULL  ,
                *for_sigpsi2		= NULL  ,
                *for_sigpsi12		= NULL  ,
                *for_muv1			= NULL  ,
                *for_muv2			= NULL  ,
                *for_sigv			= NULL  ,

                *fx_sig				= NULL  ,

                **CorrMatrix_loc	= NULL  ,

                ***cholmat			= NULL  ,

                ***save_values		= NULL  ,
                *opt_bar			= NULL  ,
                **coef_lin			= NULL;

int		stop_path;
int		evtindex;

double	dt  , sqdt;
double	dom_ft1  , dom_psi1;
double	for_ft1  , for_ft2  , for_v  , for_sqv  , for_psi1  , for_psi2  ,
for_psi12; double	fwdfx; double	*brow; double	dom_df  , for_df;
//double	domzcvol1  ,domzcvol2;
double	domzcvol1;
double	forzcvol1  , forzcvol2;

long	seed = -123456789;

        time1 = clock();


        iNumPaths = (int) (iNumPaths / 2) * 2 + 1;

        brow = dvector (0  , 3);
        matrix = dmatrix (0  , iNumPaths - 1  , 0  , 4 * (iNbTime - 1) - 1);
        res_evt = dvector (0  , iNbProduct - 1);
        path_payoff = dvector (0  , iNbProduct - 1);
        sum_payoff = dvector (0  , iNbProduct - 1);
        sum_payoff2 = dvector (0  , iNbProduct - 1);

        // for precalculations
        dom_sigf1 = dvector(0  , iNbTime - 1);
        dom_sigpsi1 = dvector(0  , iNbTime - 1);

        for_sigf1 = dvector(0  , iNbTime - 1);
        for_sigf2 = dvector(0  , iNbTime - 1);
        for_sigpsi1 = dvector(0  , iNbTime - 1);
        for_sigpsi2 = dvector(0  , iNbTime - 1);
        for_sigpsi12 = dvector(0  , iNbTime - 1);
        for_muv1 = dvector(0  , iNbTime - 1);
        for_muv2 = dvector(0  , iNbTime - 1);
        for_sigv = dvector(0  , iNbTime - 1);

        fx_sig = dvector(0  , iNbTime - 1);

        CorrMatrix_loc = dmatrix(0  , 3  , 0  , 3);

        cholmat = f3tensor(0  , iNbTime - 1  , 0  , 3  , 0  , 3);

        if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt
|| !dom_sigf1 || !dom_sigpsi1 || !for_sigf1 || !for_sigf2 || !for_sigpsi1 ||
!for_sigpsi2 || !for_sigpsi12 || !for_muv1 || !for_muv2 || !for_sigv || !fx_sig
|| !cholmat)
        {
                err = "Memory allocation failure in fxlgmsv_mc_balsam";
                goto FREE_RETURN;
        }

        if (do_optimisation)
        {
                save_values = f3tensor(0  , iNbEvent - 1  , 0  ,
params->iNbIndex  , 0  , iNumPaths - 1); opt_bar = dvector(0  , iNbEvent - 1);
                coef_lin = dmatrix(0  , iNbEvent - 1  , 0  , params->iNbIndex);

                if (!save_values || !opt_bar || !coef_lin)
                {
                        err = "Memory allocation (2) failure in
fxlgmsv_mc_balsam"; goto FREE_RETURN;
                }
        }

        memset (sum_payoff  , 0  , iNbProduct * sizeof (double));
        memset (sum_payoff2  , 0  , iNbProduct * sizeof (double));


        for (j=1; j<iNbTime; j++)
        {
                dt = dTime[j] - dTime[j-1];
                sqdt = sqrt(dt);

                dom_sigf1[j] = ddomSigma[j] * sqdt;
                dom_sigpsi1[j] = dom_sigf1[j] * dom_sigf1[j];

                for_sigf1[j] = dforSigma[j] * sqdt;
                for_sigf2[j] = for_sigf1[j] * dforAlphaLGM[j];
                for_sigpsi1[j] = for_sigf1[j] * for_sigf1[j];
                for_sigpsi2[j] = for_sigf2[j] * for_sigf2[j];
                for_sigpsi12[j] = dforRhoLGM[j] * for_sigf1[j] * for_sigf2[j];

                for_sigv[j] = dforAlpha[j] * sqdt;
                for_muv1[j] = dforLambdaEps[j] * dt;
                for_muv2[j] = dforLvlEps[j] * dt;

                fx_sig[j] = fx_vol[j] * sqdt;

                for(i=0;i<4;++i)
                {
                        for(k=0;k<4;++k)
                        {
                                CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
                        }
                }

                err = choldc (4  , CorrMatrix_loc  , cholmat[j]);
                if (err)
                {
                        goto FREE_RETURN;
                }
        }


        time2 = clock();
        smessage ("Phase 1 -BalSam generation  , time in sec: %.2f"  , (double)
(time2 - time1) / CLOCKS_PER_SEC);


        err = balsam_generation (iNumPaths  , 4 * (iNbTime - 1)  , matrix);
        if (err)
        {
                goto FREE_RETURN;
        }


        time2 = clock();
        smessage ("Phase 1 -preprocessing  , time in sec: %.2f"  , (double)
(time2 - time1) / CLOCKS_PER_SEC);



        for (i=0; i<iNumPaths; i++)
        {

                evtindex = 0;
                stop_path = 0;

                dom_ft1 = 0;
                dom_psi1 = 0;

                for_ft1 = 0;
                for_ft2 = 0;
                for_v = 1;
                for_psi1 = 0;
                for_psi2 = 0;
                for_psi12 = 0;

                fwdfx = 0;

                matrixi = matrix[i];
                memset (path_payoff  , 0  , iNbProduct * sizeof (double));

                if (func_parm_tab[0] && EvalEvent[0])
                {
                        dom_df = exp (ddomff_star[0]);
                        for_df = exp (dforff_star[0]);

                        // QUANTOLGMSV Payoff
                                err = payoff_func(	0  ,
                                                                dDate[0]  ,
                                                                dTime[0]  ,
                                                                func_parm_tab[0]
,

                                                                dom_ft1  ,
                                                                dom_psi1  ,

                                                                for_ft1  ,
                                                                for_ft2  ,
                                                                for_psi1  ,
                                                                for_psi2  ,
                                                                for_psi12  ,
                                                                for_v  ,

                                                                iNbProduct  ,
                                                                res_evt  ,
                                                                &stop_path);


                        if (err) goto FREE_RETURN;

                        for (k=0; k<iNbProduct; k++)
                        {
                                path_payoff[k] += res_evt[k] / dom_df;
                        }

                        if (do_optimisation)
                        {
                                for (i=0; i<iNumPaths-1; i++)
                                {
                                        for (k=0; k<params->iNbIndex; k++)
                                        {
                                                save_values[0][k][i] =
res_evt[params->iColBound+k];
                                        }

                                        if (!params->iKnockInCol)
                                        {
                                                save_values[0][params->iNbIndex][i]
= res_evt[params->iColPay] / dom_df;
                                        }
                                        else
                                        {
                                                save_values[0][params->iNbIndex][i]
= res_evt[params->iColPay];
                                        }
                                }
                        }

                        evtindex++;
                }

                m = 0;

                for (j=1; stop_path == 0 && j<iNbTime; j++)
                {

                        for (k=0; k<4; k++)
                        {
                                brow[k] = 0;
                                for (l=0; l<4; l++)
                                {
                                        brow[k] += cholmat[j][k][l] *
matrixi[m+l];
                                }
                        }

                        for_sqv = sqrt(for_v);

                        domzcvol1 = domzcvol_star[j] * dom_sigf1[j];
                        forzcvol1 = forzcvol1_star[j] * for_sigf1[j] * for_sqv;
                        forzcvol2 = forzcvol2_star[j] * for_sigf2[j] * for_sqv;

                        dom_ft1 += dom_sigf1[j] * brow[0];
                        dom_psi1 += dom_sigpsi1[j];

                        for_ft1 += for_sigf1[j] * for_sqv * brow[1];
                        for_ft1 += -for_sigf1[j] * for_sqv * CorrMatrix[j][1][4]
* fx_sig[j]; for_ft1 += -for_sigf1[j] * for_sqv * ( forzcvol1  +
CorrMatrix[j][1][2] * forzcvol2 ); for_ft1 += +for_sigf1[j] * for_sqv * (
CorrMatrix[j][0][1] * domzcvol1);

                        for_ft2 += for_sigf2[j] * for_sqv * brow[2];
                        for_ft2 += -for_sigf2[j] * for_sqv * CorrMatrix[j][2][4]
* fx_sig[j]; for_ft2 += -for_sigf2[j] * for_sqv * ( CorrMatrix[j][1][2] *
forzcvol1  + forzcvol2 ); for_ft2 += +for_sigf2[j] * for_sqv * (
CorrMatrix[j][2][0] * domzcvol1);

                        for_psi1 += for_sigpsi1[j] * for_v;
                        for_psi2 += for_sigpsi2[j] * for_v;
                        for_psi12 += for_sigpsi12[j] * for_v;
                        for_v += for_muv2[j] - for_muv1[j] * for_v + for_sigv[j]
* for_sqv * brow[3]; for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][3][4] *
fx_sig[j]; for_v += -for_sigv[j] * for_sqv * ( CorrMatrix[j][1][3] * forzcvol1
+ CorrMatrix[j][2][3] * forzcvol2 ); for_v += +for_sigv[j] * for_sqv * (
CorrMatrix[j][3][0] * domzcvol1); for_v = max(0.0  , for_v);

                        if (EvalEvent[j])
                        {
                                dom_df = exp (ddomff_star[evtindex] +
domgam1_star[evtindex] * dom_ft1
                                                + domgam1_2_star[evtindex] *
dom_psi1);

                                // USELESS SINCE NO FX DIFFUSION
                                //for_df = exp (dforff_star[evtindex] +
forgam1_star[evtindex] * for_ft1 + forgam2_star[evtindex] * for_ft2
                                //		+ forgam1_2_star[evtindex] *
for_psi1 + forgam2_2_star[evtindex] * for_psi2 + forgam12_star[evtindex] *
for_psi12);

                                // QUANTOLGMSV Payoff
                                err = payoff_func(	i  ,
                                                                        dDate[j]
, dTime[j]  , func_parm_tab[j]  ,

                                                                        dom_ft1
, dom_psi1  ,

                                                                        for_ft1
, for_ft2  , for_psi1  , for_psi2  , for_psi12  , for_v  ,

                                                                        iNbProduct
, res_evt  , &stop_path);


                                if (err) goto FREE_RETURN;

                                for (k=0; k<iNbProduct; k++)
                                {
                                        path_payoff[k] += res_evt[k] / dom_df;
                                }

                                if (do_optimisation)
                                {
                                        for (k=0; k<params->iNbIndex; k++)
                                        {
                                                save_values[evtindex][k][i] =
res_evt[params->iColBound+k];
                                        }

                                        if (!params->iKnockInCol)
                                        {
                                                save_values[evtindex][params->iNbIndex][i]
= res_evt[params->iColPay] / dom_df;
                                        }
                                        else
                                        {
                                                save_values[evtindex][params->iNbIndex][i]
= res_evt[params->iColPay];
                                        }
                                }

                                evtindex++;
                        }

                        m += 4;
                }

                for (k=0; k<iNbProduct; k++)
                {
                        sum_payoff[k] += path_payoff[k] / iNumPaths;
                        sum_payoff2[k] += path_payoff[k] * path_payoff[k] /
iNumPaths;
                }

                if (do_optimisation && params->iKnockInCol)
                {

                        for (j=0; j<iNbEvent; j++)
                        {
                                if (optimise[j])
                                {
                                        save_values[j][params->iNbIndex][i] =
path_payoff[(int) (save_values[j][params->iNbIndex][i] + 0.5)];
                                }
                        }
                }
        }

        for (k=0; k<iNbProduct; k++)
        {
                res[k][0] = sum_payoff[k];
                res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) /
iNumPaths;

                if (res[k][1] > 0.0)
                {
                        res[k][1] = sqrt(res[k][1]);
                }
                else
                {
                        res[k][1] = 0.0;
                }
        }


        time1 = clock();
        smessage ("Phase 2 -convolution  , time in sec: %.2f"  , (double) (time1
- time2) / CLOCKS_PER_SEC);

        if (do_optimisation)
        {

                if (matrix) free_dmatrix (matrix  , 0  , iNumPaths - 1  , 0  , 4
* (iNbTime - 1) - 1); matrix = NULL;

                time1 = clock();

                err = find_and_optimise_boundary(	save_values  ,
                                                                                        iNbEvent  ,
                                                                                        iNumPaths  ,
                                                                                        optimise  ,
                                                                                        params  ,
                                                                                        opt_bar  ,
                                                                                        coef_lin  ,
                                                                                        &(res[iNbProduct][0])  ,
                                                                                        &(res[iNbProduct][1]));

                if (err) goto FREE_RETURN;

                for (j=0; j<iNbEvent; j++)
                {
                        res[j][2] = opt_bar[j];

                        for (k=0; k<params->iNbIndex; k++)
                        {
                                res[j][3+k] = coef_lin[j][k+1];
                        }
                }

                time2 = clock();
                smessage ("Phase 3 -optimisation  , time in sec: %.2f"  ,
(double) (time2 - time1) / CLOCKS_PER_SEC);
        }

FREE_RETURN:

        if (matrix) free_dmatrix (matrix  , 0  , iNumPaths - 1  , 0  , 4 *
(iNbTime - 1) - 1); if (res_evt) free_dvector (res_evt  , 0  , iNbProduct - 1);
        if (path_payoff) free_dvector (path_payoff  , 0  , iNbProduct - 1);
        if (sum_payoff) free_dvector (sum_payoff  , 0  , iNbProduct - 1);
        if (sum_payoff2) free_dvector (sum_payoff2  , 0  , iNbProduct - 1);

        if (dom_sigf1) free_dvector (dom_sigf1  , 0  , iNbTime - 1);
        if (dom_sigpsi1) free_dvector (dom_sigpsi1  , 0  , iNbTime - 1);

        if (for_sigf1) free_dvector (for_sigf1  , 0  , iNbTime - 1);
        if (for_sigf2) free_dvector (for_sigf2  , 0  , iNbTime - 1);
        if (for_sigpsi1) free_dvector (for_sigpsi1  , 0  , iNbTime - 1);
        if (for_sigpsi2) free_dvector (for_sigpsi2  , 0  , iNbTime - 1);
        if (for_sigpsi12) free_dvector (for_sigpsi12  , 0  , iNbTime - 1);
        if (for_muv1) free_dvector (for_muv1  , 0  , iNbTime - 1);
        if (for_muv2) free_dvector (for_muv2  , 0  , iNbTime - 1);
        if (for_sigv) free_dvector (for_sigv  , 0  , iNbTime - 1);

        if (fx_sig) free_dvector (fx_sig  , 0  , iNbTime - 1);

        if (brow) free_dvector (brow  , 0  , 3);

        if (cholmat) free_f3tensor(cholmat  , 0  , iNbTime - 1  , 0  , 3  , 0  ,
3);

        if (CorrMatrix_loc) free_dmatrix(CorrMatrix_loc  , 0  , 3  , 0  , 3);

        if (save_values) free_f3tensor(save_values  , 0  , iNbEvent - 1  , 0  ,
2  , 0  , iNumPaths - 1); if (opt_bar) free_dvector(opt_bar  , 0  , iNbEvent -
1); if (coef_lin) free_dmatrix(coef_lin  , 0  , iNbEvent - 1  , 0  ,
params->iNbIndex);


        return err;
}*/

/*
// qtolgmsv2f_mc_balsam using the fxlgmsv scheme in order to check the accuracy
of the
// truncated qtolgmsv2f with respect to the fxlgmsv

Err	 qtolgmsv2f_mc_balsam(

                                        int			iNbTime  ,
                                        int			iNbEvent  ,
                                        double		*dTime  ,
                                        double		*dDate  ,

                                        int			iNumPaths  ,

                                        //Domestic

                                        double		ddomLambdaX  ,

                                        double		*ddomSigma  ,

                                        double		*domzcvol_star  ,


                                        double		*ddomff_star  ,
                                        double		*domgam1_star  ,
                                        double		*domgam1_2_star  ,

                                        //Foreign

                                        double		dforLambdaX1  ,
                                        double		dforLambdaX2  ,

                                        double		*dforSigma  ,
                                        double		*dforAlphaLGM  ,
                                        double		*dforRhoLGM  ,

                                        double		*dforAlpha  ,
                                        double		*dforLambdaEps  ,
                                        double		*dforLvlEps  ,
                                        double		*dforRho  ,
                                        double		*dforRho2  ,

                                        double		*forzcvol1_star  ,
                                        double		*forzcvol2_star  ,


                                        double		*dforff_star  ,
                                        double		*forgam1_star  ,
                                        double		*forgam2_star  ,
                                        double		*forgam1_2_star  ,
                                        double		*forgam2_2_star  ,
                                        double		*forgam12_star  ,

                                        //FX Vol
                                        double		*fx_vol  ,

                                        //	Correlation
                                        double		***CorrMatrix  ,
                                        //	0 : Dom
                                        //	1 : For1
                                        //	2 : For2
                                        //	3 : ForSV
                                        //	4 : FX


                                        void		**func_parm_tab  ,
                                        int			*EvalEvent  ,


                                        int do_optimisation  ,
                                        int			*optimise  ,
                                        MCEBPARAMS	params  ,


                                        void		(*init_func)()  ,

                                        // Payoff function of QUANTOLGMSV2F
                                        Err (*payoff_func)(	long
path_index  , double	evt_date  , double	evt_time  ,
                                                                                void	*func_parm  ,

                                                                                //Domestic

                                                                                double	domft  ,
                                                                                double	domphi  ,

                                                                                //Foreign

                                                                                double	forft1  ,
                                                                                double	forft2  ,
                                                                                double	forphi1  ,
                                                                                double	forphi2  ,
                                                                                double	forphi12  ,
                                                                                double	forv  ,


                                                                                int		nprod  ,

                                                                                double	*prod_val  ,
                                                                                int		*stop_path)  ,

                                        int			iNbProduct  ,
                                        double		**res)
{
Err		err = NULL;
clock_t	time1  , time2;
long	i  , j  , k  , l  , m;

double	*sum_payoff			= NULL  ,
                *sum_payoff2		= NULL  ,
                *path_payoff		= NULL  ,
                *res_evt			= NULL  ,
                **matrix			= NULL  ,
                *matrixi			= NULL  ,

                *dom_sigf1			= NULL  ,
                *dom_sigpsi1		= NULL  ,

                *for_sigf1			= NULL  ,
                *for_sigf2			= NULL  ,
                *for_sigpsi1		= NULL  ,
                *for_sigpsi2		= NULL  ,
                *for_sigpsi12		= NULL  ,
                *for_muv1			= NULL  ,
                *for_muv2			= NULL  ,
                *for_sigv			= NULL  ,

                *fx_sig				= NULL  ,

                **CorrMatrix_loc	= NULL  ,

                ***cholmat			= NULL  ,

                ***save_values		= NULL  ,
                *opt_bar			= NULL  ,
                **coef_lin			= NULL;

int		stop_path;
int		evtindex;

double	dt  , sqdt;
double	dom_ft1  , dom_psi1;
double	for_ft1  , for_ft2  , for_v  , for_sqv  , for_psi1  , for_psi2  ,
for_psi12; double	fwdfx; double	*brow; double	dom_df  , for_df;
//double	domzcvol1  ,domzcvol2;
double	domzcvol1;
double	forzcvol1  , forzcvol2;

long	seed = -123456789;

        time1 = clock();


        iNumPaths = (int) (iNumPaths / 2) * 2 + 1;

        brow = dvector (0  , 6);
        matrix = dmatrix (0  , iNumPaths - 1  , 0  , 7 * (iNbTime - 1) - 1);
        res_evt = dvector (0  , iNbProduct - 1);
        path_payoff = dvector (0  , iNbProduct - 1);
        sum_payoff = dvector (0  , iNbProduct - 1);
        sum_payoff2 = dvector (0  , iNbProduct - 1);

        // for precalculations
        dom_sigf1 = dvector(0  , iNbTime - 1);
        dom_sigpsi1 = dvector(0  , iNbTime - 1);

        for_sigf1 = dvector(0  , iNbTime - 1);
        for_sigf2 = dvector(0  , iNbTime - 1);
        for_sigpsi1 = dvector(0  , iNbTime - 1);
        for_sigpsi2 = dvector(0  , iNbTime - 1);
        for_sigpsi12 = dvector(0  , iNbTime - 1);
        for_muv1 = dvector(0  , iNbTime - 1);
        for_muv2 = dvector(0  , iNbTime - 1);
        for_sigv = dvector(0  , iNbTime - 1);

        fx_sig = dvector(0  , iNbTime - 1);

        CorrMatrix_loc = dmatrix(0  , 6  , 0  , 6);

        cholmat = f3tensor(0  , iNbTime - 1  , 0  , 6  , 0  , 6);

        if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt
|| !dom_sigf1 || !dom_sigpsi1 || !for_sigf1 || !for_sigf2 || !for_sigpsi1 ||
!for_sigpsi2 || !for_sigpsi12 || !for_muv1 || !for_muv2 || !for_sigv || !fx_sig
|| !cholmat)
        {
                err = "Memory allocation failure in fxlgmsv_mc_balsam";
                goto FREE_RETURN;
        }

        if (do_optimisation)
        {
                save_values = f3tensor(0  , iNbEvent - 1  , 0  ,
params->iNbIndex  , 0  , iNumPaths - 1); opt_bar = dvector(0  , iNbEvent - 1);
                coef_lin = dmatrix(0  , iNbEvent - 1  , 0  , params->iNbIndex);

                if (!save_values || !opt_bar || !coef_lin)
                {
                        err = "Memory allocation (2) failure in
fxlgmsv_mc_balsam"; goto FREE_RETURN;
                }
        }

        memset (sum_payoff  , 0  , iNbProduct * sizeof (double));
        memset (sum_payoff2  , 0  , iNbProduct * sizeof (double));


        for (j=1; j<iNbTime; j++)
        {
                dt = dTime[j] - dTime[j-1];
                sqdt = sqrt(dt);

                dom_sigf1[j] = ddomSigma[j] * sqdt;
                dom_sigpsi1[j] = dom_sigf1[j] * dom_sigf1[j];

                for_sigf1[j] = dforSigma[j] * sqdt;
                for_sigf2[j] = for_sigf1[j] * dforAlphaLGM[j];
                for_sigpsi1[j] = for_sigf1[j] * for_sigf1[j];
                for_sigpsi2[j] = for_sigf2[j] * for_sigf2[j];
                for_sigpsi12[j] = dforRhoLGM[j] * for_sigf1[j] * for_sigf2[j];

                for_sigv[j] = dforAlpha[j] * sqdt;
                for_muv1[j] = dforLambdaEps[j] * dt;
                for_muv2[j] = dforLvlEps[j] * dt;

                fx_sig[j] = fx_vol[j] * sqdt;

                for(i=0;i<7;++i)
                {
                        for(k=0;k<7;++k)
                        {
                                CorrMatrix_loc[i][k] = CorrMatrix[j][i][k];
                        }
                }

                err = choldc (7  , CorrMatrix_loc  , cholmat[j]);
                if (err)
                {
                        goto FREE_RETURN;
                }
        }


        time2 = clock();
        smessage ("Phase 1 -BalSam generation  , time in sec: %.2f"  , (double)
(time2 - time1) / CLOCKS_PER_SEC);


        err = balsam_generation (iNumPaths  , 7 * (iNbTime - 1)  , matrix);
        if (err)
        {
                goto FREE_RETURN;
        }


        time2 = clock();
        smessage ("Phase 1 -preprocessing  , time in sec: %.2f"  , (double)
(time2 - time1) / CLOCKS_PER_SEC);



        for (i=0; i<iNumPaths; i++)
        {

                evtindex = 0;
                stop_path = 0;

                dom_ft1 = 0;
                dom_psi1 = 0;

                for_ft1 = 0;
                for_ft2 = 0;
                for_v = 1;
                for_psi1 = 0;
                for_psi2 = 0;
                for_psi12 = 0;

                fwdfx = 0;

                matrixi = matrix[i];
                memset (path_payoff  , 0  , iNbProduct * sizeof (double));

                if (func_parm_tab[0] && EvalEvent[0])
                {
                        dom_df = exp (ddomff_star[0]);
                        for_df = exp (dforff_star[0]);

                        // QUANTOLGMSV Payoff
                                err = payoff_func(	0  ,
                                                                dDate[0]  ,
                                                                dTime[0]  ,
                                                                func_parm_tab[0]
,

                                                                dom_ft1  ,
                                                                dom_psi1  ,

                                                                for_ft1  ,
                                                                for_ft2  ,
                                                                for_psi1  ,
                                                                for_psi2  ,
                                                                for_psi12  ,
                                                                for_v  ,

                                                                iNbProduct  ,
                                                                res_evt  ,
                                                                &stop_path);


                        if (err) goto FREE_RETURN;

                        for (k=0; k<iNbProduct; k++)
                        {
                                path_payoff[k] += res_evt[k] / dom_df;
                        }

                        if (do_optimisation)
                        {
                                for (i=0; i<iNumPaths-1; i++)
                                {
                                        for (k=0; k<params->iNbIndex; k++)
                                        {
                                                save_values[0][k][i] =
res_evt[params->iColBound+k];
                                        }

                                        if (!params->iKnockInCol)
                                        {
                                                save_values[0][params->iNbIndex][i]
= res_evt[params->iColPay] / dom_df;
                                        }
                                        else
                                        {
                                                save_values[0][params->iNbIndex][i]
= res_evt[params->iColPay];
                                        }
                                }
                        }

                        evtindex++;
                }

                m = 0;

                for (j=1; stop_path == 0 && j<iNbTime; j++)
                {

                        for (k=0; k<7; k++)
                        {
                                brow[k] = 0;
                                for (l=0; l<7; l++)
                                {
                                        brow[k] += cholmat[j][k][l] *
matrixi[m+l];
                                }
                        }

                        for_sqv = sqrt(for_v);

                        domzcvol1 = domzcvol_star[j] * dom_sigf1[j];
                        forzcvol1 = forzcvol1_star[j] * for_sigf1[j] * for_sqv;
                        forzcvol2 = forzcvol2_star[j] * for_sigf2[j] * for_sqv;

                        dom_ft1 += dom_sigf1[j] * brow[0];
                        dom_psi1 += dom_sigpsi1[j];

                        for_ft1 += for_sigf1[j] * for_sqv * brow[3];
                        for_ft1 += -for_sigf1[j] * for_sqv * CorrMatrix[j][3][6]
* fx_sig[j]; for_ft1 += -for_sigf1[j] * for_sqv * ( forzcvol1  +
CorrMatrix[j][3][4] * forzcvol2 ); for_ft1 += +for_sigf1[j] * for_sqv * (
CorrMatrix[j][3][0] * domzcvol1);

                        for_ft2 += for_sigf2[j] * for_sqv * brow[4];
                        for_ft2 += -for_sigf2[j] * for_sqv * CorrMatrix[j][4][6]
* fx_sig[j]; for_ft2 += -for_sigf2[j] * for_sqv * ( CorrMatrix[j][4][3] *
forzcvol1  + forzcvol2 ); for_ft2 += +for_sigf2[j] * for_sqv * (
CorrMatrix[j][4][0] * domzcvol1);

                        for_psi1 += for_sigpsi1[j] * for_v;
                        for_psi2 += for_sigpsi2[j] * for_v;
                        for_psi12 += for_sigpsi12[j] * for_v;
                        for_v += for_muv2[j] - for_muv1[j] * for_v + for_sigv[j]
* for_sqv * brow[5]; for_v += -for_sigv[j] * for_sqv * CorrMatrix[j][5][6] *
fx_sig[j]; for_v += -for_sigv[j] * for_sqv * ( CorrMatrix[j][5][3] * forzcvol1
+ CorrMatrix[j][5][4] * forzcvol2 ); for_v += +for_sigv[j] * for_sqv * (
CorrMatrix[j][5][0] * domzcvol1); for_v = max(0.0  , for_v);

                        if (EvalEvent[j])
                        {
                                dom_df = exp (ddomff_star[evtindex] +
domgam1_star[evtindex] * dom_ft1
                                                + domgam1_2_star[evtindex] *
dom_psi1);

                                // USELESS SINCE NO FX DIFFUSION
                                //for_df = exp (dforff_star[evtindex] +
forgam1_star[evtindex] * for_ft1 + forgam2_star[evtindex] * for_ft2
                                //		+ forgam1_2_star[evtindex] *
for_psi1 + forgam2_2_star[evtindex] * for_psi2 + forgam12_star[evtindex] *
for_psi12);

                                // QUANTOLGMSV Payoff
                                err = payoff_func(	i  ,
                                                                        dDate[j]
, dTime[j]  , func_parm_tab[j]  ,

                                                                        dom_ft1
, dom_psi1  ,

                                                                        for_ft1
, for_ft2  , for_psi1  , for_psi2  , for_psi12  , for_v  ,

                                                                        iNbProduct
, res_evt  , &stop_path);


                                if (err) goto FREE_RETURN;

                                for (k=0; k<iNbProduct; k++)
                                {
                                        path_payoff[k] += res_evt[k] / dom_df;
                                }

                                if (do_optimisation)
                                {
                                        for (k=0; k<params->iNbIndex; k++)
                                        {
                                                save_values[evtindex][k][i] =
res_evt[params->iColBound+k];
                                        }

                                        if (!params->iKnockInCol)
                                        {
                                                save_values[evtindex][params->iNbIndex][i]
= res_evt[params->iColPay] / dom_df;
                                        }
                                        else
                                        {
                                                save_values[evtindex][params->iNbIndex][i]
= res_evt[params->iColPay];
                                        }
                                }

                                evtindex++;
                        }

                        m += 7;
                }

                for (k=0; k<iNbProduct; k++)
                {
                        sum_payoff[k] += path_payoff[k] / iNumPaths;
                        sum_payoff2[k] += path_payoff[k] * path_payoff[k] /
iNumPaths;
                }

                if (do_optimisation && params->iKnockInCol)
                {

                        for (j=0; j<iNbEvent; j++)
                        {
                                if (optimise[j])
                                {
                                        save_values[j][params->iNbIndex][i] =
path_payoff[(int) (save_values[j][params->iNbIndex][i] + 0.5)];
                                }
                        }
                }
        }

        for (k=0; k<iNbProduct; k++)
        {
                res[k][0] = sum_payoff[k];
                res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) /
iNumPaths;

                if (res[k][1] > 0.0)
                {
                        res[k][1] = sqrt(res[k][1]);
                }
                else
                {
                        res[k][1] = 0.0;
                }
        }


        time1 = clock();
        smessage ("Phase 2 -convolution  , time in sec: %.2f"  , (double) (time1
- time2) / CLOCKS_PER_SEC);

        if (do_optimisation)
        {

                if (matrix) free_dmatrix (matrix  , 0  , iNumPaths - 1  , 0  , 7
* (iNbTime - 1) - 1); matrix = NULL;

                time1 = clock();

                err = find_and_optimise_boundary(	save_values  ,
                                                                                        iNbEvent  ,
                                                                                        iNumPaths  ,
                                                                                        optimise  ,
                                                                                        params  ,
                                                                                        opt_bar  ,
                                                                                        coef_lin  ,
                                                                                        &(res[iNbProduct][0])  ,
                                                                                        &(res[iNbProduct][1]));

                if (err) goto FREE_RETURN;

                for (j=0; j<iNbEvent; j++)
                {
                        res[j][2] = opt_bar[j];

                        for (k=0; k<params->iNbIndex; k++)
                        {
                                res[j][3+k] = coef_lin[j][k+1];
                        }
                }

                time2 = clock();
                smessage ("Phase 3 -optimisation  , time in sec: %.2f"  ,
(double) (time2 - time1) / CLOCKS_PER_SEC);
        }

FREE_RETURN:

        if (matrix) free_dmatrix (matrix  , 0  , iNumPaths - 1  , 0  , 7 *
(iNbTime - 1) - 1); if (res_evt) free_dvector (res_evt  , 0  , iNbProduct - 1);
        if (path_payoff) free_dvector (path_payoff  , 0  , iNbProduct - 1);
        if (sum_payoff) free_dvector (sum_payoff  , 0  , iNbProduct - 1);
        if (sum_payoff2) free_dvector (sum_payoff2  , 0  , iNbProduct - 1);

        if (dom_sigf1) free_dvector (dom_sigf1  , 0  , iNbTime - 1);
        if (dom_sigpsi1) free_dvector (dom_sigpsi1  , 0  , iNbTime - 1);

        if (for_sigf1) free_dvector (for_sigf1  , 0  , iNbTime - 1);
        if (for_sigf2) free_dvector (for_sigf2  , 0  , iNbTime - 1);
        if (for_sigpsi1) free_dvector (for_sigpsi1  , 0  , iNbTime - 1);
        if (for_sigpsi2) free_dvector (for_sigpsi2  , 0  , iNbTime - 1);
        if (for_sigpsi12) free_dvector (for_sigpsi12  , 0  , iNbTime - 1);
        if (for_muv1) free_dvector (for_muv1  , 0  , iNbTime - 1);
        if (for_muv2) free_dvector (for_muv2  , 0  , iNbTime - 1);
        if (for_sigv) free_dvector (for_sigv  , 0  , iNbTime - 1);

        if (fx_sig) free_dvector (fx_sig  , 0  , iNbTime - 1);

        if (brow) free_dvector (brow  , 0  , 6);

        if (cholmat) free_f3tensor(cholmat  , 0  , iNbTime - 1  , 0  , 6  , 0  ,
6);

        if (CorrMatrix_loc) free_dmatrix(CorrMatrix_loc  , 0  , 6  , 0  , 6);

        if (save_values) free_f3tensor(save_values  , 0  , iNbEvent - 1  , 0  ,
2  , 0  , iNumPaths - 1); if (opt_bar) free_dvector(opt_bar  , 0  , iNbEvent -
1); if (coef_lin) free_dmatrix(coef_lin  , 0  , iNbEvent - 1  , 0  ,
params->iNbIndex);


        return err;
}*/
