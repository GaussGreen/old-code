/*--------------------------------------------------------------
        FILE: Product.cxx
        PURPOSE: Generic product routines
        AUTHOR: Dimitri Mayevski
        DATE: 21/11/2002
  --------------------------------------------------------------*/

#include "Product.h"
#include "grf_h_all.h"
#include "srt_h_all.h"

static Err GRFN_Payoff(SProductDesc *product, int idx, double time,
                       SMktData *mkt_data, double *pv) {
  Err err = NULL;
  FIRSTAllMkts *comm = (FIRSTAllMkts *)product->spec_desc;
  FIRSTMktAtT *local = comm->evt + idx;
  int i;
  double tmp;

  for (i = 0; i < product->nccy; i++)
    memcpy(local->evt->df[i], mkt_data->dfs[i],
           product->nmat[i][idx] * sizeof(double));

  if (mkt_data->type == MKT_GRFN && mkt_data->spec_desc)
    memcpy(&local->smp, mkt_data->spec_desc, sizeof(SrtSample));

  err = FIRSTEvalEvent(comm, local, product->ninst, 2, NULL, NULL, pv, &tmp);
  return err;
}

Err ProductDesc_InitGRFN(SProductDesc *product, char *und, int nevd, long *evd,
                         long ntabrows, long ntabcols, char ***tabstrs,
                         int **tabmask, long auxwidth, long *auxlen,
                         double **aux) {
  Err err = NULL;
  FIRSTAllMkts *comm = NULL;
  SrtGrfnParam defParm;
  int i, j, k, forback;

  err = srt_f_set_default_GrfnParams(&defParm);
  if (err)
    return err;

  comm = (FIRSTAllMkts *)malloc(sizeof(FIRSTAllMkts));
  if (!comm)
    return serror("Memory failure");

  err =
      FIRSTInitMktStruct(nevd, evd, ntabrows, ntabcols, tabstrs, tabmask,
                         auxwidth, auxlen, aux, und, &defParm, &forback, comm);
  if (err) {
    free(comm);
    return err;
  }

  product->type = PRODUCT_GRFN;
  product->spec_desc = comm;

  product->nex = comm->num_evt;
  product->nccy = comm->num_und;
  product->ninst = comm->num_cols;

  product->ex = comm->tms;
  product->ex_d = comm->dts;
  product->am = comm->am;

  if (product->ex[0] < 1e-5)
    return serror("Some event dates are in the past");

  product->nmat = (int **)calloc(product->nccy, sizeof(int *));
  product->mat = (double ***)calloc(product->nccy, sizeof(double **));
  product->mat_d = (long ***)calloc(product->nccy, sizeof(long **));
  if (!product->nmat || !product->mat || !product->mat_d)
    return serror("Memory failure");

  for (i = 0; i < product->nccy; i++) {
    product->nmat[i] = (int *)calloc(product->nex, sizeof(int));
    product->mat[i] = (double **)calloc(product->nex, sizeof(double *));
    product->mat_d[i] = (long **)calloc(product->nex, sizeof(long *));
    if (!product->nmat[i] || !product->mat[i] || !product->mat_d[i])
      return serror("Memory failure");

    for (j = 0; j < product->nex; j++) {
      product->nmat[i][j] = comm->evt[j].evt->dflen[i];
      product->mat[i][j] =
          (double *)calloc(product->nmat[i][j], sizeof(double));
      product->mat_d[i][j] = (long *)calloc(product->nmat[i][j], sizeof(long));
      if (!product->mat[i][j] || !product->mat_d[i][j])
        return serror("Memory failure");

      for (k = 0; k < product->nmat[i][j]; k++) {
        product->mat[i][j][k] = product->ex[j] + comm->evt[j].evt->dft[i][k];
        product->mat_d[i][j][k] = comm->evt[j].evt->dfd[i][k];
      }
    }
  }
  product->Payoff = GRFN_Payoff;

  return NULL;
}

Err ProductDesc_FreeGRFN(SProductDesc *product) {
  int i, j;
  FIRSTAllMkts *comm = (FIRSTAllMkts *)product->spec_desc;

  if (product->nmat)
    for (i = 0; i < product->nccy; i++)
      free(product->nmat[i]);
  free(product->nmat);

  if (product->mat)
    for (i = 0; i < product->nccy; i++) {
      if (product->mat[i])
        for (j = 0; j < product->nex; j++)
          free(product->mat[i][j]);
      free(product->mat[i]);
    }
  free(product->mat);

  if (product->mat_d)
    for (i = 0; i < product->nccy; i++) {
      if (product->mat_d[i])
        for (j = 0; j < product->nex; j++)
          free(product->mat_d[i][j]);
      free(product->mat_d[i]);
    }
  free(product->mat_d);

  if (comm)
    FIRSTFreeMktStruct(comm);
  free(comm);

  return NULL;
}

static Err Swaption_Payoff(SProductDesc *product, int idx, double time,
                           SMktData *mkt_data, double *pv) {
  int k;
  double iv = 0.0, **cpn = (double **)product->spec_desc;

  for (k = 0; k < product->nmat[0][idx]; k++)
    iv += mkt_data->dfs[0][k] * cpn[idx][k];

  pv[idx] = (iv > 0.0 ? iv : 0.0);
  return NULL;
}

Err ProductDesc_InitSwaptions(SProductDesc *product, SCashFlows *g) {
  product->nex = product->ninst = g->nex;
  product->nccy = 1;
  product->ex = g->ex;
  product->am = NULL;
  product->nmat = &g->nmat;
  product->mat = &g->mat;
  product->spec_desc = g->cpn;
  product->type = PRODUCT_SWAPTIONS;
  product->Payoff = Swaption_Payoff;
  return NULL;
}
