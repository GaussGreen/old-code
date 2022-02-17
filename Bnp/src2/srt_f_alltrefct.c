/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT  , Fixed Income 2020 Addins */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_ALLTREFCT.C                                     */
/*                                                                            */
/*      PURPOSE:        Functions useful to compute any trinomial tree        */
/*                                                                            */
/*      AUTHORS:        O. Van Eyseren 					      */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY:     		                                      */
/*                                                                            */
/*      DATE:           	                                              */
/*                                                                            */
/*      REASON:         		                                      */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_alltrefct.h"

#define bound(X, LUB, UPB)                                                     \
  {                                                                            \
    X = DMAX(X, LUB);                                                          \
    X = DMIN(X, UPB);                                                          \
  }

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :interp_linear
  AUTHOR          :O. Van Eyseren (from previous source code)
  DESCRIPTION     :interpolates a premium(x) linearily in x
  MODIFIES	  :
  CALL            :

<%%END---------------------------------------------------------------------*/

double interp_linear(double x, double x_lower, double prem_lower,
                     double x_upper, double prem_upper) {
  double coef, premium;

  if (x_upper != x_lower) {
    coef = (x - x_lower) / (x_upper - x_lower);
    premium = prem_lower + coef * (prem_upper - prem_lower);
  } else
    premium = prem_lower;

  return (premium);
}

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :find_closest
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :find guess such that val_array[index] is the closest to
target CALL            :

<%%END---------------------------------------------------------------------*/
void find_closest(double target, double *val_array, long min_index,
                  long max_index, long *guess) {
  long index;

  index = *guess;
  index = IMIN(max_index, index);
  index = IMAX(min_index, index);

  /* Find index just below */
  if (val_array[index] < target) {
    /* Move to the closest index strictly below the target */
    while (index < max_index && val_array[index + 1] < target)
      index++;
    /* Correct so that index really gives the closest one */
    if (index < max_index) {
      if ((val_array[index + 1] - target) < (target - val_array[index]))
        index++;
    }
  } else {
    /* Move to the next index strictly below the target */
    while (index > min_index && val_array[index] >= target)
      index--;
    /* Correct so that index really gives the closest one */
    if (index < max_index) {
      if ((val_array[index + 1] - target) < (target - val_array[index]))
        index++;
    }
  }
  *guess = index;
  return;
} /* END Err srt_f_find_closest(...) */

/* -----------------------------------------------------------------
   FUNCNAME        :find_closest_strictly_below
   AUTHOR          :O. Van Eyseren
   DESCRIPTION     :find guess such that val_array[index] is the closest to
target but always inferior to it CALL            :

<%%END---------------------------------------------------------------------*/
void find_closest_strictly_below(double target, double *val_array,
                                 long min_index, long max_index, long *guess) {
  long index;

  index = *guess;
  index = IMIN(max_index, index);
  index = IMAX(min_index, index);

  /* Find index just below */
  if (val_array[index] < target) {
    /* Move to the closest index strictly below the target */
    while (index < max_index && val_array[index + 1] < target)
      index++;
  } else {
    /* Move to the next index strictly below the target */
    while (index > min_index && val_array[index] >= target)
      index--;
  }
  *guess = index;
  return;
}

/* ---------------------------------------------------------------------------
 */

/** assume that we are probably increasing  , so we take as initial guess
    the value in node->mid_son_index.
    want index of element of prev_x closest to (drift_sam...STATEVAR);
    values in prev_x are increasing **/

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_trintrexindex
  AUTHOR          :E.Auld
  DESCRIPTION     :finds the index to center on in a trinomial tree in x
  MODIFIES	  :node->mid_son_index
  CALL            :

<%%END---------------------------------------------------------------------*/
/*---------------------------------------------------------------------------
  AMENDMENTS      :
  Reference       :
  Author          :
  Date            :
  Description     :
-----------------------------------------------------------------------------*/

/*<%%STADEC*/
void srt_f_trintrexindex(SrtStpPtr stp,  /* current step */
                         double *prev_x, /* state vars at previous step */
                         SrtTrinTreNdInf *node,  /* current node in tree */
                         SrtGrfnParam *grfparam, /* model parameter */
                         int closestmeth /* CLOSESTINX or CLOSESTINLNX */
                         )
/*<%%ENDDEC*/
{
  SrtTreStpInf *trinf;
  trinf = stp->next->trinf;

  /* Returns the index that correspond to the closest smaller value */

  find_closest(sam_get(node->drift_sam, 0, STATEVAR), prev_x,
               trinf->min_x_index, trinf->max_x_index, &node->mid_son_index);

  /* if index just above is closer  , increment */
  /*
    DO NOT NEED THAT...
    trinf = stp->trinf;
  */

  if (closestmeth == CLOSESTINX) {
    if (node->mid_son_index < trinf->max_x_index &&
        0.5 * (prev_x[node->mid_son_index] + prev_x[node->mid_son_index + 1]) <
            sam_get(node->drift_sam, 0, STATEVAR))
      node->mid_son_index++;
  } else {
    if (node->mid_son_index < trinf->max_x_index &&
        ((prev_x[node->mid_son_index] * prev_x[node->mid_son_index + 1]) <
         (sam_get(node->drift_sam, 0, STATEVAR) *
          sam_get(node->drift_sam, 0, STATEVAR))))
      node->mid_son_index++;
  }
  return;
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        : srt_f_trintredcntvec
  AUTHOR          : O. Van Eyseren from E.Auld code
  DESCRIPTION     : discounts assets to cur (cur is prob. weighted sum of three
                        vectors in node.)
  MODIFIES        :cur
  CALL            :

<%%END>---------------------------------------------------------------------*/
/*<%%STADEC>*/
void srt_f_trintredcntvec(SrtStpPtr stp, double *cur, SrtTrinTreNdInf *node,
                          double **assets, long dim)
/*<%%ENDDEC>*/
{
  SrtTreStpInf *trinf;
  int i, j, x_i;
  double tmp;

  trinf = stp->next->trinf;
  x_i = node->mid_son_index;

  if (x_i == trinf->min_x_index)
    x_i++;
  if (x_i == trinf->max_x_index)
    x_i--;

  for (i = 0; i < dim; i++) {
    tmp = 0;
    for (j = 0; j < 3; j++) {
      tmp += node->p[j] * assets[x_i - 1 + j][i];
    }
    cur[i] = tmp * node->df;
  }
  return;
}

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_trintredeltax
  AUTHOR          :O. Van Eyseren (from previous source code by E. Auld)
  DESCRIPTION     :calculates the delta_x step in a trinomial tree
                   this function assumes the the corresponding underlying
                   (specified by the und_index) is of the BS type)
  MODIFIES	  :
  CALL            :
  INPUT		  : und_index corresponds to the index of the underlying
                    used (in the order of the **TmInf)

<%%END---------------------------------------------------------------------*/
void srt_f_trintredeltax(SrtStpPtr top, double *delta_x, int und_index) {
  double tmp, max_vol_dt = 0, min_vol_dt, cum_var = 0;
  SrtStpPtr cur;
  SrtLogTmInf *tminf;

  cur = top;

  tminf = (SrtLogTmInf *)(cur->tminf[und_index]);
  min_vol_dt = tminf->int_sig2_dt;
  while (cur->next != NULL) {
    tminf = (SrtLogTmInf *)cur->tminf[und_index];
    tmp = tminf->int_sig2_dt;
    cum_var += tmp;
    max_vol_dt = DMAX(max_vol_dt, tmp);
    min_vol_dt = DMIN(min_vol_dt, tmp);
    cur = cur->next;
  }
  tmp = sqrt(cum_var / (double)cur->index) * DELTAOPTRATIO;
  /*** the order of the following two tests is vital.  It was previously in the
  opposite order and we considered that to be a bug E.Auld 17 Sep94 ***/
  max_vol_dt = sqrt(max_vol_dt);
  min_vol_dt = sqrt(min_vol_dt);

  if (tmp > DELTAMAXBOUND * min_vol_dt)
    tmp = DELTAMAXBOUND * min_vol_dt;
  if (tmp < DELTAMINBOUND * max_vol_dt)
    tmp = DELTAMINBOUND * max_vol_dt;

  *delta_x = tmp;
}

/** calculate probabilities **/
/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_trintreprob
  AUTHOR          :O. Van Eyseren (from E. Auld code)
  DESCRIPTION     :calculate probabilities  , assuming that node->drift_sam and
    node->var_at_sam are correctly set  , as well as node->mid_son_index
    If the mid_son_index hits the top (bottom)  , it is decremented
    (incremented).  Probabilities are bounded to be between PROBUPB and
    PROBLUB. They are printed if PRINFO is true.
  MODIFIES	  :node->p
  CALL            :

<%%END---------------------------------------------------------------------*/
/*---------------------------------------------------------------------------
  AMENDMENTS      :
  Reference       :
  Author          :
  Date            :
  Description     :
-----------------------------------------------------------------------------*/
/*<%%STA*/
static double PROBUPB = 2.0;
static double PROBLUB = -1.0;
static long PRINFO = 0;
/*<%%END*/

/*<%%STADEC>*/
void srt_f_trintreprob(double *prev_x, SrtTrinTreNdInf *node,
                       SrtGrfnParam *grfnparam)
/*<%%ENDDEC>*/
{
  double eta, up_x, center_x, down_x;
  double fwd = sam_get(node->drift_sam, 0, STATEVAR);
  double second_moment = node->var_at_sam + fwd * fwd;
  int bndflg = 0;

  long index = node->mid_son_index;

  up_x = prev_x[index + 1];
  center_x = prev_x[index];
  down_x = prev_x[index - 1];

  bndflg = 0;
  eta = fwd - center_x;
  if (fabs(up_x - center_x) < 1.0e-9 || fabs(down_x - center_x) < 1.0e-9 ||
      down_x > 1.0e10) {
    node->p[0] = node->p[2] = node->p[1] = 1.0 / 3.0;
  } else {
    node->p[2] =
        second_moment - center_x * center_x - eta * (center_x + down_x);
    node->p[2] /= ((up_x - center_x) * (up_x - down_x));
    node->p[0] = ((up_x - center_x) * node->p[2] - eta) / (center_x - down_x);
    node->p[1] = 1 - node->p[2] - node->p[0];
  }
  if (bndflg) {
    if (node->p[2] > 1 || node->p[2] < 0 || node->p[1] > 1 || node->p[1] < 0 ||
        node->p[0] > 1 || node->p[0] < 0)
      bound(node->p[0], PROBLUB, PROBUPB);
    bound(node->p[1], PROBLUB, PROBUPB);
    bound(node->p[2], PROBLUB, PROBUPB);
  }
  return;
}

/*-------------------------------------------------------------------------
  FUNCNAME        :srt_f_nonotreprob
  AUTHOR          :O. Van Eyseren (from E. Auld code)
  DESCRIPTION     :calculates probabilities for a nononomial tree  ,
                                   assuming that node -> moments are correctly
set  , as well as node -> mid_son_index MODIFIES		  :node ->
p[0..2][0..2]
---------------------------------------------------------------------------*/

void srt_f_nonotreprob(SrtTwoFacTreInf *trinf, SrtNonoTreeNodeInfo *node,
                       double **next_x) {
  double m, n, mn, m2, n2, m2n, mn2, m2n2, rho, e_x1, e_x2, e_x1x1, e_x2x2,
      e_x1x2, midx1, midx2;
  double e_x1x1x2, e_x1x2x2, e_x1x1x2x2;
  SrtTwoFacMmtStruct *moments;

  /* spacing */
  m = trinf->u[0];
  n = trinf->u[1];
  m2 = m * m;
  n2 = n * n;
  mn = m * n;
  m2n = m2 * n;
  mn2 = m * n2;
  m2n2 = m2 * n2;

  /* mid son node */
  midx1 = next_x[0][node->mid_son_index[0]];
  midx2 = next_x[1][node->mid_son_index[1]];

  /* moments of xi = Xi - mid son node */
  moments = &(node->moments);
  e_x1 = moments->E_X1 - midx1;
  e_x2 = moments->E_X2 - midx2;
  e_x1x1 = moments->E_X1X1 + midx1 * (midx1 - 2.0 * moments->E_X1);
  e_x2x2 = moments->E_X2X2 + midx2 * (midx2 - 2.0 * moments->E_X2);
  e_x1x2 = moments->E_X1X2 - midx2 * moments->E_X1 - midx1 * moments->E_X2 +
           midx1 * midx2;

  rho = (e_x1x2 - e_x1 * e_x2) /
        sqrt((e_x1x1 - e_x1 * e_x1) * (e_x2x2 - e_x2 * e_x2));

  e_x1x1x2 = e_x1 * e_x1 * e_x2 + 2 * e_x1 * (e_x1x2 - e_x1 * e_x2) +
             e_x2 * (e_x1x1 - e_x1 * e_x1);
  e_x1x2x2 = e_x2 * e_x2 * e_x1 + 2 * e_x2 * (e_x1x2 - e_x1 * e_x2) +
             e_x1 * (e_x2x2 - e_x2 * e_x2);
  e_x1x1x2x2 =
      e_x1 * e_x1 * e_x2 * e_x2 +
      (1.0 + 2 * rho * rho) * (e_x1x1 - e_x1 * e_x1) * (e_x2x2 - e_x2 * e_x2) +
      e_x1 * e_x1 * (e_x2x2 - e_x2 * e_x2) +
      e_x2 * e_x2 * (e_x1x1 - e_x1 * e_x1) +
      4 * e_x1 * e_x2 * (e_x1x2 - e_x1 * e_x2);

  /* compute probabilities */
  node->p[0][0] =
      e_x1x2 / mn +
      (-(e_x1x2x2 * mn) - e_x1x2 * mn2) / (2 * mn * mn2) + /* Pdd */
      (-(e_x1x1x2 * mn) - e_x1x2 * m2n) / (2 * mn * m2n) -
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (16 * mn * mn * mn * mn2 * m2n * m2n2);
  node->p[0][1] =
      -(e_x1x2 / mn) -
      (-(e_x1x2x2 * mn) - e_x1x2 * mn2) / /* Pd0 */
          (mn * mn2) +
      e_x1x1 / m2 + (-(e_x1x1 * m) - e_x1 * m2) / (2 * m * m2) -
      (-(e_x1x1x2 * mn) - e_x1x2 * m2n) / (2 * mn * m2n) +
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (8 * mn * mn * mn * mn2 * m2n * m2n2);
  node->p[0][2] =
      (-(e_x1x2x2 * mn) - e_x1x2 * mn2) / (2 * mn * mn2) - /* Pdu */
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (16 * mn * mn * mn * mn2 * m2n * m2n2);
  node->p[1][0] =
      -(e_x1x2 / mn) -
      (-(e_x1x2x2 * mn) - e_x1x2 * mn2) / (2 * mn * mn2) - /* P0d */
      (-(e_x1x1x2 * mn) - e_x1x2 * m2n) / (mn * m2n) +
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (8 * mn * mn * mn * mn2 * m2n * m2n2) +
      e_x2x2 / n2 + (-(e_x2x2 * n) - e_x2 * n2) / (2 * n * n2);
  node->p[1][1] =
      1 + e_x1x2 / mn +
      (-(e_x1x2x2 * mn) - e_x1x2 * mn2) / (mn * mn2) - /* P00 */
      e_x1x1 / m2 + (-(e_x1x1x2 * mn) - e_x1x2 * m2n) / (mn * m2n) -
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (4 * mn * mn * mn * mn2 * m2n * m2n2) -
      e_x2x2 / n2;
  node->p[1][2] =
      -(-(e_x1x2x2 * mn) - e_x1x2 * mn2) / (2 * mn * mn2) + /* P0u */
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (8 * mn * mn * mn * mn2 * m2n * m2n2) -
      (-(e_x2x2 * n) - e_x2 * n2) / (2 * n * n2);
  node->p[2][0] =
      (-(e_x1x1x2 * mn) - e_x1x2 * m2n) / (2 * mn * m2n) - /* Pud */
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (16 * mn * mn * mn * mn2 * m2n * m2n2);
  node->p[2][1] =
      -(-(e_x1x1 * m) - e_x1 * m2) / (2 * m * m2) - /* Pu0 */
      (-(e_x1x1x2 * mn) - e_x1x2 * m2n) / (2 * mn * m2n) +
      (4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 -
       2 * mn * m2n *
           (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
            2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
          (8 * mn * mn * mn * mn2 * m2n * m2n2);
  node->p[2][2] =
      -(4 * mn * mn * mn2 * (-(e_x1x1x2 * mn) - e_x1x2 * m2n) * m2n2 - /* Puu */
        2 * mn * m2n *
            (-2 * mn * (-(e_x1x2x2 * mn) - e_x1x2 * mn2) * m2n2 -
             2 * mn * mn2 * (-(e_x1x1x2x2 * mn) + e_x1x2 * m2n2))) /
      (16 * mn * mn * mn * mn2 * m2n * m2n2);
}

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_nonotregenx
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :generates the nonotree grid at the following discrete
                   date  , according to the current step parameters
                   This function is as generic as to be used in any nonotree  ,
                   whatever the underlyings  , as long as SrtTwoFacTreInf
                   has been properly initialised in the SrtStpPtr
  MODIFIES	  : **cur_x
  CALL            :
<%%END---------------------------------------------------------------------*/

/* generates the matrix of state variables at some step stp into cur_x*/
void srt_f_nonotregenx(SrtStpPtr stp, double **cur_x)

{
  int i, j;
  double delta_x[2];
  SrtTwoFacTreInf *trinf = stp->trinf;

  /* basically if it is first step
  if (trinf -> max_x_index[0] == trinf -> min_x_index[0])
  {
          cur_x[0][trinf->min_x_index[0]] = cur_x[1][trinf->min_x_index[1]] =
  0.0; return;
  }
  */
  /** generate X (the state variable) by adding delta_x over delta_x **/
  for (i = 0; i < 2; i++) {
    delta_x[i] = trinf->u[i];
    cur_x[i][trinf->min_x_index[i]] = trinf->xmin[i];
    for (j = trinf->min_x_index[i] + 1; j <= trinf->max_x_index[i]; j++)
      cur_x[i][j] = cur_x[i][j - 1] + delta_x[i];
  }
}

/*------------------------------------------------------------------------

  FUNCNAME        :srt_f_nonotredcnvec
  AUTHOR          :O. Van Eyseren

--------------------------------------------------------------------------*/

/* compute ddiscoutned expectation on some node */
void srt_f_nonotredcntvec(SrtStpPtr stp, double *cur, SrtNonoTreeNodeInfo *node,
                          double ***assets, long dim)

{
  SrtTwoFacTreInf *trinf;
  int i, j, k, x_i[2];
  double tmp;

  trinf = stp->next->trinf;
  /* Copy the (corrected if at limit) mid son index into a simple vector */
  for (i = 0; i < 2; i++) {
    x_i[i] = node->mid_son_index[i];
  }

  for (i = 0; i < dim; i++) {
    tmp = 0.0;

    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++) {
        if ((x_i[0] + j - 1 <= trinf->max_x_index[0]) &&
            (x_i[0] + j - 1 >= trinf->min_x_index[0]) &&
            (x_i[1] + k - 1 <= trinf->max_x_index[1]) &&
            (x_i[1] + k - 1 >= trinf->min_x_index[1]))
        /* standard case: point is inside tree boundaries */
        {
          tmp += node->p[j][k] * assets[x_i[0] - 1 + j][x_i[1] - 1 + k][i];
        }

        else
        /* point is outside boundaries */
        {
          tmp += node->p[j][k] * assets[x_i[0]][x_i[1]][i];
        }
      }
    }

    cur[i] = tmp * node->df;
  }
}