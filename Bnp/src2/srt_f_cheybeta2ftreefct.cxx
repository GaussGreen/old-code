/* ----------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & LLWVE

   DATE : MAI 98

   FILENAME:    srt_f_cheybeta2dtreefct.cxx

   PURPOSE:  utility functions for the new implementation  of BFL algo
                      for the Cheyette beta  2D in GRFN

   MODIFICATION:

   ----------------------------------------------------------------------------------
 */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheybeta2fdynamics.h"
#include "srt_h_cheybeta2ftreefct.h"
#include "srt_h_mc_core.h"
#include "srt_h_mc_evolve.h"
#include "srt_h_stpcorata.h"
#include <num_h_gaussj.h"

#define EPS_PHI 0.0000001

/* -------------------------------------------------------------------------------------
   FUNCTION: Err MCTree_CHEYBETA_2f_Euler_evolve(...)
   PURPOSE:  A local function to evolve the Cheyette Beta 2f from random phi's
   -------------------------------------------------------------------------------------
 */

static Err MCTree_CHEYBETA_2f_evolve_in_Grid(double **rndm_mat, SrtStpPtr step,
                                             MCTreePoint **grid,
                                             long path_index)

{
  Err err = NULL;
  SrtStpPtr top;
  SrtIRMTmInf *time_info;
  long t;
  SrtSample sample;

  long n = path_index;

  /* Start from the first time step */
  step = top = gototop(step);

  if (step != NULL) {
    /* First time step */
    t = 0;

    /* Gets the relevant time information attached to the step */
    time_info = (SrtIRMTmInf *)step->tminf[0];

    /* Initialization for the first time step (today) */
    grid[0][n].x1 = sam_get(sample, 0, X1) = 0.0;
    grid[0][n].x2 = sam_get(sample, 0, X2) = 0.0;
    grid[0][n].fw1 = 0;
    grid[0][n].fw2 = 0;
    grid[0][n].jump_x1 = 0.0;
    grid[0][n].jump_x2 = 0.0;

    grid[0][n].phi1 = sam_get(sample, 0, PHI1) = 0.0;
    grid[0][n].phi2 = sam_get(sample, 0, PHI2) = 0.0;
    grid[0][n].phi1_2 = sam_get(sample, 0, CROSSPHI) = 0.0;

    grid[0][n].short_rate = sam_get(sample, 0, SHORT_RATE) =
        sam_get(time_info->fwd_sam, 0, F_0_t);
    grid[0][n].df = exp(-grid[0][n].short_rate * step->delta_t);

    /* Loops on all the time steps */
    while (step->next != NULL) {

      /* Increment the step index by one for the grid */
      t++;

      /* The jump corresponds to the Phi for the previous time step ... */
      grid[t][n].jump_x1 =
          (sam_get(sample, 0, PHI1) + sam_get(sample, 0, CROSSPHI)) *
          step->delta_t;
      grid[t][n].jump_x2 =
          (sam_get(sample, 0, PHI2) + sam_get(sample, 0, CROSSPHI)) *
          step->delta_t;

      /* Evolve the sample for one time step ( underlying index is 0: domestic
       * in two factor...) */
      err = srt_f_CheyBeta2f_evolve_sample(step, &sample, 0, rndm_mat[0][t],
                                           rndm_mat[1][t], &sample);
      if (err)
        return err;

      /* Efectively moves to the next time step (for time info...) */
      step = step->next;
      time_info = (SrtIRMTmInf *)step->tminf[0];

      /* Stores the variables without the Non Markovian part */
      grid[t][n].x1 = sam_get(sample, 0, X1) - grid[t][n].jump_x1;
      grid[t][n].x2 = sam_get(sample, 0, X2) - grid[t][n].jump_x2;
      grid[t][n].fw1 = (1.0 - time_info->ev.twof[0].lambda * step->delta_t) *
                       sam_get(sample, 0, X1);
      grid[t][n].fw2 = (1.0 - time_info->ev.twof[1].lambda * step->delta_t) *
                       sam_get(sample, 0, X2);

      grid[t][n].phi1 = sam_get(sample, 0, PHI1);
      grid[t][n].phi2 = sam_get(sample, 0, PHI2);
      grid[t][n].phi1_2 = sam_get(sample, 0, CROSSPHI);

      grid[t][n].short_rate = sam_get(sample, 0, SHORT_RATE);
      grid[t][n].df = exp(-grid[t][n].short_rate * step->delta_t);

    } /* END while step loop */

  } /* end if != NULL */
  else {
    return "step==NULL at top of LIST";
  }

  return NULL;

} /* END MCTree_CHEYBETA_2f_evolve_in_Grid(...) */

/* ----------------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------------
   FUNCTION:     Err srt_f_cheybeta2d_make_grid()
   PURPOSE:     Generate the McGrid
   ---------------------------------------------------------------------------------
 */

Err srt_f_cheybeta2ftree_make_grid(SrtGrfnParam *grfnparam, SrtStpPtr step,
                                   SrtUndInfo *und_info, MCTreePoint **grid) {
  Err err = NULL;
  SrtSample *sam = NULL;
  SrtUndPtr und = NULL;
  long last_index;
  SrtStpPtr top;
  long n, num_time_pts;

  long seed;
  SrtMCRanSrc mcrs;

  /* -------------------- STEP 1: INITIALISATION --------------------- */

  /* Go to the first time step and store the last time point index */
  last_index = create_index(step);
  top = gototop(step);
  num_time_pts = create_index(top) + 1;

  /*	Initialises mcrs according to inputs        , makes memory allocation or
                  even generates fully the non correlated random numbers */
  seed = grfnparam->rand_seed;
  err = SrtMCRanSrc_start(&mcrs, 0, grfnparam->num_MCarlo_paths - 1, 1,
                          last_index, und_info->no_of_brownians,
                          grfnparam->sample_type, seed, top);
  if (err)
    return err;

  /* ----------------- STEP 2: PATHS SIMULATION ------------------------ */

  /* Loop on all the paths used to generate the initial draw */
  for (n = 0; n < grfnparam->num_MCarlo_paths; n++) {

    /* Generates or gets the next MC sample path in mcrs.r[mcrs.cxxur_path] */
    err = SrtMCRanSrc_next(&mcrs);
    if (err) {
      SrtMCRanSrc_end(&mcrs);
      return err;
    }

    /* Generate the current path for the Chey Beta 2f underlying in the Grid */
    err = MCTree_CHEYBETA_2f_evolve_in_Grid(mcrs.r[mcrs.cxxur_path], top, grid,
                                            n);
    if (err) {
      SrtMCRanSrc_end(&mcrs);
      return err;
    }

  } /* END of loop on all paths */

  /* Frees the random number for the MC simulation */
  SrtMCRanSrc_end(&mcrs);

  /* Return a success message */
  return NULL;

} /* END of Err  srt_f_cheybeta2d_make_grid() */

#undef EPS_PHI

/* ---------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------
   PURPOSE:  Allocation of structures for points and probabilities	of
   transition
   ------------------------------------------------------------------------------------------------
 */

#define BOUND_FACTOR 30.0
#define LAMBDA sqrt(3.0)
#define EPSICFL 0.5

static double *proba, **nodes, **realproba, **coeffs;
static double thmom[NB_NEIGHBOURS] = {1, 0, 0, 1, 1};
static TriPoint ptx[NB_NEIGHBOURS];
static double fik[NB_NEIGHBOURS][50];
static int index_pt[NB_NEIGHBOURS];

Err srt_f_alloc_computing_resources() {
  proba = dvector(1, NB_NEIGHBOURS);
  nodes = dmatrix(1, NB_NEIGHBOURS, 1, 2);
  realproba = dmatrix(1, NB_NEIGHBOURS, 1, 1);
  coeffs = dmatrix(1, NB_NEIGHBOURS, 1, NB_NEIGHBOURS);

  return NULL;
}

Err srt_f_free_computing_resources() {
  free_dvector(proba, 1, NB_NEIGHBOURS);
  free_dmatrix(nodes, 1, NB_NEIGHBOURS, 1, 2);
  free_dmatrix(realproba, 1, NB_NEIGHBOURS, 1, 1);
  free_dmatrix(coeffs, 1, NB_NEIGHBOURS, 1, NB_NEIGHBOURS);

  return NULL;
}

/*------------------------------------------------------------------------------------------------
   FUNCTION:  void transfer_info_nodes(...)
   PURPOSE:  Transfer nodes and attributes for computing the triangulation
------------------------------------------------------------------------------------------------*/
void transfer_info_nodes(MCTreePoint *mesh, double **next_assets,
                         SrtTriangulationIO *triang) {
  int n, k;
  double x, y, xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0, dx, dy;

  /* load the coordinates of points in the stucture input of triangulation */
  for (n = 0; n < triang->numberofpoints - NB_BOUND_PTS; n++) {
    x = triang->pointlist[2 * n] = mesh[n].x1;
    y = triang->pointlist[2 * n + 1] = mesh[n].x2;

    /* Determine the smallest and largest x and y coordinates. */
    if (n == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
    } else {
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
    }
  }

  /* add boundaries points in the stucture input of triangulation */
  n = triang->numberofpoints - NB_BOUND_PTS;
  dx = BOUND_FACTOR * (xmax - xmin);
  dy = BOUND_FACTOR * (ymax - ymin);
  triang->pointlist[2 * n] = xmin - dx;
  triang->pointlist[2 * n + 1] = ymin - dy;
  n++;
  triang->pointlist[2 * n] = xmin - dx;
  triang->pointlist[2 * n + 1] = ymax + dy;
  n++;
  triang->pointlist[2 * n] = xmax + dx;
  triang->pointlist[2 * n + 1] = ymin - dy;
  n++;
  triang->pointlist[2 * n] = xmax + dx;
  triang->pointlist[2 * n + 1] = ymax + dy;

  /* load the associates values for value function */
  for (n = 0; n < triang->numberofpoints; n++)
    for (k = 0; k < triang->numberofpointattributes; k++)
      triang->pointattributelist[n * triang->numberofpointattributes + k] =
          next_assets[n][k];
}

/* -------------------------------------------------------------------------------------------
   FUNCTION:  Err srt_f_cheybeta2ftree_compute_proba(...)
   PURPOSE:   compute the points and probas	of transition
 -------------------------------------------------------------------------------------------
 */
Err srt_f_cheybeta2ftree_compute_proba(SrtStpPtr step) {
  double l2 = LAMBDA * LAMBDA, sqdt = sqrt(step->delta_t);
  SrtIRMTmInf *info = (SrtIRMTmInf *)step->tminf[0];
  double v1 = sqrt(1.0 + info->correl_x), v2 = sqrt(1.0 - info->correl_x);
  double x1 = LAMBDA * sqdt * v1, x2 = LAMBDA * sqdt * v2;
  double a = 1.0 / sqrt(2.0);

  /* transition probabilities for the standard gaussian law to specified points
   */
  proba[1] = 1.0 - 2.0 / l2;
  proba[2] = 1.0 / (2.0 * l2);
  proba[3] = proba[2];
  proba[4] = proba[2];
  proba[5] = proba[2];

  /* the specified points */
  nodes[1][1] = 0.0;
  nodes[1][2] = 0.0;
  nodes[2][1] = a * x1;
  nodes[2][2] = a * x1;
  nodes[3][1] = -a * x2;
  nodes[3][2] = a * x2;
  nodes[4][1] = -a * x1;
  nodes[4][2] = -a * x1;
  nodes[5][1] = a * x2;
  nodes[5][2] = -a * x2;

  return NULL;
} /* END of Err srt_f_cheybeta2ftree_compute_A() */

/* -------------------------------------------------------------------------------------------
   FUNCTION:  Err rescal(...)
   PURPOSE:
 -------------------------------------------------------------------------------------------
 */
Err rescal(TriPoint pcur, double fw1, double fw2, double sig1, double sig2,
           int i, SrtStpPtr step) {
  double sqdt = sqrt(step->delta_t);
  SrtIRMTmInf *info = (SrtIRMTmInf *)step->tminf[0];
  double v1 = sqrt(1.0 + info->correl_x) * sqdt,
         v2 = sqrt(1.0 - info->correl_x) * sqdt;
  double a = 1.0 / sqrt(2.0);
  double xr = (pcur[0] - fw1) / sig1, yr = (pcur[1] - fw2) / sig2;
  double x = a * (xr + yr) / v1, y = a * (-xr + yr) / v2;

  /* build the matrix for transition probabilities for the standard gaussian law
   * to specified points */
  coeffs[1][i] = 1.0;
  coeffs[2][i] = x;
  coeffs[3][i] = y;
  coeffs[4][i] = x * x;
  coeffs[5][i] = y * y;

  ptx[i - 1] = pcur;

  realproba[i][1] = thmom[i - 1];

  return NULL;
} /* END of Err srt_f_cheybeta2ftree_compute_A() */

/* small routines useful for quadrature algoritm and triangulation */

static int plus1mod3[3] = {1, 2, 0};
static int minus1mod3[3] = {2, 0, 1};

#define coords_bary(x, y, a1, b1, a2, b2, a3, b3)                              \
  ((x - a2) * (b3 - b2) - (y - b2) * (a3 - a2)) /                              \
      ((a1 - a2) * (b3 - b2) - (b1 - b2) * (a3 - a2))

int not_in_set(int i, TriPoint p) {
  int k, bb = 0;
  for (k = 0; k < i - 1; k++)
    bb += (p == ptx[k]);
  return (bb == 0);
}

/* ----------------------------------------------------------------------------------------
   FUNCTION: Err srt_f_cheybeta2ftree_discount_to_node
   PURPOSE:  for one point        , fills and inverts the 6x6 symmetric matrix
   ----------------------------------------------------------------------------------------
 */
Err srt_f_cheybeta2ftree_discount_to_node(
    MCTreePoint **grid, int index_t, int index_n, SrtStpPtr step,
    double *cur_asset, /* Current assets at node */
    SrtGrfnParam *grfnparam) {
  Err err = NULL;
  int i, k, l, bb, chiure = 0;
  double sigma1, sigma2, alpha;
  SrtIRMTmInf *info = (SrtIRMTmInf *)step->tminf[0];
  TriLocateResultType location;
  TriOrientedTriangle searchtri;
  TriTriangle ptr;
  TriPoint p1, p2, p3, pt[3];
  double pt_cur[2];
  double l1, l2, l3, fxk, f1, f2, f3;

  /* The volatility terms in each direction */
  alpha = pow(fabs(grid[index_t][index_n].short_rate), info->ev.twof[0].beta);
  sigma1 = info->ev.twof[0].sig * alpha;
  sigma2 = info->ev.twof[1].sig * alpha;

  /* Find a boundary triangle. */
  searchtri.tri = getdummytri();
  searchtri.orient = 0;
  symself(searchtri);

  for (i = 1; i <= NB_NEIGHBOURS; i++)
    index_pt[i - 1] = 0;

  /* Compute the transition probas and conditionnal expectations to the current
   * node */

  /* Loop on the nodes  to which cur_node connects at next step */
  for (i = 1; i <= NB_NEIGHBOURS; i++) {
    /* Locate the position of node for quadrature formula */

    pt_cur[0] = grid[index_t][index_n].fw1 + sigma1 * nodes[i][1];
    pt_cur[1] = grid[index_t][index_n].fw2 + sigma2 * nodes[i][2];

    location = locate(pt_cur, &searchtri);

    org(searchtri, p1);
    dest(searchtri, p2);
    apex(searchtri, p3);
    l1 = coords_bary(pt_cur[0], pt_cur[1], p1[0], p1[1], p2[0], p2[1], p3[0],
                     p3[1]);
    l2 = coords_bary(pt_cur[0], pt_cur[1], p2[0], p2[1], p3[0], p3[1], p1[0],
                     p1[1]);
    l3 = 1.0 - l1 - l2;

    /* Determine the closest point and transfer its rescaled coordinates in the
       weight matrix as well as the value fonction */
    if (l1 >= l2) {
      if (l1 >= l3) {
        pt[0] = p1;
        if (l2 >= l3) {
          pt[1] = p2;
          pt[2] = p3;
        } else {
          pt[1] = p3;
          pt[2] = p2;
        }
      } else {
        pt[0] = p3;
        pt[1] = p1;
        pt[2] = p2;
      }
    } else {
      if (l2 >= l3) {
        pt[0] = p2;
        if (l1 >= l3) {
          pt[1] = p1;
          pt[2] = p3;
        } else {
          pt[1] = p3;
          pt[2] = p1;
        }
      } else {
        pt[0] = p3;
        pt[1] = p2;
        pt[2] = p1;
      }
    }

    if ((not_in_set(i, pt[0])) && (pt[0][2] != DBL_MAX))
      rescal(pt[0], grid[index_t][index_n].fw1, grid[index_t][index_n].fw2,
             sigma1, sigma2, i, step);
    else if ((not_in_set(i, pt[1])) && (pt[1][2] != DBL_MAX))
      rescal(pt[1], grid[index_t][index_n].fw1, grid[index_t][index_n].fw2,
             sigma1, sigma2, i, step);
    else if ((not_in_set(i, pt[2])) && (pt[2][2] != DBL_MAX))
      rescal(pt[2], grid[index_t][index_n].fw1, grid[index_t][index_n].fw2,
             sigma1, sigma2, i, step);
    else {
      rescal(pt_cur, grid[index_t][index_n].fw1, grid[index_t][index_n].fw2,
             sigma1, sigma2, i, step);
      index_pt[i - 1] = 1;
    }

    /* Stock values fk */
    for (k = 0; k < grfnparam->node_dim; k++) {
      /* Linear interpolation - first treatment of the boundaries */
      f1 = p1[2 + k];
      f2 = p2[2 + k];
      f3 = p3[2 + k];

      /* Checks for values outside the convex hull */
      if (f1 == DBL_MAX) {
        if (f2 == DBL_MAX) {
          if (f3 == DBL_MAX) {
            chiure = 1;
          } else {
            f1 = f3;
            f2 = f3;
          }
        } else {
          f1 = f2;
          if (f3 == DBL_MAX)
            f3 = f2;
        }
      } else if (f2 == DBL_MAX) {
        f2 = f1;
        if (f3 == DBL_MAX)
          f3 = f1;
      } else if (f3 == DBL_MAX)
        f3 = f2;

      fik[i - 1][k] = l1 * f1 + l2 * f2 + l3 * f3;
    }
  }

  /* Inverts matrix and compute transition probabilities */
  err = gaussj(coeffs, NB_NEIGHBOURS, realproba, 1);
  if (err)
    return (err);

  bb = 0;
  for (i = 1; i <= NB_NEIGHBOURS; i++)
    bb += ((realproba[i][1] <= -EPSICFL) || (realproba[i][1] >= 1.0 + EPSICFL));

  if (chiure == 1) {
    l = 0;
    while (fik[l][0] == DBL_MAX)
      l++;

    for (i = 1; i <= NB_NEIGHBOURS; i++)
      for (k = 0; k < grfnparam->node_dim; k++)
        fik[i - 1][k] = fik[l][k];
  }

  /* Loop on all the cash flows */
  for (k = 0; k < grfnparam->node_dim; k++)
    cur_asset[k] = 0.0;

  /* Loop on the nodes  to which cur_node connects at next step */
  for (i = 1; i <= NB_NEIGHBOURS; i++) {
    /* Discount the assets */
    for (k = 0; k < grfnparam->node_dim; k++) {
      if (bb > 0) /* linear interpolation */
        cur_asset[k] += grid[index_t][index_n].df * proba[i] * fik[i - 1][k];
      else {
        /* order 2 quadrature        , fits the moments */
        if (index_pt[i - 1] == 1)
          fxk = fik[i - 1][k];
        else
          fxk = ptx[i - 1][2 + k];
        cur_asset[k] += grid[index_t][index_n].df * realproba[i][1] * fxk;
      }
    }

  } /* END of the i LOOP */

  return NULL;
} /* END of Err srt_f_cheybeta2ftree_discount_to_node() */

#undef coords_bary
#undef LAMBDA
#undef BOUND_FACTOR
#undef EPSICFL

/*--------------------------------- End of File
 * -------------------------------------*/
