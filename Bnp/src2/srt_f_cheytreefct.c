/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT
 */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_CHEYTREEFCT.C                                     */
/*                                                                            */
/*      PURPOSE:        Functions to compute Cheyette tree                    */
/*                                                                            */
/*      AUTHORS:        E. Auld, K.L. Chau, J. Malhi, O. Van Eyseren          */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheydynamics.h"
#include "srt_h_cheytreefct.h"

/** assume that we are probably increasing, so we take as a guess
    the value in node->mid_son_index.
    want index of element of prev_r closest to drift_sam.short_rate;
    values in prev_r are increasing **/

/*------------------------------------------------------------------------
  FUNCNAME        :srt_f_chetrerindex
  AUTHOR          :E.Auld
  DESCRIPTION     :finds the index to center on in cheyette tree
  MODIFIES		  :node->mid_son_index
  ---------------------------------------------------------------------------*/

void srt_f_chetrerindex(
    SrtStpPtr        stp,        /* current step */
    double*          prev_r,     /* rates at previous step */
    SrtTrinTreNdInf* node,       /* current node in tree */
    SrtGrfnParam*    grfnparam,  /* model parameter */
    int              closestmeth /* CLOSESTINX or CLOSESTINLNX */
)

{
    SrtCheTreInf* trinf = stp->next->trinf;

    find_closest_strictly_below(
        sam_get(node->drift_sam, 0, SHORT_RATE),
        prev_r,
        trinf->min_r_index,
        trinf->max_r_index,
        &node->mid_son_index);

    /* if index just above is closer, increment */
    trinf = stp->trinf;

    if (closestmeth == CLOSESTINX)
    {
        if (node->mid_son_index < trinf->max_r_index &&
            0.5 * (prev_r[node->mid_son_index] + prev_r[node->mid_son_index + 1]) <
                sam_get(node->drift_sam, 0, SHORT_RATE))
        {
            node->mid_son_index++;
        }
    }
    else
    {
        if (node->mid_son_index < trinf->max_r_index &&
            ((prev_r[node->mid_son_index] * prev_r[node->mid_son_index + 1]) <
             (sam_get(node->drift_sam, 0, SHORT_RATE) * sam_get(node->drift_sam, 0, SHORT_RATE))))
        {
            node->mid_son_index++;
        }
    }
    return;
}

#define bound(X, LUB, UPB) \
    {                      \
        X = DMAX(X, LUB);  \
        X = DMIN(X, UPB);  \
    }

/** calculate probabilities **/
/*-----------------------------------------------------------------------
  FUNCNAME        :srt_f_chetreprob
  AUTHOR          :E.Auld
  DESCRIPTION     :calculate probabilities,
                                   assuming that node->mid_son_inde is correctly set,
                                   and also node->drift_sam and node->var_at_sam.
                                   If the mid_son_index hits the top (bottom),
                                   it is decremented (incremented).
                                   Probabilities are bounded to be between PROBUPB and
                                   PROBLUB
  MODIFIES		  :node->p
  CALL            :

---------------------------------------------------------------------------*/

static double PROBUPB = 2.0;
static double PROBLUB = -1.0;
static long   PRINFO  = 0;

void srt_f_chetreprob(SrtStpPtr stp, double* prev_r, SrtTrinTreNdInf* node, SrtGrfnParam* grfnparam)

{
    SrtCheTreInf* trinf;
    double        eta, up_x, center_x, down_x;
    double        fwd           = sam_get(node->drift_sam, 0, SHORT_RATE);
    double        second_moment = node->var_at_sam + fwd * fwd;
    int           bndflg        = 0;

    long index = node->mid_son_index;

    trinf = stp->next->trinf;

    if (index == trinf->max_r_index)
    {
        index--;
    }

    if (index == trinf->min_r_index)
    {
        index++;
        bndflg = 1;
    }

    up_x     = prev_r[index + 1];
    center_x = prev_r[index];
    down_x   = prev_r[index - 1];

    bndflg = 0;
    eta    = fwd - center_x;
    if (fabs(up_x - center_x) < 1.0e-9 || fabs(down_x - center_x) < 1.0e-9 || down_x > 1.0e10)
    {
        node->p[0] = node->p[2] = 1.0 / 3.0;
        node->p[1]              = 1 - node->p[0] - node->p[2];
    }
    else
    {
        node->p[2] = second_moment - center_x * center_x - eta * (center_x + down_x);
        node->p[2] /= ((up_x - center_x) * (up_x - down_x));
        node->p[0] = ((up_x - center_x) * node->p[2] - eta) / (center_x - down_x);
        node->p[1] = 1 - node->p[2] - node->p[0];
    }
    if (node->p[2] > 1 || node->p[2] < 0 || node->p[1] > 1 || node->p[1] < 0 || node->p[0] > 1 ||
        node->p[0] < 0)
    {
        if (bndflg)
        {
            bound(node->p[0], PROBLUB, PROBUPB);
            bound(node->p[1], PROBLUB, PROBUPB);
            bound(node->p[2], PROBLUB, PROBUPB);
        }
    }
    return;
}

/*-----------------------------------------------------------------------
  FUNCNAME			: srt_f_chetredcntvec
  AUTHOR			: E.Auld
  DESCRIPTION		: discounts three vectors in che tree
                                          i.e. cur will be discounted probability weighted sum
                                          of three of the elements of assets
  MODIFIES			: cur
---------------------------------------------------------------------------*/
void srt_f_chetredcntvec(
    SrtStpPtr        stp,
    double*          cur,
    SrtTrinTreNdInf* node,
    double***        assets,
    int              node_dim,
    double**         prev_phi)
{
    int           i, j, r_i;
    double        cur_val;
    double        asset_phi1, asset_phi2;
    static long   phiind[3];
    SrtCheTreInf* trinf;

    trinf = stp->next->trinf;

    r_i = node->mid_son_index;

    for (j = 0; j < 3; j++)
    {
        find_closest_strictly_below(
            sam_get(node->drift_sam, 0, PHI),
            prev_phi[r_i - 1 + j],
            trinf->min_phi_index,
            trinf->max_phi_index,
            &phiind[j]);
        if (phiind[j] == trinf->max_phi_index)
        {
            phiind[j] -= 1;
        }
    }

    for (i = 0; i < node_dim; i++)
    {
        cur_val = 0;

        for (j = 0; j < 3; j++)
        {
            /* added AS: better */
            if (r_i - 1 + j >= trinf->max_r_index)
            {
                asset_phi1 = assets[trinf->max_r_index][phiind[j]][i];
                asset_phi2 = assets[trinf->max_r_index][phiind[j] + 1][i];
            }
            else if (r_i - 1 + j <= trinf->min_r_index)
            {
                asset_phi1 = assets[trinf->min_r_index][phiind[j]][i];
                asset_phi2 = assets[trinf->min_r_index][phiind[j] + 1][i];
            }
            else
            {
                asset_phi1 = assets[r_i - 1 + j][phiind[j]][i];
                asset_phi2 = assets[r_i - 1 + j][phiind[j] + 1][i];
            }
            /* end added AS */

            cur_val += node->p[j] * interp_linear(
                                        sam_get(node->drift_sam, 0, PHI),
                                        prev_phi[r_i - 1 + j][phiind[j]],
                                        asset_phi1,
                                        prev_phi[r_i - 1 + j][phiind[j] + 1],
                                        asset_phi2);
        }
        cur[i] = cur_val * node->df;
    }
} /* void srt_f_chetredcntvec(...) */

/* ----------------------------------------------------------------------------- */

/* cutting */
#define CUTAFTERSTDEV 4
#define CUTYESNO 1

static void det_tre_lims(
    SrtStpPtr     top,
    SrtStpPtr     stp,
    SrtCheTreInf* trinf,
    SrtIRMTmInf*  tminf,
    SrtCheTreInf* prvtrinf,
    SrtIRMTmInf*  prvtminf,
    SrtGrfnParam* grfnparam);

/* --------------------------------------------------------------------------------
   FUNCNAME        :srt_f_cheytreelim
   DESCRIPTION     :allocates and populates SrtCheTreInf struct at each step
                                   population costs in calling det_tre_lims at each step
                                   maxtreeinf will contain the biggest maxr and maxphi,
                                   and the smallest minr and minphi
   -------------------------------------------------------------------------------- */
Err srt_f_cheytreelim(SrtStpPtr stp, SrtGrfnParam* grfnparam, SrtCheTreInf* maxtrinf)
{
    Err           err;
    SrtStpPtr     top;
    SrtCheTreInf *trinf, *prvtrinf;
    SrtIRMTmInf * tminf, *prvtminf;

    top = stp = gototop(stp);
    memset(maxtrinf, 0, sizeof(SrtCheTreInf));

    /* allocate space */
    err = srtstpalloc(top, sizeof(SrtCheTreInf), 1, 0);
    if (err)
    {
        return err;
    }

    prvtminf = NULL;
    prvtrinf = NULL;

    /* loop on steps */
    while (stp)
    {
        /* get time info and tree info */
        tminf = (SrtIRMTmInf*)stp->tminf[0];
        trinf = (SrtCheTreInf*)stp->trinf;

        /* compute tree limits */
        det_tre_lims(top, stp, trinf, tminf, prvtrinf, prvtminf, grfnparam);

        maxtrinf->max_r_index   = IMAX(maxtrinf->max_r_index, trinf->max_r_index);
        maxtrinf->max_phi_index = IMAX(maxtrinf->max_phi_index, trinf->max_phi_index);
        maxtrinf->min_r_index   = IMIN(maxtrinf->min_r_index, trinf->min_r_index);
        maxtrinf->min_phi_index = IMIN(maxtrinf->min_phi_index, trinf->min_phi_index);

        prvtrinf = trinf;
        prvtminf = tminf;
        stp      = stp->next;
    }

    maxtrinf->max_r_index++;
    maxtrinf->min_r_index--;

    return NULL;

} /* END Err  srt_f_cheytreelim(...) */

/* --------------------------------------------------------------------------------- */
/* this function must be called from the top in order to work correctly */
static void det_tre_lims(
    SrtStpPtr     top,      /* first step */
    SrtStpPtr     stp,      /* current step */
    SrtCheTreInf* trinf,    /* current tree info (on output) */
    SrtIRMTmInf*  tminf,    /* current time info (on input) */
    SrtCheTreInf* prvtrinf, /* previous tree info (on input) */
    SrtIRMTmInf*  prvtminf, /* previous time info (on input) */
    SrtGrfnParam* grfnparam /* model parameters*/
)
{
    static SrtSample top_sam, top_drf_sam, bot_sam, bot_drf_sam;
    static double    r_mid, r_max, r_min;
    static double    cum_var_logr;
    double           stdev_logr;
    double           logu;
    long             topind;

    /* static variables are initialised when called from first step,
       they hold results from previous call (previous step) */

    /* CAREFUL:
       top_sam = node where r is maximum
       top_drf_sam = expectations at this node
       bot_... = where r is minimum */

    /* first step: init */
    if (stp == top)
    {
        trinf->max_r_index   = 0;
        trinf->min_r_index   = 0;
        trinf->max_phi_index = 0;
        trinf->min_phi_index = 0;
        trinf->u             = 1.0;
        r_mid                = sam_get(tminf->fwd_sam, 0, F_0_t);
        trinf->logrmin       = log(r_mid);
        top_sam = bot_sam = tminf->fwd_sam;

        if (CUTYESNO)
        {
            cum_var_logr = 0.0;
        }
    }

    /* i-th step: populate */
    else
    {
        /* U_FACTOR = 2, logu = spacing = 2 * stdev */
        trinf->u             = exp(U_FACTOR * prvtminf->ev.onef.sig * stp->prev->sqrt_delta_t);
        logu                 = log(trinf->u);
        trinf->max_phi_index = IMIN(grfnparam->width_phi, stp->index);
        trinf->min_phi_index = 0;
        r_mid                = sam_get(tminf->fwd_sam, 0, F_0_t);

        /* cut the tree */
        if (CUTYESNO)
        {
            cum_var_logr += prvtminf->ev.onef.sig2 * stp->prev->delta_t;
            stdev_logr = sqrt(cum_var_logr);
            topind     = DTOL(CUTAFTERSTDEV * stdev_logr / logu) + 1;
            topind     = IMIN(topind, stp->index);
        }
        else
        {
            topind = stp->index;
        }

        trinf->max_r_index = topind;
        trinf->min_r_index = -topind;
        trinf->logrmin     = log(r_mid) + logu * trinf->min_r_index;

        /* determine max drift at top and bottom of tree */
        srt_f_chedrfatsam(stp->prev, &top_sam, &top_drf_sam, 0);
        srt_f_chedrfatsam(stp->prev, &bot_sam, &bot_drf_sam, 0);

        /** populate top_sam and bot_sam for next time this function is called **/
        /* r at top sam */
        sam_get(top_sam, 0, SHORT_RATE) =
            exp(trinf->logrmin + logu * (trinf->max_r_index - trinf->min_r_index + 1));
        /* maxphi and minphi and midphi at top sam */
        srt_f_chephilim(
            top,
            stp,
            tminf,
            &top_sam,
            &(sam_get(bot_drf_sam, 0, PHI)),
            &(sam_get(top_drf_sam, 0, PHI)));
        sam_get(top_sam, 0, PHI) = sam_get(top_drf_sam, 0, PHI);

        /* same at bot sam */
        sam_get(bot_sam, 0, SHORT_RATE) = exp(trinf->logrmin);
        srt_f_chephilim(
            top,
            stp,
            tminf,
            &bot_sam,
            &(sam_get(bot_drf_sam, 0, PHI)),
            &(sam_get(top_drf_sam, 0, PHI)));
        sam_get(bot_sam, 0, PHI) = sam_get(bot_drf_sam, 0, PHI);
    }
    return;
}

/* ---------------------------------------------------------------------------------- */