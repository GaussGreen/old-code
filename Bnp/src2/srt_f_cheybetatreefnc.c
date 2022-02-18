/* ============================================================================

   FILENAME:			srt_f_cheybetatreefnc.c

   DESCRIPTION:			utility functions for local vol Cheyette tree

   ============================================================================ */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheybetatreefnc.h"

/* ---------------------------------------------------------------------------- */

/* A numerically big number to initialise minvar */
#define SRT_BIG 1.0e+20

/* ---------------------------------------------------------------------------- */

/*	Function:		srt_f_locvol_init_trinf
    Object:			initialise the tree information in step ( geomtry and limits )
                        allocate and populate
*/

Err srt_f_CheyBetaTree_init_trinf(
    SrtStpPtr     step,      /*	Linked list of time steps*/
    SrtGrfnParam* grfnparam, /*	Mdl discretisation parameters */
    SrtUndPtr     und,       /*	Interest Rate Underlying */
    SrtLocTreInf* maxtrinf)  /*	Will contain the smallest min_r_index
                                                     and the greates max_r_index
                                                     (for allocation) */
{
    Err           err;
    SrtStpPtr     top;
    SrtLocTreInf* trinf;
    SrtLocTreInf* prvtrinf = NULL;
    SrtIRMTmInf*  tminf;
    SrtIRMTmInf*  prvtminf = NULL;
    TermStruct*   ts;
    SrtMdlType    mdl_type;
    /* double          convexified_forward; */
    double minvar;
    double cur_cum_var;
    double prev_cum_var;
    double delta_var;
    double delta_state_var;
    long   nstep;

    /* Get mdl_type */
    err = get_underlying_mdltype(und, &mdl_type);
    if (err)
    {
        return err;
    }

    /* Get Term Struct */
    err = get_underlying_ts(und, &ts);
    if (err)
    {
        return err;
    }

    /* Allocate space for the tree information pointer */
    top = gototop(step);
    err = srtstpalloc(top, sizeof(SrtLocTreInf), 1, 0);
    if (err)
    {
        return err;
    }

    /* --------------------- Compute an average variance for the tree ---------------- */
    /* Go to the first time step */
    top = gototop(step);

    /* Reset minvar and prev_cum_var (variance until previous step) */
    minvar       = SRT_BIG;
    prev_cum_var = 0, 0;

    /* Move on to the next time step */
    step = step->next;

    /* Set nstep to 0 */
    nstep = 0;

    /* Main loop on steps */
    while (step)
    {
        /* Get tminf and trinf at the current time step */
        tminf = (SrtIRMTmInf*)(step->tminf[0]);
        trinf = (SrtLocTreInf*)(step->trinf);

        /*  Get tminf and trinf at the previous time step */
        prvtrinf = (SrtLocTreInf*)(step->prev->trinf);
        prvtminf = (SrtIRMTmInf*)(step->prev->tminf[0]);

        /* Get lgm equivalent vol at previous step */
        cur_cum_var = sam_get(tminf->fwd_sam, 0, PHI);
        delta_var   = cur_cum_var - prev_cum_var;

        /* Set the "average" phi (spot rate follows forward) in the tree information */
        trinf->mid_phi = sam_get(tminf->fwd_sam, 0, PHI);

        /* Update variance */
        minvar = DMIN(minvar, delta_var);

        /* Set max r and min r in the tree information attached to the step */
        /*		convexified_forward = sam_get(tminf->fwd_sam, 0, SHORT_RATE)
                                + sam_get(tminf->fwd_sam, 0, PHI);
        */
        trinf->min_state_var = -CUTTREEATSTDEV * sqrt(cur_cum_var);
        trinf->max_state_var = +CUTTREEATSTDEV * sqrt(cur_cum_var);
        /*		trinf->min_r = convexified_forward + trinf->min_state_var;
                        trinf->max_r = convexified_forward + trinf->max_state_var;
        */
        /* Increment nstep */
        nstep++;

        /* Move on to next time step */
        step = step->next;

        /* Update var */
        prev_cum_var = cur_cum_var;

    } /* END of loop on steps */

    /* Set the r spacing delta_r to sqrt of minimum of (average_stdev, min_stdev * 2) */
    delta_state_var           = sqrt(DMIN(cur_cum_var / nstep, minvar * 2));
    maxtrinf->delta_state_var = delta_state_var;

    /* --------------------- Set up the tree limits according to delta_r ---------------- */

    /* Initialise the limits for the first step */
    step                       = top;
    trinf                      = (SrtLocTreInf*)(step->trinf);
    trinf->delta_state_var     = delta_state_var;
    trinf->delta_phi           = 0.0;
    trinf->max_state_var_index = trinf->min_state_var_index = maxtrinf->min_state_var_index =
        maxtrinf->max_state_var_index                       = 0;
    trinf->max_phi_index = trinf->min_phi_index = maxtrinf->min_phi_index =
        maxtrinf->max_phi_index                 = 0;
    /* Initialise the limits for the following steps */
    step = step->next;
    while (step)
    {
        /* Get the tree information attached to the step */
        trinf = (SrtLocTreInf*)(step->trinf);

        /* Set up the relevant limits for the r discretisation*/
        trinf->delta_state_var = delta_state_var;
        trinf->max_state_var_index =
            IMIN(step->index, DTOL(trinf->max_state_var / delta_state_var) + 1);
        trinf->max_state_var = trinf->max_state_var_index * delta_state_var;
        maxtrinf->max_state_var_index =
            IMAX(maxtrinf->max_state_var_index, trinf->max_state_var_index);
        trinf->min_state_var_index =
            IMAX(-step->index, DTOL(trinf->min_state_var / delta_state_var) - 1);
        trinf->min_state_var = trinf->min_state_var_index * delta_state_var;
        maxtrinf->min_state_var_index =
            IMIN(maxtrinf->min_state_var_index, trinf->min_state_var_index);

        /* Sets up the relevant limits for the phi dicretisation */
        trinf->max_phi_index    = IMIN(grfnparam->width_phi, step->index);
        maxtrinf->max_phi_index = IMAX(maxtrinf->max_phi_index, trinf->max_phi_index);
        trinf->min_phi_index    = 0;
        trinf->max_phi          = DMAX(2.0, step->time) * trinf->mid_phi;
        trinf->min_phi          = DMIN(0.5, 1.0 / step->time) * trinf->mid_phi;
        trinf->delta_phi =
            (trinf->max_phi - trinf->min_phi) / (trinf->max_phi_index - trinf->min_phi_index);

        /* Move on to the next step */
        step = step->next;

    } /* END of loop on steps */

    return NULL;

} /* END of function  Err srt_f_CheyBetaTree_init_trinf(...) */

#undef SRT_BIG

/* ----------------------------------------------------------------------------- */

/*	Function:		srt_f_locvol_make_r_grid
        Object:			make r grid, once for all */

Err srt_f_CheyBetaTree_make_state_var_grid(SrtLocTreInf maxtrinf, double* r_grid)
{
    long i;

    for (i = maxtrinf.min_state_var_index; i <= maxtrinf.max_state_var_index; i++)
    {
        r_grid[i] = i * maxtrinf.delta_state_var;
    }

    return NULL;
}

/* ----------------------------------------------------------------------------- */

/*	Function:		srt_f_locvol_make_phi_grid
        Object:			make phi grid, for a given time step */

Err srt_f_CheyBetaTree_make_phi_grid(SrtLocTreInf* trinf, double* cur_phi_grid)
{
    long i;

    for (i = trinf->min_phi_index; i <= trinf->max_phi_index; i++)
    {
        cur_phi_grid[i] = trinf->min_phi + i * trinf->delta_phi;
    }

    return NULL;
}

/* ----------------------------------------------------------------------------- */

/*	Function:		srt_f_locvol_recombine_in_r
        Object:			Compute recombination of the tree in r */

Err srt_f_CheyBetaTree_recombine_in_r(SrtTrinTreNdInf* node, double* r_grid, SrtLocTreInf* nxttrinf)
{
    double stdev;
    double target;

    stdev = sqrt(node->state_var_variance);

    /* Find the mid son index the closest to the forward of the SHORT RATE */
    target = sam_get(node->drift_sam, 0, STATEVAR);
    find_closest(
        target,
        r_grid,
        nxttrinf->min_state_var_index + 1,
        nxttrinf->max_state_var_index - 1,
        &(node->son_index[TREE_MID]));

    node->son_state_var_value[TREE_MID] = r_grid[node->son_index[TREE_MID]];

    /* Find down son indexes depending on the local variance (from guess = node -1 ) */
    target                     = node->state_var_expectation - NUMSTDEVINSPACING * stdev;
    node->son_index[TREE_DOWN] = node->son_index[TREE_MID] - 1;
    find_closest(
        target,
        r_grid,
        nxttrinf->min_state_var_index,
        node->son_index[TREE_MID] - 1,
        &(node->son_index[TREE_DOWN]));
    node->son_state_var_value[TREE_DOWN] = r_grid[node->son_index[TREE_DOWN]];

    /* Find up son indexes depending on the local variance (from guess = node + 1) */
    target                   = node->state_var_expectation + NUMSTDEVINSPACING * stdev;
    node->son_index[TREE_UP] = node->son_index[TREE_MID] + 1;
    find_closest(
        target,
        r_grid,
        node->son_index[TREE_MID] + 1,
        nxttrinf->max_state_var_index,
        &(node->son_index[TREE_UP]));
    node->son_state_var_value[TREE_UP] = r_grid[node->son_index[TREE_UP]];

    return NULL;
}

/* ----------------------------------------------------------------------------- */

/*	Function:		srt_f_locvol_calc_prob
        Object:			Compute probabilities */

Err srt_f_CheyBetaTree_calc_prob(SrtTrinTreNdInf* node)
{
    double Su;
    double Sm;
    double Sd;
    double E;
    double V;

    Su = node->son_state_var_value[TREE_UP], Sm = node->son_state_var_value[TREE_MID],
    Sd = node->son_state_var_value[TREE_DOWN], E = node->state_var_expectation,
    V = node->state_var_variance;

    node->p[TREE_UP] =
        1.0 + (-E + Su) / (-Su + Sm) +
        (-((-E + Su) * (-Su * Su + Sm * Sm)) + (-Su + Sm) * (-E * E + Su * Su - V)) /
            ((Sd * Sd - Su * Su) * (-Su + Sm) - (Sd - Su) * (-Su * Su + Sm * Sm)) -
        ((Sd - Su) * (-((-E + Su) * (-Su * Su + Sm * Sm)) + (-Su + Sm) * (-E * E + Su * Su - V))) /
            ((-Su + Sm) * ((Sd * Sd - Su * Su) * (-Su + Sm) - (Sd - Su) * (-Su * Su + Sm * Sm)));

    node->p[TREE_MID] =
        -((-E + Su) / (-Su + Sm)) +
        ((Sd - Su) * (-((-E + Su) * (-Su * Su + Sm * Sm)) + (-Su + Sm) * (-E * E + Su * Su - V))) /
            ((-Su + Sm) * ((Sd * Sd - Su * Su) * (-Su + Sm) - (Sd - Su) * (-Su * Su + Sm * Sm)));

    node->p[TREE_DOWN] = 1.0 - node->p[TREE_MID] - node->p[TREE_UP];

    if (node->p[TREE_UP] < 0.0 || node->p[TREE_MID] < 0.0 || node->p[TREE_DOWN] < 0.0)
    {
        /*	smessage ("Negative probabilities in local vol tree"); */
    }

    return NULL;
}

/* ----------------------------------------------------------------------------- */

/*	Function:		srt_f_locvol_discount_assets
        Object:			Discount asset prices (as many assets as Grfn cols ) */

Err srt_f_CheyBetaTree_discount_assets(
    double*          cur_assets,
    double***        next_assets,
    SrtTrinTreNdInf* node,
    long             nassets,
    SrtLocTreInf*    nxttrinf,
    double*          next_r_grid,
    double           next_phi,
    double*          next_phi_grid)
{
    int         i, j;
    double      asset_here_and_now;
    long        next_phi_index_high;
    static long next_phi_index_low;
    double      asset_at_phi_low;
    double      asset_at_phi_high;

    /* Find the closest Phi node to the value of Phi for the next step (below it) */
    find_closest_strictly_below(
        next_phi,
        next_phi_grid,
        nxttrinf->min_phi_index,
        nxttrinf->max_phi_index - 1,
        &next_phi_index_low);

    /* The second phi that brackets it is the following one */
    next_phi_index_high = next_phi_index_low + 1;

    /* Compute the expected value for each phi, than interpolate linearily */
    for (i = 0; i < nassets; i++)
    {
        asset_here_and_now = 0.0;

        for (j = 0; j < 3; j++)
        {
            asset_at_phi_low  = next_assets[node->son_index[j]][next_phi_index_low][i];
            asset_at_phi_high = next_assets[node->son_index[j]][next_phi_index_high][i];

            asset_here_and_now += node->p[j] * interp_linear(
                                                   next_phi,
                                                   next_phi_grid[next_phi_index_low],
                                                   asset_at_phi_low,
                                                   next_phi_grid[next_phi_index_high],
                                                   asset_at_phi_high);
        }
        cur_assets[i] = asset_here_and_now * node->df;
    }

    /* Return a success message */
    return NULL;
}

/* ----------------------------------------------------------------------------- */
