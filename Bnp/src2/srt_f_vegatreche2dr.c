/* SRT_F_VEGATRECHE2DR.c */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheydynamics.h"
#include "srt_h_cheytreefct.h"

/* -------------------------------------------------------------------------- */

#define PSWAP(X, Y)          \
    {                        \
        void* _VPTR = X;     \
        X           = Y;     \
        Y           = _VPTR; \
    }

/* -------------------------------------------------------------------------- */

/* returns pointer to step number index */
SrtStpPtr gotoindex(SrtStpPtr step, long index)
{
    SrtStpPtr cur;

    for (cur = gototop(step); cur != NULL; cur = cur->next)
    {
        if (cur->index == index)
            return cur;
    }
    return NULL;
}

/* -------------------------------------------------------------------------- */

/*	This function must answer the following question:
        Can the function compute the request with type 'request_type' ? */
SRT_Boolean srt_f_treche_requestOK(int request_type)
{
    switch (request_type)
    {
    case IO_PREMIUM:
        return SRT_YES;
        break;
    case IO_STDEV:
        return SRT_NO;
        break;
    case IO_SIGMA_SHIFT:
        return SRT_NO;
        break; /* set to yes when use it */
    case IO_TAU_SHIFT:
        return SRT_NO;
        break; /* set to yes when use it */
    case IO_RATE_SHIFT:
        return SRT_NO;
        break;
    default:
        return SRT_NO;
        break;
    }
}

/*---------------------------------------------------------------------------
  FUNCTION        :srt_f_treche2dr
  AUTHOR          :E.Auld from code by G.Amblard
  MODIFIED		  :A.Savine
  DESCRIPTION     :values a structure according to a 2d tree discretization of
                                   the cheyette interest rate model, by which we mean
                                   the one where sigma(t,r) = sigma(t)*r
-----------------------------------------------------------------------------*/

Err srt_f_vegatreche2dr(
    SrtUndPtr     und,       /* market parameters */
    SrtGrfnParam* grfnparam, /* model parameters */
    SrtStpPtr     stp,       /* step pointer */
    GrfnDeal*     gd,        /* deal description */
    EvalEventFct  evalcf,    /* cashflow evalution function*/
    SrtIOStruct*  iolist,    /* list of requests */
    SrtUndInfo*   und_info)    /* underlying info */

{
    /* standard */
    Err       err = NULL;
    int       i, j, max_r_index = 0;
    double    zc_yield;
    double    cashflow;
    SrtStpPtr top, bot;

    /* information about the  particular  node within the tree we are at */
    SrtTrinTreNdInf node;

    /* information about the dimensions of this tree */
    SrtCheTreInf maxtrinf, *trinf;

    /* information needed for the model at a particular time */
    SrtIRMTmInf* tminf;

    /* val of assets being computed through the tree at cur and prev time step
    4 dimensions:	- request - phi	- rate	- assets */

    double ****prev_assets, ****cur_assets;

    /* state variables at current time step and previous time step */
    double *cur_r, *prev_r, **cur_phi, **prev_phi;

    /* amount of space for assets */
    long sz;

    /* variables used for computing greeks */
    double sigma_shift, sigma_bucket_start, sigma_bucket_end;
    int    sigma_shift_type;
    double tau_shift, tau_bucket_start, tau_bucket_end;
    int    tau_shift_type;
    int    xtra_calc_num = 0;
    int    xtra_calc_index;

    /* dealing with requests */
    SrtStpPtr    my_pointer;
    SrtIOVal*    io_request_val_p;
    SrtListAtom* io_request;

    /* --------------- end of declarations --------------------------------------- */

    /* set top (1st step), bot (last step) and stp (set to top) */
    top = stp = gototop(stp);
    bot       = gotobot(stp);

    /* for each request */
    for (io_request = iolist->head; io_request != NULL; io_request = io_request->next)
    {
        /* whithin request */
        io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);

        /* can it be computed? */
        if (srt_f_treche_requestOK(io_request_val_p->type) == SRT_NO)
        {
            continue;
        }
        /* if yes */
        else
        {
            /* number of computations increases */
            xtra_calc_num++;

            /* initialise shifts for computation of greeks */
            sigma_shift        = 0;
            sigma_bucket_start = -1;
            sigma_bucket_end   = -1;
            sigma_shift_type   = SH_NONE;
            tau_shift          = 0;
            tau_bucket_start   = -1;
            tau_bucket_end     = -1;
            tau_shift_type     = SH_NONE;
            switch (io_request_val_p->type)
            {
            case IO_PREMIUM:
                break;
            case IO_SIGMA_SHIFT:
                sigma_shift        = io_request_val_p->shift_value;
                sigma_shift_type   = io_request_val_p->shift_type;
                sigma_bucket_start = io_request_val_p->bucket_start;
                sigma_bucket_end   = io_request_val_p->bucket_end;
                break;
            case IO_TAU_SHIFT:
                tau_shift        = io_request_val_p->shift_value;
                tau_shift_type   = io_request_val_p->shift_type;
                tau_bucket_start = io_request_val_p->bucket_start;
                tau_bucket_end   = io_request_val_p->bucket_end;
                break;
            default:
                return serror("invalid request in srt_f_vegatreche2dr");
                break;
            }

            /* compute premium */
            if (io_request_val_p->type == IO_PREMIUM)
            {
                /* init step list, tminf etc. */
                if (err = srt_f_vegashiftirministp(
                        stp,
                        und,
                        und_info,
                        sigma_bucket_start,
                        sigma_bucket_end,
                        sigma_shift,
                        sigma_shift_type,
                        tau_bucket_start,
                        tau_bucket_end,
                        tau_shift,
                        tau_shift_type))
                {
                    return err;
                }

                /* INITIALISES TREE INFO */
                /* compute tree limits and geometry */
                if (err = srt_f_cheytreelim(stp, grfnparam, &maxtrinf))
                {
                    return err;
                }

                ((SrtIOVal*)(io_request->element->val.pval))->pval = (void*)stp;
            } /* END of || (io_request_val_p->type == IO_PREMIUM) */
            /* vega shifts */
            else
            {
                /* my_pointer points to a copy of the SrtStp pointed by stp */
                if (err = srt_f_stpdup(stp, &my_pointer))
                {
                    return err;
                }

                /*------------------------------------------------------------------------------
                NOTE:
                The trinf field will be the same for all the request since the tree
                is FROZEN. Thus the copy will have a trinf field pointing to NULL.
                If you pass a copy of the Stp structure to a function that needs to
                access the trinf field, it will crash.........
                -------------------------------------------------------------------------------*/

                if (err = srt_f_vegashiftirministp(
                        my_pointer,
                        und,
                        und_info,
                        sigma_bucket_start,
                        sigma_bucket_end,
                        sigma_shift,
                        sigma_shift_type,
                        tau_bucket_start,
                        tau_bucket_end,
                        tau_shift,
                        tau_shift_type))
                {
                    return err;
                }

                /* each request of the SrtIOStruct has now a field that points to a SrtStp
                   structure. For this structure, the parameters have been computed with the shift
                   corresponding to the request */
                ((SrtIOVal*)(io_request->element->val.pval))->pval = (void*)my_pointer;

            } /* END if (io_request_val_p->type != IO_PREMIUM) */

        } /* END if is request OK */

    } /* END of the for loop on requests */

    /** Memory Allocation **/
    sz = xtra_calc_num * (maxtrinf.max_r_index - maxtrinf.min_r_index + 1) *
         (maxtrinf.max_phi_index - maxtrinf.min_phi_index + 1) * grfnparam->node_dim;

    prev_assets = f4tensor(
        0,
        xtra_calc_num - 1,
        maxtrinf.min_r_index,
        maxtrinf.max_r_index,
        maxtrinf.min_phi_index,
        maxtrinf.max_phi_index,
        0,
        grfnparam->node_dim - 1);
    cur_assets = f4tensor(
        0,
        xtra_calc_num - 1,
        maxtrinf.min_r_index,
        maxtrinf.max_r_index,
        maxtrinf.min_phi_index,
        maxtrinf.max_phi_index,
        0,
        grfnparam->node_dim - 1);

    cur_r   = dvector(maxtrinf.min_r_index, maxtrinf.max_r_index);
    prev_r  = dvector(maxtrinf.min_r_index, maxtrinf.max_r_index);
    cur_phi = dmatrix(
        maxtrinf.min_r_index, maxtrinf.max_r_index, maxtrinf.min_phi_index, maxtrinf.max_phi_index);
    prev_phi = dmatrix(
        maxtrinf.min_r_index, maxtrinf.max_r_index, maxtrinf.min_phi_index, maxtrinf.max_phi_index);

    if (!prev_assets || !cur_assets || !cur_r || !prev_r || !cur_phi || !prev_phi)
        return serror("memory allocation failure in srt_f_treche2dr");

    memset(&cur_assets[0][maxtrinf.min_r_index][maxtrinf.min_phi_index][0], 0, sz * sizeof(double));

    /***** TREE *****/
    /* backward induction */

    /** For each time stp **/
    for (stp = bot; stp != NULL; stp = stp->prev)
    {
        /* get trinf and tminf */
        trinf = (SrtCheTreInf*)stp->trinf;
        tminf = (SrtIRMTmInf*)stp->tminf[0];

        /** create grid of state variables **/
        srt_f_chesammatrix(top, stp, cur_r, cur_phi, grfnparam);

        /** For each request **/
        for (io_request = iolist->head; io_request != NULL; io_request = io_request->next)
        {
            io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);

            /* IF stdev or rubbish requested, do not do anything */
            if (srt_f_treche_requestOK(io_request_val_p->type) == SRT_NO)
                continue;

            xtra_calc_index = 0;

            my_pointer = gotoindex((SrtStpPtr)(io_request_val_p->pval), stp->index);
            tminf      = (SrtIRMTmInf*)my_pointer->tminf[0];

            /* For each r */
            for (i = trinf->min_r_index; i <= trinf->max_r_index; i++)
            {
                sam_get(node.cur_sam, 0, SHORT_RATE) = cur_r[i];

                /* To ensure compatibility with the Y_T_at_t functions */
                sam_get(node.cur_sam, 0, STATEVAR) = cur_r[i] - sam_get(tminf->fwd_sam, 0, F_0_t);

                /** For each phi **/
                for (j = trinf->min_phi_index; j <= trinf->max_phi_index; j++)
                {
                    /*  Now the time info is taken from the specific shifted price structure.
                    The tree being frozen, the tree info keeps being taken from stp.	*/
                    sam_get(node.cur_sam, 0, PHI) = cur_phi[i][j];

                    /** if not at the end: **/
                    if (stp->next)
                    {
                        /** compute discount factor for that stp **/
                        Y_T_at_t_compute(1, &node.cur_sam, &tminf->yp, &zc_yield, 0, ONE_FAC, CHEY);
                        node.df = exp(-zc_yield);

                        /** calculate expectations at this node **/
                        srt_f_chedrfatsam(my_pointer, &node.cur_sam, &node.drift_sam, 0);

                        /** calculate r index to center on in previous level **/
                        /* The function uses inside only field trinf => stp */
                        srt_f_chetrerindex(stp, prev_r, &node, grfnparam, CLOSESTINLNX);

                        /** calculate variance at this node **/
                        /* The function uses inside only field tminf -> my_pointer instead of stp */
                        srt_f_chevaratsam(my_pointer, &node.cur_sam, &node.var_at_sam, 0);

                        /** calculate probabilities **/
                        /* the function uses inside only field trinf -> stp */
                        srt_f_chetreprob(stp, prev_r, &node, grfnparam);

                        /** calculate discounted prob weighted sum of previous cashflows **/
                        /* The function uses inside only field trinf -> stp */
                        srt_f_chetredcntvec(
                            stp,
                            cur_assets[xtra_calc_index][i][j],
                            &node,
                            prev_assets[xtra_calc_index],
                            grfnparam->node_dim,
                            prev_phi);
                    }

                    /** call cashflow function **/
                    err = evalcf(
                        (GrfnEvent*)stp->e,
                        &node.cur_sam,
                        gd,
                        (double*)cur_assets[xtra_calc_index][i][j],
                        (EvalEventDfsFct)srt_f_calc_grfn_event_dfs,
                        und_info,
                        &cashflow);
                    if (err)
                    {
                        return err;
                    }
                } /* phi loop*/
            }     /* r loop */

            xtra_calc_index++;
        } /* request loop */

        PSWAP(cur_assets, prev_assets);
        PSWAP(cur_r, prev_r);
        PSWAP(cur_phi, prev_phi);

    } /* END of loop on time steps */

    /* Stores all the results for all the requests in the IO list */
    xtra_calc_index = 0;
    for (io_request = iolist->head; io_request != NULL; io_request = io_request->next)
    {
        io_request_val_p = (SrtIOVal*)(io_request->element->val.pval);

        if (srt_f_treche_requestOK(io_request_val_p->type) == SRT_NO)
            continue;
        strncpy(io_request_val_p->computation_origin, "Tree_CHE_2d", strlen("Tree_CHE_2d"));
        io_request_val_p->dval = prev_assets[xtra_calc_index][0][0][grfnparam->node_dim - 1];

        /*	free the SrtStp list duplicated from stp except for
                the PREMIUM which points to the original Stp    */
        if (io_request_val_p->type != IO_PREMIUM)
        {
            free_list(io_request_val_p->pval);
        }
        xtra_calc_index++;
    }

    /* Stores all the columns PV in the I/O list */
    err = srt_f_IOstructsetcolpvs(
        (SrtIOStruct*)(iolist), SRT_NO, (double*)prev_assets[0][0][0], grfnparam->node_dim, "");
    if (err)
        return err;

    /** free the world **/
    free_f4tensor(
        prev_assets,
        0,
        xtra_calc_num - 1,
        maxtrinf.min_r_index,
        maxtrinf.max_r_index,
        maxtrinf.min_phi_index,
        maxtrinf.max_phi_index,
        0,
        grfnparam->node_dim - 1);
    free_f4tensor(
        cur_assets,
        0,
        xtra_calc_num - 1,
        maxtrinf.min_r_index,
        maxtrinf.max_r_index,
        maxtrinf.min_phi_index,
        maxtrinf.max_phi_index,
        0,
        grfnparam->node_dim - 1);

    free_dvector(cur_r, maxtrinf.min_r_index, maxtrinf.max_r_index);
    free_dvector(prev_r, maxtrinf.min_r_index, maxtrinf.max_r_index);
    free_dmatrix(
        cur_phi,
        maxtrinf.min_r_index,
        maxtrinf.max_r_index,
        maxtrinf.min_phi_index,
        maxtrinf.max_phi_index);
    free_dmatrix(
        prev_phi,
        maxtrinf.min_r_index,
        maxtrinf.max_r_index,
        maxtrinf.min_phi_index,
        maxtrinf.max_phi_index);

    return NULL;
}
