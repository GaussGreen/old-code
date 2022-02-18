/* ----------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & VE

   DATE : JULY 98

   FILENAME:  srt_f_lgm2ftree.c

   PURPOSE:  main implementation of the 2 factor LGM tree in GRFN

   MODIFICATION:

   ----------------------------------------------------------------------------------*/
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgm2ftreefct.h"
//#include "C:\WORKAREA\Development\CHART\inc\ChartTools.h"

#define SWAP(a, b)         \
    {                      \
        void* temp = a;    \
        a          = b;    \
        b          = temp; \
    }

#define ALLOC_BIGMEM_ASSETS         \
    {                               \
        *assets = f4tensor(         \
            0,                      \
            num_time_pts - 1,       \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1],  \
            0,                      \
            nbi);                   \
        *prob = f4tensor(           \
            0,                      \
            num_time_pts - 1,       \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1],  \
            0,                      \
            4);                     \
        *df = f3tensor(             \
            0,                      \
            num_time_pts - 1,       \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1]); \
        *dom = f3tensor(            \
            0,                      \
            num_time_pts - 1,       \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1]); \
    }

#define FREE_BIGMEM_ASSETS                                                            \
    {                                                                                 \
        free_f4tensor(assets, 0, num_time_pts - 1, -maxx, maxx, -maxy, maxy, 0, nbi); \
        free_f4tensor(prob, 0, num_time_pts - 1, -maxx, maxx, -maxy, maxy, 0, 4);     \
        free_f3tensor(df, 0, num_time_pts - 1, -maxx, maxx, -maxy, maxy);             \
        free_f3tensor(dom, 0, num_time_pts - 1, -maxx, maxx, -maxy, maxy);            \
    }

/* --------------------------------------------------------------------------------------------- */
/* --------------------------------------------------------------------------------------------- */

Err srt_f_lgm2dtree_proj(
    SrtUndPtr     und,       /* LGM2F underlying */
    SrtGrfnParam* grfnparam, /* Discretisation parameters */
    SrtStpPtr     stp,       /* Pointer on step list */
    GrfnDeal*     gd,        /* Deal description */
    EvalEventFct  evalcf,    /* DF evaluation function */
    SrtIOStruct*  iolist,    /* Requests */
    SrtUndInfo*   und_info   /* Underlying info (used in evalcf) */
)
{
    Err        err = NULL;
    Err        srt_f_lgm2dtree_forminim();
    Err        simul_and_regress_mcstatic(), simul_and_regress_pdestatic();
    int        nbi_inf = grfnparam->colmininf, nbi_sup = grfnparam->colminsup;
    int        nbi = grfnparam->colminsup - grfnparam->colmininf + 1;
    double **  vec, **cov, *alpha;
    double     residu;
    double ****assets, ****prob, ***df, ***dom;
    int        k, num_time_pts, n;
    long       maxx, maxy;
    char       buffer[200];

    num_time_pts = create_index(stp) + 1;

    // allocation for projection
    cov   = dmatrix(1, nbi, 1, nbi);
    vec   = dmatrix(1, nbi, 1, 1);
    alpha = dvector(0, nbi);

    /* call to tree lgm2d */
    err = srt_f_lgm2dtree_forminim(
        und, grfnparam, stp, gd, evalcf, iolist, und_info, &assets, &prob, &df, &dom, &maxx, &maxy);
    if (err)
    {
        FREE_BIGMEM_ASSETS
        return err;
    }

    /* simulation et construction de la matrice */
    /* err = simul_and_regress_pdestatic(grfnparam, stp, assets, prob,
                                            vec, cov, nbi, &residu, maxx, maxy, num_time_pts); */
    err = simul_and_regress_mcstatic(grfnparam, stp, assets, prob, df, dom, vec, cov, nbi, &residu);
    if (err)
    {
        FREE_BIGMEM_ASSETS
        return err;
    }

    alpha[0] = residu;
    for (k = 1; k <= nbi; k++)
    {
        alpha[k] = vec[k][1];
        n        = sprintf(buffer, "%lf\n", alpha[k]);
        smessage(buffer);
    }

    /* Stores the premium in the Input/Output list */
    err = srt_f_IOstructsetpremium((iolist), SRT_NO, residu, "Tree_LGM_2f");

    /* Stores all the columns PV in the I/O list */
    err = srt_f_IOstructsetcolpvs((iolist), SRT_NO, (double*)alpha, nbi + 1, "");

    // desallocation for projection
    free_dmatrix(cov, 1, nbi, 1, nbi);
    free_dmatrix(vec, 1, nbi, 1, 1);
    free_dvector(alpha, 0, nbi);

    /* Free the tensors previously allocated */
    FREE_BIGMEM_ASSETS

    return err;
}

/* ----------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------*/

Err simul_and_regress_mcstatic(
    SrtGrfnParam* grfnparam,
    SrtStpPtr     step,
    double****    assets,
    double****    prob,
    double***     df,
    double***     dom,
    double**      vec,
    double**      cov,
    int           nbi,
    double*       residu)
{
    Err                    err = NULL;
    double                 rand, port, dt, res, r2num, r2den, dfcur, time;
    long                   idum;
    long                   k, n, l, i, j, t;
    SrtDiagTwoFacTreeInfo* trinf;
    SrtStpPtr              top = gototop(step);

    /*--- Initialise generator ---*/
    idum = -1;
    uniform(&idum);

    /*--- initialisation ---*/
    for (k = 1; k <= nbi; k++)
    {
        for (l = 1; l <= nbi; l++)
            cov[k][l] = 0.0;
        vec[k][1] = 0.0;
    }

    /* ----------------- PATHS SIMULATION ------------------------ */
    /* Loop on all the paths used to generate the initial draw */
    for (n = 0; n < grfnparam->num_MCarlo_paths; n++)
    {
        /*--- Start from the first time step ---*/
        step = top;
        t = i = j = 0;
        time      = 0.0;
        dfcur     = 1.0;
        while ((step != NULL) && (time <= grfnparam->minmaxtime))
        {
            trinf = (SrtDiagTwoFacTreeInfo*)(step->trinf);
            dt    = step->delta_t;

            i = min(max(i, -trinf->max_index[0]), trinf->max_index[0]);
            j = min(max(j, -trinf->max_index[1]), trinf->max_index[1]);

            /* matrix-vector */
            for (k = 1; k <= nbi; k++)
            {
                vec[k][1] += (dfcur * assets[t][i][j][k - 1] * assets[t][i][j][nbi]);
                for (l = 1; l <= nbi; l++)
                    cov[k][l] += (dfcur * assets[t][i][j][k - 1] * assets[t][i][j][l - 1]);
            }

            if (dom[t][i][j] > 0.0)
                break;  // exit to free boundary

            dfcur *= df[t][i][j];

            /* next point simulated */
            if (step->next != NULL)
            {
                rand = uniform(&idum);
                if (rand < prob[t][i][j][1])
                    i++;
                else if (rand < (prob[t][i][j][1] + prob[t][i][j][2]))
                    j++;
                else if (rand < (prob[t][i][j][1] + prob[t][i][j][2] + prob[t][i][j][3]))
                    i--;
                else if (
                    rand <
                    (prob[t][i][j][1] + prob[t][i][j][2] + prob[t][i][j][3] + prob[t][i][j][4]))
                    j--;
            }

            time += dt;
            step = step->next;
            t++;

        } /* END of loop on t time step */

    } /* END of loop on all paths */

    /* solve the quadratic minimization  */
    gaussj(cov, nbi, vec, 1);

    /* compute the residu with the same trajectories */
    idum = -1;
    uniform(&idum);

    /* ----------------- PATHS SIMULATION ------------------------ */
    /* Loop on all the paths used to generate the initial draw */
    res = 0.0;
    for (n = 0; n < grfnparam->num_MCarlo_paths; n++)
    {
        /*--- Start from the first time step ---*/
        step = top;
        t = i = j = 0;
        r2num = r2den = 0.0;
        time          = 0.0;
        dfcur         = 1.0;
        while ((step != NULL) && (time <= grfnparam->minmaxtime))
        {
            trinf = (SrtDiagTwoFacTreeInfo*)(step->trinf);
            dt    = step->delta_t;

            i = min(max(i, -trinf->max_index[0]), trinf->max_index[0]);
            j = min(max(j, -trinf->max_index[1]), trinf->max_index[1]);

            /* sum residu */
            /*
            port  = 0.0;
            for (k = 1 ; k <= nbi; k++)
                    port += (vec[k][1] * assets[t][i][j][k-1]);
            res += (dfcur * (port - assets[t][i][j][nbi]) * (port - assets[t][i][j][nbi]) * dt);

            port  = 0.0;
            for (k = 1 ; k <= nbi; k++)
                    port += (vec[k][1] * assets[t][i][j][k-1]);
            r2num += dfcur * port * port * dt;
            r2den += dfcur * assets[t][i][j][nbi] * assets[t][i][j][nbi] * dt;
            */

            port = 0.0;
            for (k = 1; k <= nbi; k++)
                port += (vec[k][1] * assets[t][i][j][k - 1]);
            res += (dfcur * (port - assets[t][i][j][nbi]) * (port - assets[t][i][j][nbi]) * dt);

            if (dom[t][i][j] > 0.0)
                break;  // exit to free boundary

            dfcur *= df[t][i][j];

            /* next point simulated */
            if (step->next != NULL)
            {
                rand = uniform(&idum);
                if (rand < prob[t][i][j][1])
                    i++;
                else if (rand < (prob[t][i][j][1] + prob[t][i][j][2]))
                    j++;
                else if (rand < (prob[t][i][j][1] + prob[t][i][j][2] + prob[t][i][j][3]))
                    i--;
                else if (
                    rand <
                    (prob[t][i][j][1] + prob[t][i][j][2] + prob[t][i][j][3] + prob[t][i][j][4]))
                    j--;
            }

            time += dt;
            step = step->next;
            t++;

        } /* END of loop on t time step */

        /* res += r2num/r2den; */

    } /* END of loop on all paths */

    *residu = res / grfnparam->num_MCarlo_paths;

    /* Return a success message */
    return NULL;
}

/* ----------------------------------------------------------------------------------*/
/* ----------------------------------------------------------------------------------*/

#define ALLOC_MEM_ASSETS            \
    {                               \
        cur_assets = f3tensor(      \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1],  \
            0,                      \
            node_dim - 1);          \
        next_assets = f3tensor(     \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1],  \
            0,                      \
            node_dim - 1);          \
    }

#define FREE_MEM_ASSETS             \
    {                               \
        free_f3tensor(              \
            cur_assets,             \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1],  \
            0,                      \
            node_dim - 1);          \
        free_f3tensor(              \
            next_assets,            \
            -maxtrinf.max_index[0], \
            maxtrinf.max_index[0],  \
            -maxtrinf.max_index[1], \
            maxtrinf.max_index[1],  \
            0,                      \
            node_dim - 1);          \
    }

/* --------------------------------------------------------------------------------------------- */

Err srt_f_lgm2dtree_forminim(
    SrtUndPtr     und,       /* LGM2F underlying */
    SrtGrfnParam* grfnparam, /* Discretisation parameters */
    SrtStpPtr     stp,       /* Pointer on step list */
    GrfnDeal*     gd,        /* Deal description */
    EvalEventFct  evalcf,    /* DF evaluation function */
    SrtIOStruct*  iolist,    /* Requests */
    SrtUndInfo*   und_info,  /* Underlying info (used in evalcf) */
    double*****   assets,
    double*****   prob,
    double****    df,
    double****    dom,
    long*         maxx,
    long*         maxy)
{
    Err       err = NULL;
    int       i, j;
    double    cashflow;
    double    zc_yield;
    SrtStpPtr top, bot;
    /* Information about the  tree */
    SrtPentoTreeNodeInfo  node;
    SrtDiagTwoFacTreeInfo maxtrinf, *trinf;
    /* Information needed for the model at a particular time*/
    SrtIRMTmInf* tminf;
    /* Functions values */
    double ***next_assets, ***cur_assets;
    /* Dimensions */
    int node_dim = grfnparam->node_dim;
    /* For basis calculus */
    double y1, y2;
    /* for the storage */
    int num_time_pts, nt, k;
    int nbi = grfnparam->colminsup - grfnparam->colmininf + 1;

    // not initialise : default node_dim
    if (grfnparam->colmincible == 0)
        grfnparam->colmincible = grfnparam->node_dim - 1;

    /*--------- END of declarations -------------*/
    top = stp    = gototop(stp);
    bot          = gotobot(stp);
    num_time_pts = create_index(top) + 1;

    /* Allocates time info pointer at each step of the stp structure */
    if (err = srtstptminfalloc(stp, 1))
        return err;

    /* Initialises the steps for an LGM2F model: allocates and populates time info */
    if (err = srt_f_irministp(stp, und, 0, und, und_info))
        return err;

    /* Allocates and populates tree info */
    if (err = lgm2f_trelim(stp, &maxtrinf))
        return err;
    *maxx = maxtrinf.max_index[0];
    *maxy = maxtrinf.max_index[1];

    /* Memory Allocation */
    ALLOC_MEM_ASSETS
    ALLOC_BIGMEM_ASSETS
    if (!next_assets || !cur_assets)
    {
        FREE_MEM_ASSETS
        return serror("allocation failure srt_f_lgm2dtree");
    }

    /* Backward induction ---- Back for the futures */
    for (nt = num_time_pts - 1, stp = bot; stp != NULL; stp = stp->prev, nt--)
    {
        /* Get trinf and tminf */
        trinf = (SrtDiagTwoFacTreeInfo*)(stp->trinf);
        tminf = (SrtIRMTmInf*)(stp->tminf[0]);

        /* Stores PHIs independently of the node in LGM (depend only on time) */
        sam_get(node.cur_sam, 0, PHI1)     = sam_get(tminf->fwd_sam, 0, PHI1);
        sam_get(node.cur_sam, 0, PHI2)     = sam_get(tminf->fwd_sam, 0, PHI2);
        sam_get(node.cur_sam, 0, CROSSPHI) = sam_get(tminf->fwd_sam, 0, CROSSPHI);

        /* Loop on XiYj  in dual basis */
        for (i = -trinf->max_index[0]; i <= trinf->max_index[0]; i++)
            for (j = -trinf->max_index[1]; j <= trinf->max_index[1]; j++)
            {
                /* Transfer statevars back into original basis */
                y1 = i * trinf->spacing[0];
                y2 = j * trinf->spacing[1];
                sam_get(node.cur_sam, 0, X1) =
                    trinf->dual_basis[0][0] * y1 + trinf->dual_basis[0][1] * y2;
                sam_get(node.cur_sam, 0, X2) =
                    trinf->dual_basis[1][0] * y1 + trinf->dual_basis[1][1] * y2;

                /* If it is not the last step, compute forward, connections and probas */
                if (stp->next)
                    populate_lgm2f_tree_node(tminf, &node, stp->delta_t, stp->next->trinf);

                /* r1 = X1 + H11 + H12, r2 = X2 + H22 + H21,  r = F(0,t) + r1 + r2 */
                sam_get(node.cur_sam, 0, X1) += tminf->rf.twof[0][0].H + tminf->rf.twof[0][1].H;
                sam_get(node.cur_sam, 0, X2) += tminf->rf.twof[1][0].H + tminf->rf.twof[1][1].H;
                sam_get(node.cur_sam, 0, SHORT_RATE) = sam_get(tminf->fwd_sam, 0, F_0_t) +
                                                       sam_get(node.cur_sam, 0, X1) +
                                                       sam_get(node.cur_sam, 0, X2);
                sam_get(node.cur_sam, 0, STATEVAR) = sam_get(node.cur_sam, 0, X1);

                /* If it is not the last step, discount to node */
                if (stp->next)
                {
                    /* Compute discount factor for node */
                    Y_T_at_t_compute(1, &node.cur_sam, &tminf->yp, &zc_yield, 0, TWO_FAC, LGM);
                    node.df = exp(-zc_yield);

                    /* Compute discounted asset prices at current node */
                    lgm2f_tree_expectation(&node, next_assets, cur_assets[i][j], node_dim);
                }

                /* Evaluate cash flows in Grfn tableau */
                if (err = evalcf(
                        (GrfnEvent*)stp->e,
                        &node.cur_sam,
                        gd,
                        (double*)cur_assets[i][j],
                        (EvalEventDfsFct)srt_f_calc_grfn_event_dfs,
                        und_info,
                        &cashflow))
                {
                    FREE_MEM_ASSETS
                    return err;
                }

                /* ----- store the data, */
                for (k = 0; k <= 4; k++)
                    (*prob)[nt][i][j][k] = node.p[k];

                (*df)[nt][i][j] = node.df;

                if (grfnparam->minfreedom != 0)
                {
                    if (cur_assets[i][j][grfnparam->minfreedom] < 1.0)
                        (*dom)[nt][i][j] = 0.0;
                    else
                        (*dom)[nt][i][j] = 1.0;
                }

                for (k = 0; k < nbi; k++)
                    (*assets)[nt][i][j][k] = cur_assets[i][j][grfnparam->colmininf + k];
                (*assets)[nt][i][j][nbi] = cur_assets[i][j][grfnparam->colmincible];

            } /* End of loop on i, j */

        SWAP(cur_assets, next_assets);
    } /* End of backward induction */

    /* Free the tensors previously allocated */
    FREE_MEM_ASSETS

    return err;

} /* END Err srt_f_lgm2dtree(...) */

/* -------------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------------- */

Err simul_and_regress_pdestatic(
    SrtGrfnParam* grfnparam,
    SrtStpPtr     step,
    double****    assets,
    double****    prob,
    double***     df,
    double**      vec,
    double**      cov,
    int           nbi,
    double*       residu,
    long          maxx,
    long          maxy,
    long          maxt)
{
    Err                    err = NULL;
    double                 rand, port, dt, res;
    long                   idum;
    long                   k, n, l, i, j, t, im, ip, jm, jp;
    SrtDiagTwoFacTreeInfo* trinf;
    double **              Rt_cur, **Rt_next;

    /*--- Initialise generator ---*/
    idum = -1;
    uniform(&idum);

    /*--- initialisation ---*/
    for (k = 1; k <= nbi; k++)
    {
        for (l = 1; l <= nbi; l++)
            cov[k][l] = 0.0;
        vec[k][1] = 0.0;
    }

    /* ---------------- PATHS SIMULATION AND QUAD MINIMIZATION -------------- */
    for (n = 0; n < grfnparam->num_MCarlo_paths; n++)
    {
        /*--- Start from the first time step ---*/
        step = gototop(step);
        t    = 0;
        i    = 0;
        j    = 0;
        while (step->next != NULL)
        {
            trinf = (SrtDiagTwoFacTreeInfo*)(step->trinf);
            dt    = step->delta_t;
            /* matrix-vector */
            for (k = 1; k <= nbi; k++)
            {
                vec[k][1] += (assets[t][i][j][k - 1] * assets[t][i][j][nbi] * dt);
                for (l = 1; l <= nbi; l++)
                    cov[k][l] += (assets[t][i][j][k - 1] * assets[t][i][j][l - 1] * dt);
            }

            /* next point simulated */
            rand = uniform(&idum);
            if (rand < prob[t][i][j][1])
            {
                if (i < trinf->max_index[0])
                    i++;
            }
            else if (rand < (prob[t][i][j][1] + prob[t][i][j][2]))
            {
                if (j < trinf->max_index[1])
                    j++;
            }
            else if (rand < (prob[t][i][j][1] + prob[t][i][j][2] + prob[t][i][j][3]))
            {
                if (i > -trinf->max_index[0])
                    i--;
            }
            else if (
                rand < (prob[t][i][j][1] + prob[t][i][j][2] + prob[t][i][j][3] + prob[t][i][j][4]))
            {
                if (j > -trinf->max_index[1])
                    j--;
            }

            step = step->next;
            t++;
        } /* END of loop on t time step */
    }     /* END of loop on all paths */

    /* solve the quadratic minimization  */
    gaussj(cov, nbi, vec, 1);

    /*----------- compute the residu by PDE ----------------*/
    Rt_cur  = dmatrix(-maxx, maxx, -maxy, maxy);
    Rt_next = dmatrix(-maxx, maxx, -maxy, maxy);
    for (t = maxt - 1, step = gotobot(step)->prev; step != NULL; step = step->prev, t--)
    {
        trinf = (SrtDiagTwoFacTreeInfo*)(step->trinf);
        dt    = step->delta_t;
        /* Loop on XiYj */
        for (i = -trinf->max_index[0]; i <= trinf->max_index[0]; i++)
            for (j = -trinf->max_index[1]; j <= trinf->max_index[1]; j++)
            {
                /* residu */
                port = 0.0;
                for (k = 1; k <= nbi; k++)
                    port += (vec[k][1] * assets[t][i][j][k - 1]);
                res = ((port - assets[t][i][j][nbi]) * (port - assets[t][i][j][nbi]) * dt);

                im = (i == -trinf->max_index[0]) ? i : i - 1;
                ip = (i == trinf->max_index[0]) ? i : i + 1;
                jm = (j == -trinf->max_index[1]) ? j : j - 1;
                jp = (j == trinf->max_index[1]) ? j : j + 1;

                Rt_cur[i][j] =
                    res + prob[t][i][j][0] * Rt_next[i][j] + prob[t][i][j][1] * Rt_next[ip][j] +
                    prob[t][i][j][2] * Rt_next[i][jp] + prob[t][i][j][3] * Rt_next[im][j] +
                    prob[t][i][j][4] * Rt_next[i][jm];
            } /* End of loop on i, j */

        SWAP(Rt_cur, Rt_next);
    } /* End of backward induction */

    *residu = Rt_next[0][0];

    free_dmatrix(Rt_cur, -maxx, maxx, -maxy, maxy);
    free_dmatrix(Rt_next, -maxx, maxx, -maxy, maxy);

    /* Return a success message */
    return NULL;
}

#undef FREE_MEM_ASSETS
#undef ALLOC_MEM_ASSETS
#undef SWAP
#undef FREE_BIGMEM_ASSETS
#undef ALLOC_BIGMEM_ASSETS

/* ========= END OF FILE =================================================== */