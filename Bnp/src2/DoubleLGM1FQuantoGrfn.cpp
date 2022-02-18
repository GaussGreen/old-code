#include "DoubleLGM1FQuantoGrfn.h"

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

Err payoff_doublelgm1fquanto_pde(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* dom_yc,
    void* for_yc,

    /* Model data	*/
    double domlam,
    double forlam,
    double domphi,
    double forphi,

    /* Grid data	*/
    int      l1,
    int      u1,
    int      l2,
    int      u2,
    double*  r1,
    double** r2,

    /* Vector of results to be updated */
    int       nprod,
    double*** prod_val)
{
    GRFNPARMTREE   total;
    GRFNCOMMSTRUCT global;
    FIRSTMktAtT*   local;
    double *       dom_dff = NULL, *dom_gam = NULL, *fx_dff = NULL, *fx_gam = NULL, *fx_gam2 = NULL,
           *dom_gam2 = NULL, *for_dff = NULL, *for_gam = NULL, *for_gam2 = NULL;
    int            i, j, l;
    double         temp;
    unsigned short do_dom, do_for, do_fx;
    Err            err = NULL;

    /* Get the event */
    total  = (GRFNPARMTREE)func_parm;
    global = total->global;
    local  = total->local;

    /* Precalculate DFF, gamma and 0.5 * gamma * gamma */

    if (total->num_dom_df > 0 && total->dom_idx != -1)
    {
        dom_dff  = dvector(0, total->num_dom_df - 1);
        dom_gam  = dvector(0, total->num_dom_df - 1);
        dom_gam2 = dvector(0, total->num_dom_df - 1);

        if (!dom_dff || !dom_gam || !dom_gam2)
        {
            err = "Memory allocation error (1) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < total->num_dom_df; i++)
        {
            dom_dff[i]  = swp_f_df(evt_date, total->dom_df_dts[i], (char*)dom_yc);
            dom_gam[i]  = (1.0 - exp(-domlam * total->dom_df_tms[i])) / domlam;
            dom_gam2[i] = 0.5 * dom_gam[i] * dom_gam[i] * domphi;
        }

        do_dom = 1;
    }
    else
    {
        do_dom = 0;
    }

    if (total->num_fx_df > 0 && total->fx_idx != -1)
    {
        fx_dff  = dvector(0, total->num_fx_df - 1);
        fx_gam  = dvector(0, total->num_fx_df - 1);
        fx_gam2 = dvector(0, total->num_fx_df - 1);

        if (!fx_dff || !fx_gam || !fx_gam2)
        {
            err = "Memory allocation error (2) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < total->num_fx_df; i++)
        {
            fx_dff[i]  = swp_f_df(evt_date, total->fx_df_dts[i], (char*)dom_yc);
            fx_gam[i]  = (1.0 - exp(-domlam * total->fx_df_tms[i])) / domlam;
            fx_gam2[i] = 0.5 * fx_gam[i] * fx_gam[i] * domphi;
        }

        do_fx = 1;
    }
    else
    {
        do_fx = 0;
    }

    if (total->num_for_df > 0 && total->for_idx != -1)
    {
        for_dff  = dvector(0, total->num_for_df - 1);
        for_gam  = dvector(0, total->num_for_df - 1);
        for_gam2 = dvector(0, total->num_for_df - 1);

        if (!for_dff || !for_gam || !for_gam2)
        {
            err = "Memory allocation error (3) in payoff_doublelgm1fquanto_pde";
            goto FREE_RETURN;
        }

        for (i = 0; i < total->num_for_df; i++)
        {
            for_dff[i]  = swp_f_df(evt_date, total->for_df_dts[i], (char*)for_yc);
            for_gam[i]  = (1.0 - exp(-forlam * total->for_df_tms[i])) / forlam;
            for_gam2[i] = 0.5 * for_gam[i] * for_gam[i] * forphi;
        }

        do_for = 1;
    }
    else
    {
        do_for = 0;
    }

    /* Eval payoff */
    for (i = l1; i < u1; i++)
    {
        for (j = l2; j < u2; j++)
        {
            local->smp.und[total->fx_idx].sv[SPOT] = -999999;

            if (do_dom)
            {
                for (l = 0; l < total->num_dom_df; l++)
                {
                    local->evt->df[total->dom_idx][l] =
                        dom_dff[l] * exp(-dom_gam[l] * r1[i] - dom_gam2[l]);
                }
            }
            if (do_fx)
            {
                for (l = 0; l < total->num_fx_df; l++)
                {
                    local->evt->df[total->fx_idx][l] =
                        fx_dff[l] * exp(-fx_gam[l] * r1[i] - fx_gam2[l]);
                }
            }
            if (do_for)
            {
                for (l = 0; l < total->num_for_df; l++)
                {
                    local->evt->df[total->for_idx][l] =
                        for_dff[l] * exp(-for_gam[l] * r2[i][j] - for_gam2[l]);
                }
            }
            err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val[i][j], &temp);

            if (err)
            {
                goto FREE_RETURN;
            }
        }
    }

FREE_RETURN:

    if (dom_dff)
    {
        free_dvector(dom_dff, 0, total->num_dom_df - 1);
    }

    if (dom_gam)
    {
        free_dvector(dom_gam, 0, total->num_dom_df - 1);
    }

    if (dom_gam2)
    {
        free_dvector(dom_gam2, 0, total->num_dom_df - 1);
    }

    if (fx_dff)
    {
        free_dvector(fx_dff, 0, total->num_fx_df - 1);
    }

    if (fx_gam)
    {
        free_dvector(fx_gam, 0, total->num_fx_df - 1);
    }

    if (fx_gam2)
    {
        free_dvector(fx_gam2, 0, total->num_fx_df - 1);
    }

    if (for_dff)
    {
        free_dvector(for_dff, 0, total->num_for_df - 1);
    }

    if (for_gam)
    {
        free_dvector(for_gam, 0, total->num_for_df - 1);
    }

    if (for_gam2)
    {
        free_dvector(for_gam2, 0, total->num_for_df - 1);
    }

    return err;
}
