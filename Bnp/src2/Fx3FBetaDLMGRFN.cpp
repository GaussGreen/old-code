
#include "Fx3FBetaDLMGRFN.h"

#include "math.h"
#include "srt_h_all.h"

Err grfn_payoff_4_3dfxBetaDLM_mc(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double Xdom,
    double Yfor,
    double Zfx,
    /* Results */
    int     num_col,
    double* res,
    int*    stop_path)
{
    GRFNPARMBETAMCDLM total;
    GRFNCOMMSTRUCT    global;
    FIRSTMktAtT*      local;
    int               l;
    double            temp;
    Err               err;

    if (func_parm == NULL)
    {
        memset(res, 0, num_col * sizeof(double));
        *stop_path = 0;
        return NULL;
    }
    else
    {
        err = NULL;
    }

    memset(res, 0, num_col * sizeof(double));

    /* Get the event */
    total  = (GRFNPARMBETAMCDLM)func_parm;
    global = total->global;
    local  = total->local;

    /* Calc market data */
    local->smp.und[total->fx_idx].sv[SPOT] =
        exp(total->fx_coef_const + (total->fx_coef_lin + total->fx_coef_quad * Zfx) * Zfx +
            total->fx_coef_dom * Xdom + total->fx_coef_for * Yfor);

    if (total->do_dom)
    {
        for (l = 0; l < total->num_dom_df; l++)
        {
            local->evt->df[total->dom_idx][l] =
                exp(total->df_dom_coef_const[l] + total->df_dom_coef_lin[l] * Xdom);
        }
    }
    if (total->do_fx)
    {
        for (l = 0; l < total->num_fx_df; l++)
        {
            local->evt->df[total->fx_idx][l] =
                exp(total->df_fx_coef_const[l] + total->df_fx_coef_lin[l] * Xdom);
        }
    }
    if (total->do_for)
    {
        for (l = 0; l < total->num_for_df; l++)
        {
            local->evt->df[total->for_idx][l] =
                exp(total->df_for_coef_const[l] + total->df_for_coef_lin[l] * Yfor);
        }
    }

    err = FIRSTEvalEvent(global, local, num_col, 2, NULL, NULL, res, &temp);

    *stop_path = 0;

    return err;
}

Err grfn_payoff_4_3dfxBetaDLMQBeta_Tree(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Market data */
    double spot_fx,
    char*  dom_yc,
    char*  for_yc,
    /* Nodes data */
    long n1,
    long n2,
    long n3,
    /* i: d1, j: d2, k: d3, l = {0: xDom, 1: xFor, 2: log (Fx/Fx0)} */
    double**** sv,
    /* Vector of results to be updated */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl)
{
    GRFNPARMBETATREEDLM total;
    GRFNCOMMSTRUCT      global;
    FIRSTMktAtT*        local;
    double *dom_const = NULL, *dom_lin = NULL, *fx_const = NULL, *fx_lin = NULL, *for_const = NULL,
           *for_lin = NULL;

    double ***     svi, **svij, *svijk;
    int            i, j, k, l;
    double         temp;
    double         dom_phi, dom_lam, for_phi, for_lam;
    unsigned short do_dom, do_for, do_fx;
    double         coef_A, coef_B, coef_C;
    double         coef1_S, coef2_S;
    double         temp_X, const_X, coef1_X, coef2_X, fwd_fx;

    Err err = NULL;

    /* Get the event */
    total  = (GRFNPARMBETATREEDLM)func_parm;
    global = total->global;
    local  = total->local;

    /* Precalculate DFF, gamma and 0.5 * gamma * gamma */

    dom_phi = total->phi_dom;
    for_phi = total->phi_for;
    dom_lam = total->lam_dom;
    for_lam = total->lam_for;

    if (total->num_dom_df > 0 && total->dom_idx != -1)
    {
        dom_const = dvector(0, total->num_dom_df - 1);
        dom_lin   = dvector(0, total->num_dom_df - 1);

        if (!dom_const || !dom_lin)
        {
            err = "Memory allocation error (1) in grfn_payoff_4_3dfx_tree";
            goto FREE_RETURN;
        }

        for (i = 0; i < total->num_dom_df; i++)
        {
            dom_const[i] = log(swp_f_df(evt_date, total->dom_df_dts[i], (char*)dom_yc));
            dom_lin[i]   = -(1.0 - exp(-dom_lam * total->dom_df_tms[i])) / dom_lam;
            dom_const[i] -= 0.5 * dom_lin[i] * dom_lin[i] * dom_phi;
        }

        do_dom = 1;
    }
    else
    {
        do_dom = 0;
    }

    if (total->num_fx_df > 0 && total->fx_idx != -1)
    {
        fx_const = dvector(0, total->num_fx_df - 1);
        fx_lin   = dvector(0, total->num_fx_df - 1);

        if (!fx_const || !fx_lin)
        {
            err = "Memory allocation error (2) in grfn_payoff_4_3dfx_tree";
            goto FREE_RETURN;
        }

        for (i = 0; i < total->num_fx_df; i++)
        {
            fx_const[i] = log(swp_f_df(evt_date, total->fx_df_dts[i], (char*)dom_yc));
            fx_lin[i]   = -(1.0 - exp(-dom_lam * total->fx_df_tms[i])) / dom_lam;
            fx_const[i] -= 0.5 * fx_lin[i] * fx_lin[i] * dom_phi;
        }

        do_fx = 1;
    }
    else
    {
        do_fx = 0;
    }

    if (total->num_for_df > 0 && total->for_idx != -1)
    {
        for_const = dvector(0, total->num_for_df - 1);
        for_lin   = dvector(0, total->num_for_df - 1);

        if (!for_const || !for_lin)
        {
            err = "Memory allocation error (3) in grfn_payoff_4_3dfx_tree";
            goto FREE_RETURN;
        }

        for (i = 0; i < total->num_for_df; i++)
        {
            for_const[i] = log(swp_f_df(evt_date, total->for_df_dts[i], (char*)for_yc));
            for_lin[i]   = -(1.0 - exp(-for_lam * total->for_df_tms[i])) / for_lam;
            for_const[i] -= 0.5 * for_lin[i] * for_lin[i] * for_phi;
        }

        do_for = 1;
    }
    else
    {
        do_for = 0;
    }

    /* spot Fx reconstruction */
    coef_B = 1.0 / (1.0 + 2.0 * total->C0 * total->varFFX);
    coef_C = total->C0 * coef_B;
    coef_A = -0.5 * (total->varFFX * coef_B * total->B0 * total->B0 - log(coef_B));
    coef_B *= total->B0;

    const_X = 0.5 * (total->varFFX + total->beta_dom * total->beta_dom * total->phi_dom -
                     total->beta_for * total->beta_for * total->phi_for);
    coef1_X = total->beta_dom;
    coef2_X = -total->beta_for;

    fwd_fx = spot_fx * swp_f_df(total->today, evt_date, for_yc) /
             swp_f_df(total->today, evt_date, dom_yc);

    coef_A += log(fwd_fx);
    coef_A += -0.5 * (total->beta_dom * total->beta_dom * total->phi_dom -
                      total->beta_for * total->beta_for * total->phi_for);

    coef1_S = -coef1_X;
    coef2_S = -coef2_X;

    /* Eval payoff */
    for (i = 0; i < n1; i++)
    {
        svi = sv[i];
        for (j = 0; j < n2; j++)
        {
            svij = svi[j];
            for (k = 0; k < n3; k++)
            {
                svijk = svij[k];

                temp_X = const_X + coef1_X * svijk[0] + coef2_X * svijk[1] + svijk[2];

                local->smp.und[total->fx_idx].sv[SPOT] =
                    exp(coef_A + temp_X * (coef_B + temp_X * coef_C) + coef1_S * svijk[0] +
                        coef2_S * svijk[1]);

                if (do_dom)
                {
                    for (l = 0; l < total->num_dom_df; l++)
                    {
                        local->evt->df[total->dom_idx][l] =
                            exp(dom_const[l] + dom_lin[l] * svijk[0]);
                    }
                }
                if (do_fx)
                {
                    for (l = 0; l < total->num_fx_df; l++)
                    {
                        local->evt->df[total->fx_idx][l] = exp(fx_const[l] + fx_lin[l] * svijk[0]);
                    }
                }
                if (do_for)
                {
                    for (l = 0; l < total->num_for_df; l++)
                    {
                        local->evt->df[total->for_idx][l] =
                            exp(for_const[l] + for_lin[l] * svijk[1]);
                    }
                }

                err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val[i][j][k], &temp);

                /* temp */
                /*
                prod_val[i][j][k][0] = svijk[0];
                prod_val[i][j][k][1] = svijk[0] * svijk[0];
                prod_val[i][j][k][2] = svijk[1];
                prod_val[i][j][k][3] = svijk[1] * svijk[1];
                prod_val[i][j][k][4] = svijk[2];
                prod_val[i][j][k][5] = svijk[2] * svijk[2];
                */

                if (err)
                {
                    goto FREE_RETURN;
                }
            }
        }
    }

FREE_RETURN:

    if (dom_const)
        free_dvector(dom_const, 0, total->num_dom_df - 1);
    if (dom_lin)
        free_dvector(dom_lin, 0, total->num_dom_df - 1);
    if (fx_const)
        free_dvector(fx_const, 0, total->num_fx_df - 1);
    if (fx_lin)
        free_dvector(fx_lin, 0, total->num_fx_df - 1);
    if (for_const)
        free_dvector(for_const, 0, total->num_for_df - 1);
    if (for_lin)
        free_dvector(for_lin, 0, total->num_for_df - 1);

    return err;
}

Err grfn_payoff_4_3dfxBetaDLMQTStar_Tree(
    /* Event */
    double evt_date,
    double evt_time,
    void*  func_parm,
    /* Nodes data */
    long       n1,
    long       n2,
    long       n3,
    double**** sv,
    /* Results */
    long       nprod,
    double**** prod_val,
    /* Barrier details */
    int    is_bar,
    int    bar_k,
    int**  bar_idx,
    int    bar_col,
    double bar_lvl)
{
    GRFNPARMBETAMCDLM total;
    GRFNCOMMSTRUCT    global;
    FIRSTMktAtT*      local;
    int               i, j, k, l;
    double ***        svi, **svij, *svijk;
    double ***        prod_vali, **prod_valij, *prod_valijk;
    double            Xdom, Yfor, Zfx, df;
    double            temp;
    Err               err;

    /* Get the event */
    total  = (GRFNPARMBETAMCDLM)func_parm;
    global = total->global;
    local  = total->local;

    for (i = 0; i < n1; i++)
    {
        svi       = sv[i];
        prod_vali = prod_val[i];
        for (j = 0; j < n2; j++)
        {
            svij       = svi[j];
            prod_valij = prod_vali[j];
            for (k = 0; k < n3; k++)
            {
                svijk       = svij[k];
                prod_valijk = prod_valij[k];

                Xdom = svijk[0];
                Yfor = svijk[1];
                Zfx  = svijk[2];

                /* Calc Bond */
                df = exp(total->bond_pay_const + total->bond_pay_lin * Xdom);

                /* Calc Spot Fx */
                local->smp.und[total->fx_idx].sv[SPOT] = exp(
                    total->fx_coef_const + (total->fx_coef_lin + total->fx_coef_quad * Zfx) * Zfx +
                    total->fx_coef_dom * Xdom + total->fx_coef_for * Yfor);

                /* Calc all the DF */
                if (total->do_dom)
                {
                    for (l = 0; l < total->num_dom_df; l++)
                    {
                        local->evt->df[total->dom_idx][l] =
                            exp(total->df_dom_coef_const[l] + total->df_dom_coef_lin[l] * Xdom);
                    }
                }
                if (total->do_fx)
                {
                    for (l = 0; l < total->num_fx_df; l++)
                    {
                        local->evt->df[total->fx_idx][l] =
                            exp(total->df_fx_coef_const[l] + total->df_fx_coef_lin[l] * Xdom);
                    }
                }
                if (total->do_for)
                {
                    for (l = 0; l < total->num_for_df; l++)
                    {
                        local->evt->df[total->for_idx][l] =
                            exp(total->df_for_coef_const[l] + total->df_for_coef_lin[l] * Yfor);
                    }
                }

                /* Pretreatment */
                for (l = 0; l < nprod; l++)
                {
                    prod_valijk[l] *= df;
                }

                err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_valijk, &temp);

                /* Postreatment */
                for (l = 0; l < nprod; l++)
                {
                    prod_valijk[l] /= df;
                }
            }
        }
    }

    return err;
}