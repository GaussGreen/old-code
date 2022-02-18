/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_ts_fx_fct.c

   PURPOSE:     Function to work with an FX TermStruct:
   ------------------------------------------------------------------------ */

#define SWAP(a, b)     \
    {                  \
        double tempr;  \
        tempr = (a);   \
        (a)   = (b);   \
        (b)   = tempr; \
    }

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ts.h"
#include "srt_h_ts_fx.h"

/* ------------------------------------------------------------------------ */

/* Finds the Value of the loval volatility at one given point in time */

double find_fx_sig(double time, TermStruct* l)
{
    SrtLst* ls;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;

    if (ls != NULL)
    {
        return ((FxTermStructVal*)ls->element->val.pval)->sigx;
    }
    else
    {
        return ((FxTermStructVal*)l->tail->element->val.pval)->sigx;
    }
}

/* ------------------------------------------------------------------------- */

double fx_cum_vol_func(double time, SRT_Boolean sigsq, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal* tsval;
    double           val, cum_vol_time = 0.0, prev_time;

    ls = l->head;

    if (l->head == l->tail)
    {
        val = ((FxTermStructVal*)ls->element->val.pval)->sigx;
        if (sigsq == SRT_YES)
            val *= val;
        return (val * time);
    }
    /* first case : time is before the first sigma bucket time */
    if (((FxTermStructVal*)ls->element->val.pval)->time >= time)
    {
        prev_time = 0.0;
        val       = ((FxTermStructVal*)ls->element->val.pval)->sigx;
    }
    /* general case */
    else
    {
        while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (FxTermStructVal*)l->tail->element->val.pval;
            val   = tsval->sigx;
        }
        else
        {
            tsval = (FxTermStructVal*)ls->previous->element->val.pval;
            val   = ((FxTermStructVal*)ls->element->val.pval)->sigx;
        }

        prev_time    = tsval->time;
        cum_vol_time = (sigsq == SRT_YES ? tsval->int_sig2_dt : tsval->int_sig_dt);
    }
    if (sigsq == SRT_YES)
        val *= val;
    cum_vol_time += val * (time - prev_time);

    return cum_vol_time;
}

/* V_dx def: int(s = 0, s = t) [pdx(s)*s_x(s)*s_dom(s)/F_dom(s)] ds */
double V_dx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           temp;
    double           dom_sig, dom_tau, dom_F;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    temp    = (exp(time / dom_tau) - 1) * dom_tau;

    result = tsval->V_dx + tsval->rhodx * tsval->sigx * (dom_sig / dom_F) * temp;

    return result;
}

/* W_dx def: int(s = 0, s = t) [pdx(s)*s_x(s)*s_dom(s)*psi_dom(s)/F_dom(s)] ds */
double W_dx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           dom_sig, dom_tau, dom_F, dom_Psi;
    double           temp;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    dom_Psi = tsval->StochRates_Ts[0].Psi;

    temp = (exp(time / dom_tau) - 1) * dom_tau;

    result = tsval->W_dx + tsval->rhodx * tsval->sigx * dom_sig * (dom_Psi / dom_F) * temp +
             tsval->rhodx * tsval->sigx * dom_sig * dom_tau * (temp - time);

    return result;
}

/*O_fd def: int(s = 0, s = t) [pfd(s)*s_dom(s)*s_for(s)/(F_dom(s)*F_for(s))] ds */
double O_fd_func(double time, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           dom_sig, dom_tau, dom_F, for_sig, for_tau, for_F;
    double           temp, temp1, temp2, temp3, lambda;

    result = 0.0;
    ls     = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_F   = tsval->StochRates_Ts[0].F;
    dom_sig = tsval->StochRates_Ts[0].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    for_F   = tsval->StochRates_Ts[1].F;
    for_sig = tsval->StochRates_Ts[1].sig;
    for_tau = tsval->StochRates_Ts[1].tau;
    temp    = dom_sig / dom_F;
    temp1   = for_sig / for_F;
    temp2   = temp1 * temp;
    lambda  = 1.0 / dom_tau + 1.0 / for_tau;
    temp3   = (exp(time * lambda) - 1) / lambda;

    result = tsval->O_fd + tsval->rhofd * temp2 * temp3;

    return result;
}

/* --------------------------------------------------------------------------------------- */

double Phi_func(double time, TermStruct* l)
{
    double dValG, dValH;

    G_H_func(time, l, &dValG, &dValH);

    return dValG * F_func(time, l) * F_func(time, l);
}

/* int ( s = 0, s = t) [pfd(s)*s_ifr_dom(s,t)*s_ifr_for(s,t)] ds = F_dom(t)*F_for(t)*O_fd(t) */
double Phi_fd_func(double time, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p; /* -p for previous */
    double           result;
    double           temp, temp1;

    result = 0.0;
    temp   = 0.0;
    temp1  = 0.0;
    ls     = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    temp = tsval->rhofd * tsval->StochRates_Ts[1].sig * tsval->StochRates_Ts[0].sig *
           (exp((tsval->StochRates_Ts[0].lambda + tsval->StochRates_Ts[1].lambda) * time) - 1) *
           1.0 / (tsval->StochRates_Ts[0].lambda + tsval->StochRates_Ts[1].lambda);

    temp1 = exp(-time * (tsval->StochRates_Ts[0].lambda + tsval->StochRates_Ts[1].lambda));

    result = tsval->Phi_fd * temp1 + temp1 * temp;

    return result;
}

/* --------------------------------------------------------------------------------------- */

double H_fd_func(double time, TermStruct* domts, TermStruct* forts, TermStruct* l)
{
    double result;
    result = -Psi_func(time, forts) * F_func(time, forts) * O_fd_func(time, l) +
             P_fd_func(time, l) * F_func(time, forts);

    return result;
}

/* M_fx def: int(s = 0, s = time) [pfx*sx(sf(s)/Ff(s))ds] */
double M_fx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           temp;
    double           result;
    double           sig, F, tau;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    sig  = tsval->StochRates_Ts[1].sig;
    F    = tsval->StochRates_Ts[1].F;
    tau  = tsval->StochRates_Ts[1].tau;
    temp = (exp(time / tau) - 1) * tau;

    result = tsval->M_fx + tsval->rhofx * tsval->sigx * (sig / F) * temp;
    return result;
}

double N_fx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           temp;
    double           for_F, for_tau, for_sig;
    double           result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    for_tau = tsval->StochRates_Ts[1].tau;
    for_sig = tsval->StochRates_Ts[1].sig;
    for_F   = tsval->StochRates_Ts[1].F;
    temp    = (1 - exp(-time / for_tau)) * for_tau;

    result = tsval->N_fx + tsval->M_fx * for_F * temp +
             tsval->rhofx * tsval->sigx * for_sig * for_tau * (time - temp);

    return result;
}

double P_fd_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           temp, temp1;
    double           temp2;
    double           dom_sig, dom_tau, dom_F, forPsi, for_sig, for_tau, for_F, lambda;
    double           result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_F   = tsval->StochRates_Ts[1].F;
    forPsi  = tsval->StochRates_Ts[1].Psi;
    temp    = (dom_sig / dom_F) * (for_sig / for_F);
    lambda  = (1.0 / dom_tau) + (1.0 / for_tau);
    temp1   = (exp(time * lambda) - 1) / lambda;
    temp2   = (exp(time / dom_tau) - 1) * dom_tau;

    result = tsval->P_fd + tsval->rhofd * forPsi * temp * temp1 +
             tsval->rhofd * for_F * temp * (temp1 - temp2) * for_tau;
    return result;
}

/* R_fd_func: int (s = 0, s = t) [rho_fd(s)*s_dom(s)*s_for(s)*Psi_dom(s)/(F_for(s)*F_dom(s))] ds */
double R_fd_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           temp, temp1;
    double           temp2;
    double           dom_sig, dom_tau, dom_F, dom_Psi, for_sig, for_tau, for_F, lambda;
    double           result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_F   = tsval->StochRates_Ts[1].F;
    dom_Psi = tsval->StochRates_Ts[0].Psi;
    temp    = (for_sig / for_F) * (dom_sig / dom_F);
    lambda  = (1.0 / dom_tau) + (1.0 / for_tau);
    temp1   = (exp(time * lambda) - 1) / lambda;
    temp2   = (exp(time / for_tau) - 1) * for_tau;

    result = tsval->R_fd + tsval->rhofd * dom_Psi * temp * temp1 +
             tsval->rhofd * dom_F * temp * (temp1 - temp2) * dom_tau;
    return result;
}

double M_fd_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           dom_sig, for_sig, dom_tau, for_tau, dom_F, for_F, dom_Psi, forPsi, lambda;
    double           temp, temp1, temp2, temp3;
    double           result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_F   = tsval->StochRates_Ts[1].F;
    dom_Psi = tsval->StochRates_Ts[0].Psi;
    forPsi  = tsval->StochRates_Ts[1].Psi;
    temp    = dom_sig * for_sig;
    lambda  = (1.0 / dom_tau) + (1.0 / for_tau);
    temp1   = (exp(time / dom_tau) - 1) * dom_tau;
    temp2   = (exp(time / for_tau) - 1) * for_tau;
    temp3   = (exp(time * lambda) - 1) / lambda;

    result = tsval->M_fd + tsval->rhofd * (temp / (dom_F * for_F)) * dom_Psi * forPsi * temp3 +
             tsval->rhofd * (temp / dom_F) * for_tau * dom_Psi * (temp3 - temp1) +
             tsval->rhofd * (temp / for_F) * dom_tau * forPsi * (temp3 - temp2) +
             tsval->rhofd * temp * dom_tau * for_tau * (temp3 + time - temp1 - temp2);

    return result;
}

double Q_fd_func(double time, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           lambda;
    double           temp, temp1, temp2;
    double           dom_tau, for_tau, dom_F, for_F, dom_Psi, forPsi, dom_sig, for_sig;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_F   = tsval->StochRates_Ts[1].F;
    dom_Psi = tsval->StochRates_Ts[0].Psi;
    forPsi  = tsval->StochRates_Ts[1].Psi;
    lambda  = 1 / dom_tau + 1 / for_tau;
    temp    = (1 - exp(-time / for_tau)) * for_tau;
    temp1   = (1 - exp(-time * lambda)) / lambda;
    temp2   = (exp(time / dom_tau) - 1.0) * dom_tau;

    result = tsval->Q_fd - for_F * tsval->O_fd * dom_Psi * temp -
             for_F * dom_F * tsval->O_fd * dom_tau * (temp - temp1) -
             tsval->rhofd * dom_sig * for_sig * dom_Psi / (dom_F * lambda) * (temp2 - temp) -
             tsval->rhofd * dom_sig * for_sig * dom_tau / lambda * (temp2 - temp - time + temp1) +
             for_F * tsval->P_fd * temp +
             tsval->rhofd * dom_Psi * dom_sig * for_sig / (lambda * dom_F) * (temp2 - temp) +
             tsval->rhofd * dom_sig * for_sig * dom_tau *
                 (temp2 / lambda - temp / lambda - time * for_tau + temp * for_tau);

    return result;
}

double S_fd_func(double time, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           lambda;
    double           temp, temp1, temp2;
    double           dom_sig, for_sig, dom_tau, for_tau, dom_F, for_F;
    double           result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_sig = tsval->StochRates_Ts[1].sig;
    for_F   = tsval->StochRates_Ts[1].F;
    temp    = dom_sig / dom_F;
    temp1   = for_sig / for_F;
    lambda  = (1 / dom_tau + 1 / for_tau);
    temp2   = (exp(lambda * time) - 1) / lambda;

    result = tsval->S_fd + tsval->rhofd * temp * temp1 * temp2;
    return result;
}

double T_fd_func(double time, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           temp, temp1, temp2, temp3;
    double           dom_sig, for_sig, dom_tau, for_tau, dom_F, for_F;
    double           lambda;
    double           result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }
    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_F   = tsval->StochRates_Ts[1].F;
    temp    = (1 - exp(-time / for_tau)) * for_tau;
    lambda  = 1 / dom_tau + 1 / for_tau;
    temp1   = (dom_sig * for_sig) / (dom_F * for_F);

    temp2 = (exp(time / dom_tau) - 1) * dom_tau / lambda;
    temp3 = (1 - exp(-time / for_tau)) * for_tau / lambda;

    result =
        tsval->T_fd + for_F * tsval->S_fd * temp + tsval->rhofd * for_F * temp1 * (temp2 - temp3);

    return result;
}

/* U_fd_func: int (s = 0, s = t) [psi_dom(s)*phi_df(s)/F_dom(s)] ds */
double U_fd_func(double time, TermStruct* l)
{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double dom_tau, for_tau, dom_sig, for_sig, dom_F, for_F, dom_Psi, temp, temp1, temp2, lambda;
    double result;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_tau = tsval->StochRates_Ts[0].tau;
    for_tau = tsval->StochRates_Ts[1].tau;
    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_F   = tsval->StochRates_Ts[0].F;
    for_F   = tsval->StochRates_Ts[1].F;
    dom_Psi = tsval->StochRates_Ts[0].Psi;
    temp    = (1 - exp(-time / for_tau)) * for_tau;
    temp1   = (exp(time / dom_tau) - 1) * dom_tau;
    lambda  = 1 / dom_tau + 1 / for_tau;
    temp2   = (1 - exp(-lambda * time)) / lambda;

    result = tsval->U_fd + for_F * dom_Psi * tsval->O_fd * temp +
             tsval->O_fd * dom_F * for_F * dom_tau * (temp - temp2) +
             dom_Psi * (tsval->rhofd / lambda) * (dom_sig * for_sig / (dom_F)) * (temp1 - temp) +
             tsval->rhofd * (dom_sig * for_sig * dom_tau / lambda) * (temp1 - temp - time + temp2);

    return result;
}

/* V_fd def: int(s=0, s = t) [ p_fd(s)*s_dom(s)*s_for(s)/F_dom(s)] ds */
double V_fd_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result, dom_tmp;
    double           dom_sig, for_sig, dom_tau, dom_F;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    dom_tmp = (exp(time / dom_tau) - 1) * dom_tau;

    result = tsval->V_fd + (tsval->rhofd * dom_sig * for_sig / dom_F) * dom_tmp;

    return result;
}

/* W_fd def: int(s=0, s = t) [ p_fd(s)*s_dom(s)*s_for(s)*Psi_dom(s)/F_dom(s)] ds */
double W_fd_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result, dTemp;
    double           dom_sig, for_sig, dom_tau, dom_F, dom_psi;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    for_sig = tsval->StochRates_Ts[1].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    dom_psi = tsval->StochRates_Ts[0].Psi;
    dTemp   = (exp(time / dom_tau) - 1) * dom_tau;

    result = tsval->W_fd + (tsval->rhofd * dom_sig * for_sig / dom_F) * dom_psi * dTemp +
             (tsval->rhofd * dom_sig * for_sig * dom_tau) * (dTemp - time);

    return result;
}

double X_dx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           temp;
    double           dom_sig, dom_tau, dom_F, SqdomF;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    SqdomF  = dom_F * dom_F;
    temp    = 0.5 * (exp(2 * time / dom_tau) - 1) * dom_tau;

    result = tsval->X_dx + tsval->rhodx * tsval->sigx * (dom_sig / SqdomF) * temp;

    return result;
}

double Y_dx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           dom_sig, dom_tau, dom_F, dom_Psi, SqdomF;
    double           temp, temp1;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    dom_sig = tsval->StochRates_Ts[0].sig;
    dom_tau = tsval->StochRates_Ts[0].tau;
    dom_F   = tsval->StochRates_Ts[0].F;
    SqdomF  = dom_F * dom_F;
    dom_Psi = tsval->StochRates_Ts[0].Psi;

    temp  = (exp(time / dom_tau) - 1) * dom_tau;
    temp1 = 0.5 * (exp(2 * time / dom_tau) - 1) * dom_tau;

    result = tsval->Y_dx + tsval->rhodx * tsval->sigx * dom_sig * (dom_Psi / SqdomF) * temp1 +
             tsval->rhodx * tsval->sigx * dom_sig * dom_tau * (temp1 - temp);

    return result;
}

/*---------------------------------------------------------------------------*/

double X_fx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           temp;
    double           for_sig, for_tau, for_F, sqforF;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    for_sig = tsval->StochRates_Ts[1].sig;
    for_tau = tsval->StochRates_Ts[1].tau;
    for_F   = tsval->StochRates_Ts[1].F;
    sqforF  = for_F * for_F;
    temp    = 0.5 * (exp(2 * time / for_tau) - 1) * for_tau;

    result = tsval->X_fx + tsval->rhofx * tsval->sigx * (for_sig / sqforF) * temp;

    return result;
}

/*----------------------------------------------------------------------------*/

double Y_fx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           for_sig, for_tau, for_F, forPsi, sqforF;
    double           temp, temp1;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    for_sig = tsval->StochRates_Ts[1].sig;
    for_tau = tsval->StochRates_Ts[1].tau;
    for_F   = tsval->StochRates_Ts[1].F;
    sqforF  = for_F * for_F;
    forPsi  = tsval->StochRates_Ts[1].Psi;

    temp  = (exp(time / for_tau) - 1) * for_tau;
    temp1 = 0.5 * (exp(2 * time / for_tau) - 1) * for_tau;

    result = tsval->Y_fx + tsval->rhofx * tsval->sigx * for_sig * (forPsi / sqforF) * temp1 +
             tsval->rhofx * tsval->sigx * for_sig * for_tau * (temp1 - temp);

    return result;
}

/* --------------------------------------------------------------------------------------- */

double V_fx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           temp;
    double           for_sig, for_tau, for_F;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    for_sig = tsval->StochRates_Ts[1].sig;
    for_tau = tsval->StochRates_Ts[1].tau;
    for_F   = tsval->StochRates_Ts[1].F;
    temp    = (exp(time / for_tau) - 1) * for_tau;

    result = tsval->V_fx + tsval->rhofx * tsval->sigx * (for_sig / for_F) * temp;

    return result;
}

/* W_fx_func def: int(s = 0, s = t) [pfx(s)*sx(s)*sf(s)*psif(s)/F(s)] ds */
double W_fx_func(double time, TermStruct* l)

{
    SrtLst*          ls;
    FxTermStructVal *tsval, *tsval_p;
    double           result;
    double           for_sig, for_tau, for_F, forPsi;
    double           temp;
    result = 0.0;

    ls = l->head;

    while ((ls != NULL) && (((FxTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (FxTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (FxTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    for_sig = tsval->StochRates_Ts[1].sig;
    for_tau = tsval->StochRates_Ts[1].tau;
    for_F   = tsval->StochRates_Ts[1].F;
    forPsi  = tsval->StochRates_Ts[1].Psi;

    temp = (exp(time / for_tau) - 1) * for_tau;

    result = tsval->W_fx + tsval->rhofx * tsval->sigx * for_sig * (forPsi / for_F) * temp +
             tsval->rhofx * tsval->sigx * for_sig * for_tau * (temp - time);

    return result;
}

/* -------------------------------------------------------------------------------------------
   This Function computes the BS Implied Volatility of an FX Option in a stochastic world
   framework. The Domestic Underlying and the Foreign Underlying has to follow a LGM dynamic.

   The domestic, the foreign and the FX Underlying must be first initialised
   ------------------------------------------------------------------------------------------- */

Err srt_f_get_fx_implied_vol(double yr_to_exp, String fx_und_names, double* fx_bs_vol)
{
    double            dom_bd_vol, for_bd_vol;
    double            dom_G, dom_H, dom_L, dom_O, dom_psi;
    double            sq_dom_psi;
    double            for_G, for_H, for_L, for_O, for_psi;
    double            sq_for_psi;
    double            Vdx, Wdx, Vfx, Wfx;
    double            Ofd, Rfd, Pfd, Mfd;
    double            spot_fx_cum_vol, sq_fx_bs_vol;
    double            cov_dom_int_fx, cov_for_int_fx, cov_dom_int_for_int;
    TermStruct *      fx_ts, *dom_ts, *for_ts;
    Err               err;
    String            dom_ir_und_name;
    String            for_ir_und_name;
    SrtUndPtr         fx_und;
    SrtUndPtr         dom_und;
    SrtUndPtr         for_und;
    SrtModelType      model_type;
    SrtModelType      dom_model_type;
    SrtModelType      for_model_type;
    SrtUnderlyingType und_type;

    /* Gets the Fx Underlying */
    fx_und = lookup_und(fx_und_names);
    if (fx_und == NULL)
        return serror("Undefined Underlying %s", fx_und_names);

    /* Gets the Underlying Type and Checks it is an FX */
    und_type = get_underlying_type(fx_und);
    if (und_type != FOREX_UND)
        return serror("Need an Fx Underlying in srt_f_get_fx_implied_vol !");

    /* Gets the Model Name Type */
    model_type = get_mdltype_from_fxund(fx_und);

    /* Gets the Fx term Structure */
    err = get_underlying_ts(fx_und, &fx_ts);
    if (err)
        return serror("Can not get the fx term structure in srt_f_get_fx_implied_vol ");

    /* If only a Standard Fx Stoch Rates model, return the Sqrt (cumVol / T) */
    if (model_type != FX_STOCH_RATES)
    {
        spot_fx_cum_vol = fx_cum_vol_func(yr_to_exp, SRT_YES, fx_ts);
        if (yr_to_exp > 0.0)
            *fx_bs_vol = sqrt(spot_fx_cum_vol / yr_to_exp);
        else
            *fx_bs_vol = find_fx_sig(0.0, fx_ts);

        return NULL;
    }

    /* We NOW have an FX STOCH RATES model */

    /* Gets the Domestic Interest Rate Underlying */
    dom_ir_und_name = get_domname_from_fxund(fx_und);
    dom_und         = lookup_und(dom_ir_und_name);
    if (!dom_und)
        return serror("Underlying %s not defined", dom_ir_und_name);

    /* Gets the Foreign Interest Rate Underlying */
    for_ir_und_name = get_forname_from_fxund(fx_und);
    for_und         = lookup_und(for_ir_und_name);
    if (!for_und)
        return serror("Underlying %s not defined", for_ir_und_name);

    /* verify that we are in the lgm/lgm(vasicek) configuration */
    dom_model_type = get_mdltype_from_irund(dom_und);
    for_model_type = get_mdltype_from_irund(for_und);

    if ((dom_model_type != LGM) && (dom_model_type != VASICEK))
        return serror("srt_f_get_fx_implied_vol will only work for lgm/lgm ");

    if ((for_model_type != LGM) && (for_model_type != VASICEK))
        return serror("srt_f_get_fx_implied_vol will only work for lgm/lgm ");

    /* Get the dom & for & fx ts */
    err = get_underlying_ts(dom_und, &dom_ts);
    if (err)
        return serror("Can not get the domestic term structure in srt_f_get_fx_implied_vol");

    err = get_underlying_ts(for_und, &for_ts);
    if (err)
        return serror("Can not get the foreign term structure in srt_f_get_fx_implied_vol");

    /* Check on Maturity first not to devide by zero */
    if (yr_to_exp <= 0.0)
    {
        *fx_bs_vol = find_fx_sig(0.0, fx_ts);
        return NULL;
    }

    /* Compute the Domestic BS Implied Volatility */
    G_H_func(yr_to_exp, dom_ts, &dom_G, &dom_H);
    dom_psi    = Psi_func(yr_to_exp, dom_ts);
    sq_dom_psi = dom_psi * dom_psi;
    dom_L      = L_func(yr_to_exp, dom_ts);
    dom_O      = O_func(yr_to_exp, dom_ts);

    dom_bd_vol =
        (-sq_dom_psi * dom_G + 2 * dom_psi * dom_L + dom_G * sq_dom_psi - 2 * dom_O) / yr_to_exp;

    /* Compute the Foreign BS Implied Volatility */
    G_H_func(yr_to_exp, for_ts, &for_G, &for_H);
    for_psi    = Psi_func(yr_to_exp, for_ts);
    sq_for_psi = for_psi * for_psi;
    for_L      = L_func(yr_to_exp, for_ts);
    for_O      = O_func(yr_to_exp, for_ts);

    for_bd_vol =
        (-sq_for_psi * for_G + 2 * for_psi * for_L + for_G * sq_for_psi - 2 * for_O) / yr_to_exp;

    /* Compute the BS Implied Volatility of the Spot FX */
    spot_fx_cum_vol = 1.0 / yr_to_exp * fx_cum_vol_func(yr_to_exp, SRT_YES, fx_ts);

    /* Compute the Covariance between the Spot Fx and the Domestic zero coupon bond */
    Vdx = V_dx_func(yr_to_exp, fx_ts);
    Wdx = W_dx_func(yr_to_exp, fx_ts);

    cov_dom_int_fx = -dom_psi * Vdx + Wdx;

    /* Compute the Covariance between the Spot Fx and the Foreign zero coupon bond */
    Vfx = V_fx_func(yr_to_exp, fx_ts);
    Wfx = W_fx_func(yr_to_exp, fx_ts);

    /* Compute the Covariance between the  dom. zc bond and the for zc bond */
    Ofd = O_fd_func(yr_to_exp, fx_ts);
    Rfd = R_fd_func(yr_to_exp, fx_ts);
    Pfd = P_fd_func(yr_to_exp, fx_ts);
    Mfd = M_fd_func(yr_to_exp, fx_ts);

    cov_dom_int_for_int = dom_psi * for_psi * Ofd - for_psi * Rfd - dom_psi * Pfd + Mfd;

    cov_for_int_fx = -for_psi * Vfx + Wfx;

    sq_fx_bs_vol = spot_fx_cum_vol + dom_bd_vol + for_bd_vol + (2 / yr_to_exp) * cov_for_int_fx -
                   (2 / yr_to_exp) * cov_dom_int_fx - (2 / yr_to_exp) * cov_dom_int_for_int;

    *fx_bs_vol = sqrt(sq_fx_bs_vol);

    return err;

} /* Err srt_f_get_fx_implied_vol(...) */

/*==============================================================================*/

/* This routines returns the correlation between the following
                three underlyings : log(FXt/FXt-1)
                                                        rdom(t)
                                                        rfor(t)

                as well as the variance of log(FXt/FXt-1)
          in a FX StochRates (LGM\LGM) Framework.

*/

Err srt_f_get_fx_stoch_rates_correl(
    double   start_date,
    double   end_date,
    String   fx_und_names,
    double** corr_matrix,
    double*  var_ln_fx_out)
{
    TermStruct *      fx_ts, *dom_ts, *for_ts;
    Err               err = NULL;
    String            dom_ir_und_name;
    String            for_ir_und_name;
    SrtUndPtr         fx_und;
    SrtUndPtr         dom_und;
    SrtUndPtr         for_und;
    SrtModelType      model_type;
    SrtModelType      dom_model_type;
    SrtModelType      for_model_type;
    SrtUnderlyingType und_type;
    double            var_dom_sr, var_for_sr, std_dom_sr, std_for_sr;
    double            cov_ln_fx_dom_sr, cov_ln_fx_for_sr, cov_dom_for_sr;
    double            sq_dom_psi, prev_dom_psi, sq_prev_dom_psi;
    double            sq_for_psi, prev_for_psi, sq_prev_psi;
    double            dom_F, dom_G, prev_dom_G, dom_H, prev_dom_H, dom_psi, dom_S_dt;
    double            for_F, for_G, prev_for_G, for_H, prev_for_H, for_psi, for_S_dt;
    double            dom_L_dt, dom_O_dt;
    double            for_L_dt, for_O_dt;
    double            Pfd_dt, Rfd_dt, Ofd_dt, Mfd_dt;
    double            Vdx_dt, Wdx_dt;
    double            Vfx_dt, Wfx_dt;
    double            var_int_dom_sr, var_int_for_sr, var_int_fx;
    double            cov_int_dom_for, cov_int_dom_fx, cov_int_for_fx, var_ln_fx, std_ln_fx;

    /* Gets the Fx Underlying */
    fx_und = lookup_und(fx_und_names);

    /* Gets the Underlying Type and Checks it is an FX */
    und_type = get_underlying_type(fx_und);
    if (und_type != FOREX_UND)
        return serror("Need an Fx Underlying in srt_f_get_fx_stoch_rates_correl !");

    /* Gets the Model Name Type */
    model_type = get_mdltype_from_fxund(fx_und);

    /* Gets the Fx term Structure */
    err = get_underlying_ts(fx_und, &fx_ts);
    if (err)
        return serror("Can not get the FX term structure in srt_f_get_fx_stoch_rates_correl");

    /* If a Standard Fx model, return NULL */
    if (model_type != FX_STOCH_RATES)
    {
        return NULL;
    }

    /* We NOW have an FX STOCH RATES model */

    /* Gets the Domestic Interest Rate Underlying */
    dom_ir_und_name = get_domname_from_fxund(fx_und);
    dom_und         = lookup_und(dom_ir_und_name);
    if (!dom_und)
        return serror("Underlying %s not defined", dom_ir_und_name);

    /* Gets the Foreign Interest Rate Underlying */
    for_ir_und_name = get_forname_from_fxund(fx_und);
    for_und         = lookup_und(for_ir_und_name);
    if (!for_und)
        return serror("Underlying %s not defined", for_ir_und_name);

    /* Checks that we are in the LGM/LGM configuration */

    dom_model_type = get_mdltype_from_irund(dom_und);
    for_model_type = get_mdltype_from_irund(for_und);
    if (dom_model_type != LGM || for_model_type != LGM)
        return serror("srt_f_get_fx_stoch_rates_correl(...)  will only work for LGM/LGM .... ");

    /* Get the dom & for & Fx ts */
    err = get_underlying_ts(dom_und, &dom_ts);
    if (err)
        return serror("Can not get the domestic term structure in srt_f_get_fx_stoch_rates_correl");

    err = get_underlying_ts(for_und, &for_ts);
    if (err)
        return serror("Can not get the foreign term structure in srt_f_get_fx_stoch_rates_correl");

    /* Call to domestic functions */

    dom_F = F_func(end_date, dom_ts);
    G_H_func(end_date, dom_ts, &dom_G, &dom_H);
    G_H_func(start_date, dom_ts, &prev_dom_G, &prev_dom_H);
    prev_dom_psi    = Psi_func(start_date, dom_ts);
    sq_prev_dom_psi = prev_dom_psi * prev_dom_psi;
    dom_psi         = Psi_func(end_date, dom_ts);
    sq_dom_psi      = dom_psi * dom_psi;
    dom_L_dt        = L_func(end_date, dom_ts) - L_func(start_date, dom_ts);
    dom_O_dt        = O_func(end_date, dom_ts) - O_func(start_date, dom_ts);
    dom_S_dt        = S_func(end_date, dom_ts) - S_func(start_date, dom_ts);

    /* Call to foreign functions */

    for_F = F_func(end_date, for_ts);
    G_H_func(end_date, for_ts, &for_G, &for_H);
    G_H_func(start_date, for_ts, &prev_for_G, &prev_for_H);
    prev_for_psi = Psi_func(start_date, for_ts);
    sq_prev_psi  = prev_for_psi * prev_for_psi;
    for_psi      = Psi_func(end_date, for_ts);
    sq_for_psi   = for_psi * for_psi;
    for_L_dt     = L_func(end_date, for_ts) - L_func(start_date, for_ts);
    for_O_dt     = O_func(end_date, for_ts) - O_func(start_date, for_ts);
    for_S_dt     = S_func(end_date, for_ts) - S_func(start_date, for_ts);

    /* Call to mixed functions */

    Pfd_dt = P_fd_func(end_date, fx_ts) - P_fd_func(start_date, fx_ts);
    Rfd_dt = R_fd_func(end_date, fx_ts) - R_fd_func(start_date, fx_ts);
    Ofd_dt = O_fd_func(end_date, fx_ts) - O_fd_func(start_date, fx_ts);
    Vdx_dt = V_dx_func(end_date, fx_ts) - V_dx_func(start_date, fx_ts);
    Wdx_dt = W_dx_func(end_date, fx_ts) - W_dx_func(start_date, fx_ts);
    Vfx_dt = V_fx_func(end_date, fx_ts) - V_fx_func(start_date, fx_ts);
    Wfx_dt = W_fx_func(end_date, fx_ts) - W_fx_func(start_date, fx_ts);
    Mfd_dt = M_fd_func(end_date, fx_ts) - M_fd_func(start_date, fx_ts);

    /* Get the variance of log(Xt) */

    var_int_dom_sr = sq_dom_psi * (dom_G - prev_dom_G) -
                     2 * dom_psi * (dom_psi * dom_G - prev_dom_psi * prev_dom_G - dom_L_dt) +
                     dom_G * sq_dom_psi - prev_dom_G * sq_prev_dom_psi - 2 * dom_O_dt;

    /* Computation of the variance of the integral of the foreign short rate */
    var_int_for_sr = sq_for_psi * (for_G - prev_for_G) -
                     2 * for_psi * (for_psi * for_G - prev_for_psi * prev_for_G - for_L_dt) +
                     for_G * sq_for_psi - prev_for_G * sq_prev_psi - 2 * for_O_dt;

    /* Computation of the variance of the FX integral */
    var_int_fx = fx_cum_vol_func(end_date - start_date, SRT_YES, fx_ts);

    /* Computation of the Variance the ln of the FX  */

    cov_int_dom_for = dom_psi * for_psi * Ofd_dt - for_psi * Rfd_dt - dom_psi * Pfd_dt + Mfd_dt;

    cov_int_dom_fx = dom_psi * Vdx_dt - Wdx_dt;
    cov_int_for_fx = for_psi * Vfx_dt - Wfx_dt;

    var_ln_fx = var_int_dom_sr + var_int_for_sr + var_int_fx - 2 * cov_int_dom_for +
                2 * cov_int_dom_fx - 2 * cov_int_for_fx;
    std_ln_fx = sqrt(var_ln_fx);

    (*var_ln_fx_out) = var_ln_fx;

    /* Get the correlations between the three underlyings */

    var_dom_sr = dom_F * dom_F * (dom_G - prev_dom_G);
    var_for_sr = for_F * for_F * (for_G - prev_for_G);
    std_dom_sr = sqrt(var_dom_sr);
    std_for_sr = sqrt(var_for_sr);

    cov_ln_fx_dom_sr = dom_F * (dom_psi * (dom_G - prev_dom_G) - dom_S_dt) +
                       dom_F * (Pfd_dt - for_psi * Ofd_dt) + dom_F * Vdx_dt;

    cov_ln_fx_for_sr = -for_F * (for_psi * (for_G - prev_for_G) - for_S_dt) -
                       for_F * (Rfd_dt - dom_psi * Ofd_dt) + for_F * Vfx_dt;

    cov_dom_for_sr = dom_F * for_F * Ofd_dt;

    (*corr_matrix)[0] = cov_ln_fx_dom_sr / (std_dom_sr * std_ln_fx);
    (*corr_matrix)[1] = cov_ln_fx_for_sr / (std_for_sr * std_ln_fx);
    (*corr_matrix)[2] = cov_dom_for_sr / (std_for_sr * std_dom_sr);

    return err;
}

Err srt_f_get_fx_stoch_rates_ind_correl(
    double  start_date,
    double  end_date,
    String  first_und_name,
    String  second_und_name,
    String  fx_und_name,
    double* corr_or_var)
{
    double *  corr_matrix, var_ln_fx;
    Err       err;
    SrtUndPtr fx_und;
    String    dom_und_name, for_und_name;
    long      today;

    corr_matrix = dvector(0, 2);

    fx_und = lookup_und(fx_und_name);

    today = get_today_from_underlying(fx_und);

    err = srt_f_get_fx_stoch_rates_correl(
        (start_date - today) * YEARS_IN_DAY,
        (end_date - today) * YEARS_IN_DAY,
        fx_und_name,
        &corr_matrix,
        &var_ln_fx);
    if (err)
        return err;

    dom_und_name = get_domname_from_fxund(fx_und);
    for_und_name = get_forname_from_fxund(fx_und);

    if (((!strcmp(first_und_name, fx_und_name)) && (!strcmp(second_und_name, dom_und_name))) ||
        ((!strcmp(first_und_name, dom_und_name)) && (!strcmp(second_und_name, fx_und_name))))
        *corr_or_var = corr_matrix[0];
    else if (
        ((!strcmp(first_und_name, fx_und_name)) && (!strcmp(second_und_name, for_und_name))) ||
        ((!strcmp(first_und_name, for_und_name)) && (!strcmp(second_und_name, fx_und_name))))
        *corr_or_var = corr_matrix[1];
    else if (
        ((!strcmp(first_und_name, dom_und_name)) && (!strcmp(second_und_name, for_und_name))) ||
        ((!strcmp(first_und_name, for_und_name)) && (!strcmp(second_und_name, dom_und_name))))
        *corr_or_var = corr_matrix[2];
    else if ((!strcmp(first_und_name, fx_und_name)) && (!strcmp(second_und_name, fx_und_name)))
        *corr_or_var = var_ln_fx;

    if (corr_matrix)
        free_dvector(corr_matrix, 0, 2);
    corr_matrix = NULL;

    return err;
}
