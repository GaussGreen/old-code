#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

static double Z_Func(double x, double t1, double t2)
{
    return ((exp(x * t2) - exp(x * t1)) / x);
}

/* integral(s, t, exp(-x * (T - u))) */
static double Phi_Func(double x, double T, double s, double t)
{
    double result;

    result = (exp(-x * (T - t)) - exp(-x * (T - s))) / x;

    return result;
}

/* integral(s, t, sigB(u, T)) */
static double Etha_Func(double x, double T, double s, double t)
{
    double result;

    result = (t - s - Phi_Func(x, T, s, t)) / x;

    return result;
}

/* integral of Etha */
static double B_Func(double x, double T, double s, double t)
{
    double result;

    result = -(t * t - s * s) / 2.0 +
             1.0 / x * ((t - 1.0 / x) * exp(-x * (T - t)) - (s - 1.0 / x) * exp(-x * (T - s)));

    return result;
}

/* integral(s, t, sigB(u, T, x) * sigB(u, T, y)) */
static double Psi_Func(double x, double y, double T, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, T, s, t) - Phi_Func(y, T, s, t) + Phi_Func(x + y, T, s, t));

    return result;
}

/* integral(s, t, sigB(u, Tx, x) * sigB(u, Ty, y)) */
static double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t)
{
    double result;

    result = 1.0 / (x * y) *
             (t - s - Phi_Func(x, Tx, s, t) - Phi_Func(y, Ty, s, t) +
              exp(-x * (Tx - Ty)) * Phi_Func(x + y, Ty, s, t));

    return result;
}

Err fill_mc_init(
    long    pay_date,
    double  pay_time,
    double* date,
    double* time,
    long    nb_dates,
    double* sig_dates,
    long    nb_sig_dates,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* sig_curve_fx,
    double  correl_dom_for,
    double  correl_dom_fx,
    double  correl_for_fx,
    char*   dom_yc,
    char*   for_yc,
    double* dom_ifr,
    double* dom_fwd,
    double* dom_std,
    double* dom_phi,
    double* dom_beta,
    double* dom_bond_pay,
    double* dom_beta_pay,
    double* for_ifr,
    double* for_fwd,
    double* for_std,
    double* for_phi,
    double* for_beta,
    double* fx_fwd,
    double* fx_std,
    double* dom_for_cov,
    double* dom_fx_cov,
    double* for_fx_cov)
{
    double sig_dom, sig_dom2;
    double sig_for, sig_for2;
    double sig_fx, sig_domfor;
    double T1, T2, start_date, end_date, start_mat, end_mat;
    double var_dom, expect_dom;
    double var_for, expect_for;
    double expect_fx, QTexpect, fx_vol;
    double adj_fx_pay, adj_quanto, adj_pay_dom, adj_pay_for;
    double x_dom, y_dom, x_domfor;
    double x_for, y_for;
    double phi_dom, phi_for, zc_dom, zc_for, zc_pay;
    double mat, pay_mat, mat_pay;
    double correl12, correl13_1, correl13_2, correl13_3, correl23_1, correl23_2, correl23_3;
    int    i, k;
    long   StartIndex, EndIndex;
    Err    err = NULL;

    dom_fwd[0] = 0;
    dom_std[0] = 0;
    dom_phi[0] = 0;
    phi_dom    = 0;

    for_fwd[0] = 0;
    for_std[0] = 0;
    for_phi[0] = 0;
    phi_for    = 0;

    mat_pay = (pay_date - date[0]) / 365.0;

    start_date = date[0];
    start_mat  = time[0];
    StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);

    for (k = 0; k < nb_dates - 1; k++)
    {
        end_date = date[k + 1];
        end_mat  = time[k + 1];
        EndIndex = Get_Index(end_mat, sig_dates, nb_sig_dates);
        mat      = (end_mat - start_mat);

        var_dom     = 0;
        expect_dom  = 0;
        adj_pay_dom = 0;

        var_for     = 0;
        expect_for  = 0;
        adj_quanto  = 0;
        adj_pay_for = 0;

        dom_ifr[k] = swp_f_zr(start_date, start_date + 1, dom_yc);
        for_ifr[k] = swp_f_zr(start_date, start_date + 1, for_yc);

        zc_dom = swp_f_zr(start_date, end_date, dom_yc);
        zc_for = swp_f_zr(start_date, end_date, for_yc);

        dom_beta[k] = -1.0 / lda_dom * (1.0 - exp(-lda_dom * mat));
        for_beta[k] = -1.0 / lda_for * (1.0 - exp(-lda_for * mat));

        /* QTexpect of the log Fx */
        QTexpect = -mat * (zc_for - zc_dom) - 0.5 * (for_beta[k] * for_beta[k] * for_phi[k] -
                                                     dom_beta[k] * dom_beta[k] * dom_phi[k]);

        zc_pay          = swp_f_zr(start_date, pay_date, dom_yc);
        pay_mat         = (pay_date - start_date) / 365.0;
        dom_bond_pay[k] = exp(-zc_pay * pay_mat);
        dom_beta_pay[k] = -1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));

        /*	Implied Fx Vol */
        err = Fx3DtsImpliedVol(
            end_mat,
            start_mat,
            end_mat,
            sig_dates,
            nb_sig_dates,
            sig_curve_dom,
            lda_dom,
            sig_curve_for,
            lda_for,
            sig_dates,
            sig_curve_fx,
            nb_sig_dates,
            correl_dom_for,
            correl_dom_fx,
            correl_for_fx,
            &fx_vol);

        /*	Calculate expectation of the log Fx under Q-Tfix */
        expect_fx = QTexpect - 0.5 * fx_vol * fx_vol * mat;

        adj_fx_pay = 0;

        correl12   = 0;
        correl13_1 = correl13_2 = correl13_3 = 0;
        correl23_1 = correl23_2 = correl23_3 = 0;

        for (i = StartIndex; i < EndIndex + 1; i++)
        {
            if (i > StartIndex)
            {
                T1 = sig_dates[i - 1];
            }
            else
            {
                /* First part */
                T1 = start_mat;
            }

            if (i == EndIndex || StartIndex == EndIndex)
            {
                /* Last part */
                T2 = end_mat;
            }
            else
            {
                T2 = sig_dates[i];
            }

            sig_dom = sig_dom2 = sig_domfor = sig_curve_dom[i];
            sig_dom2 *= sig_dom2;
            sig_for = sig_for2 = sig_curve_for[i];
            sig_for2 *= sig_for2;
            sig_domfor *= sig_for;
            sig_fx = sig_curve_fx[i];

            x_dom = Z_Func(lda_dom, T1, T2);
            y_dom = Z_Func(2 * lda_dom, T1, T2);

            x_for = Z_Func(lda_for, T1, T2);
            y_for = Z_Func(2 * lda_for, T1, T2);

            x_domfor = Z_Func(lda_dom + lda_for, T1, T2);

            expect_dom += sig_dom2 * (Phi_Func(lda_dom, end_mat, T1, T2) -
                                      Phi_Func(2 * lda_dom, end_mat, T1, T2));

            /* Q-Tpay adjustment */
            adj_pay_dom += sig_dom2 * (x_dom - exp(-lda_dom * mat_pay) * y_dom);

            var_dom += sig_dom2 * y_dom;

            expect_for += sig_for2 * (Phi_Func(lda_for, end_mat, T1, T2) -
                                      Phi_Func(2 * lda_for, end_mat, T1, T2));

            /* Quanto and Q-T pay adjustment */
            adj_quanto += sig_for * sig_fx * x_for;
            adj_pay_for += sig_domfor * (x_for - exp(-lda_dom * mat_pay) * x_domfor);

            var_for += sig_for2 * y_for;

            /* Fx adjustment to go to Qpay */
            adj_fx_pay +=
                correl_dom_fx * sig_dom * sig_fx *
                    (Etha_Func(lda_dom, end_mat, T1, T2) - Etha_Func(lda_dom, mat_pay, T1, T2)) -
                correl_dom_for * sig_domfor *
                    (Psi_Func(lda_for, lda_dom, end_mat, T1, T2) -
                     Psi2_Func(lda_for, lda_dom, end_mat, mat_pay, T1, T2)) +
                sig_dom2 * (Psi_Func(lda_dom, lda_dom, end_mat, T1, T2) -
                            Psi2_Func(lda_dom, lda_dom, end_mat, mat_pay, T1, T2));

            /* Correlations */

            /* domestic / foreign */
            correl12 += sig_dom * sig_for * x_domfor;

            /* domestic / fx */
            correl13_1 += sig_dom2 * (x_dom - exp(-lda_dom * end_mat) * y_dom);
            correl13_2 += sig_domfor * (x_dom - exp(-lda_for * end_mat) * x_domfor);
            correl13_3 += correl_dom_fx * sig_dom * sig_fx * x_dom;

            /* foreign fx */
            correl23_1 += sig_for2 * (x_for - exp(-lda_for * end_mat) * y_for);
            correl23_2 += sig_domfor * (x_for - exp(-lda_dom * end_mat) * x_domfor);
            correl23_3 += sig_for * sig_fx * x_for;
        }

        phi_dom += var_dom;
        phi_for += var_for;

        /* Forward and Standard deviation */
        dom_fwd[k + 1] =
            1.0 / lda_dom * (expect_dom - exp(-lda_dom * end_mat) * adj_pay_dom) +
            dom_phi[k] * exp(-lda_dom * mat) * Phi_Func(-lda_dom, start_mat, start_mat, end_mat);

        dom_std[k + 1] = exp(-lda_dom * end_mat) * sqrt(var_dom);
        for_fwd[k + 1] =
            1.0 / lda_for * expect_for +
            exp(-lda_for * end_mat) *
                (-correl_dom_for * adj_pay_for / lda_dom - correl_for_fx * adj_quanto) +
            for_phi[k] * exp(-lda_for * mat) * Phi_Func(-lda_for, start_mat, start_mat, end_mat);

        for_std[k + 1] = exp(-lda_for * end_mat) * sqrt(var_for);
        fx_fwd[k + 1]  = expect_fx + adj_fx_pay;
        fx_std[k + 1]  = fx_vol * sqrt(mat);

        /* Covariance*/
        dom_for_cov[k + 1] = correl_dom_for * correl12 * exp(-(lda_dom + lda_for) * end_mat);

        dom_fx_cov[k + 1] =
            exp(-lda_dom * end_mat) *
            (correl13_1 / lda_dom - correl_dom_for * correl13_2 / lda_for + correl13_3);

        for_fx_cov[k + 1] = exp(-lda_for * end_mat) *
                            (-correl23_1 / lda_for + correl_dom_for * correl23_2 / lda_dom +
                             correl_for_fx * correl23_3);

        /* Phi */
        dom_phi[k + 1] = phi_dom * exp(-2 * lda_dom * end_mat);
        for_phi[k + 1] = phi_for * exp(-2 * lda_for * end_mat);

        start_date = end_date;
        start_mat  = end_mat;
        StartIndex = EndIndex;
    }

    dom_ifr[k]      = swp_f_zr(date[k], date[k] + 1, dom_yc);
    for_ifr[k]      = swp_f_zr(date[k], date[k] + 1, for_yc);
    zc_pay          = swp_f_zr(date[k], pay_date, dom_yc);
    pay_mat         = (pay_date - date[k]) / 365.0;
    dom_bond_pay[k] = exp(-zc_pay * pay_mat);
    dom_beta_pay[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));

    return err;
}

Err fill_mc_init_corr(
    long    pay_date,
    double  pay_time,
    double* date,
    double* time,
    long    nb_dates,
    double* sig_dates,
    long    nb_sig_dates,
    double* sig_curve_dom,
    double  lda_dom,
    double* sig_curve_for,
    double  lda_for,
    double* sig_curve_fx,
    double* correl_dom_for_ts,
    double* correl_dom_fx_ts,
    double* correl_for_fx_ts,
    char*   dom_yc,
    char*   for_yc,
    double* dom_ifr,
    double* dom_fwd,
    double* dom_std,
    double* dom_phi,
    double* dom_beta,
    double* dom_bond_pay,
    double* dom_beta_pay,
    double* for_ifr,
    double* for_fwd,
    double* for_std,
    double* for_phi,
    double* for_beta,
    double* fx_fwd,
    double* fx_std,
    double* dom_for_cov,
    double* dom_fx_cov,
    double* for_fx_cov)
{
    double sig_dom, sig_dom2;
    double sig_for, sig_for2;
    double sig_fx, sig_domfor;
    double T1, T2, start_date, end_date, start_mat, end_mat;
    double var_dom, expect_dom;
    double var_for, expect_for;
    double expect_fx, QTexpect, fx_vol;
    double adj_fx_pay, adj_quanto, adj_pay_dom, adj_pay_for;
    double x_dom, y_dom, x_domfor;
    double x_for, y_for;
    double phi_dom, phi_for, zc_dom, zc_for, zc_pay;
    double mat, pay_mat, mat_pay;
    double correl_dom_for, correl_dom_fx, correl_for_fx;
    double correl12, correl13_1, correl13_2, correl13_3, correl23_1, correl23_2, correl23_3;
    int    i, k;
    long   StartIndex, EndIndex;
    Err    err = NULL;

    dom_fwd[0] = 0;
    dom_std[0] = 0;
    dom_phi[0] = 0;
    phi_dom    = 0;

    for_fwd[0] = 0;
    for_std[0] = 0;
    for_phi[0] = 0;
    phi_for    = 0;

    mat_pay = (pay_date - date[0]) / 365.0;

    start_date = date[0];
    start_mat  = time[0];
    StartIndex = Get_Index(start_mat, sig_dates, nb_sig_dates);

    for (k = 0; k < nb_dates - 1; k++)
    {
        end_date = date[k + 1];
        end_mat  = time[k + 1];
        EndIndex = Get_Index(end_mat, sig_dates, nb_sig_dates);
        mat      = (end_mat - start_mat);

        var_dom     = 0;
        expect_dom  = 0;
        adj_pay_dom = 0;

        var_for     = 0;
        expect_for  = 0;
        adj_quanto  = 0;
        adj_pay_for = 0;

        dom_ifr[k] = swp_f_zr(start_date, start_date + 1, dom_yc);
        for_ifr[k] = swp_f_zr(start_date, start_date + 1, for_yc);

        zc_dom = swp_f_zr(start_date, end_date, dom_yc);
        zc_for = swp_f_zr(start_date, end_date, for_yc);

        dom_beta[k] = -1.0 / lda_dom * (1.0 - exp(-lda_dom * mat));
        for_beta[k] = -1.0 / lda_for * (1.0 - exp(-lda_for * mat));

        /* QTexpect of the log Fx */
        QTexpect = -mat * (zc_for - zc_dom) - 0.5 * (for_beta[k] * for_beta[k] * for_phi[k] -
                                                     dom_beta[k] * dom_beta[k] * dom_phi[k]);

        zc_pay          = swp_f_zr(start_date, pay_date, dom_yc);
        pay_mat         = (pay_date - start_date) / 365.0;
        dom_bond_pay[k] = exp(-zc_pay * pay_mat);
        dom_beta_pay[k] = -1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));

        /*	Implied Fx Vol */
        err = Fx3DtsImpliedVol_corr(
            end_mat,
            start_mat,
            end_mat,
            sig_dates,
            nb_sig_dates,
            sig_curve_dom,
            lda_dom,
            sig_curve_for,
            lda_for,
            sig_dates,
            sig_curve_fx,
            nb_sig_dates,
            sig_dates,
            correl_dom_for_ts,
            correl_dom_fx_ts,
            correl_for_fx_ts,
            nb_sig_dates,
            &fx_vol);

        /*	Calculate expectation of the log Fx under Q-Tfix */
        expect_fx = QTexpect - 0.5 * fx_vol * fx_vol * mat;

        adj_fx_pay = 0;

        correl12   = 0;
        correl13_1 = correl13_2 = correl13_3 = 0;
        correl23_1 = correl23_2 = correl23_3 = 0;

        for (i = StartIndex; i < EndIndex + 1; i++)
        {
            if (i > StartIndex)
            {
                T1 = sig_dates[i - 1];
            }
            else
            {
                /* First part */
                T1 = start_mat;
            }

            if (i == EndIndex || StartIndex == EndIndex)
            {
                /* Last part */
                T2 = end_mat;
            }
            else
            {
                T2 = sig_dates[i];
            }

            sig_dom = sig_dom2 = sig_domfor = sig_curve_dom[i];
            sig_dom2 *= sig_dom2;
            sig_for = sig_for2 = sig_curve_for[i];
            sig_for2 *= sig_for2;
            sig_domfor *= sig_for;
            sig_fx         = sig_curve_fx[i];
            correl_dom_for = correl_dom_for_ts[i];
            correl_dom_fx  = correl_dom_fx_ts[i];
            correl_for_fx  = correl_for_fx_ts[i];

            x_dom = Z_Func(lda_dom, T1, T2);
            y_dom = Z_Func(2 * lda_dom, T1, T2);

            x_for = Z_Func(lda_for, T1, T2);
            y_for = Z_Func(2 * lda_for, T1, T2);

            x_domfor = Z_Func(lda_dom + lda_for, T1, T2);

            expect_dom += sig_dom2 * (Phi_Func(lda_dom, end_mat, T1, T2) -
                                      Phi_Func(2 * lda_dom, end_mat, T1, T2));

            /* Q-Tpay adjustment */
            adj_pay_dom += sig_dom2 * (x_dom - exp(-lda_dom * mat_pay) * y_dom);

            var_dom += sig_dom2 * y_dom;

            expect_for += sig_for2 * (Phi_Func(lda_for, end_mat, T1, T2) -
                                      Phi_Func(2 * lda_for, end_mat, T1, T2));

            /* Quanto and Q-T pay adjustment */
            adj_quanto += correl_for_fx * sig_for * sig_fx * x_for;
            adj_pay_for +=
                correl_dom_for * sig_domfor * (x_for - exp(-lda_dom * mat_pay) * x_domfor);

            var_for += sig_for2 * y_for;

            /* Fx adjustment to go to Qpay */
            adj_fx_pay +=
                correl_dom_fx * sig_dom * sig_fx *
                    (Etha_Func(lda_dom, end_mat, T1, T2) - Etha_Func(lda_dom, mat_pay, T1, T2)) -
                correl_dom_for * sig_domfor *
                    (Psi_Func(lda_for, lda_dom, end_mat, T1, T2) -
                     Psi2_Func(lda_for, lda_dom, end_mat, mat_pay, T1, T2)) +
                sig_dom2 * (Psi_Func(lda_dom, lda_dom, end_mat, T1, T2) -
                            Psi2_Func(lda_dom, lda_dom, end_mat, mat_pay, T1, T2));

            /* Correlations */

            /* domestic / foreign */
            correl12 += correl_dom_for * sig_dom * sig_for * x_domfor;

            /* domestic / fx */
            correl13_1 += sig_dom2 * (x_dom - exp(-lda_dom * end_mat) * y_dom);
            correl13_2 +=
                correl_dom_for * sig_domfor * (x_dom - exp(-lda_for * end_mat) * x_domfor);
            correl13_3 += correl_dom_fx * sig_dom * sig_fx * x_dom;

            /* foreign fx */
            correl23_1 += sig_for2 * (x_for - exp(-lda_for * end_mat) * y_for);
            correl23_2 +=
                correl_dom_for * sig_domfor * (x_for - exp(-lda_dom * end_mat) * x_domfor);
            correl23_3 += correl_for_fx * sig_for * sig_fx * x_for;
        }

        phi_dom += var_dom;
        phi_for += var_for;

        /* Forward and Standard deviation */
        dom_fwd[k + 1] =
            1.0 / lda_dom * (expect_dom - exp(-lda_dom * end_mat) * adj_pay_dom) +
            dom_phi[k] * exp(-lda_dom * mat) * Phi_Func(-lda_dom, start_mat, start_mat, end_mat);

        dom_std[k + 1] = exp(-lda_dom * end_mat) * sqrt(var_dom);
        for_fwd[k + 1] =
            1.0 / lda_for * expect_for +
            exp(-lda_for * end_mat) * (-adj_pay_for / lda_dom - adj_quanto) +
            for_phi[k] * exp(-lda_for * mat) * Phi_Func(-lda_for, start_mat, start_mat, end_mat);

        for_std[k + 1] = exp(-lda_for * end_mat) * sqrt(var_for);
        fx_fwd[k + 1]  = expect_fx + adj_fx_pay;
        fx_std[k + 1]  = fx_vol * sqrt(mat);

        /* Covariance*/
        dom_for_cov[k + 1] = correl12 * exp(-(lda_dom + lda_for) * end_mat);

        dom_fx_cov[k + 1] =
            exp(-lda_dom * end_mat) * (correl13_1 / lda_dom - correl13_2 / lda_for + correl13_3);

        for_fx_cov[k + 1] =
            exp(-lda_for * end_mat) * (-correl23_1 / lda_for + correl23_2 / lda_dom + correl23_3);

        /* Phi */
        dom_phi[k + 1] = phi_dom * exp(-2 * lda_dom * end_mat);
        for_phi[k + 1] = phi_for * exp(-2 * lda_for * end_mat);

        start_date = end_date;
        start_mat  = end_mat;
        StartIndex = EndIndex;
    }

    dom_ifr[k]      = swp_f_zr(date[k], date[k] + 1, dom_yc);
    for_ifr[k]      = swp_f_zr(date[k], date[k] + 1, for_yc);
    zc_pay          = swp_f_zr(date[k], pay_date, dom_yc);
    pay_mat         = (pay_date - date[k]) / 365.0;
    dom_bond_pay[k] = exp(-zc_pay * pay_mat);
    dom_beta_pay[k] = 1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));

    return err;
}
