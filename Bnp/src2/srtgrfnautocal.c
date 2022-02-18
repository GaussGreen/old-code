
#include "srtgrfnautocal.h"

#include "SrtAccess.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtgrfnmainQuanto.h"

void fx3datc_set_default_param(FX3DATC_PARAM param)
{
    param->num_stp   = 50;
    param->num_paths = 25000;
    param->num_stpx  = 50; /*-----------For Quanto Product Only----------*/
    param->do_pecs   = 0;
}

#define MAX_CAL_DATES 512
Err GrfnAutocalCaller(
    /*	Today */
    long   today,
    double spot_fx,
    /*	Get Cash Vol function */
    Err (*get_cash_vol)(/*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    /*	Domestic market */
    char* dom_und_name,
    char* dom_ccy,
    char* dom_yc_name,        /*	Name of the yield curve */
    char* dom_vol_curve_name, /*	Name of the market vol curve */
    char* dom_ref_rate_name,  /*	Name of the reference rate */
    char* dom_instr_freq,     /*	Frequency and basis of instruments */
    char* dom_instr_basis,
    /*	Domestic Model */
    double dom_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Foreign market */
    char* for_und_name,
    char* for_ccy,
    char* for_yc_name,        /*	Name of the yield curve */
    char* for_vol_curve_name, /*	Name of the market vol curve */
    char* for_ref_rate_name,  /*	Name of the reference rate */
    char* for_instr_freq,     /*	Frequency and basis of instruments */
    char* for_instr_basis,
    /*	Domestic Model */
    double for_lambda, /*	Lambda: may NOT be changed in the process */
    /*	Fx market */
    long*   fx_mkt_vol_date, /*	Option maturity dates */
    double* fx_mkt_vol,      /*	Option BS vol */
    int     num_fx_mkt_vol,  /*	Number of options */
    /*	Fx model */
    double* corr_times,
    double* correl_dom_for,
    double* correl_dom_fx,
    double* correl_for_fx,
    long    corr_n_times,
    /*	Structure */
    /*	If ex_date is NULL,
    exercise dates will be generated 2bd before start */
    int num_ex_dates, /*	Exercise dates,
                                                      all supposed to be on or after today */
    long* ex_date,    /*	Supposed to be sorted */
    int*  cal_date,   /*	1: use ex_date as calibration date, 0: don't */
    /*	If calibration dates are overwritten */
    int     owr_cal_dates, /*	1: use overwritten dates, 0: use ex_date */
    int     num_owr_dates, /*	Number of overwritten dates */
    long*   owr_date,      /*	Overwritten dates */
    char**  dom_end_tenor, /*	Tenors of the underlying instruments or "DIAG" */
    char**  for_end_tenor,
    long    end_date,   /*	End date for diagonal */
    double* dom_strike, /*	Domestic strikes 0: ATM */
    double* for_strike, /*	Foreign strikes 0: ATM */
    /*	Calibration parameters */
    CPD_DIAG_CALIB_PARAM calib_param,
    /*	GRFN Tableau */
    long     tableauRows,
    long     tableauCols,
    char***  tableauStrings,
    int**    tableauMask,
    long     auxWidth,
    long*    auxLen,
    double** aux,
    /* End of Day Flags */
    int is_end_of_day_fixing,
    int is_end_of_day_payment,
    /*	Pricing method */
    int mth, /*	0: Tree, 1: Monte-Carlo */
    /*	Pricing parameters */
    FX3DATC_PARAM pricing_param,
    /*	Calibration output */
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** dom_sig,
    double** for_sig,
    int*     num_fx_vol,
    double** fx_vol_time,
    double** fx_vol,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA dom_inst_data, /*	NULL = don't save calibration instrument data */
    CPD_CALIB_INST_DATA for_inst_data, /*	NULL = don't save calibration instrument data */
    /*	Pricing output */
    int*     num_prod,
    double** prod_val,
    /*	Whether to keep the data */
    int keep_calib_data,
    int keep_und,
    /*	Quanto or Hybrid ?	*/
    int quantoflag)
{
    Err      err            = NULL;
    int      free_inst_data = 0, free_und = 0;
    double** temp_prod_val = NULL;
    int      i;
    char     fx_und_name[256];
    double*  barrier;
    int*     bar_col;
    int      use_ndates;
    long     use_dates[MAX_CAL_DATES];
    int      use_cal[MAX_CAL_DATES];

    /*	Init */
    *prod_val    = NULL;
    *sig_time    = NULL;
    *dom_sig     = NULL;
    *for_sig     = NULL;
    *fx_vol_time = NULL;
    *fx_vol      = NULL;
    *prod_val    = NULL;

    if (owr_cal_dates)
    {
        use_ndates = num_owr_dates;
        for (i = 0; i < use_ndates; i++)
        {
            use_dates[i] = owr_date[i];
            use_cal[i]   = 1;
        }
    }
    else
    {
        use_ndates = num_ex_dates;
        for (i = 0; i < use_ndates; i++)
        {
            use_dates[i] = ex_date[i];
            use_cal[i]   = cal_date[i];
        }
    }

    /*	Calibrate all */
    err = cpd_calib_all(
        today,
        get_cash_vol,
        dom_yc_name,
        dom_vol_curve_name,
        dom_ref_rate_name,
        dom_instr_freq,
        dom_instr_basis,
        dom_lambda,
        for_yc_name,
        for_vol_curve_name,
        for_ref_rate_name,
        for_instr_freq,
        for_instr_basis,
        for_lambda,
        fx_mkt_vol_date,
        fx_mkt_vol,
        num_fx_mkt_vol,
        corr_times,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        corr_n_times,
        use_ndates,
        use_dates,
        use_cal,
        dom_end_tenor,
        for_end_tenor,
        end_date,
        dom_strike,
        for_strike,
        num_sig,
        sig_time,
        dom_sig,
        for_sig,
        num_fx_vol,
        fx_vol_time,
        fx_vol,
        calib_param,
        dom_inst_data,
        for_inst_data);
    if (err)
    {
        goto FREE_RETURN;
    }
    free_inst_data = 1;

    /*	Create underlyings */
    err = cpd_calib_all_makeund(
        today,
        dom_ccy,
        dom_und_name,
        dom_yc_name,
        *dom_sig,
        dom_lambda,
        for_ccy,
        for_und_name,
        for_yc_name,
        *for_sig,
        for_lambda,
        corr_times,
        correl_dom_for,
        correl_dom_fx,
        correl_for_fx,
        corr_n_times,
        *num_sig,
        *sig_time,
        *num_fx_vol,
        spot_fx,
        *fx_vol_time,
        *fx_vol,
        fx_und_name);
    if (err)
    {
        goto FREE_RETURN;
    }
    free_und = 1;

    /*	Price */
    if (mth == 0)
    {
        if (quantoflag == 1) /*	Quanto Product	*/
        {
            /*	PDE */
            err = SrtGrfnQuantoPDEWithCorrelTS2(
                fx_und_name,
                num_ex_dates,
                ex_date,
                tableauRows,
                tableauCols,
                tableauStrings,
                tableauMask,
                auxWidth,
                auxLen,
                aux,
                pricing_param->num_stp,
                pricing_param->num_stpx,
                num_prod,
                prod_val);

            if (err)
            {
                goto FREE_RETURN;
            }
        }

        else /*	Hybrid Product	*/
        {
            barrier = (double*)calloc(num_ex_dates, sizeof(double));
            bar_col = (int*)calloc(num_ex_dates, sizeof(int));
            if (!barrier || !bar_col)
            {
                goto FREE_RETURN;
            }

            /*	Tree */
            err = SrtGrfn3DFXTree(
                fx_und_name,
                num_ex_dates,
                ex_date,
                tableauRows,
                tableauCols,
                tableauStrings,
                tableauMask,
                auxWidth,
                auxLen,
                aux,
                is_end_of_day_fixing,
                is_end_of_day_payment,
                barrier,
                bar_col,
                pricing_param->num_stp,
                num_prod,
                1,
                prod_val);
            free(barrier);
            free(bar_col);
            if (err)
            {
                goto FREE_RETURN;
            }
        }
    }
    else
    {
        /*	MC */
        err = SrtGrfn3DFXMc(
            fx_und_name,
            num_ex_dates,
            ex_date,
            tableauRows,
            &tableauCols,
            tableauStrings,
            tableauMask,
            auxWidth,
            auxLen,
            aux,
            is_end_of_day_fixing,
            is_end_of_day_payment,
            pricing_param->num_paths,
            pricing_param->do_pecs,
            &temp_prod_val);
        if (err)
        {
            goto FREE_RETURN;
        }

        *num_prod = tableauCols;
        *prod_val = (double*)calloc(*num_prod, sizeof(double));
        if (!*prod_val)
        {
            err = "Memory allocation failure in GrfnAutocalCaller";
            goto FREE_RETURN;
        }
        for (i = 0; i < *num_prod; i++)
        {
            (*prod_val)[i] = temp_prod_val[i][0];
        }
        free_dmatrix(temp_prod_val, 0, *num_prod - 1, 0, 1);

        temp_prod_val = NULL;
    }

FREE_RETURN:

    if (temp_prod_val)
    {
        free_dmatrix(temp_prod_val, 0, *num_prod - 1, 0, 1);
        temp_prod_val = NULL;
    }

    if (err)
    {
        if (*sig_time)
        {
            free(*sig_time);
            *sig_time = NULL;
        }
        if (*dom_sig)
        {
            free(*dom_sig);
            *dom_sig = NULL;
        }
        if (*for_sig)
        {
            free(*for_sig);
            *for_sig = NULL;
        }
        if (*fx_vol_time)
        {
            free(*fx_vol_time);
            *fx_vol_time = NULL;
        }
        if (*fx_vol)
        {
            free(*fx_vol);
            *fx_vol = NULL;
        }
        if (free_inst_data)
        {
            cpd_free_calib_inst_data(dom_inst_data);
            cpd_free_calib_inst_data(for_inst_data);
        }

        if (free_und)
        {
            srt_f_destroy_und(dom_und_name);
            srt_f_destroy_und(for_und_name);
            srt_f_destroy_und(fx_und_name);
            destroy_correlation_list();
        }

        if (*prod_val)
        {
            free(*prod_val);
        }
    }

    if (!keep_calib_data)
    {
        if (*sig_time)
        {
            free(*sig_time);
            *sig_time = NULL;
        }
        if (*dom_sig)
        {
            free(*dom_sig);
            *dom_sig = NULL;
        }
        if (*for_sig)
        {
            free(*for_sig);
            *for_sig = NULL;
        }
        if (*fx_vol_time)
        {
            free(*fx_vol_time);
            *fx_vol_time = NULL;
        }
        if (*fx_vol)
        {
            free(*fx_vol);
            *fx_vol = NULL;
        }
        if (free_inst_data)
        {
            cpd_free_calib_inst_data(dom_inst_data);
            cpd_free_calib_inst_data(for_inst_data);
        }
    }

    if (!keep_und)
    {
        if (free_und)
        {
            srt_f_destroy_und(dom_und_name);
            srt_f_destroy_und(for_und_name);
            srt_f_destroy_und(fx_und_name);
            destroy_correlation_list();
        }
    }

    return err;
}
