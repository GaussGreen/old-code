#ifndef SRT_H_RNGFLTR_H
#define SRT_H_RNGFLTR_H

SrtErr srt_f_rngfltr(Date today, Date dcnt_start_date, Date start_date,
                     Date end_date, char *cap_floor_code, char *win_lose,
                     int lag, int end_lag, double strike, double float_fwd,
                     Date float_set, double correlation, double float_tenor,
                     Date *bus_day_array, double *rate_array, int ra_len,
                     double *vol_date_array, double *vol_array,
                     double *strikes_vec, double **mkt_vol_2darray,
                     long m_mat_dates, long m_strikes, int float_flg,
                     int smile_flg, int FHLB_flg, int reset_flg,
                     double *answer);

SrtErr srt_f_rngfltr_resettable(
    Date today, Date dcnt_start_date, Date start_date, Date end_date,
    char *win_lose_cap, char *win_lose_floor, int lag, int end_lag,
    double strike, double band_width, double float_fwd_fwd, double float_fwd,
    double spr_or_cpn, Date float_set, double correlation, double float_tenor,
    Date *bus_day_array, double *rate_array, int ra_len, double *vol_date_array,
    double *vol_array, double *strikes_vec, double **mkt_vol_2darray,
    long m_mat_dates, long m_strikes, int float_flg, int smile_flg,
    int FHLB_flg, double *answer);

#endif
