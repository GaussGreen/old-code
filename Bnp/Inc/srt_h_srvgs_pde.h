
#include "SRT_H_ALL.H>

Err srt_f_stoch_rate_and_vol_gamma_smile_pde(
    Date event_date,

    double max_time, long min_node, long min_num_mesh,

    char *und_name, char *dom_und_name, char *yc_name,

    double equity_strike, char *rec_pay_str, char *greek_str, double *greeks);
