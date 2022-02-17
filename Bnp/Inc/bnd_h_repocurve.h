/********************************************************************************
  File name: bnd_f_repocurve.h
********************************************************************************/

Err srt_f_GCRepoCurve(long TradeDate, int N, long *TermDates, double *TermRates,
                      long *FwdDates, double *DailyFwdRates,
                      double *DailyTermRates);

Err srt_f_SpecialRepoCurve(long TradeDate, int N, long *TermDates,
                           double *TermRates, long *FwdDates,
                           double *DailyFwdRates, double *DailyTermRates);