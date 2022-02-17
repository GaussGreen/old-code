int princ_calc(int argc, UF_ARG *argv, UF_ARG *ret);

int princ_eq(int argc, UF_ARG *argv, UF_ARG *ret);

int tenor(int argc, UF_ARG *argv, UF_ARG *ret);
int agg_sensit(int argc, UF_ARG *argv, UF_ARG *ret);

int clean_sy(int argc, UF_ARG *argv, UF_ARG *ret);
int sy_yield(int argc, UF_ARG *argv, UF_ARG *ret);
int yield_sy(int argc, UF_ARG *argv, UF_ARG *ret);
int sy_clean(int argc, UF_ARG *argv, UF_ARG *ret);

#define SENSIT_UFDESC                                                          \
  {"princ_calc", "CGGVGVI", "D",                                               \
   "@princ_calc(yc  , dates range  , flow range  , [aggregation dates  , "     \
   "method(1/0)])",                                                            \
   princ_calc},                                                                \
      {"princ_eq", "I", "D", "@princ_eq(index)", princ_eq},                    \
      {"agg_sensit", "I", "D", "@agg_sensit(index)", agg_sensit},              \
      {"tenor", "I", "D", "@tenor(index)", tenor},                             \
      {"sy_clean", "IIDD", "D", "@sy_clean(start  ,end  ,coupon  ,sy)",        \
       sy_clean},                                                              \
      {"clean_sy", "IIDD", "D", "@clean_sy(start  ,end  ,coupon  ,clean)",     \
       clean_sy},                                                              \
      {"sy_yield", "IIWWDDZD", "D",                                            \
       "@sy_yield(start  ,end  ,compd  ,basis  ,cp  ,sy  ,rdmptn  "            \
       ",first_coupon)",                                                       \
       sy_yield},                                                              \
  {                                                                            \
    "yield_sy", "IIWWDDZD", "D",                                               \
        "@yield_sy(start  ,end  ,compd  ,basis  ,cp  ,yld  ,rdmptn  "          \
        ",first_coupon)",                                                      \
        yield_sy                                                               \
  }
