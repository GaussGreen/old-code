/****************************************************************************/
/*      C interface for DR prepay model                              .      */
/****************************************************************************/
/*      DRIntfc_c.h                                                        */
/****************************************************************************/

#include <ppdrfrm.h> /* SRM3 DR prepay header */

#ifdef __cplusplus
extern "C" {
#endif

void* DR_create_handle(
    PP_POOL*,       /* underlying pool */
    PP_DRFRM*,      /* DR fixed rate model parameters */
    char* ratePath, /* path of file that contains DRFRM historical data */
    long today);    /* the date today */

/* calculate CPR; calcs are as of BOM(simDate)*/
int DR_Cpr(
    void* handle,   /* handle created by DR_create_handle(...) */
    long simDate,   /* the simulation date YYYYMMDD */
    PP_POOL* ,      /* the current pool */
    int sz,         /* # rates. should equal #mos from today to simDate */
    double* r,      /* the diffused mortgage rates */
    double* cpr);   /* (O) */

void DR_destroy_handle(void*);

void DR_get_error( char*, int sz);

#ifdef __cplusplus
}
#endif
