/* ==================================================================================

   FILENAME :  BGMEval.h

   PURPOSE :  declares the functions in BGMEval.c for the evaluation
                                of a Deal with a BGM model

   ==================================================================================
 */
#ifndef BGMEval_H
#define BGMEval_H
#include "grf_h_mdlcomm.h"

Err FirstValueDealInBGMModel(
    /*	GRFN Tableau	*/
    int init_num_evt_dts, long *init_evt_dts, long tab_rw, long tab_cl,
    char ***tab_str, int **mask, long aux_width, long *aux_len, double **aux,
    /*	Domestic underlying name	*/
    char *dom_nme,
    /*	GRFN param	*/
    SrtGrfnParam *prm,
    /*	Premium	*/
    double *price, double *stdev);

Err get_bgm_TenorDatesFromUnd(SrtUndPtr und_ptr, int numtenors,
                              long *Tenordates);

Err get_bgm_VolMatrixFromUnd(SrtUndPtr und_ptr, int numtenors, int numfactors,
                             double **VolMatrix);

Err get_bgm_HistCorrelFromUnd(SrtUndPtr und_ptr, int numtenors,
                              double **histcorrel);

Err get_bgm_FIRSTGetEventInfoFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	As output from FIRSTGetEvtDatesFromDeal	*/
    int num_evt,
    /*	As output from FIRSTGetUndFromDeal	*/
    int num_und,
    /*	As output from FIRSTGetMaxNumDfFromDeal	*/
    int max_num_df,

    /*	All pointers are allocated inside      ,
            and must be freed by call to get_bgm_FIRSTFreeEventInfoFromDeal
     */

    /*	Events themselves	*/
    FIRSTMktAtT **evt,

    /*	Wether American	*/
    int **am,

    /*	Information relative to df required for event evaluation	*/

    /*	num_df_mat[i][j] = number of df required for event i from underlying j
     */
    int ***num_df_mat,
    /*	df_mat_dts[i][j][k] and df_mat_tms[i][j][k]
            = maturity of the df number k required from underlying j at event i
     */
    long ****df_mat_dts, double ****df_mat_tms);

Err FirstValueBermudaDealOptimizationInBGMModel(
    /*	GRFN Tableau	*/
    int init_num_evt_dts, long *init_evt_dts, long tab_rw, long tab_cl,
    char ***tab_str, int **mask, long aux_width, long *aux_len, double **aux,
    /*	Domestic underlying name	*/
    char *dom_nme,
    /*	GRFN param	*/
    SrtGrfnParam *prm,
    /*	Premium	*/
    double *price, double *stdev, double *exfrontier);

double BermudaNewton(double epsmin, double epsmax, double tol,
                     int init_num_evt_dts, long *init_evt_dts, long tab_rw,
                     long tab_cl, char ***tab_str, int **mask, long aux_width,
                     long *aux_len, double **aux, char *dom_nme,
                     SrtGrfnParam *prm, double *stdev, int n1, int n2);

Err FirstValueBermudaDealDerivativesInBGMModel(
    /*	GRFN Tableau	*/
    int init_num_evt_dts, long *init_evt_dts, long tab_rw, long tab_cl,
    char ***tab_str, int **mask, long aux_width, long *aux_len, double **aux,
    char *dom_nme, SrtGrfnParam *prm, double *stdev, int nex, double *price,
    double *deriv, double *secondderiv);

double FirstValueBermudaDealInBGMModel(
    /*	GRFN Tableau	*/
    int init_num_evt_dts, long *init_evt_dts, long tab_rw, long tab_cl,
    char ***tab_str, int **mask, long aux_width, long *aux_len, double **aux,
    char *dom_nme, SrtGrfnParam *prm, double *stdev, int numex);

Err srt_f_BGMMidat(int num_evt, long *evt_dts, SrtUndPtr und_ptr, long numpaths,
                   int **num_df_mat, long ***df_mat_dts, double ***df_mat_tms,
                   double *dfs, int numpoints, double *exfrontier,
                   double Strike, SrtReceiverType PayRec, SrtBasisCode Basis,
                   double *pv);

Err BGMMidatGetCube(int num_evt, long *evt_dts, SrtUndPtr und_ptr,
                    long nummcpaths, int **num_df_mat, long ***df_mat_dts,
                    double ***df_mat_tms, double Strike, SrtReceiverType PayRec,
                    SrtBasisCode Basis, double ***SwapLevelCube);

Err OptimizeFrontierBGMMidat(double epsmin, double epsmax, int num_evt,
                             long numpaths, double *exfrontier, double *dfs,
                             double ***SwapLevelCube, int n1, int n2);

double ComputeBGMMidat(int num_evt, long numpaths, double *exfrontier,
                       double *dfs, double ***SwapLevelCube);

#endif
