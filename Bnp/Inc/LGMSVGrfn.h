#ifndef LGMSVGrfnH
#define LGMSVGrfnH

#include "grf_h_mdlcomm.h"

typedef struct
{
    GRFNCOMMSTRUCT global;
    FIRSTMktAtT*   local;

    long    num_df;
    double* df_tms;
    long*   df_dts;

    /* For DF reconstruction */
    /* B(t,T) = B(0,T) / B(0, t) * Exp[dff - gam * f(t,T*) - gam2 * psi */
    double* dff;
    double* gam;
    double* gam2;

} grfn_parm_lgm_sv, *GRFNPARMLGMSV;

typedef struct
{
    GRFNCOMMSTRUCT global;
    FIRSTMktAtT*   local;

    long    num_df;
    double* df_tms;
    long*   df_dts;

    /* For DF reconstruction */
    /* B(t,T) = B(0,T) / B(0, t) * Exp[dff - gam1 * f1(t,T*) - gam2 * f2(t,T*)
                                                                                    - gam1_2  * ps1
       - gam2_2  * psi1 - gam12 * psi12 */
    double* dff;
    double* gam1;
    double* gam2;
    double* gam1_2;
    double* gam2_2;
    double* gam12;

} grfn_parm_lgm_sv_2F, *GRFNPARMLGMSV2F;

/****************************************************************
 *			MODEL UNDER QTstar g(t)	piecewise constant		    *
 ****************************************************************/

Err payoff_lgmsv_pde(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,
    double dTstar, /* In time */

    /* Grid data	*/
    int     lphi,
    int     uphi,
    int     lx,
    int     ux,
    int     leps,
    int     ueps,
    double* phi,
    double* ftTstar,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit,Xt,Epst,Product	*/
    double**** prod_val);

Err payoff_lgmsv_rec_swaption(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,
    double dTstar, /* In time */

    /* Grid data	*/
    int     lphi,
    int     uphi,
    int     lx,
    int     ux,
    int     leps,
    int     ueps,
    double* phi,
    double* x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit,Xt,Epst,Product	*/
    double**** prod_val);

/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

/************************************************************
 *						MODEL UNDER QBeta
 **
 ************************************************************/
Err payoff_lgmsv_pde_QBeta(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,

    /* Grid data	*/
    int     lphi,
    int     uphi,
    int     lx,
    int     ux,
    int     leps,
    int     ueps,
    double* phi,
    double* x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit,Xt,Epst,Product	*/
    double**** prod_val);

Err payoff_lgmsv_rec_swaption_QBeta(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,

    /* Grid data	*/
    int     lphi,
    int     uphi,
    int     lx,
    int     ux,
    int     leps,
    int     ueps,
    double* phi,
    double* x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit,Xt,Epst,Product	*/
    double**** prod_val);

/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/
/************************************************************************************************************************/

/********************************************************************************************
 *			MODEL UNDER QTstar u(t)=g(t)exp(-lambdaX*(Tstar-t))	piecewise constant
 **
 ********************************************************************************************/
Err payoff_lgmsv_pde_UtPieceWise(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,
    double dTstar, /* In time */
    double alpha,
    double rho,
    double ut,

    /* Grid data	*/
    int     lpsi,
    int     upsi,
    int     lx,
    int     ux,
    int     lz,
    int     uz,
    double* psi,
    double* x,
    double* z,

    /* Vector of results to be updated */
    int* nprod,
    /* 4 dimensions : Psit,Xt,Zt,Product	*/
    double**** prod_val);

Err payoff_lgmsv_rec_swaption_UtPieceWise(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,
    double dTstar, /* In time */

    /* Grid data	*/
    int     lpsi,
    int     upsi,
    int     lf,
    int     uf,
    int     lz,
    int     uz,
    double* psi,
    double* x,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Psit,Xt,Epst,Product	*/
    double**** prod_val);

/****************************************************************
 *						Monte Carlo
 **
 ****************************************************************/
Err payoff_lgmsv_mc(
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,

    double ft,
    double phi,
    double v,

    /* Vector of results to be updated */
    int nprod,
    /* Result	*/
    double* prod_val,
    int*    stop_path);

Err payoff_lgmsv_FFT(
    double evt_date,
    double evt_time,
    void*  func_parm,

    /* Market data	*/
    void* yc,

    /* Model data	*/
    double lamx,
    double dTstar, /* In time */

    /* Grid data	*/
    int iNbPhi,
    int iNbft,

    int iIndexPhiMean,
    int iIndexft0,

    double dPhiStep,
    double dftStep,
    double dPhitMean,

    /* Vector of results to be updated */
    int nprod,
    /* 4 dimensions : Phit,Xt,Epst,Product	*/
    double*** prod_val);

Err payoff_lgmsv2F_mc(
    long   path_index,
    double evt_date,
    double evt_time,
    void*  func_parm,

    double ft1,
    double ft2,
    double psi1,
    double psi2,
    double psi12,
    double v,

    /* Vector of results to be updated */
    int nprod,
    /* Result	*/
    double* prod_val,
    int*    stop_path);

#endif