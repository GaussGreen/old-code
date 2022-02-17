// FxSVCalib.h : FXSV underlying description and calibration

#ifndef __FXSVCALIB_H__
#define __FXSVCALIB_H__

typedef struct _SFxSVUndDesc {
  long today;
  char dom_und[256], for_und[256];      // IR LGM1f underlyings with const tau
  double spot;                          // Spot FX rate
  long *dates;                          // Term
  double *times;                        // structure
  int ntimes;                           // of
  double *beta, *alpha, *gamma, ***rho; // model parameters
} SFxSVUndDesc;

// Memory allocation and copy construction:

Err InitFxSVUnd(SFxSVUndDesc *und, char *dom_und, char *for_und, double spot,
                long *dates, int ndates, double *beta, double *alpha,
                double *gamma, double ***rho);

// Destructor:
char *FreeFxSVUnd(void *ptr);

// Full underlying description with all the times merged:
typedef struct _SFxSVUndFull {
  double *times, *sigtms_d, *sigtms_f, *sig_d, *sig_f;
  double lam_d, lam_f;
  int ntimes, nsig_d, nsig_f;
  SFxSVUndDesc *pmdl;
} SFxSVUndFull;

Err FxSVFullFromUnd(SFxSVUndFull *o, SrtUndPtr und, long last_date);
Err FxSVFullFree(SFxSVUndFull *o);

// Structure for communication with NAG Runge-Kutta integrator:
typedef struct _SFxSVComm_ODE {
  double sig_d, sig_f, lam_d, lam_f;
  double beta, alpha, gamma, **rho;
  double T_pay, u_re, u_im;
  int want[4]; // Flags "want derivatives" : { beta  , alpha  , gamma  , rho }
  int is_cur_idx; // Flag showing whether the current piece is being calibrated
} SFxSVComm_ODE;

// Cache to store calculated integrand values:
typedef struct _SFxSVCache SFxSVCache;
struct _SFxSVCache {
  int maxpts;       // Max number of stored points in cache
  int cur_pt;       // Currently requested point
  double *v;        // Stored points
  double ***fn;     // Stored integrand values: per v  , per k  , per fn
  SFxSVCache *next; // Linked list
};

Err FxSVInitCache(SFxSVCache *cache, int maxpts, int nk);
Err FxSVFreeCache(SFxSVCache *cache, int nk);

// Structure for communication with a quadrature integrator
typedef struct _SFxSVComm_InvFT {
  SFxSVUndFull *o;
  double T_pay; // Option pay time
  int idx_from; // Index of the first piece being calibrated
  int idx_to;   // Index of the last piece being calibrated
  int want[4];  // Flags "want derivatives" : { beta  , alpha  , gamma  , rho }
  double v_max; // Upper limit of integration

  int nk;      // Number of strikes
  int want_k;  // Current strike
  int want_fn; // Current function (integrand or its derivative)
  double *k;   // Strikes ( log(K/FFX0) )

  // Integrand cache:
  SFxSVCache *c_head, *c_cur;

  // Counter for function calls:
  int count;

} SFxSVComm_InvFT;

Err FxSVDensityFT(char *fxundname, long fix_date, long pay_date, double u_re,
                  double u_im, double *h_re, double *h_im, int *want, int idx);

Err FxSVOptions(char *fxundname, long fix_date, long pay_date, int nK,
                double *K, char **rec_pay_str, int *want, int idx,
                double **res);

typedef struct _SFxSVComm_Calib {
  SFxSVUndFull *o;
  double *T_pay; // Option pay times (per ex date)
  int idx_from;  // Index of the first piece smile params of which are being
                 // calibrated
  int idx_to;    // Index of the last piece smile params of which are being
                 // calibrated
  int **calib; // Flags "calibrate parameters" : { beta  , alpha  , gamma  , rho
               // } (per ex date)

  int *nk;          // Number of strikes (per ex date)
  double **k;       // Strikes ( log(K/FFX0) ) (per ex date)
  double **mkt_val; // Market prices (per ex date)

  double bl[4],
      bu[4]; // Lower and upper bounds { beta  , alpha  , gamma  , rho }

  int ntimes_last;

} SFxSVComm_Calib;

Err FxSVCalibrate(SrtUndPtr und, // Starting points and params not calibrated
                                 // must be initialized
                  int **calib,   // Flags calibrate { beta  , alpha  , gamma  ,
                                 // rho } (per ex date)
                  int *nK,       // Number of strikes per exercise date
                  double **K,    // Strikes per exercise date
                  double **vol); // Vols per exercise date per strike

#endif // #ifndef __FXSVCALIB_H__
