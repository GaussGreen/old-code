#ifndef LGMSVPDEH
#define LGMSVPDEH

#include "math.h"
#include "srt_h_all.h"

#define LGMSV_TRUE 1
#define LGMSV_FALSE 0

/* Type of interpolation */

#define LGMSV_LINEAR 0
#define LGMSV_QUADRATIC 1
#define LGMSV_MIXQUADRALINEAR 2

/* ---------------------------------
         Default Value of LGMSVParam
   --------------------------------- */
#define LGMSV_IsQTstarModel 3

/* For X Grid */
#define LGMSV_IsNbSigmaMaxXGridRightAuto LGMSV_TRUE
#define LGMSV_iNbSigmaMaxXGridRight 15.0
#define LGMSV_iNbSigmaXGridRight 3.5

#define LGMSV_IsNbSigmaMaxXGridLeftAuto LGMSV_FALSE
#define LGMSV_iNbSigmaMaxXGridLeft 5.0
#define LGMSV_iNbSigmaXGridLeft 3.5

#define LGMSV_dPercentExtremePointXGrid 0.1
#define LGMSV_IsEqSpacedExtremePointXGrid LGMSV_FALSE

/* Use only for old Version */
#define LGMSV_IsNbSigmaXGridRightAuto 0

/* For Z Grid */
#define LGMSV_IsNbSigmaMaxZGridRightAuto LGMSV_FALSE
#define LGMSV_iNbSigmaMaxZGridRight 7.0
#define LGMSV_iNbSigmaZGridRight 5.0
#define LGMSV_iNbSigmaZGridLeft 5.0

#define LGMSV_dPercentExtremePointZGrid 0.1
#define LGMSV_IsEqSpacedExtremePointZGrid LGMSV_FALSE

#define LGMSV_TypeInterpolationPhi LGMSV_MIXQUADRALINEAR
#define LGMSV_IsExtrapolationFlat LGMSV_TRUE

/* For Phi Grid */
#define LGMSV_iNbSigmaPhiGrid 3
#define LGMSV_IsExtremePoints LGMSV_TRUE
#define LGMSV_iNbSigmaExtremePoints 4
#define LGMSV_IsPhiBoundExp LGMSV_TRUE

/* For Jump Treatment */
#define LGMSV_TypeInterpolationJump LGMSV_MIXQUADRALINEAR
#define LGMSV_IsExtrapolationFlatJump LGMSV_TRUE

#define LGMSV_VerifExpectation LGMSV_FALSE

/* For Debug */
#define LGMSV_ModeExport LGMSV_FALSE
#define LGMSV_SaveTensor LGMSV_FALSE

/* Use only for old Version */
#define LGMSV_VFloor 1e-4

/* For Calculation Optimisation */
#define LGMSV_GridXADIOptim LGMSV_TRUE
#define LGMSV_GridZADIOptim LGMSV_TRUE
#define LGMSV_GridPhiOptim LGMSV_TRUE

#define LGMSV_Tstar 10 /* In years */

/* -----------------------------------------------------------------------------------------------
        Definition of the structure LGMSVParam
   ------------------------------------------------------------------------------------------------
 */
typedef struct {
  int IsQTstarModel;

  /* For X Grid */
  int IsNbSigmaMaxXGridRightAuto;
  double iNbSigmaMaxXGridRight;
  double iNbSigmaXGridRight;

  double IsNbSigmaMaxXGridLeftAuto;
  double iNbSigmaMaxXGridLeft;
  double iNbSigmaXGridLeft;

  double dPercentExtremePointXGrid;
  int IsEqSpacedExtremePointXGrid;

  /* Use only for old Version */
  int IsNbSigmaXGridRightAuto;

  /* For Z Grid */
  int IsNbSigmaMaxZGridRightAuto;
  double iNbSigmaMaxZGridRight;
  double iNbSigmaZGridRight;
  double iNbSigmaZGridLeft;
  double dPercentExtremePointZGrid;
  int IsEqSpacedExtremePointZGrid;

  /* For Phi Grid */
  int TypeInterpolationPhi;
  int IsExtrapolationFlat;

  double iNbSigmaPhiGrid;
  int IsExtremePoints;
  double iNbSigmaExtremePoints;
  int IsPhiBoundExp;

  /* For Jump Treatment */
  int TypeInterpolationJump;
  int IsExtrapolationFlatJump;

  int VerifExpectation;

  /* For Debug */
  int ModeExport;
  int SaveTensor;

  /* Use only for old Version */
  double VFloor;

  /* For Calculation Optimisation */
  int GridXADIOptim;
  int GridZADIOptim;
  int GridPhiOptim;

  /*	For Monte Carlo */
  int UseBalsamGen;
  int UseReverseMC;
  int UseNewTStarMC;
  int iSchemeOrder;

  /*	Numerical param of Multiple Integral */
  double dMultiIntegMinTime;
  int iUseOldExpect;

  double Tstar;
} LGMSVParam, *LGMSVPARAM;

/* -----------------------------------------------------------------------------------------------
        Main function
        -------------


        Modelisation
        -------------

        Xt= rt - f(0      ,t)
        dXt =  [Phit-lambda_X*Xt]dt + sigma(Xt      ,t      ,omega) dWt
        dPhit = [sigma*sigma - 2*lambda_X*Phit] dt

        sigma : time      , spot dependant and stochastic of the following form
                sigma(Xt      ,t      ,omega)  = sig(t)*sigma* Epst

                Where
                        The stochastic part of the sigma : Epst follows :
                                dEpst = -lambda_Eps*[Epst-1]dt + alpha * Epst *
   dZt (mean reverting to 1      , log normal vol)

   ------------------------------------------------------------------------------------------------
 */

Err lgmSV_adi(
    /*	Time Information  */
    int iNbTime, double *dTime, double *dDate,

    /*	Space Discretisation	*/
    int iNbPhi, int iNbX, int iNbEps,

    /*	Model data Information	*/
    double dLambdaX, double dSigma,

    double *SigTime, double *Sig, int iNbSigTime,

    double dLambdaEps, double dAlpha, double dRho,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void **func_parm_tab, int *EvalEvent,

    /*	Market data */
    double *Ifr, char *cYieldCurve,

    /*	Payoff function */
    Err (*payoff_func)(/* Time Information */
                       double dDate, double dTime, void *func_parm,

                       /* Market data	*/
                       void *cYieldCurve,

                       /*	Model data Information	*/
                       double dLambdaX,

                       /* Martingale Mesure QTstar information */
                       double dTstar, /* In time */

                       /* Grid data Information	*/
                       int iIndPhitMin, int iIndPhitMax, int iIndXtMin,
                       int iIndXtMax, int iIndEpstMin, int iIndEpstMax,

                       double *GridPhi, double *GridftTstar,

                       /* Tensor of results to be updated		*/
                       /* 4 dimensions : Phit      ,Xt      ,Epst      ,Product
                        */
                       int iNbProduct, double ****PayoffTensor),
    /*	Result */
    int iNbProduct, double *dProductArrayPv);

Err lgmSV_adi_QBeta(
    /*	Time Information  */
    int iNbTime, double *dTime, double *dDate,

    /*	Space Discretisation	*/
    int iNbPhi, int iNbX, int iNbEps,

    /*	Model data Information	*/
    double dLambdaX, double dSigma,

    double *SigTime, double *Sig, int iNbSigTime,

    double dLambdaEps, double dAlpha, double dRho,

    /* Parameters */
    LGMSVParam Params,

    /*	Product data */
    void **func_parm_tab, int *EvalEvent,

    /*	Market data */
    double *Ifr, char *cYieldCurve,

    /*	Payoff function */
    Err (*payoff_func)(/* Time Information */
                       double dDate, double dTime, void *func_parm,

                       /* Market data	*/
                       void *cYieldCurve,

                       /*	Model data Information	*/
                       double dLambdaX,

                       /* Grid data Information	*/
                       int iIndPhitMin, int iIndPhitMax, int iIndXtMin,
                       int iIndXtMax, int iIndEpstMin, int iIndEpstMax,

                       double *GridPhi, double *GridX,

                       /* Tensor of results to be updated		*/
                       /* 4 dimensions : Phit      ,Xt      ,Epst      ,Product
                        */
                       int iNbProduct, double ****PayoffTensor),
    /*	Result */
    int iNbProduct, double *dProductArrayPv);

Err Fill_lgmSV_defaultParam(LGMSVPARAM Params);

Err lgmSV_adi_UtPieceWise(
    /*	Time Information  */
    int iNbTime, double *dTime, double *dDate,

    /*	Space Discretisation	*/
    int iNbPsi, int iNbX, int iNbZ,

    /*	Model data Information	*/
    double dLambdaX, double dSigma, /* No More Use */

    double *SigTime, double *Sig, int iNbSigTime,

    double *dLambdaEps, double *dLvlEps, double *dAlphaEps, double *dRho,

    /* Parameters */
    LGMSVPARAM Params,

    /*	Product data */
    void **func_parm_tab, int *EvalEvent,

    /*	Market data */
    double *Ifr, char *cYieldCurve,

    /*	Payoff function */
    Err (*payoff_func)(/* Time Information */
                       double dDate, double dTime, void *func_parm,

                       /* Market data	*/
                       void *cYieldCurve,

                       /*	Model data Information	*/
                       double dLambdaX,

                       /* Martingale Mesure QTstar information */
                       double dTstar, /* In time */

                       double dAlphaV, double dRho, double dSigt,

                       /* Grid data Information	*/
                       int iIndPsitMin, int iIndPsitMax, int iIndXtMin,
                       int iIndXtMax, int iIndZtMin, int iIndZtMax,

                       double *GridPsi, double *GridX, double *GridZ,

                       /* Tensor of results to be updated		*/
                       /* 4 dimensions : Psit      ,Xt      ,Zt      ,Product
                        */
                       int *iNbProduct, double ****PayoffTensor),
    /*	Result */
    int iNbProduct, int iMaxNbProduct, double *dProductArrayPv);

#endif