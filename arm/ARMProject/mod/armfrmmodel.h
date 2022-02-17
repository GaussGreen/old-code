/*
 * $Log: armfrmmodel.h,v $
 * Revision 1.39  2004/05/12 14:02:50  emezzine
 *  Added a new parameter to control vega perturbation
 *
 * Revision 1.38  2004/04/07 15:45:09  emezzine
 * Added new param Presicion.
 *
 * Revision 1.37  2004/03/31 15:27:33  rguillemot
 * CMT Bug Fix
 *
 * Revision 1.36  2004/03/19 15:49:38  emezzine
 * Added new function to pre initialize beta curve.
 *
 * Revision 1.35  2004/02/17 10:56:34  emezzine
 * Added Bootsrap2D
 *
 * Revision 1.34  2004/02/16 14:03:07  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.33  2004/02/09 18:12:46  emezzine
 * Modif calibrate().
 *
 * Revision 1.32  2004/02/09 08:54:08  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.31  2004/02/04 15:17:27  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.30  2004/01/26 13:35:17  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.29  2004/01/12 07:13:19  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.28  2003/12/16 16:11:04  arm
 * Added : class ARM_Swaption
 *
 * Revision 1.27  2003/12/16 15:44:29  emezzine
 * added computeStubSpot funcion
 *
 * Revision 1.26  2003/10/31 08:11:22  emezzine
 * New version to calibrate LogProba Curve.
 *
 * Revision 1.25  2003/10/13 11:22:58  emezzine
 * improvement
 *
 * Revision 1.40  2003/10/10 11:44:09  sgasquet
 * Ajout fonction EuroDigitale
 *
 * Revision 1.39  2003/10/07 13:18:36  emezzine
 *  Correct a bug in calibration(initialise sigma if beta non null)
 *
 * Revision 1.23  2003/10/06 15:51:33  rguillemot
 * Ajout d'un spread dans le modele FRM (Modele Mixture)
 *
 * Revision 1.22  2003/09/29 08:12:29  emezzine
 * Modif  name
 *
 * Revision 1.20  2003/09/18 05:28:54  ebenhamou
 * added the correct function
 *
 * Revision 1.19  2003/09/17 13:10:36  emezzine
 * Added Beta Calibration " Beta is a probability of lognormality of SFRM ".
 *
 * Revision 1.18  2003/07/25 14:58:21  emezzine
 *  Add the 3th portfolio and its vectors
 *
 * Revision 1.17  2003/07/11 12:00:29  emezzine
 * Cut builVolatility()
 *
 * Revision 1.16  2003/07/07 09:21:05  emezzine
 * Add itsError2
 *
 * Revision 1.15  2003/06/23 09:55:51  emezzine
 * create a new face for constructor and remove itsYieldCurve
 *
 * Revision 1.14  2003/06/18 09:12:16  emezzine
 * Add CalibParams
 *
 * Revision 1.13  2003/06/12 09:08:51  jmprie
 * ajout du calage de la mean reversion
 *
 * Revision 1.12  2003/06/10 17:18:09  emezzine
 * merge armfrmstructure, armfrmanalyticmodel and armfrmmodel in one class
 *
 * Revision 1.11  2003/05/26 09:42:36  emezzine
 * modif constrictor
 *
 * Revision 1.10  2003/04/28 14:58:17  emezzine
 * Ajout de view()
 *
 * Revision 1.9  2003/03/28 13:14:25  emezzine
 * Use EuroCaplet et EuroSwaption to calibrate
 *
 * Revision 1.8  2003/03/26 15:31:43  mab
 * Recup. suite incident bizarre
 *
 * Revision 1.6  2003/03/24 16:05:13  mab
 * #define finite
 *
 * Revision 1.5  2003/03/24 14:39:50  emezzine
 * Added EuroSwaption:
 *
 * Revision 1.4  2003/03/19 14:49:01  jmprie
 * nouveau constructeur + EuroCaplet() et ExpectedFwdYield()
 *
 */
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* armfrmmodel.h: interface for the ARM_FRMModel class.                     */
/*                                                                          */
/*--------------------------------------------------------------------------*/
#ifndef _ARMFRMMODEL_H
#define _ARMFRMMODEL_H

#include "model.h"

#define SIGMACALONLY 0
#define SIGMAANDSHIFTCALONLY 1
#define SIGMAANDMEANREVCALONLY 2
#define ALLPARAMETERSCAL 3


class ARM_FRMVol; 
class ARM_Security;
class ARM_VolCurve;
class ARM_Portfolio;
class ARM_Swaption;
class ARM_QuasiNewton;
class ARM_Model;

class ARM_FRMModel;

typedef void (*pFunctionToCalibrate)( ARM_FRMModel* const myObject, double& lambda );
void NAGFunctionToCalibrate( ARM_FRMModel* const myObject, double& lambda );


class ARM_FRMModel : public ARM_Model  
{
	private:
          
            ARM_Portfolio* itsPfToCalibrateSigma;
            ARM_Portfolio* itsPfToCalibrateShift; 
            ARM_Portfolio* itsPfToCalibrateMRS; 
            ARM_Vector* itscalibParams;
            ARM_QuasiNewton* itsOptimizer;

            ARM_FRMVol* itsFRMVol;

            ARM_Vector* itsOneMu;
            ARM_Vector** itsMu;
            int       itsNbMu;                  

            int itsCalibrationMode;
            int itsSizeToNag;
            int itsIsCalibrated;
            int itsVolFraSwap;	
            int itsMFlag; 
            int itsIndexProduct;
            double itsPresicion;  
            double itsVegaLevel;

            ARM_Vector* itsErrorSigma;
            ARM_Vector* itsErrorShift;
            ARM_Vector* itsErrorMRS;

            double localVolatility(double StartDate,double EndDate,ARM_Security* Data,ARM_Vector* Mu);

	public:

            //Constructeur, Destructeur 
            ARM_FRMModel(void);
            ARM_FRMModel(const ARM_FRMModel& Data);
            ARM_FRMModel(ARM_Date& asOfDate,
                        ARM_ZeroCurve* zc,
                        ARM_FRMVol*FRMVol,
                        ARM_Portfolio* Pf1,
                        ARM_Portfolio* Pf2 = NULL,
                        ARM_Portfolio* Pf3 = NULL,
                        ARM_Vector* calibParams = NULL,
                        long paramToCal = 0,
                        ARM_Model* model = NULL,
                        size_t WithPreInitialise=0,
                        double Presicion = 1.0e-3,
                        double VegaLevel = 1.0e-10);

            ~ARM_FRMModel(void);
            void FreeVectors(void);
            inline ARM_FRMVol* GetFRMVol() const {return itsFRMVol;}
            inline void SetFRMVol(ARM_FRMVol* frmvol) {itsFRMVol = frmvol;}

            void PreInitialiseBetaCurve(ARM_Model* model);
            ARM_Portfolio* CreatePortfolioWithShiftedStrike(const ARM_Portfolio* Pf,
                                                            ARM_Model* model,
                                                            double shift = 0.1);
            void CalibrateBetaAndSigma(ARM_Vector* var, 
                                         double min , 
                                         double max);  
            double TargetFunction(double* x);
            double RootFinding();
            void Bootstrap2D();
     
            ARM_FRMVol* GetFRMVol(void);			  
          
            double SimpleFuncToMinimize(double x);
            void Calibrate(void);
            void Bootstrap1D();
            void CalibrateMeanReversion();
            void CalibrateLogProbaCurve(void);
            
            void SetMu(ARM_Portfolio* FRMPortfolio) ;
            ARM_Vector** GetMu(void);
            void SetNbMu(int nbMu);
            ARM_Vector* GetMu(size_t i);
            void SetMu(ARM_Vector* Mu) ;
            void SetOneMu(ARM_Vector* Mu);
			void SetGlobalMu(ARM_Vector** Mu, int NbMu);
            ARM_Vector* ComputeMu(ARM_Security* Data);
            
            double Mdec(ARM_Security* Data);
			double Spread(ARM_Security* Data);

            void SetVolFraSwap(int flag);
			void SetFlagToCalibrate(int flag);

            double Price(ARM_Security* ProductToPrice);
            double Volatility(ARM_Security* Data,ARM_Vector* Mu);  
            double ComputeStubSpot(ARM_Swaption* Data);


            double ExpectedFwdYield(double startDate,double endDate,
                            double payDate,int compMeth,int dayCount,
                            int domOrFrg,int discYC,
							int YieldDecomp = 1,
							double Margin = 0.0,
							StoreFwdRateInfo* StoreInfo = NULL,
							int IndexFreq = -1);

            double EuroCaplet(double settlement,double resetMaturity,
                              double startMaturity,double endMaturity,
							  double payMaturity,
                              int optionType,double fwd,double strike,
                              int DomOrFrg = 1,
                              int compMeth = 0, 
                              int dayCount = KACTUAL_365,
							  bool IsCMS = false,
                              bool isTreasury = false,
							  int FixFrequency = 1,
							  double UnderlyingTenor = 0.0,
							  int YieldDecomp = K_COMP_PROP,
							  double Margin = 0.0,
							  ARM_Currency* ccy = NULL,
							  StoreFwdRateAndCapletInfo* StoreInfo = NULL);

			double EuroDigitale(double settlement,double resetMaturity, 
								double startMaturity,double endMaturity,
								int optionType,double fwd,double strike,
								int DomOrFrg = 1,
                                int compMeth = 0,
                                int dayCount = KACTUAL_365);

            double EuroSwaption(ARM_Security* sec, double startswap,
                                double matswaption, int optionType,
                                double swapfwd, double Strike,
                                ARM_Vector* Forwds, ARM_Vector* poids,
                                int DomOrFrg=1);
            // Standard functions
            void Init(void);
            void BitwiseCopy(const ARM_Object *ARM_FRMModel);
            void Copy(const ARM_Object *ARM_FRMModel);
            ARM_Object* Clone(void);
            void View(char* id, FILE* ficOut);
            void GlobalCalibrate(double Tol, size_t MaxIter, size_t GradCal){};

            ARM_Vector* GetResult(void);
            inline long GetParamToCal() const { return itsMFlag;}
            inline int GetIndexProduct() const { return itsIndexProduct;}
            
            ARM_Portfolio *Get2ndPortfolio();
            double EvaluatePortfolioError(double x);
            double EvaluateFunction(int n,double* x);

            void UpDateErrors(void);

			/// for a more generic calibration

			/// accessors for the calibration of the mean reversion portfolio
            ARM_Portfolio* GetPfToCalibrateMRS() const { return itsPfToCalibrateMRS; }
            void callRPObj(int *nfct, int *sizex, double *x, int *spadim, double *grid, int *liudat, int *iudata, int *ldudat, double *dudata, double *funct, double *jacob, int *info) ;
            ARM_Vector* GetCalibParams() const { return itscalibParams; }

			/// the static function to calibrate
			static pFunctionToCalibrate FUNCTIONTOCALIBRATE;
};


#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
