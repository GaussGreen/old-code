/*
 * $Log: bsmodel.h,v $
 * Revision 1.35  2004/04/07 07:17:08  rguillemot
 * Smile interpole
 *
 * Revision 1.34  2004/04/05 13:20:03  mab
 * in : void SetZeroCurve(ARM_ZeroCurve* zc)
 * corrected : ARM_Model::SetZeroCurve(zc);
 *
 * Revision 1.33  2004/04/05 12:25:44  mab
 * Added : void SetZeroCurve(ARM_ZeroCurve* zc)
 *
 * Revision 1.32  2004/04/02 10:02:29  mcampet
 * MC remove itsDiscounteRate
 *
 * Revision 1.31  2004/03/31 15:27:51  rguillemot
 * CMT Bug Fix
 *
 * Revision 1.30  2004/03/22 08:36:16  mab
 * Added: void SetCvxAdjVolatility(ARM_VolCurve* newVol);
 *
 * Revision 1.29  2004/03/05 17:25:04  mab
 * Added : void SetDiscPricingMode(int discPricingMode)
 *
 * Revision 1.28  2004/02/26 14:49:06  rguillemot
 * Discount Curve in BSModel
 *
 * Revision 1.27  2004/02/26 10:05:09  mab
 * Formatting
 *
 * Revision 1.26  2004/02/17 14:14:14  mcampet
 * MC add View Method
 *
 * Revision 1.25  2004/02/16 13:58:55  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.24  2004/02/13 14:48:55  mcampet
 *  MC add DomorFgn  in ExpectedSwapRate
 *
 * Revision 1.23  2004/02/09 08:54:32  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.22  2004/02/04 15:19:14  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.21  2004/02/02 10:59:13  mab
 * A comment added
 *
 * Revision 1.20  2004/02/02 10:52:28  mab
 * Just formatting!
 *
 * Revision 1.19  2004/01/26 13:43:28  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.18  2004/01/21 11:41:03  emezzine
 *  new version after merge with extendedBS
 *
 * Revision 1.16  2004/01/05 10:38:15  arm
 * MA: Just formatting
 *
 * Revision 1.15  2003/12/23 16:59:05  rguillemot
 * Conv Adj Merge
 *
 * Revision 1.14  2003/12/22 09:51:13  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.13  2003/12/09 08:37:25  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.12  2003/10/13 10:28:04  emezzine
 *  A new virtual function "ExpectedSwapRate"
 *
 * Revision 1.11  2003/10/10 15:38:39  mab
 * Added: void SetDiscountCurve(..)
 *
 * Revision 1.10  2003/09/22 15:32:50  mab
 * Added : double ComputeVol(..)
 *
 * Revision 1.9  2003/03/19 16:47:42  mab
 * Added : SetType(int aType)
 *
 * Revision 1.8  2003/02/11 16:49:32  mab
 * a parameter : double volTenor = 0 Added in bsOption
 *
 * Revision 1.7  2002/11/15 10:26:27  mab
 * after SetName(ARM_BSMODEL); SetBSCategory(1); added
 *
 * Revision 1.6  2002/10/09 09:33:46  mab
 * added key word const in some constructors
 *
 * Revision 1.5  2002/02/25 13:55:08  mab
 * added : virtual double bsOption(double spot, double strike, ....)
 *
 * Revision 1.4  2001/07/18 11:44:12  abizid
 * Ajouts et Corrections pour le calcul des Expected
 *
 * Revision 1.3  2000/02/01 17:14:37  nicolasm
 * Ajout ExpectedForwardYield
 *
 * Revision 1.2  1999/02/25 09:38:27  vberger
 * Rajout d'un constructeur ARM_BSModel(ARM_ZeroCurve* zerocurve, ARM_VolCurve* volatility,int volType = K_INDEX_RATE);
 *
 */

/*----------------------------------------------------------------------------*
 bsmodel.h
 
 Header for the ARM_BSModel class, an class for the Black & Scholes framework
*----------------------------------------------------------------------------*/
#ifndef _BSMODEL_H
#define _BSMODEL_H



// Includes
#include "bsxtic.h"
#include "model.h"
#include "ycmodel.h"

#include "foncstd.h"


class ARM_Date;
class ARM_Vector;
class ARM_ConvAdjustManager;
class ARM_CorrelManager;
class ARM_VolCurve;
class CapletStoreInfo;
class ARM_BSSmiledModel;



class ARM_BSModel : public ARM_Model
{
    private:

        double itsSpot; // initial (at start date) spot price 
        int    itsType; // (Yield or Price) used for option valuation with B&S


        //---  Not-Shared  Objects : 

        ARM_ZeroCurve*     itsDividendYield;    // Continuous dividend yield

        ARM_VolCurve*      itsVolatility;       // Volatility in % (12 for 12%)
        ARM_VolCurve*      itsCvxAdjVolatility; // Volatility in % used for 
                                                // convexity adjustment
        ARM_VolCurve*      itsSpreadLock;       // Spread lock in BP
        ARM_YCModel*       itsYCModel;          // Used when bsmodel is propagated

		ARM_BSSmiledModel* itsSabrModel;        // SABR model if needed


        //-- Variables for spreadOption ---

		ARM_VolCurve*  itsCorrelationMatrix; // Correlation
        ARM_VolCurve*  itsCashVolMatrix;     // Leg CMT
        ARM_VolCurve*  itsSpreadVolCurve; 

        int            itsModelType;         // 2LOG, LOG, NOR
		int            itsSpreadVolType;     // K_INPUTED,K_COMPUTED

		ARM_Vector*    itsCalibInfos;        // for on the fly marginal calibration
		int            itsNumSteps;          // for numerical formulae

        void Init(void);
        
    public:

       // Constructors and destructor
       ARM_BSModel(void);
	   
	   ARM_BSModel(ARM_Date& AsOfDate, 
                   ARM_ZeroCurve* YielCurve,
                   ARM_VolCurve* SpreadLockCurve,          // Spread vol  
                   ARM_VolCurve* ConvAdjVolatility = NULL, // IRG
                   ARM_VolCurve* Volatility = NULL,        // Pricing vol : IRG: 
                                                           // Or Swopt : case CMS
                   ARM_CorrelManager* CorrelManager = NULL,
                   ARM_ConvAdjustManager* ConvAdjustManager = NULL,
                   ARM_ZeroCurve* DiscountCurve = NULL);

       // BSModelGEN constructor
	   ARM_BSModel(ARM_Date& AsOfDate,
                   ARM_ZeroCurve* YieldCurve,
                   ARM_VolCurve* SpreadLockCurve,
                   ARM_VolCurve* convAdjVolatility,
                   ARM_VolCurve* priceVolatility,                         
                   ARM_CorrelManager* CorrelManager,
                   ARM_ConvAdjustManager* ConvAdjustManager,
                   ARM_ZeroCurve* DiscountCurve,
				   ARM_VolCurve* Correl,          // Correlations
				   ARM_VolCurve* CashVol,         // Cash Vol (particulièrement pour CMT)
				   ARM_VolCurve* SpreadVolCurve,
				   int ModelType = K_2LOG,		  // 2LOG, LOG, NOR
				   int SpreadVolType= K_COMPUTED, // K_INPUTED,K_COMPUTED
				   ARM_BSSmiledModel* sabrMod = NULL,
				   ARM_Vector* calibInfos = NULL,
				   int NumSteps = 0); 


       ARM_BSModel(ARM_Date& startDate,
                   double spot, 
                   ARM_ZeroCurve* dividend,
                   ARM_ZeroCurve* discountRate,					// forecast rates curve
                   ARM_VolCurve* volatility,
                   int volType = K_PRICE,
				   ARM_ZeroCurve* realDiscountRate = NULL);		// discount curve
 
       ARM_BSModel(ARM_Date& startDate, 
                   double spot, 
                   double flatdiv,
                   double flatdiscrate,
                   double flatvol,
                   int volType = K_PRICE);

        ARM_BSModel(ARM_ZeroCurve* zerocurve,
                    ARM_VolCurve* volatility,
                    int volType = K_INDEX_RATE);

        virtual bool IsNormal(void)
		{
            return(false);
		}

        ARM_BSModel(const ARM_BSModel& bs);    
        ARM_BSModel& operator = (const ARM_BSModel& bs);    

        virtual ~ARM_BSModel(void);

        void SetDiscPricingMode(int discPricingMode)
        {
            ARM_Model::SetDiscPricingMode(discPricingMode);

            if (itsYCModel)
               itsYCModel->SetDiscPricingMode(discPricingMode);
        }
              
        void SetZeroCurve(ARM_ZeroCurve* zc)
        {
            ARM_Model::SetZeroCurve(zc);

            if (itsYCModel)
            {
               ARM_ZeroCurve* clonedZC = NULL;

               if (zc)
                  clonedZC = (ARM_ZeroCurve *) zc->Clone();

               delete itsYCModel->GetZeroCurve();
               itsYCModel->SetZeroCurve(clonedZC);
            }

        }

        void BitwiseCopy(const ARM_Object* srcBSModel);   
        void Copy(const ARM_Object* srcBSModel);
        virtual ARM_Object* Clone(void);
        virtual void View(char* id = NULL, FILE* ficOut = NULL);

        // return the convexity adjustment of a CMS, including convexity adj + timelag
        double ConvexAdjustCMS(double cmsvalue,
                               double Tr,
                               double Ts, double Te, double Tp,
                               long frequency, int dayCount,
                               int YieldDecompFreq,
                               double Margin,
                               double& volSwap,
                               StoreFwdRateInfo* StoreInfo = NULL);
		
		// return the convexity adjustment of CMS + CMS Quanto , but not include the timelag
		double ConvAdjCMSWithQuanto(double cmsvalue, 
                                    double ResetJul,
								    double StartJul, double EndJul,
								    long frequency, int dayCount,
								    int YieldDecompFreq,
								    double Margin,
								    double paymentlag,
								    double& volSwap,
								    StoreFwdRateInfo* StoreInfo = NULL,
								    int DomOrFrg = 1,
  								    ARM_VolCurve* fFxCorr = NULL,
						            ARM_VolCurve* FxVol = NULL);

		// return the timelag adjustment
		double PaymentLagAdjCMS(double cmsvalue, 
                                double ResetDateJul,
								double StartDateJul, double EndDateJul, double PayDateJul,
								long frequency, int dayCount,
								int YieldDecompFreq,
								double Margin,
								StoreFwdRateInfo* StoreInfo = NULL);

		double PaymentLagAdjCMSQuanto(double volswapF, 
									  double CorrDF, 
									  double ResetDateJul, 
									  double StartDateJul, 
									  double EndDateJul, 
									  double PayDateJul,
									  int frequency, 
									  int FwdRule, 
									  double& volswapD,
									  int TypeStub = K_LONGSTART,
									  int intRule = K_ADJUSTED, 
									  int adjFirstdate = 1,
									  int dayCount = KACTUAL_365,
									  bool isTreasury = false);

		// Compute decapitalized Forward
		double CptForwardDecap(double fwd, 
							   double spread, 
							   double vol, 
							   double decompFreq, 
							   double timelag);

        // To calculate the swap rate in model non in security
        virtual double SwapRate(ARM_Date& StartDate,
                                ARM_Date& EndDate,
                                int frequency, 
                                int FwdRule,
                                ARM_Currency* ccy,
                                int TypeStub = K_LONGSTART, 
                                int intRule = K_ADJUSTED,
                                int adjFirstdate = 1,
                                int dayCount = KACTUAL_365,
                                bool isTreasury = false);

        virtual double ExpectedSwapRate(ARM_Date& resetDate,
                                ARM_Date& StartDate,
                                ARM_Date& EndDate,
                                ARM_Date& PayDate,
                                int frequency,
                                int FwdRule,
                                int YieldDecompFreq,
                                double Margin,
                                ARM_Currency* ccy,
                                int TypeStub = K_SHORTSTART, 
                                int intRule = K_ADJUSTED,                          
                                int adjFirstdate = 1,
                                int dayCount=KACTUAL_365,
                                int DomOrFrgRate = 1,
                                bool isTreasury = false,
                                StoreFwdRateInfo* StoreInfo = NULL);

        void Set(double spot, 
                 ARM_ZeroCurve* dividend, 
                 ARM_ZeroCurve* discountRate,			// forecast rates curve
                 ARM_VolCurve* volatility,
                 int type,
				 ARM_ZeroCurve* realDiscountRate = NULL);		// discount curve

		void Set(ARM_ZeroCurve* zeroCurve,
                 ARM_ZeroCurve* discountCurve,
                 ARM_VolCurve* spreadLock, 
                 ARM_VolCurve* capVol, 
                 ARM_VolCurve* indexVol = NULL,
				 ARM_VolCurve* Correlations = NULL,  
				 ARM_VolCurve* CashVol = NULL,
				 ARM_VolCurve* SpreadVolCurve= NULL,
				 int ModelType = K_2LOG, /* K_2LOG */
				 int spreadVolType = K_COMPUTED);
		        
		inline double GetSpot(void) const 
        { 
            return itsSpot;
        }

		inline void SetSpot(double spot) 
        { 
            itsSpot=spot;
        }

        inline int GetType(void) const
        {
            return itsType;
        }

        virtual ARM_VolCurve* GetVolatility(int mode = -1) const
        { 
            return(itsVolatility);
        }

        inline ARM_VolCurve* GetCvxAdjVolatility(void) const 
        { 
            return itsCvxAdjVolatility;
        }

        inline ARM_VolCurve* GetSpreadLock(void) const 
        { 
            return itsSpreadLock;
        }

        inline ARM_ZeroCurve* GetDividend(void) const 
        { 
            return itsDividendYield;
        }

        inline ARM_YCModel* GetYCModel(void) const 
        {
            return(itsYCModel);
        }
        
		inline ARM_VolCurve* GetCorrelationMatrix(void) const
		{
			return (itsCorrelationMatrix);
		}

        inline ARM_VolCurve* GetCashVolMatrix(void) const
		{
			return(itsCashVolMatrix);
		}

		inline ARM_VolCurve* GetSpreadVolCurve(void) const
		{
			return(itsSpreadVolCurve);
		}
		inline int GetModelType(void) const
		{
			return(itsModelType);
		}

		inline int GetSpreadVolType(void) const
		{
			return(itsSpreadVolType);
		}

		inline ARM_BSSmiledModel* GetSabrModel(void) const
		{
			return(itsSabrModel);
		}

		inline ARM_Vector* GetCalibInfos(void) const
		{
			return(itsCalibInfos);
		}

        void SetType(int aType) 
        {
            itsType = aType;
        } 

        void SetCvxAdjVolatility(ARM_VolCurve* newVol);

        void SetVolatility(ARM_VolCurve* newVol);

		void Set_XXX_Volatility(ARM_VolCurve* vol);

		void SetCorrelationMatrix(ARM_VolCurve* newCorrelationMatrix);

        void SetCashVolMatrix(ARM_VolCurve* newCashVolMatrix); 
        
		void SetSpreadVolCurve(ARM_VolCurve* newSpreadVolCurve);

		void SetModelType(int aModelType)
		{
			itsModelType = aModelType ; // 2LOG, LOG, NOR
		}

		void SetSpreadVolType(int aSpreadVolType)
		{
			itsSpreadVolType = aSpreadVolType ; // 2LOG, LOG, NOR
		}		

		void SetSabrModel(ARM_BSSmiledModel* SabrMod);

		void SetCalibInfos(ARM_Vector* calibInfos)
        {
            if (itsCalibInfos)
               delete itsCalibInfos;

            if (calibInfos)    
               itsCalibInfos = (ARM_Vector *) calibInfos->Clone();
        }

		void SetNumSteps(int NumSteps)
		{
			itsNumSteps = NumSteps;
		}

		int GetNumSteps(void)
		{
			return(itsNumSteps);
		}

        ARM_Model* Shift(int param, ARM_BucketShift* h);

virtual double ComputeVol(double matu,
                          double tenor,
                          double fwd,
                          double strike,
						  int mode = -1); // mode could be K_LIBOR, K_CMS, K_IRG, K_SWOPT, 
										  // but in this model, only two cases: K_LIBOR and all the others!!
		                                  // K_LIBOR implies using SABR if available otherwise use standard BS
        
        double ComputeCvxAdjVol(double matu,
                                double tenor,
                                double fwd, 
                                double strike);

        double ZeroPrice(double calcDate,
                         double zMaturity,
                         int DomOrFrg = 1);

        double ZeroFwdPrice(double calcDate,
                            double fwdDate, 
                            double zMaturity, 
                            int discYC = 1); 

        double CFPrice(double calcDate,
                       ARM_Vector* cfTerms, 
                       ARM_Vector *cfValues,
                       int discYC = 1);

        double CFFwdPrice(double calcDate,
                          double fwdDate, 
                          ARM_Vector* cfTerms,
                          ARM_Vector* cfValues, 
                          int discYC = 1);
        
        virtual double ForwardYield(double calcDate,
                                    double resetDate, 
                                    double maturityDate,
                                    double yieldMaturity, 
                                    int compMeth = 0, 
                                    int dayCount = KACTUAL_365, 
                                    int DomOrFrgRate = 1); 

        virtual double ExpectedFwdYield(ARM_Date& fwdDate,
                                        ARM_Date& maturity, 
                                        ARM_Date& payDate,
                                        int compMeth = 0, 
                                        int dayCount = KACTUAL_365,
                                        int DomOrFrgRate = 1, 
                                        int discYC = 1,
                                        int YieldDecomp = K_COMP_PROP,
                                        double Margin = 0.0,
                                        StoreFwdRateInfo* StoreInfo = NULL,
										int indexFreq = -1);

        virtual double ExpectedFwdYield(double fwdJulDate,
                                        double maturity,
                                        double payJulDate,
                                        int compMeth = 0,
                                        int dayCount = KACTUAL_365,
                                        int DomOrFrgRate = 1, 
                                        int discYC = 1,
                                        int YieldDecomp = K_COMP_PROP,
                                        double Margin = 0.0,
                                        StoreFwdRateInfo* StoreInfo = NULL,
										int indexFreq = -1);

		virtual double ExpectedFwdYieldNonAdj(ARM_Date& fwdDate,
											  ARM_Date& maturity, 
											  ARM_Date& payDate,
											  int compMeth = 0, 
											  int dayCount = KACTUAL_365,
											  int DomOrFrgRate = 1, 
											  int discYC = 1,
											  int YieldDecomp = K_COMP_PROP,
											  double Margin = 0.0,
											  StoreFwdRateInfo* StoreInfo = NULL);

		virtual double ExpectedFwdYieldNonAdj(double fwdJulDate,
											  double maturity,
											  double payJulDate,
											  int compMeth = 0,
											  int dayCount = KACTUAL_365,
											  int DomOrFrgRate = 1, 
											  int discYC = 1,
											  int YieldDecomp = K_COMP_PROP,
											  double Margin = 0.0,
											  StoreFwdRateInfo* StoreInfo = NULL);
       
        // Compute Caplet price, only use for replication at this moment
        double EuroCaplet(double Settlement, 
                          double ResetMaturity,
                          double StartMaturity, 
                          double EndMaturity,
                          double PayMaturity,
                          int OptionType, 
                          double Fwd, 
                          double Strike, 
                          int DomOrFrg,
                          int CompMeth, 
                          int DayCount,
                          bool IsCMS = false,
                          bool isTreasury = false,
                          int FixFrequency = 1,
                          double UnderlyingTenor = 0.0,
                          int YieldDecomp = K_COMP_PROP,
                          double Margin = 0.0,
                          ARM_Currency* ccy = NULL,
                          StoreFwdRateAndCapletInfo* StoreInfo = NULL);

        // Compute BS parameters (forward, strike and volatility adjustment) to put in the
        // BS formula to compute CMS caplet price
        void ComputeCapletBSParameters(double ResetDate,
                                       double StartDate,
                                       double EndDate,
                                       double PayMaturity,
                                       int OptionType,
                                       double Fwd, 
                                       double OptionStrike,
                                       int DomOrFrg,
                                       int CompMeth,
                                       int DayCount,
                                       bool IsCMS,
                                       bool isTreasury,
                                       int FixFrequency,
                                       double UnderlyingTenor,
                                       int YieldDecomp,
                                       double Margin,
                                       ARM_Currency* ccy,
                                       // Outputs
                                       double& FwdAdj,
                                       double& FwdAdjNoDecap,
                                       double& StrikeAdj,
                                       double& VolAdj,
                                       StoreFwdRateAndCapletInfo* StoreInfo = NULL);

		double SpreadOptionPrice(double fwd1, double fwd2, 
								 double vol1, 
								 double vol2, double Correl, 
								 double strike, 
								 double optMat, int optType,
                                 double w1 = 1, double w2 = 1,
				                 int ComputedFormula = 1/* 0 formule d'Olivier, 1 Formule Classique*/);

		double GenericSpreadOptionPrice(double fwd1, double fwd2, 
                                        marginal_descr* m1, 
                                        marginal_descr* m2, 
                                        jointdistrib_descr* jd, 
                                        double strike, 
                                        double optMat, int optType,
                                        double w1, double w2,
									    int ComputedFormula = 1,/* 0: O.C Formulae, 1: Classic Formulae */
									    int NumSteps = 100);

        // Compute Swaption price, only use for replication at this moment
        virtual double EuroSwaption(double ResetSwap, 
                                    double StartSwap,
                                    double Endswap, 
                                    int OptionType,
                                    double SwapFwd, 
                                    double Strike,
                                    double Annuity, 
                                    int DomOrFrg = 1,
                                    ARM_Vector* fixCfPrices = NULL,
                                    ARM_Vector* fixCfTerms  = NULL,
                                    double dt = 0);

        virtual double bsOption(double spot,
                                double strike,
                                double volatility,
                                double dividend,
                                double discountRate,
                                double maturity,
                                double CallPut,
                                double volTenor = 0)
        {
            double bsOpt = ::bsOption(spot, strike, volatility,
                                       dividend, discountRate,
                                       maturity, CallPut);
            return(bsOpt);
        } 

        virtual double bsDelta(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut)
        {
            double delta = ::bsDelta(spot, strike, volatility, dividend, discountRate,
                                     maturity, CallPut);
            return(delta);
        }

        virtual double bsGamma(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut)
        {
            double gamma = ::bsGamma(spot, strike, volatility, dividend, discountRate, 
                                     maturity, CallPut);
            return(gamma);
        }

        virtual double bsVega(double spot,
                              double strike,
                              double volatility,
                              double dividend,
                              double discountRate,
                              double maturity,
                              double CallPut)
        {
            double vega = ::bsVega(spot, strike, volatility, dividend, 
                                   discountRate, maturity, CallPut);
            return(vega);
        }

        virtual double bsTheta(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut)
        {
            double theta = ::bsTheta(spot, strike, volatility, dividend, discountRate, 
                                     maturity, CallPut);
            return(theta);
        }

        virtual double bsKappa(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut)
        {
            double kappa = ::bsKappa(spot, strike, volatility, dividend, 
                                     discountRate, maturity, CallPut);
            return(kappa);
        }

        virtual double bsRho(double spot,
                             double strike,
                             double volatility,
                             double dividend,
                             double discountRate,
                             double maturity,
                             double CallPut)
        {
            double rho = ::bsRho(spot,strike, volatility, dividend, 
                                 discountRate, maturity, CallPut);
            return(rho);
        }

        bool HasCorrelationStructure()
        {
            if (itsCorrelationMatrix)
			   return(true);
			else
			   return(false);
        }

};



/*=======================================*
    ARM_BSNorModel Class
 *=======================================*/
class ARM_BSNorModel : public ARM_BSModel
{
    public:

        ARM_BSNorModel(void) : ARM_BSModel() 
        {}

	    ARM_BSNorModel(ARM_Date& AsOfDate, 
                       ARM_ZeroCurve* YielCurve,
                       ARM_VolCurve* SpreadLockCurve,          // Spread vol  
                       ARM_VolCurve* ConvAdjVolatility = NULL, // IRG
                       ARM_VolCurve* Volatility = NULL,        // Pricing vol : IRG: 
                                                               // Or Swopt : case CMS
                       ARM_CorrelManager* CorrelManager = NULL,
                       ARM_ConvAdjustManager* ConvAdjustManager = NULL,
                       ARM_ZeroCurve* DiscountCurve = NULL)
              : ARM_BSModel(AsOfDate, 
                            YielCurve,
                            SpreadLockCurve,
                            ConvAdjVolatility,
                            Volatility,
                            CorrelManager,
                            ConvAdjustManager,
                            DiscountCurve) 
        {}


	    ARM_BSNorModel(ARM_Date& AsOfDate,
                       ARM_ZeroCurve* YieldCurve,
                       ARM_VolCurve* SpreadLockCurve,
                       ARM_VolCurve* convAdjVolatility,
                       ARM_VolCurve* priceVolatility,                         
                       ARM_CorrelManager* CorrelManager,
                       ARM_ConvAdjustManager* ConvAdjustManager,
                       ARM_ZeroCurve* DiscountCurve,
				       ARM_VolCurve* Correl,  // Correlations
				       ARM_VolCurve* CashVol, // Cash Vol (particulièrement pour CMT)
				       ARM_VolCurve* SpreadVolCurve,
				       int ModelType = K_2LOG,		  // 2LOG, LOG, NOR
				       int SpreadVolType= K_COMPUTED, // K_INPUTED,K_COMPUTED
				       ARM_BSSmiledModel* sabrMod = NULL)
          : ARM_BSModel(AsOfDate,
                        YieldCurve,
                        SpreadLockCurve,
                        convAdjVolatility,
                        priceVolatility,           
                        CorrelManager,
                        ConvAdjustManager,
                        DiscountCurve,
                        Correl,
                        CashVol,
                        SpreadVolCurve,
                        ModelType,
                        SpreadVolType,
                        sabrMod) 
        {}


        ARM_BSNorModel(ARM_Date& startDate,
                       double spot, 
                       ARM_ZeroCurve* dividend,
                       ARM_ZeroCurve* discountRate, 
                       ARM_VolCurve* volatility,
                       int volType = K_PRICE)
                      : ARM_BSModel(startDate,
                                    spot, 
                                    dividend,
                                    discountRate, 
                                    volatility,
                                    volType) 
        {}


        ARM_BSNorModel(ARM_Date& startDate, 
                       double spot, 
                       double flatdiv,
                       double flatdiscrate,
                       double flatvol,
                       int volType = K_PRICE)
                   : ARM_BSModel(startDate, 
                                 spot, 
                                 flatdiv,
                                 flatdiscrate,
                                 flatvol,
                                 volType) 
        {}

        ARM_BSNorModel(ARM_ZeroCurve* zerocurve,
                       ARM_VolCurve* volatility,
                       int volType = K_INDEX_RATE)
                      : ARM_BSModel(zerocurve,
                                    volatility,
                                    volType) 
        {}
        
        virtual ~ARM_BSNorModel(void)
        {}


        ARM_BSNorModel(const ARM_BSNorModel& rhs) : ARM_BSModel(rhs) 
        {};
        
        ARM_BSNorModel& operator = (const ARM_BSNorModel& rhs)
        {
	        if( this != & rhs )
	        {
		        ARM_BSModel::operator=(rhs);
	        }
	        return *this;
        }

		virtual bool IsNormal(void)
		{
            return(true);
		}

        virtual ARM_Object* Clone(void) 
        { 
             return(new ARM_BSNorModel(*this)); 
        }


        double ComputeVol(double matu,
                          double tenor,
                          double fwd,
                          double strike,
						  int mode = -1)
        {
            double res = ARM_BSModel::ComputeVol(matu,
                                                 tenor,
                                                 fwd,
                                                 strike,
						                         mode);
            
            res *= 100.0; // Because this is the dimension used generally in the pricing!

            return(res);
        }

        virtual double bsOption(double spot,
                                double strike,
                                double volatility,
                                double dividend,
                                double discountRate,
                                double maturity,
                                double CallPut,
                                double volTenor = 0);

        virtual double bsDelta(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut);

        virtual double bsGamma(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut);

        virtual double bsVega(double spot,
                              double strike,
                              double volatility,
                              double dividend,
                              double discountRate,
                              double maturity,
                              double CallPut);

        virtual double bsTheta(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut) 
        { 
            return(0.0); 
        }

        virtual double bsKappa(double spot,
                               double strike,
                               double volatility,
                               double dividend,
                               double discountRate,
                               double maturity,
                               double CallPut) 
        { 
            return(0.0); 
        }

        virtual double bsRho(double spot,
                             double strike,
                             double volatility,
                             double dividend,
                             double discountRate,
                             double maturity,
                             double CallPut) 
        { 
            return(0.0); 
        }
};



                          /********************************************/
                          /*                                          */
                          /*   The IRG/SWOPT BS like pricing model    */
                          /*                                          */
                          /********************************************/



class ARM_BSPricingModel : public ARM_BSModel
{
    private:

        ARM_BSModel*   itsIRGModel;
        ARM_BSModel*   itsSWOPTModel;

        
        ARM_ZeroCurve* itsForcastCurve;
        ARM_ZeroCurve* itsDiscountCurve;

        ARM_VolCurve*  itsSpreadVol;
        ARM_VolCurve*  itsAdjConvVol;


        int            itsModelType;

        ARM_Vector*    itsCalibInfos;
        ARM_Vector*    itsNumInfos;


        // shared objects: inherited from ARM_Model Class
        /* 
           ARM_CorrelManager* itsCorrelMgr;

           ARM_ConvAdjustManager* itsConvAdjustManager;
         */

        // This is a GENERATED BS model which in fact will be used
        // effectively in the pricing process

        ARM_BSModel* itsRealBSPricingModel; 


        void Init(void)
        {
            SetName(ARM_BSPRICINGMODEL);

            itsIRGModel   = NULL;
            itsSWOPTModel = NULL;

            itsForcastCurve  = NULL;
            itsDiscountCurve = NULL;

            itsSpreadVol     = NULL;
            itsAdjConvVol    = NULL;

            itsModelType     = K_2LOG;

            itsCalibInfos    = NULL;
            itsNumInfos      = NULL;

            SetBSCategory(1);

            // The generated BS model

            itsRealBSPricingModel = NULL;
        }
 
    public:

        ARM_BSPricingModel(void);  

        ARM_BSPricingModel(ARM_BSModel* IRGModel, ARM_BSModel* SWOPTModel = NULL);

        ARM_BSPricingModel(ARM_Date AsOfDate,
                           ARM_ZeroCurve* ForcastCurve,
                           ARM_BSModel*   IRGModel      = NULL, 
                           ARM_BSModel*   SWOPTModel    = NULL,
                           ARM_ZeroCurve* DiscountCurve = NULL,
                           ARM_VolCurve*  SpreadVol     = NULL,
                           ARM_VolCurve*  AdjConvVol    = NULL,
                           int            ModelType     = K_2LOG,
                           ARM_Vector*    CalibInfos    = NULL,
                           ARM_Vector*    NumInfos      = NULL,
                           ARM_CorrelManager* CorrelMgr = NULL,
                           ARM_ConvAdjustManager* ConvAdjustManager = NULL);

        ARM_BSPricingModel(const ARM_BSPricingModel& bs);    
        ARM_BSPricingModel& operator = (const ARM_BSPricingModel& bs); 

        virtual ~ARM_BSPricingModel(void);  
        
        void BitwiseCopy(const ARM_Object* srcBSModel);   
        
        void Copy(const ARM_Object* srcBSModel);
        
        virtual ARM_Object* Clone(void);

        ARM_BSModel* GetRealBSPricingModel(ARM_Security* Security)
        {
            return(itsRealBSPricingModel);
        }

        ARM_BSModel* GetGeneratedRealBSPricingModel(ARM_Security* Security);
        
        void View(char* id = NULL, FILE* ficOut = NULL);
};













#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/