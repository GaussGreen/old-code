#include "DKMaille.h"
#include "DKMaille2D.h"
#include "Tree.h"


extern FILE* fic;

// FIXMEFRED: mig.vc8 (25/05/2007 11:12:06):missing return types
class Swaption 
{
    public:
        
        Swaption();
	    
        virtual ~Swaption();
        
        void Init(   double dSpotDate,
                DKMaille<double> dNoticeDates,
                double dNumTimeLinesBeforeFirstNotice,
                double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                double dNumTimeLinesPerYear,
                double dOptimal,
                double LastDate,
                double dMeanReversion,
                double smileParameter,
                DKMaille<double> dFloatResetDates,
                DKMaille<double> dFixedPaymentDates,
                DKMaille<double> dFloatPaymentDates,
                DKMaille<double> dIndexStartDates,
                DKMaille<double> dIndexEndDates,
                DKMaille<double> dFixedAccrualBasis,
                DKMaille<double> dFloatAccrualBasis,
                DKMaille<double> dFixedCoupon,
                DKMaille<double> dBasisSpread);
        
        void ModifyBoosterData(  DKMaille<double> dFloatResetDates,
                            DKMaille<double> dFixedPaymentDates,
                            DKMaille<double> dFloatPaymentDates,
                            DKMaille<double> dIndexStartDates,
                            DKMaille<double> dIndexEndDates,
                            DKMaille<double> dFixedAccrualBasis,
                            DKMaille<double> dFloatAccrualBasis,
                            DKMaille<double> dFixedCoupon,
                            DKMaille<double> dBasisSpread);
        
        void SetVol( DKMaille<double> dSigmaDates,
                DKMaille<double> dSigma);
        
        void CreateTree();

        void CleanTree();
        
        void SetMeanReversionForTree();
        
        void CalibrateZCinTree();
        
        DKMaille<double>& Price(int PriceOnWhatNotice);

        DKMaille<double>& StraddlePrice(int PriceOnWhatNotice);
		
        DKMaille<double>& PriceInTree(int PriceOnWhatNotice);
		
        DKMaille<double>& PriceInTreeWithAnalytics(int PriceOnWhatNotice);
        
        int GetFirstResetDateAfterNotice(int uiNotice);

        void BootstrapVolForQModelInTree(DKMaille<double>& Sigma,
                                    DKMaille<double> InputVol,
                                    double dNoticePeriodInDays);

        void Print(FILE* file);
        
        void SetK();
		
        void Setq(double SmileParameter);

        double SpotDate;
        double MeanReversion;
		DKMaille<double> MeanReversionInTree;
        int NbPasTotal;
        int NbPasBeforeLastNotice;
        DKMaille<double> T;
        DKMaille<double> dT;
        double NbNotices;
        double NbCoupons;
        double CallPut;
        bool PriceZCinTree;
        bool QModel;
        DKMaille<double> q;
        DKMaille<double> K;
        DKMaille<double> NoticeDates;
        DKMaille<double> FloatResetDates;
        DKMaille<double> FixedPaymentDates;
        DKMaille<double> FloatPaymentDates;
        DKMaille<double> IndexStartDates;
        DKMaille<double> IndexEndDates;
        DKMaille<double> FixedAccrualBasis;
        DKMaille<double> FloatAccrualBasis;
        DKMaille<double> FixedCoupon;
        DKMaille<double> BasisSpread;
        DKMaille2D<double> A_FixedPayment;
        DKMaille2D<double> A_FloatPayment;
        DKMaille2D<double> A_IndexStart;
        DKMaille2D<double> A_IndexEnd;
        DKMaille2D<double> B_FixedPayment;
        DKMaille2D<double> B_FloatPayment;
        DKMaille2D<double> B_IndexStart;
        DKMaille2D<double> B_IndexEnd;
        

        CcyMarketData ccyMarketData;
        Tree tree;
};




DKMaille<double>& SwaptionPrice(double dSpotDate,
                                DKMaille<double> &dDates,
                                DKMaille<double> &dRates,
                                DKMaille<double> &dRatesNoBasis,
                                DKMaille<double> &dStdDevX,
                                DKMaille<double> &dStdDevY,
                                DKMaille2D<double> &dStdDevZ,
                                double dMeanReversion,
                                double smileParameter1,
                                double smileParameter2,
                                double dIsSwaptionCalibrationWithBasis,
                                DKMaille2D<double> &dBoosterData,
                                double dNumTimeLinesBeforeFirstNotice,
                                double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                double dNumTimeLinesPerYear,
                                double dOptimal);



void GetBootstrappedVolForQModelInTree( DKMaille<double>& SigmaDates,
                                        DKMaille<double>& Sigma,
                                        DKMaille<double> dNoticeDates,
                                        DKMaille<double> ZC_Dates,
                                        DKMaille<double> ZC_Rates,
                                        DKMaille<double> Basis_ZC_Rates,
                                        DKMaille<double> dStdDevX,
                                        DKMaille<double> dStdDevY,
                                        DKMaille2D<double> dStdDevZ,
                                        double dSpotDate,
                                        double dFinalMaturity,
                                        double dNumTimeLinesBeforeFirstNotice,
                                        double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                        double dNumTimeLinesPerYear,
                                        double dOptimal,
                                        double dMeanReversion,
                                        double dSmileParameter);

DKMaille<double> Bootstrapping_DK3F_Numerical_SP(  DKMaille2D<double>& dStdDevZ,  
                                                    DKMaille<double>& dStdDevX,
                                                    DKMaille<double>& dStdDevY,  
                                                    DKMaille<double>& ZC_Dates,
                                                    DKMaille<double>& ZC_Rates,
                                                    DKMaille<double>& Basis_ZC_Rates,
                                                    DKMaille<double>& dNoticeDates,
                                                    DKMaille<double>& dSwapStartDate,
                                                    DKMaille<double>& dSwapEndDate,
                                                    DKMaille<double>& dModelParameter,                                                            
                                                    double dJulianObservationDate,
                                                    double dNumTimeLinesBeforeFirstNotice,
                                                    double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                                    double dNumTimeLinesPerYear,
                                                    double dOptimal,                                                            
                                                    double dSmileParameter);
