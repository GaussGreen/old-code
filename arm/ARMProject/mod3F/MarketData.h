#include "DKMaille.h"
#include "DKMaille2D.h"

// FIXMEFRED: mig.vc8 (25/05/2007 11:09:50): missing return types
class CcyMarketData
{
    public :
        
        CcyMarketData();
        
        virtual ~CcyMarketData();
        
        void Init(	double SpotDate, 
				DKMaille<double> dDates, 
				DKMaille<double> dRates, 
				DKMaille<double> dDatesNoBasis, 
				DKMaille<double> dRatesNoBasis);

        double TauxZC(double T);
        
        double ZC(double T);

        double BasisTauxZC(double T);
        
        double BasisZC(double T);

        void GetConstantSigma();
        
        void SetForwardRates(DKMaille<double> T);
        
        double GetInstantaneousRate(double T);
        
        double GetBforBondPrice(double MeanReversion, 
                                double t, 
                                double T,
                                DKMaille<double> q,
                                DKMaille<double> K);
        
        void Print(FILE* file);
        

        DKMaille<double> ZC_Dates;
        DKMaille<double> ZC_YearFraction;
        DKMaille<double> ZC_Rates;
            
        DKMaille<double> Basis_ZC_Dates;
        DKMaille<double> Basis_ZC_YearFraction;
        DKMaille<double> Basis_ZC_Rates;
        
        DKMaille<double> ForwardRates;
        DKMaille<double> ForwardDates;

        DKMaille<double> SigmaDates;
        DKMaille<double> Sigma;
		DKMaille<double> ConstSigma;
};


class MarketData
{
    public :
        MarketData();
        virtual ~MarketData();
        
        // Attributs

        double AsOf;
        
        // IR
        CcyMarketData DomCcy;
        CcyMarketData ForCcy;
    
        // FX
        DKMaille<double> FX_BS_Vol_Dates;
        DKMaille<double> FX_BS_Vol;
        DKMaille<double> FX_Spot_Vol_Dates;
        DKMaille<double> FX_Spot_Vol;
               
};



