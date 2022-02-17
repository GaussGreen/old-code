 
#include "Magnet/Magnet.h"
#include "General/General.h"

#include "Convolution.h"
#include "Calibrator.h"
#include "Quadrature.h"
#include "Tranche.h"
#include "CorrSkew.h"

#include "imsl.h"
#include "imsls.h"

#include <vector>

#undef  _MAGNET_WITH_GENERIC
#define _MAGNET_WITHOUT_GENERIC

using namespace CM;


MAGNET_PREFIX("ModelLoss_");
MAGNET_CATEGORY("ModelLoss");

MAGNET_DESCRIPTION("Code Interface for the Random Beta- Random Recovery Model")
MAGNET_RESULT("Loss Distributions, Tranche Priches and Base Correlation Skew")


// Calculates the entire Correlation Skew for a homogenous portfolio 
MAGNET_X_FUNCTION7(
	   Array<double>,  CorrelationSkewSimple,
	   int					NoNames,				MAGNET_MANDATORY,
	   int					NotperName,				MAGNET_MANDATORY,
	   double				Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   double				DefaultProb,			MAGNET_MANDATORY,
	   double				recovery,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{	
	double Loss=(1-recovery)*NotperName;
	return CorrelationSkewSimple(NoNames,NotperName,MktLossTranche,Beta, Attachment, Loss, DefaultProb);  
}


// Calculates the entire Correlation Skew for a  portfolio with homogenous beta
// and heterogenous losses and Default Probability 
MAGNET_X_FUNCTION6(
	   Array<double>,  CorrelationSkew,
	   double				Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{
	return CorrelationSkew(Notional,MktLossTranche,Beta, Attachment, Loss, DefaultProb);
}


// Calculates the entire Correlation Skew for a  portfolio with heterogenous beta
// losses and Default Probability 
MAGNET_X_FUNCTION6(
    Array<double>,  CorrelationSkewGeneral,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{
	return CorrelationSkewGeneral(Notional,MktLossTranche,Beta, Attachment, Loss, DefaultProb);
}


// Calculates the entire Correlation Skew for a  portfolio with homogenous beta
// and heterogenous losses and Default Probability 
MAGNET_X_FUNCTION6(
	   Array<double>,  CorrelationSkew2,
	   double				Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Matrix<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{
	return CorrelationSkew2(Notional,MktLossTranche,Beta, Attachment, Loss, DefaultProb);
}


// Calculates the entire Correlation Skew for a  portfolio with heterogenous beta
// losses and Default Probability 
MAGNET_X_FUNCTION6(
    Array<double>,  CorrelationSkewGeneral2,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Matrix<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{
	return CorrelationSkewGeneral2(Notional,MktLossTranche,Beta, Attachment, Loss, DefaultProb);
}




// Calculates the entire Correlation Skew for a  portfolio with heterogenous beta
// losses and Default Probability using the Secant False Position Method
MAGNET_X_FUNCTION6(
    Array<double>,  CorrelationSkewGeneralSec,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{
	return CorrelationSkewGeneralSecant(Notional,MktLossTranche,Beta, Attachment, Loss, DefaultProb);  
}


// Calculates the Base Correlation for a single Attachment Point for a  portfolio with heterogenous beta
// losses and Default Probability using the Secant False Position Method
MAGNET_X_FUNCTION6(
    Array<double>,  CorrelationSkewGeneralSecSingleK,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   double				Attachment,				MAGNET_MANDATORY)        // Wrapper type and name
{
	Array<double> Attach(3);
	Attach[0]=0; Attach[1]=Attachment; Attach[2]=1;

	return CorrelationSkewGeneralSecant(Notional,MktLossTranche,Beta, Attach, Loss, DefaultProb);  
}



// Calculates the Base Correlation for a single Attachment Point for a  portfolio with heterogenous beta
// losses and Default Probability using a simple search method
MAGNET_X_FUNCTION6(
    Array<double>,  CorrelationSkewGeneralSingleK,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&MktLossTranche,		MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   double				Attachment,				MAGNET_MANDATORY)        // Wrapper type and name
{
	Array<double> AttachmentVect(3);
	AttachmentVect[0]=0;
	AttachmentVect[1]=Attachment;
	AttachmentVect[2]=1;
	
	return CorrelationSkewGeneral(Notional,MktLossTranche,Beta, AttachmentVect, Loss, DefaultProb);  
}


// Calc expected tranche loss for heterogenous portfolio using creditMetric model
MAGNET_X_FUNCTION6(
    Array<double>,  TrancheLossCreditmetrics,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name
{
	Array<int> LossParam=CalcLossUnit(n,Loss);
    Array<double> lossDist = LossDist(n, Beta, Loss, DefaultProb);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,LossParam[0],Attachment,lossDist);

	return TrancheLoss;  
}


// Calc expected loss for a single tranche given a base correlation
MAGNET_X_FUNCTION7(
    Array<double>,  BaseCorrelationLossSingleK_pdf,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   double				BaseCorrelation,		MAGNET_MANDATORY,			
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   double				singleK,				MAGNET_MANDATORY)        // Wrapper type and name
{
	Array<double> Attachment(3);
	Attachment[0]=0; Attachment[1]=singleK;	Attachment[2]=1;

	Array<int> LossParam=CalcLossUnit(n,Loss);

    // calc average beta
	int m=Beta.size();
	double betaavg=0;
	for (int i=0; i<m; i++)
	{
		betaavg=betaavg+Beta[i];
	}
	betaavg=betaavg/n;

    // scale beta to match base correlation, average of new betas is bc
	double tweak=(BaseCorrelation-betaavg)/(1-betaavg);
	BetaTweak(Beta,tweak);

    // calc expected tranche loss using tweaked beta
    Array<double> lossDist = LossDist(n, Beta, Loss, DefaultProb);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,LossParam[0],Attachment,lossDist);

	return TrancheLoss;  
}



// Calc tranche loss
MAGNET_X_FUNCTION10(
    Array<double>,  TrancheLossPolynomialExp,
	   int					n,					    MAGNET_MANDATORY,
	   double				a0,						MAGNET_MANDATORY,
	   double				a1,						MAGNET_MANDATORY,
	   double				a2,						MAGNET_MANDATORY,
	   double				aexp,					MAGNET_MANDATORY,
	   double				kexp,					MAGNET_MANDATORY, 
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY)        // Wrapper type and name

{
	Array<int> LossParam=CalcLossUnit(n,Loss);
	Array<double> lossDist = LossDistRandomFactorPolExp(n,a0,a1,a2,aexp,kexp,Loss,DefaultProb);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,LossParam[0],Attachment,lossDist);

	return TrancheLoss;  
}


// Calc tranche loss
MAGNET_X_FUNCTION8(
    Array<double>,  TrancheLossNFactDiscrete,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&theta,					MAGNET_MANDATORY,
	   Array<double>		&probMarket,			MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				recovery,				MAGNET_MANDATORY)        // Wrapper type and name

{
	Array<int> LossParam=CalcLossUnit(n,Loss);
    Array<double>  lossDist = LossDistRandomDiscrete(n, probMarket,theta, Loss, DefaultProb);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,LossParam[0],Attachment,lossDist);

    return TrancheLoss;  
}

// Calibrate two betas and one theta to target spread
MAGNET_X_FUNCTION13(
    Array<double>, Calibrate2Factor,
        CalibrationMode     mode,                   MAGNET_MANDATORY,
        int                 n,                      MAGNET_MANDATORY,
        Array<double>&      seeds,                  MAGNET_MANDATORY, // beta1, beta2, theta
        Array<int>&         masks,                  MAGNET_MANDATORY, // which seed to tweak
        Array<double>&      targetSpreads,          MAGNET_MANDATORY, 
        Array<int>&         loss,                   MAGNET_MANDATORY,
        Array<double>&      defaultProb,            MAGNET_MANDATORY,
        int                 notional,               MAGNET_MANDATORY,
        Array<double>&      attachment,             MAGNET_MANDATORY,
        double              noIntervals,            MAGNET_MANDATORY,
        double              lbound,                 MAGNET_MANDATORY,
        double              ubound,                 MAGNET_MANDATORY,
        double              mbound,                 MAGNET_MANDATORY)
{
    return Calibrate2Fact(mode, n, seeds, masks, targetSpreads,
                          loss, defaultProb, notional, attachment,
                          noIntervals, lbound, ubound, mbound);
}


// Calc tranche loss
MAGNET_X_FUNCTION11(
    Array<double>,  TrancheLossNFact,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&Theta,					MAGNET_MANDATORY,
	   Array<int>			&Loss,					MAGNET_MANDATORY,
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   int					Notional,				MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				NoIntervals,			MAGNET_MANDATORY,
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY)
{
	Array<int> LossParam=CalcLossUnit(n,Loss);
    Array<double> lossDist = LossDistRandomNFactorsSimp(n,Beta,Theta, Loss, DefaultProb, NoIntervals, lbound, ubound, mbound);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,LossParam[0],Attachment,lossDist);

	return TrancheLoss;  
}


// Calc tranche loss
MAGNET_X_FUNCTION12(
    Array<double>,  TrancheLossNFactJumpHet,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&Theta,					MAGNET_MANDATORY,  
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY,
	   Array<int>			&Notional,				MAGNET_MANDATORY,
	   Array<double>		&Recovery,				MAGNET_MANDATORY,  
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				NoLossUnits,			MAGNET_MANDATORY,	
	   double				NoIntervals,			MAGNET_MANDATORY,
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY)        // Wrapper type and name

{
	Matrix<double>  lossDist = LossDistRandomNFactorsSimpsBuck(n, 0,0,0,0,0,Beta,Theta, Notional,Recovery,DefaultProb, NoIntervals, lbound, ubound, mbound,NoLossUnits);

	double TotalNotional=0; double TotalLosses=0;
	for (int l=0; l<n; l++)
	{ 
        TotalNotional=TotalNotional+Notional[l];
        TotalLosses=TotalNotional+(1-Recovery[l])*Notional[l];  
	}

	double LossUnit=TotalLosses/NoLossUnits;

	int MaxNoUnits=ceil(NoLossUnits)+1;
	
	Array<double> TrancheLoss = ExpectedLossTrancheBuck(TotalNotional,LossUnit,MaxNoUnits,Attachment,lossDist);

	return TrancheLoss;  
}


// Calc tranche loss
MAGNET_X_FUNCTION11(
    Array<double>,  TrancheLossNFactHomog,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&Beta,					MAGNET_MANDATORY,
	   Array<double>		&Theta,					MAGNET_MANDATORY,
	   int					indNotional,			MAGNET_MANDATORY,
	   double				recovery,				MAGNET_MANDATORY,
	   double				DefaultProb,			MAGNET_MANDATORY,
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				NoIntervals,			MAGNET_MANDATORY,
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY)        // Wrapper type and name

{
	double Notional=indNotional*n;
 
	Array<double>  lossDist = LossDistRandomNFactorsHomogenSimp(n, Beta, Theta, DefaultProb, NoIntervals, lbound, ubound, mbound);

	double lossUnit = (1-recovery)*indNotional;

	Array<double> TrancheLoss = ExpectedLossTranche(Notional,lossUnit,Attachment,lossDist);
	
	return TrancheLoss;  
}


// Calibrate 3 factor, 2 theta, and 1 r1 model 
MAGNET_X_FUNCTION14(
    Array<double>,  Calibrate3FactorHomogRR,
       CalibrationMode      mode,                   MAGNET_MANDATORY,
	   int					n,					    MAGNET_MANDATORY,
       Array<double>        &seeds,                 MAGNET_MANDATORY, // beta1, 2, 3, theta1, 2, r1
       Array<int>           &masks,                 MAGNET_MANDATORY,
       Array<double>        &targetSpreads,         MAGNET_MANDATORY,
       double               r2,                     MAGNET_MANDATORY,
	   Array<double>		&Recovery,				MAGNET_MANDATORY, // 3 point recovery average recovery is middle value
	   int					indNotional,			MAGNET_MANDATORY, // notional per name
	   double				DefaultProb,			MAGNET_MANDATORY, // individual name default prob
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				NoIntervals,			MAGNET_MANDATORY, // numerical integration params
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY)
{
    return Calibrate3FactHomogRR(mode, n, seeds, masks, targetSpreads, r2,
                                 Recovery, indNotional, DefaultProb, Attachment,
                                 NoIntervals, lbound, ubound, mbound);
}



// Calc tranche loss for combined N random factor and 3 random recovery regime
MAGNET_X_FUNCTION12(
    Array<double>,  TrancheLossNFactHomogRR,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&Beta,					MAGNET_MANDATORY, // N betas
	   Array<double>		&Theta,					MAGNET_MANDATORY, // N-1 thetas
	   Array<double>		&Recovery,				MAGNET_MANDATORY, // 3 point recovery average recovery is middle value
	   Array<double>		&ThresholdR,			MAGNET_MANDATORY, // 2 r boundaries
	   int					indNotional,			MAGNET_MANDATORY, // notional per name
	   double				DefaultProb,			MAGNET_MANDATORY, // individual name default prob
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				NoIntervals,			MAGNET_MANDATORY, // numerical integration params
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY)
{
	double loss0=100*(1-Recovery[0]);
	double loss1=100*(1-Recovery[1]);
	double loss2=100*(1-Recovery[2]);

	int minloss=GCD(int(loss0),GCD(int(loss1),int(loss2)));

	int rat0=loss0/minloss;
	int rat1=loss1/minloss;
	int rat2=loss2/minloss;

    double b1 = Beta[0];
    double b2 = Beta[1];
    double b3 = Beta[2];
	
	double Notional=indNotional*n;
	double lossUnit=(1-Recovery[0])*indNotional/(rat0);
	
	Array<double> lossDist = 
        LossDistRandomNFactorsHomogenSimpRR(n, Beta, Theta, 0, DefaultProb, NoIntervals, 
                                            Recovery, ThresholdR, lbound, ubound, mbound);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,lossUnit,Attachment,lossDist);

	return TrancheLoss;  
}

class LossDistribution : public CM::Object{
public:
    int lossUnit;
    Array<double> lossDist;
};

MAGNET_X_FUNCTION14(
    SharedPointer<LossDistribution>, LossDistNFactRR,
       int					n,					    MAGNET_MANDATORY,
	   Array<double>		&BetaCo,    			MAGNET_MANDATORY, // N Beta coefficients
	   Array<double>		&Theta,					MAGNET_MANDATORY, // N-1 thetas
       Array<double>		&RecoveryCo,			MAGNET_MANDATORY, // 3 point recovery coefficients
	   Array<double>		&ThresholdR,			MAGNET_MANDATORY, // 2 r boundaries
       Array<double>        &Beta,                  MAGNET_MANDATORY, // individual name betas
       Array<double>        &ExpectedRecovery,      MAGNET_MANDATORY, // individual expected recovery
	   Array<int>   		&Loss,					MAGNET_MANDATORY, // individual expected loss
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY, // individual default probability)
       double				NoIntervals,			MAGNET_MANDATORY,
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY,
       bool                 useLinearBeta,          = false)
{
    SharedPointer<LossDistribution> rslt = new LossDistribution;
    rslt->lossDist = LossDistRandomNFactorsSimpRR(n,BetaCo,Theta, RecoveryCo, ThresholdR, 
                                                          Beta, ExpectedRecovery, Loss, DefaultProb, 
                                                          NoIntervals, lbound, ubound, mbound, rslt->lossUnit,
                                                          useLinearBeta);
    return rslt;
}

MAGNET_X_FUNCTION2(
    Matrix<double>, GetLossDist,
        SharedPointer<LossDistribution> PortLossDist,   MAGNET_MANDATORY,
        double                          threshold,      MAGNET_MANDATORY)
{
    std::vector<double> loss;
    std::vector<double> dist;

    Array<double>& data = PortLossDist->lossDist;
    for (int i = 0; i < data.size(); ++i) {
        if (data[i] >= threshold) {
            loss.push_back(i);
            dist.push_back(data[i]);
        }
    }

    Matrix<double> rslt(dist.size() + 1, 2);
    rslt[0][0] = PortLossDist->lossUnit;
    rslt[0][1] = dist.size();
    for (i = 1; i <= rslt[0][1]; ++i) {
        rslt[i][0] = loss[i-1];
        rslt[i][1] = dist[i-1];
    }
    return rslt;
}

MAGNET_X_FUNCTION3(
    Array<double>, TrancheExpectedLossNFactRR,
        int					Notional,				MAGNET_MANDATORY, // total notional
        Array<double>		&Attachment,			MAGNET_MANDATORY,
        SharedPointer<LossDistribution> PortLossDist,   MAGNET_MANDATORY)
{
  	Array<double> TrancheLoss = ExpectedLossTranche(Notional,PortLossDist->lossUnit,
                                                    Attachment,PortLossDist->lossDist);
    return TrancheLoss;
}

MAGNET_X_FUNCTION3(
    Array<double>, TrancheletExpectedLossNFactRR,
        int					Notional,				MAGNET_MANDATORY, // total notional
        Matrix<double>		&Attachment,			MAGNET_MANDATORY,
        SharedPointer<LossDistribution> PortLossDist,   MAGNET_MANDATORY)
{
  	Array<double> TrancheLoss = ExpectedLossTranchelet(Notional,PortLossDist->lossUnit,
                                                    Attachment,PortLossDist->lossDist);
    return TrancheLoss;
}


// Calc tranche loss using N+3 regime for heterogenous portfolio
MAGNET_X_FUNCTION16(
    Array<double>,  TrancheLossNFactRR,
	   int					n,					    MAGNET_MANDATORY,
	   Array<double>		&BetaCo,    			MAGNET_MANDATORY, // N Beta coefficients
	   Array<double>		&Theta,					MAGNET_MANDATORY, // N-1 thetas
       Array<double>		&RecoveryCo,			MAGNET_MANDATORY, // 3 point recovery coefficients
	   Array<double>		&ThresholdR,			MAGNET_MANDATORY, // 2 r boundaries
       Array<double>        &Beta,                  MAGNET_MANDATORY, // individual name betas
       Array<double>        &ExpectedRecovery,      MAGNET_MANDATORY, // individual expected recovery
	   Array<int>   		&Loss,					MAGNET_MANDATORY, // individual expected loss
	   Array<double>		&DefaultProb,			MAGNET_MANDATORY, // individual default probability
	   int					Notional,				MAGNET_MANDATORY, // total notional
	   Array<double>		&Attachment,			MAGNET_MANDATORY,
	   double				NoIntervals,			MAGNET_MANDATORY,
	   double				lbound,					MAGNET_MANDATORY,							
	   double				ubound,					MAGNET_MANDATORY,
	   double				mbound,					MAGNET_MANDATORY,
       bool                 useLinearBeta,          = false)
{
    int lossunit;
    Array<double> lossDist = LossDistRandomNFactorsSimpRR(n,BetaCo,Theta, RecoveryCo, ThresholdR, 
                                                          Beta, ExpectedRecovery, Loss, DefaultProb, 
                                                          NoIntervals, lbound, ubound, mbound, lossunit,
                                                          useLinearBeta);
	Array<double> TrancheLoss = ExpectedLossTranche(Notional,lossunit,Attachment,lossDist);
    return TrancheLoss;
}



