//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVQmcPortfolioLoss.hpp
//
//   Description : QMC interface of Portfolio Loss State variable
//
//----------------------------------------------------------------------------

#ifndef SVQmcPortfolioLoss_HPP
#define SVQmcPortfolioLoss_HPP

#include "edginc/IQMCStateVariableBase.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/PortfolioName.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"

//#include "edginc/SVGenPortfolioLoss.hpp"


DRLIB_BEGIN_NAMESPACE

/************************************************************************/
/* Portfolio loss derived SV:                                           */
/************************************************************************/
class SVQmcPortfolioLoss : public PortfolioName::SVGen::SV //, public IQMCStateVariableBase
{
public:
	virtual ~SVQmcPortfolioLoss() {}

    SVQmcPortfolioLoss(
	    const string& name,
	    double notional,
	    const CreditAssetConstSP& asset,
	    double beta,
	    const CDoubleConstSP& recoveryOverride,
	    const DateTime& protectionStartDate,
	    const DateTime& maturityCutOff,
	    const CtgLegLossPerDefaultArrayConstSP& histContLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeRecovered,
	    const DateTimeLiteVectorConstSP& timeline
        ) :
    SV(name, notional, asset, beta, recoveryOverride, protectionStartDate,
        maturityCutOff, histContLosses, feeLosses, feeRecovered, timeline)
    {}

protected:
	SVQmcPortfolioLoss();
};

DECLARE(SVQmcPortfolioLoss);

class SVQmcPortfolioLossFullMC : public SVQmcPortfolioLoss
{
public:
	SVQmcPortfolioLossFullMC(
		const vector<SVDateOfDefaultSP>& _ddSVs, // specific for this object (Full MC of QMC)
	    const string& name,
	    double notional,
	    const CreditAssetConstSP& asset,
	    double beta,
	    const CDoubleConstSP& recoveryOverride,
	    const DateTime& protectionStartDate,
	    const DateTime& maturityCutOff,
	    const CtgLegLossPerDefaultArrayConstSP& histContLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeRecovered,
	    const DateTimeLiteVectorConstSP& timeline
        ) :
    SVQmcPortfolioLoss(name, notional, asset, beta, recoveryOverride, protectionStartDate,
        maturityCutOff, histContLosses, feeLosses, feeRecovered, timeline),
        ddSVs(_ddSVs),
        mySubResults(new PortfolioNameSVSubResultVector(_ddSVs.size()))
    {    }

private:

   	virtual void calculateLossEvents() const;
    virtual void setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const;

    PortfolioNameSVSubResultVectorSP mySubResults;
    SVQmcPortfolioLossFullMC();
    vector<SVDateOfDefaultSP> ddSVs;
};


// ********
// Note, that Fast MC is incompatible with the SVGen::SV, but only compatible with IndexedSVGen::SV
// ********

class SVQmcPortfolioLossFastMC : public SVQmcPortfolioLoss
{
public:
	SVQmcPortfolioLossFastMC(
		const vector<SVSurvivalDiscFactorSP>& _sdfSVs,// specific for this object (Full MC of QMC)
	    const string& name,
	    double notional,
	    const CreditAssetConstSP& asset,
	    double beta,
	    const CDoubleConstSP& recoveryOverride,
	    const DateTime& protectionStartDate,
	    const DateTime& maturityCutOff,
	    const CtgLegLossPerDefaultArrayConstSP& histContLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeRecovered,
	    const DateTimeLiteVectorConstSP& timeline
        ) :
    SVQmcPortfolioLoss(name, notional, asset, beta, recoveryOverride, protectionStartDate,
        maturityCutOff, histContLosses, feeLosses, feeRecovered, timeline), sdfSVs(_sdfSVs)
    {}

	//virtual void prepare(bool mm);
	//virtual DoubleArrayConstSP getPortfolioLoss();
private:

   	virtual void calculateLossEvents() const;
    virtual void setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const;

	SVQmcPortfolioLossFastMC();
	vector<SVSurvivalDiscFactorSP> sdfSVs;
};


/***    Wonderful world of indexed SV... ***/


class SVQmcIndexedPortfolioLoss : public PortfolioName::IndexedSVGen::SV //, public IQMCStateVariableBase
{
public:
	virtual ~SVQmcIndexedPortfolioLoss() {}

    SVQmcIndexedPortfolioLoss(
	    const string& name,
	    double notional,
	    const CreditAssetConstSP& asset,
	    double beta,
	    bool hasDefaulted,
	    double lossGivenDefault,
	    const DateTimeConstSP& defaultSettDate,
	    const CDoubleConstSP& recoveryOverride,
	    const DateTime& protectionStartDate,
	    const DateTime& maturityCutOff,
	    const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeRecovered
        ) :
/*
	Anatoly: here, instead of protectionStartDate, pass protectionStartDateIndex
	instead of maturityCutOff, pass maturityCutOffIndex - on the product timeline
	here, I have used "1" just to ensure compilation
	also, the first argument needs to be productTimeline
*/
    SV(DateTimeArrayConstSP(),name, notional, asset, beta, hasDefaulted, lossGivenDefault, defaultSettDate, recoveryOverride,
        1, 1, feeLosses, feeRecovered)
    {}

protected:
	SVQmcIndexedPortfolioLoss();
};

DECLARE(SVQmcIndexedPortfolioLoss);

class SVQmcIndexedPortfolioLossFullMC : public SVQmcIndexedPortfolioLoss
{
public:
	SVQmcIndexedPortfolioLossFullMC(
		const vector<SVDateOfDefaultSP>& _ddSVs,// specific for this object (Full MC of QMC)
	    const string& name,
	    double notional,
	    const CreditAssetConstSP& asset,
	    double beta,
	    bool hasDefaulted,
	    double lossGivenDefault,
	    const DateTimeConstSP& defaultSettDate,
	    const CDoubleConstSP& recoveryOverride,
	    const DateTime& protectionStartDate,
	    const DateTime& maturityCutOff,
	    const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeRecovered
        ) :
    SVQmcIndexedPortfolioLoss(name, notional, asset, beta, hasDefaulted,
        lossGivenDefault, defaultSettDate, recoveryOverride,
        protectionStartDate, maturityCutOff, feeLosses, feeRecovered),
    ddSVs(_ddSVs)
    {}

private:

   	virtual void calculateLossEvents() const;
    virtual void setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const;

	SVQmcIndexedPortfolioLossFullMC();

    PortfolioNameSVSubResultVectorSP mySubResults;
    vector<SVDateOfDefaultSP> ddSVs;
};


class SVQmcIndexedPortfolioLossFastMC : public SVQmcIndexedPortfolioLoss
{
public:
	SVQmcIndexedPortfolioLossFastMC(
		const vector<SVSurvivalDiscFactorSP>& _sdfSVs,// specific for this object (Full MC of QMC)
	    const string& name,
	    double notional,
	    const CreditAssetConstSP& asset,
	    double beta,
	    bool hasDefaulted,
	    double lossGivenDefault,
	    const DateTimeConstSP& defaultSettDate,
	    const CDoubleConstSP& recoveryOverride,
	    const DateTime& protectionStartDate,
	    const DateTime& maturityCutOff,
	    const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	    const FeeLegReductionPerDefaultArrayConstSP& feeRecovered
        ) :
    SVQmcIndexedPortfolioLoss(name, notional, asset, beta, hasDefaulted,
        lossGivenDefault, defaultSettDate, recoveryOverride,
        protectionStartDate, maturityCutOff, feeLosses, feeRecovered),
        sdfSVs(_sdfSVs)
    {}

private:

   	virtual void calculateLossEvents() const;
    virtual void setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const;

	SVQmcIndexedPortfolioLossFastMC();
	vector<SVSurvivalDiscFactorSP> sdfSVs;
};


DRLIB_END_NAMESPACE
#endif // EDR_IQMCSTATEVARIABLE_HPP

