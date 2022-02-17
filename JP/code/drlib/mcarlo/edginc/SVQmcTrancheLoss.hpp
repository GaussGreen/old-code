//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVQmcTrancheLoss.hpp
//
//   Description : QMC interface of Tranche Loss State variable
//
//----------------------------------------------------------------------------

#ifndef SVQmcTrancheLoss_HPP
#define SVQmcTrancheLoss_HPP

#include "edginc/IQMCStateVariableBase.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SVGenTrancheLoss.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"


DRLIB_BEGIN_NAMESPACE

/************************************************************************/
/* Tranche loss derived SV:                                             */
/************************************************************************/
class SVQmcTrancheLoss : public IQMCStateVariableBase, public SVTrancheLoss
{
public:
	virtual ~SVQmcTrancheLoss() {}
	SVQmcTrancheLoss(
		const DoubleArraySP _notionals,
		const DateTimeArraySP _portfolioLossDates,
		const DoubleArraySP _lowerPct,
		const DoubleArraySP _upperPct);

	virtual double element(int dtIndex) const
	{
		throw ModelException("SVQmcPortfolioLoss::element",
			"incorrect way to access this information, as it is inefficient. Use getPortfolioLoss() instead");
		return 0.0;
	}
	virtual void prepare(bool mm);

protected:
	SVQmcTrancheLoss();

	DoubleArraySP notionals;
	DateTimeArraySP portfolioLossDates;
	DoubleArray lowerLoss;
	DoubleArray rangeOfLoss;

	double        totalNotional;
	CDoubleMatrixSP trancheLoss;

	DoubleArray   portfolioLoss; // internal, used inside calculations
	CreditTrancheLossConfigArraySP trancheLossConfig;

};

DECLARE(SVQmcTrancheLoss);

class SVQmcTrancheLossFullMC : public SVQmcTrancheLoss
{
public:
	SVQmcTrancheLossFullMC(
		const vector<SVDateOfDefaultSP>& _ddSVs,
		const DoubleArraySP _notionals,
		const DateTimeArraySP _portfolioLossDates,
		const DoubleArraySP _lowerPct,
		const DoubleArraySP _upperPct);

	virtual CDoubleMatrixConstSP getTrancheLosses();
private:
	SVQmcTrancheLossFullMC();
	vector<SVDateOfDefaultSP> ddSVs;
};

class SVQmcTrancheLossFastMC : public SVQmcTrancheLoss
{
public:
	SVQmcTrancheLossFastMC(
		const vector<SVSurvivalDiscFactorSP>& _sdfSVs,
		const DoubleArraySP _notionals,
		const DateTimeArraySP _portfolioLossDates,
		const DoubleArraySP _lowerPct,
		const DoubleArraySP _upperPct,
		const IConvolutorSP _convolutor);

	virtual void prepare(bool mm);
	virtual CDoubleMatrixConstSP getTrancheLosses();
private:
	SVQmcTrancheLossFastMC();
	vector<SVSurvivalDiscFactorSP> sdfSVs;
	IConvolutorSP convolutor;
};


DRLIB_END_NAMESPACE
#endif // SVQmcTrancheLoss_HPP

