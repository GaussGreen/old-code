
#ifndef _INGPMODELS_HWxSVMMSpread_MODEL
#define _INGPMODELS_HWxSVMMSpread_MODEL

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpmodels/MultiAssets.h"
#include "gpinfra/curvemodelparam.h"

#include "gpinfra/modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_SVMMSpread;
class ARM_HullWhite1F;
class ARM_HullWhite2F;

class ARM_HWxSVMMSpread : public ARM_MultiAssetsModel
{
private:
	ARM_SVMMSpread *	itsSpreadModel;
	ARM_HullWhite1F *	itsHWModel;

	std::vector<double>		itsResetTimes;
	int					itsNbEffectiveReset;
	ARM_VectorVector	itsEigenValues;

	ARM_VectorVector	itsxCorrel;
	std::vector<double>		itsCorrEndTimes;
	
	mutable bool		itsConvexAdjDone;
	double				itsLastTHW;

	std::vector<double>		itsSpreadGVol;
	ARM_CurveModelParam	itsFakeParam;

public:

	ARM_HWxSVMMSpread(const ARM_ModelNameMap& modelNameMap, ARM_HullWhite2F * correlModel, const std::vector<double>& corrEndTimes, double ConstantCrossCorrel);

	ARM_HWxSVMMSpread(const ARM_HWxSVMMSpread& rhs);

	ASSIGN_OPERATOR(ARM_HWxSVMMSpread)

	virtual ~ARM_HWxSVMMSpread();

	virtual ARM_Object*	Clone() const	{ return new ARM_HWxSVMMSpread(*this);};

	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void	NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances) const;
	virtual void	ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances) const;
	virtual void	ModelStateLocalVariancesAndStdDev(const std::vector<double>& timeSteps);
	virtual void	NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances) const;
	virtual void	NumMethodStateLocalVariancesAndStdDev(const std::vector<double>& timeSteps);
	
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual ARM_BoolVector NeedMCIntegProcess() const;

	virtual string	toString(const string& indent, const string& nextIndent) const;

protected:

	void	computeCrossCorrel(ARM_HullWhite2F * correlModel, double ConstantCrossCorrel);

	void	computeConvexAdj() const;

	void	computeConvexAdj(ARM_PricingStatesPtr& states, double fromTime, double toTime, int iFirst, int iLast) const;

	void	HWMCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const;
};

CC_END_NAMESPACE()

#endif