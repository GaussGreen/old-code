#ifndef _INGPMODELS_HWSVBGMQto_MODEL
#define _INGPMODELS_HWSVBGMQto_MODEL

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpmodels/MultiAssets.h"

#include "gpinfra/modelnamemap.h"
#include "gpinfra/curvemodelparam.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_HullWhite1F;
class ARM_BGMSV1F;
class ARM_LN_Fx;

class ARM_HWSVBGMQtoModel : public ARM_MultiAssetsModel 
{
private:
	ARM_BGMSV1F *			itsSVMMModel;
	ARM_HullWhite1F *		itsHWModel;
	ARM_LN_Fx *				itsFXModel;
	ARM_VectorVector		itsEigenValues;
	std::vector<double>			itsResetTimes;
	int						itsNbEffectiveReset;
	ARM_CurveModelParam		itsFakeParam;

public:

	ARM_HWSVBGMQtoModel(const ARM_ModelNameMap& modelNameMap, const ARM_CurveMatrix& correlCurveMatrix);

	ARM_HWSVBGMQtoModel(const ARM_HWSVBGMQtoModel& rhs);

	ASSIGN_OPERATOR(ARM_HWSVBGMQtoModel)

	virtual ~ARM_HWSVBGMQtoModel();

	virtual ARM_Object*	Clone() const	{ return new ARM_HWSVBGMQtoModel(*this);};

	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void	NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances) const;
	virtual void	NumMethodStateLocalVariancesAndStdDev(const std::vector<double>& timeSteps);
	virtual void	NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances) const;
	virtual void	ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;			
	virtual void	ModelStateLocalVariancesAndStdDev(const std::vector<double>& timeSteps);

	
	virtual ARM_BoolVector NeedMCIntegProcess() const;

	virtual string	toString(const string& indent, const string& nextIndent) const;

private:
	void	computeConvexAdj(ARM_PricingStatesPtr& states, double fromTime, double toTime, int iFirst, int iLast) const;

	void	HWMCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const;
};	

CC_END_NAMESPACE()

#endif