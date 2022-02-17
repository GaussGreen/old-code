
#ifndef _INGPMODELS_BISVMMMODEL_H
#define _INGPMODELS_BISVMMMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpinfra/modelnamemap.h"

#include "MultiAssets.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsBGMSV1F;
class ARM_BGMSV1F;

class ARM_BiSVMM : public ARM_MultiAssetsModel
{
protected:
	
	ARM_BGMSV1F *			itsModel1;
	ARM_BGMSV1F *			itsModel2;
	std::vector<double>			itsResetTimes;
	int						itsNbEffectiveReset;
	ARM_VectorVector		itsEigenValues;

public:
	ARM_BiSVMM(const ARM_ModelNameMap& modelNameMap, const ARM_CurveMatrix& correlationMatrix);
	ARM_BiSVMM(const ARM_BiSVMM& rhs);

	ASSIGN_OPERATOR(ARM_BiSVMM)

	virtual ~ARM_BiSVMM();

	virtual ARM_Object*	Clone() const	{ return new ARM_BiSVMM(*this);};
	
	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void	NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances) const;
	virtual void	ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances) const;
	virtual void	ModelStateLocalVariancesAndStdDev(const std::vector<double>& timeSteps);
	virtual void	NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances) const;

	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;

	virtual string	toString(const string& indent, const string& nextIndent) const;
};

CC_END_NAMESPACE()

#endif
