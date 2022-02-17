#ifndef _INGPMODELS_HWSBGMQto_MODEL
#define _INGPMODELS_HWSBGMQto_MODEL

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpmodels/MultiAssets.h"

#include "gpinfra/modelnamemap.h"
#include "gpinfra/curvemodelparam.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_HullWhite1F;
class ARM_SmiledFRM;
class ARM_LN_Fx;

class ARM_HWSBGMQtoModel : public ARM_MultiAssetsModel 
{
private:
	ARM_SmiledFRM *			itsSBGMModel;
	ARM_HullWhite1F *		itsHWModel;
	ARM_LN_Fx *				itsFXModel;
	ARM_CurveModelParam		itsFakeParam;
	int						itsSBGMoffsetIdx;
	ARM_VectorVector		itsEigenValues;

public:

	ARM_HWSBGMQtoModel(const ARM_ModelNameMap& modelNameMap, const ARM_CurveMatrix& correlCurveMatrix);

	ARM_HWSBGMQtoModel(const ARM_HWSBGMQtoModel& rhs);

	ASSIGN_OPERATOR(ARM_HWSBGMQtoModel)

	virtual ~ARM_HWSBGMQtoModel();

	virtual ARM_Object*	Clone() const	{ return new ARM_HWSBGMQtoModel(*this);};

	virtual void	MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void	NumMethodStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const;
	virtual void	NumMethodStateLocalVariancesAndStdDev(const ARM_GP_Vector& timeSteps);
	virtual void	NumMethodStateGlobalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& globalVariances) const;
	virtual void	ModelStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances ) const;			
	virtual void	ModelStateLocalVariancesAndStdDev(const ARM_GP_Vector& timeSteps);

	
	virtual ARM_BoolVector NeedMCIntegProcess() const;

	virtual string	toString(const string& indent, const string& nextIndent) const;

private:
	ARM_GP_Vector	GetQuantoConvexAdj(double t, double T, int iFirst, int iLast) const;

	double			SBGMIntegratedCovariance(double s, double t, int i, int j) const;
	double			SBGMIntegratedCovariance(double t, int i, int j) const;
};	

CC_END_NAMESPACE()

#endif