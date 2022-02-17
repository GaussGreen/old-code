/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file MultiAssetsReverting.h
 *
 *  \brief
 *
 *  \brief multi assets mean reverting model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */



#ifndef _INGPMODELS_MULTIASSETSMEANREVERTING_H
#define _INGPMODELS_MULTIASSETSMEANREVERTING_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssets.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_MultiAssetsMeanReverting :	public ARM_MultiAssetsModel
{
public:
	ARM_MultiAssetsMeanReverting(
		const ARM_ModelNameMap*	modelNameMap	= NULL,
		const ARM_GP_Matrix* correlationMatrix	= NULL );
	
	ARM_MultiAssetsMeanReverting(
		const ARM_ModelNameMap*	modelNameMap,
		const ARM_CurveMatrix* correlationMatrix);
	
	ARM_MultiAssetsMeanReverting(const ARM_MultiAssetsMeanReverting& rhs);
	ASSIGN_OPERATOR(ARM_MultiAssetsMeanReverting)	
	virtual ~ARM_MultiAssetsMeanReverting();	

	/// function to tells whether all the base models are mean reverting compatible
	static bool MeanRevertingCompatible(const ARM_ModelNameMap& modelNameMap);

	/// function to cache in stdDev (for performance reason!)
	virtual void ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);

	/// function to compute the integrated local covariances for Ornstein Uhlenbeck models!
    virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	/// function to compute the integrated global covariances for Ornstein Uhlenbeck models!
	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(FactorCount(), true); };

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_MultiAssetsMeanReverting(*this);}
	virtual string ExportShortName() const { return "LMAMR";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
