/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file NP1IRNFXModel.h
 *
 *  \brief
 *
 *  \brief N+1 interest rate + N fx multi assets model
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date Jully 2006
 */

#ifndef _INGPMODELS_NP1IRNFXMODEL_H
#define _INGPMODELS_NP1IRNFXMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssetsMeanReverting.h"

#include "gpinfra/modelnamemap.h"

#include "gpcalib/densityfunctors.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_NP1IRNFXModel :	public ARM_MultiAssetsMeanReverting
{
private:
	/// to initialise the sub models
	void Validate();

	int itsNbCcy;
public:
	ARM_NP1IRNFXModel(
		const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlationMatrix);

	ARM_NP1IRNFXModel(const ARM_NP1IRNFXModel& rhs);
	ASSIGN_OPERATOR(ARM_NP1IRNFXModel)
	virtual ~ARM_NP1IRNFXModel();

    const ARM_PricingModelPtr& GetModel(int modelIdx) const { return (* GetModelMap())[modelIdx]->Model(); }

    /// No Cholesky decomposition needed at brownian level
    virtual bool NeedsToCholeskyDecomposeFactors() const { return false; }

	virtual size_t FactorCount() const;
	virtual ARM_BoolVector NeedMCIntegProcess() const;

	ARM_PricingModelPtr DomModel() const { return (*(GetModelMap()))[0]->Model(); }
	ARM_PricingModelPtr DomBasisModel() const { return (*(GetModelMap()))[2*itsNbCcy-1]->Model(); }

	ARM_PricingModelPtr ForModel(int i) const;
	ARM_PricingModelPtr FxModel(int i) const;

    virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

    virtual void	NumMethodStateLocalVariances(	
						const std::vector<double>& timeSteps,
						ARM_MatrixVector& localVariances ) const;		

	virtual void	NumMethodStateGlobalVariances( 
						const std::vector<double>& timeSteps,
						ARM_MatrixVector& globalVariances ) const;

    virtual void ModelStateLocalVariances( 
		const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	void ModelStateLocalVariancesAndStdDev( 
		const std::vector<double>& timeSteps );

    void NumMethodStateLocalVariancesAndStdDev( 
		const std::vector<double>& timeSteps ) {};

	void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;
    virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_NP1IRNFXModel(*this);}
	virtual string ExportShortName() const { return "LNIFX";}

	void CalibrateFunctional (
		const std::vector<double>& ResetTimes,
		vector<ARM_DensityFunctor*> densities,
		int nbRows,
		int nbCols,
		int sizeGrid = 501,
		double nbStdDev = 6.,
		bool rescaling = true);
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

