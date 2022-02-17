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
 *	\date June 2007
 */

#ifndef _INGPMODELS_2IRNFXSV_H
#define _INGPMODELS_2IRNFXSV_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssetsMeanReverting.h"

#include "gpinfra/modelnamemap.h"

#include "gpcalib/densityfunctors.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_2IRFXSV :	public ARM_MultiAssetsMeanReverting
{
private:
	/// to initialise the sub models
	void Validate();

public:
	enum modelsAlias
    {
        DomModel	=0,         /// the domestic stochastic IR model
        ForModel,			    /// the foreign stochastic IR model
        FxModel,			    /// the stochastic FX model
		DomBasisModel,		    /// pure basis model
		ForBasisModel,		    /// pure basis model
        NbModels
    };

	enum hestonModelAlias
	{
		SpotHeston	=3,
		VarHeston	=4
	};

	ARM_2IRFXSV(
		const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlationMatrix);

	ARM_2IRFXSV(const ARM_2IRFXSV& rhs);
	ASSIGN_OPERATOR(ARM_2IRFXSV)
	virtual ~ARM_2IRFXSV();

    const ARM_PricingModelPtr& GetModel(int modelIdx) const { return (* GetModelMap())[modelIdx]->Model(); }

    /// No Cholesky decomposition needed at brownian level
    virtual bool NeedsToCholeskyDecomposeFactors() const { return false; }

	virtual size_t FactorCount() const;
	virtual ARM_BoolVector NeedMCIntegProcess() const;
	bool SupportAnalyticMarginal() const { return false; };

	ARM_PricingModelPtr GetDomModel() const { return (*(GetModelMap()))[DomModel]->Model(); }
	ARM_PricingModelPtr GetForModel() const { return (*(GetModelMap()))[ForModel]->Model(); }
	ARM_PricingModelPtr GetFxModel() const { return (*(GetModelMap()))[FxModel]->Model(); }
	ARM_PricingModelPtr GetDomBasisModel() const { return (*(GetModelMap()))[DomBasisModel]->Model(); }
	ARM_PricingModelPtr GetForBasisModel() const { return (*(GetModelMap()))[ForBasisModel]->Model(); }

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
	virtual ARM_Object* Clone() const { return new ARM_2IRFXSV(*this);}
	virtual string ExportShortName() const { return "LFXSV";}

	void CalibrateFunctional (
		const std::vector<double>& ResetTimes,
		vector<ARM_DensityFunctor*> densities,
		int sizeGrid = 501,
		double nbStdDev = 6.);
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

