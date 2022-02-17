/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file 1IRFXModel.h
 *
 *  \brief
 *
 *  \brief 1 interest + fx multi assets model
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date Jully 2006
 */



#ifndef _INGPMODELS_1IRFXMODEL_H
#define _INGPMODELS_1IRFXMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssetsMeanReverting.h"
#include "gpmodels/typedef.h"

#include "gpinfra/modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsQ1F;
class ARM_ModelParamsSmiled;
class ARM_2IRFXModel;

////////////////////////////////////////////////////////////////////////////////
//
// This class defines a multi asset with 2 models:
// _ a 1 factor Hull & White model for the domestic model associated to a 
// basis model
// _ a p factors Smiled FX model for the N forward FX rates
//
////////////////////////////////////////////////////////////////////////////////


class ARM_1IRFXModel :	public ARM_MultiAssetsMeanReverting
{
private:
	mutable ARM_GP_Vector			itsTimeSteps;
	mutable ARM_VectorVector		itsEigenValues;
	mutable ARM_MatrixVector		itsEigenVectors;
	mutable ARM_GP_Matrix			itsFwdFxDrifts;
	mutable std::vector<double>			itsHWRelDrifts;
	mutable ARM_GP_Matrix					itsIntAbsDrifts;
	// Correlation calculated based on correlation matrix
	ARM_2IRFXModelPtr				its2IRFXModel;
	ARM_MatrixVector				its2IRFXCorrelMatrix;
	mutable std::vector<double>			itsACPErrors;

	/// to initialise the sub models
	void Validate();

    bool IsBasisRefModel() const { return GetRefModel() == &*((*(GetModelMap()))[DomBasisModel]->Model()); }

    void ComputeIntegratedVCV(
			double step,
			double nextStep,
			const std::vector<double>& resetTimes,
			size_t fxModelNb,
            const ARM_ModelParamsQ1F* const domModelParams,
            const ARM_ModelParamsSmiled* const fxModelParams,
            const ARM_GP_Matrix& correlMatrix,
            ARM_GP_Matrix& variances) const;

	void ComputeDrifts(
			const std::vector<double>& timeSteps,
			const std::vector<double>& resetTimes,
			const std::vector<double>& settlementTimes,
			const ARM_ModelParamsQ1F* const domModelParams,
			const ARM_ModelParamsSmiled* const fxModelParams,
			const ARM_MatrixVector& localVariances  );


	void ComputeCorrelMatrix(
			const std::vector<double>& resetTimes,
			const std::vector<double>& settlementTimes,
			ARM_MatrixVector& correlMatrix, 
			bool withDomesticIR);
	
public:
    static const size_t NbLocalModels;

    enum modelsAlias
    {
        DomModel	=0,         /// the domestic stochastic IR model
        FxModel,			    /// the stochastic FX model
		DomBasisModel,		    /// pure basis model
        NbModels
    };

	ARM_1IRFXModel(const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlationMatrix,
		ARM_2IRFXModelPtr Model2IRFX  = ARM_2IRFXModelPtr(NULL));

	ARM_1IRFXModel(const ARM_1IRFXModel& rhs);
	ASSIGN_OPERATOR(ARM_1IRFXModel)
	virtual ~ARM_1IRFXModel();

	void CalibrateCorrelation( ARM_2IRFXModelPtr Model2IRFX );

    const ARM_PricingModelPtr& GetModel(modelsAlias modelIdx) const { return (* GetModelMap())[modelIdx]->Model(); }

    /// No Cholesky decomposition needed at brownian level
    virtual bool NeedsToCholeskyDecomposeFactors() const { return false; }

	virtual size_t FactorCount() const;
	virtual ARM_BoolVector NeedMCIntegProcess() const;

    virtual bool SupportAnalyticMarginal() const { return false; }

	ARM_2IRFXModelPtr Get2IRFXModel() { return its2IRFXModel; }

    virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
	void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr);
	void SetNumeraire(const ARM_NumerairePtr& numerairePtr);

    virtual void	NumMethodStateLocalVariances(	const std::vector<double>& timeSteps,
												ARM_MatrixVector& localVariances ) const;		

	virtual void	NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
												ARM_MatrixVector& globalVariances ) const;

    virtual void ModelStateLocalVariances( 
		const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;
	void ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps );
    void NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps ) {};
	void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;
	
	void From2to1IRFXModelStates(ARM_PricingStatesPtr& states, int timeIndex);

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;
    virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_1IRFXModel(*this);}
	virtual string ExportShortName() const { return "L1IFX";}

	virtual string toString(const string& indent="",const string& nextIndent="") const;

protected:
	/// Specialised init method for tree because of 1D tree calibrations
	virtual ARM_PricingStatesPtr BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
