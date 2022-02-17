/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file HWHW2FQtoModel.h
 *
 *  \brief
 *
 *  \brief a little less simpler Quanto Model
 *
 *	\author  L. du Bois
 *	\version 1.0
 *	\date May 2007
 */



#ifndef _INGPMODELS_HWHW2FQTOMODEL_H
#define _INGPMODELS_HWHW2FQTOMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssetsMeanReverting.h"

#include "gpinfra/modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingContext;
class ARM_ModelParamsHW1FStd;
class ARM_ModelParamsHW2F;
class ARM_CurveModelParam;


//! Model designed to price callable instrument on quanto spread.
/*!	This class defines a multi-asset with 2 models and half:
	\l a 1 factor Hull & White model for domestic currency
	\l a 2-factor Hull & White model for foreign currency
	\l a fake factor for FX
 */
class ARM_HWHW2FQtoModel :	public ARM_MultiAssetsMeanReverting
{
private:
	bool itsForLocalModel;

    ARM_GP_MatrixPtr itsDriftCorrections;

    void AddIntegratedLocalCorrections( const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// to initialise the sub models
	void Validate(bool fxFlag);

	typedef void (ARM_HWHW2FQtoModel::*InitCalibrationFunc)(
		double evalTime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

    void ComputeIntegratedFwdFxVCV(double step, double nextStep,
            const ARM_ModelParamsHW1FStd* const domModelParams,
            const ARM_ModelParamsHW2F* const forModelParams,
            const ARM_ModelParamsHW1FStd* const fxModelParams,
            const ARM_CurveModelParam& fxVol,
            const ARM_GP_Matrix& correlMatrix,
            ARM_GP_Matrix& variances) const;

protected:
	/// Specialised init method for tree because of 1D tree calibrations
	virtual ARM_PricingStatesPtr BackwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);

	/// Temporary specialised init method for MC waiting for actual use of sampler
	virtual ARM_PricingStatesPtr ForwardLookingInit(ARM_NumMethodPtr& numMethod, double firstInductTime);
	
public:
    static const size_t NbLocalModels;

    enum modelsAlias
    {
        DomModel	=0,         /// the domestic stochastic IR model
        ForModel,			    /// the foreign stochastic IR model first factor
        FxModel,			    /// the stochastic FX model
		DomBasisModel,		    /// pure basis model
		ForBasisModel,
		ForLocalModel,
		NbModels
    };

	ARM_HWHW2FQtoModel(
		const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlCurveMatrix,
		bool fxFlag);

	ARM_HWHW2FQtoModel(const ARM_HWHW2FQtoModel& rhs);
	ASSIGN_OPERATOR(ARM_HWHW2FQtoModel)
	virtual ~ARM_HWHW2FQtoModel(){}

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_HWHW2FQtoModel(*this);}
	virtual string ExportShortName() const { return "LHW2Q";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;


    const ARM_PricingModelPtr& GetModel(modelsAlias modelIdx) const { return (* GetModelMap())[modelIdx]->Model(); }


	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const;


	virtual bool SupportAnalyticMarginal() const;

	virtual void IntegratedLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	void NumMethodStateLocalGlobalVariances( const std::vector<double>& timeSteps,
	    ARM_MatrixVector& localVariances,
	    ARM_MatrixVector& variances ) const;

    virtual ARM_BoolVector NeedMCIntegProcess() const;

    virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	/// No Cholesky decomposition needed at brownian level
    virtual bool NeedsToCholeskyDecomposeFactors() const { throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "coucou"); return false; }
    virtual bool NeedArrowDebreuPrices() const { return false; }

	/// No more used because variance initialisation is done through the numerical method sampler
    void ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps) {}
    void NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps ) {}

    virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;


    virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

    virtual ARM_VectorPtr DiscountFactor( 
		const string& modelName,
        double evalTime, 
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

    virtual ARM_VectorPtr Forward(
	    const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
