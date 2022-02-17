/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file HWHWQtoModel.h
 *
 *  \brief
 *
 *  \brief simple Quanto Model
 *
 *	\author  J. Messines
 *	\version 1.0
 *	\date February 2007
 */



#ifndef _INGPMODELS_HWHWQTOMODEL_H
#define _INGPMODELS_HWHWQTOMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssetsMeanReverting.h"

#include "gpinfra/modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingContext;
class ARM_ModelParamsHW1FStd;
class ARM_CurveModelParam;

////////////////////////////////////////////////////////////////////////////////
//
// This class defines a multi asset with 2 models et demi :
// _ a 1 factor Hull & White model for each currency
// _ a fake factor for FX
//
////////////////////////////////////////////////////////////////////////////////


class ARM_HWHWQtoModel :	public ARM_MultiAssetsMeanReverting
{
private:
	bool itsForLocalModel;

    ARM_GP_MatrixPtr itsDriftCorrections;

    void AddIntegratedLocalCorrections( const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// to initialise the sub models
	void Validate(bool fxFlag);

	typedef void (ARM_HWHWQtoModel::*InitCalibrationFunc)(
		double evalTime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const;

    void ComputeIntegratedFwdFxVCV(double step, double nextStep,
            const ARM_ModelParamsHW1FStd* const domModelParams,
            const ARM_ModelParamsHW1FStd* const forModelParams,
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
        ForModel,			    /// the foreign stochastic IR model
        FxModel,			    /// the stochastic FX model
		DomBasisModel,		    /// pure basis model
		ForBasisModel,
		ForLocalModel,
		NbModels
    };

	ARM_HWHWQtoModel(
		const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlCurveMatrix,
		bool fxFlag);

	ARM_HWHWQtoModel(const ARM_HWHWQtoModel& rhs);
	ASSIGN_OPERATOR(ARM_HWHWQtoModel)
	virtual ~ARM_HWHWQtoModel(){}

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_HWHWQtoModel(*this);}
	virtual string ExportShortName() const { return "LHWQT";}
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

	void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& variances ) const;

	void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

    virtual ARM_BoolVector NeedMCIntegProcess() const;

    virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	/// No Cholesky decomposition needed at brownian level
    virtual bool NeedsToCholeskyDecomposeFactors() const { return false; }
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
