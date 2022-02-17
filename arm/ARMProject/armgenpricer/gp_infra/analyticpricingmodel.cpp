/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file AnalyticPricingModel.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date April 2004
 */

#include "gpinfra/analyticpricingmodel.h"
#include "gpinfra/pricingstates.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_AnalyticPricingModel
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_AnalyticPricingModel::ARM_AnalyticPricingModel(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params)
:	ARM_PricingModel(zc,params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticPricingModel
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_AnalyticPricingModel::ARM_AnalyticPricingModel(const ARM_AnalyticPricingModel& rhs)
:	ARM_PricingModel(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticPricingModel
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_AnalyticPricingModel::~ARM_AnalyticPricingModel()
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticPricingModel
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_AnalyticPricingModel& ARM_AnalyticPricingModel::operator=(const ARM_AnalyticPricingModel& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModel::operator=(rhs);
        // Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticPricingModel
///	Routine: all no relevant virtual fcts
///	Returns: error
///	Action : unimplemented method that throws an error
////////////////////////////////////////////////////
void ARM_AnalyticPricingModel::Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "Re_InitialiseCalibParams : unimplemented function !"); }

void ARM_AnalyticPricingModel::PreProcessing(ARM_ModelFitter& modelFitter)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "PreProcessing : unimplemented function !"); }

void ARM_AnalyticPricingModel::PostProcessing(const ARM_ModelFitter& modelFitter)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "PostProcessing : unimplemented function !"); }

void ARM_AnalyticPricingModel::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "AdviseCurrentCalibSecIndex : unimplemented function !"); }

void ARM_AnalyticPricingModel::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "AdviseCurrentCalib : unimplemented function !"); }

void ARM_AnalyticPricingModel::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio,
								ARM_ModelParam* modelParam, 
								size_t factorNb)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "AdviseBreakPointTimes : unimplemented function !"); }

ARM_PricingStatesPtr ARM_AnalyticPricingModel::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "Init : unimplemented function !"); return ARM_PricingStatesPtr(NULL); }


ARM_PricingStatesPtr ARM_AnalyticPricingModel::FirstPricingStates( size_t bucketSize ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "FirstPricingStates : unimplemented function !"); return ARM_PricingStatesPtr(NULL); }

std::vector<double>* ARM_AnalyticPricingModel::ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos )
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "ComputeModelTimes : unimplemented function !"); 
	return nullptr; 
}

bool ARM_AnalyticPricingModel::SupportBackwardInduction() const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "SupportBackwardInduction : unimplemented function !"); return false; }

bool ARM_AnalyticPricingModel::SupportForwardInduction() const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "SupportForwardInduction : unimplemented function !"); return false; }

bool ARM_AnalyticPricingModel::SupportAnalyticMarginal() const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "SupportAnalyticMarginal : unimplemented function !"); return false; }

void ARM_AnalyticPricingModel::IntegratedLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "IntegratedLocalDrifts : unimplemented function !"); }

void ARM_AnalyticPricingModel::ModelStateLocalVariances( const std::vector<double>& timeSteps,ARM_MatrixVector& localVariances ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "LocalVariances : unimplemented function !"); }

void ARM_AnalyticPricingModel::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,ARM_MatrixVector& localVariances ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "LocalVariances : unimplemented function !"); }


bool ARM_AnalyticPricingModel::NeedsToCholeskyDecomposeFactors( ) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "NeedsToCholeskyDecomposeFactors : unimplemented function !"); return false; }

double ARM_AnalyticPricingModel::VarianceToTime(double var,double minTime,double maxTime) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "VarianceToTime : unimplemented function !"); return 0.0; }

void ARM_AnalyticPricingModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "MCModelStatesFromToNextTime : unimplemented function !");}

void ARM_AnalyticPricingModel::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "TreeStatesToModelStates : unimplemented function !"); }

void ARM_AnalyticPricingModel::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr )
{}



CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
