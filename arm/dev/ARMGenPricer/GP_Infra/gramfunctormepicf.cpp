/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctormepicf.cpp
 *
 *  \brief gramfunctorcf is the file for closed form related grammar functor
 *	the functor design allows to have a context and a unified command like
 *	interface
 *	\author  O. Croissant
 *	\version 1.0
 *	\date June 2005
 */




/// gpbase
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/checkinputs.h"
#include "gpbase/stringmanip.h"
#include "gpbase/singleton.h"
#include "gpbase/datemanip.h"

//#include <libCCTools++\CCString.h>

/// gpinfra
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramfunctormepicf.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/curvemodelparam.h"


/// gpmodel
#include "gpmodels/MultiAssets.h"
#include "gpmodels/SABR_Eq.h"
#include "gpmodels/SABR_ModelParams.h"
#include "gpmodels/typedef.h"

/// gpclosedforms
#include "gpclosedforms/cppi_options.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/sabrimpliedvol.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"

/// gpnummethods
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/finummethod.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/scheduler.h"

CC_BEGIN_NAMESPACE( ARM )


// MC default pricing values
const double MC_NB_INTER_STEPS		= 1;

/// the corresponding static string
string ARM_CF_MepiCallFctor::itsFuncName		= "CFMepiCall";
string ARM_CF_MepiGreekFctor::itsFuncName		= "CFMepiCallGreek";

/*
ARM_VectorPtr ARM_GramFctorArg_ConvertToVector( const ARM_GramFctorArg& arg, size_t vecSize )
{
	if( GFAT_VECTOR_TYPE == arg.GetType() )
		return arg.GetVector();
	else if( GFAT_DOUBLE_TYPE == arg.GetType() )
		return ARM_VectorPtr( new ARM_GP_Vector( vecSize, arg.GetDouble() ) );
	else
       ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : can convert only double to vector!");
}
*/

////////////////////////////////////////////////////
///	Class  : ARM_CF_MepiCallFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the closed form function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_CF_MepiCallFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 22, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[0],  GFAT_VECTOR_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[1],  GFAT_VECTOR_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[2],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[3],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[4],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[5],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[6],  GFAT_VECTOR_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[7],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[8],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[9],  GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[10], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[11], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[12], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[13], GFAT_VECTOR_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[14], GFAT_STRING_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[15], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[16], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[17], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[18], GFAT_DOUBLE_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[19], GFAT_STRING_TYPE, ARM_CF_MepiCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[20], GFAT_STRING_TYPE, ARM_CF_MepiCallFctor::itsFuncName );

	ARM_CF_MepiDispatcher::DistributionType distType = (ARM_CF_MepiDispatcher::DistributionType) ARM_ArgConv_cfMepiDistType.GetNumber( arg[21].GetString() );

	ARM_VectorPtr result;
	bool vectorialFormula = false;
	size_t vecSize = 0;

	for( size_t i=0; i<13; ++i )
	{
		if(	GFAT_VECTOR_TYPE == arg[i].GetType() )
		{
			vectorialFormula = true;
			vecSize = arg[i].GetVector()->size();
		}
	}

	// fix for the unability ofthe parser to see that the mot clé is stochastic : its depends on the states via the model
		vectorialFormula = true;
		ARM_PricingStates* pt=&(*states);
		vecSize = pt->size();
		// fix

		ARM_CountedPtr<ARM_ZeroCurve> zerocurvePtr;
		ARM_PricingModelPtr EquityModelPtr;
		ARM_PricingModelPtr  DiscountingModelPtr;
		ARM_MultiAssetsModel* multiassetPtr2;
		
				// MRGK5 Random Generator
		ARM_RandomGeneratorPtr  pBaseRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::MRGK5,
				ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

		/// antithetic box muller
		ARM_RandomGeneratorPtr normRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::BoxMuller,
		pBaseRandomGen ) );

		/// antithetic variates!
		ARM_RandomGeneratorPtr numRandGen = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::AntitheticOne,
		normRandGen ) );
		
		if ( vectorialFormula )
	{
		ARM_VectorPtr spot				= ARM_GramFctorArg_ConvertToVector( arg[0], vecSize );
		ARM_VectorPtr strike			= ARM_GramFctorArg_ConvertToVector( arg[1], vecSize );
		double relative_maturity		= arg[2].GetDouble();
		double Risk_Factor				= arg[3].GetDouble();
		double Min_Exposure				= arg[4].GetDouble();
		double Max_Exposure				= arg[5].GetDouble();
		ARM_VectorPtr Initial_Protection		= ARM_GramFctorArg_ConvertToVector( arg[6], vecSize );
		double Final_Protection			= arg[7].GetDouble();
		double borrowing_spread			= arg[8].GetDouble();
		double Yearly_Fees				= arg[9].GetDouble();
		double CashSpread				= arg[10].GetDouble();
		double periods_number			= arg[11].GetDouble();		/// in fact integer coded as a double
		double avg_period_number		= arg[12].GetDouble();		/// in fact integer coded as a double
		ARM_VectorPtr AlreadyAsianReset	= ARM_GramFctorArg_ConvertToVector( arg[13], vecSize );
		double callOrPut				= ARM_ArgConv_CallPut.GetNumber( arg[14].GetString() );
		double PortMin					= arg[15].GetDouble();
		double PortMax					= arg[16].GetDouble();
		double port_discret_size		= arg[17].GetDouble();		/// in fact integer coded as a double
		double hermite_or_MC_size		= arg[18].GetDouble();		/// in fact integer coded as a double
		
		string CurveName				= arg[19].GetString();		
		string EquityName				= arg[20].GetString();		
		
		

		size_t resultSize			= spot->size();
		result = ARM_VectorPtr( new ARM_GP_Vector( resultSize ) );

		ARM_MultiAssetsModel* multiassetPtr=dynamic_cast<ARM_MultiAssetsModel*>(mod);
		if (multiassetPtr)
		{
			ARM_ModelNameMap* modelMapPtr=multiassetPtr->GetModelMap();
			ARM_ModelNameMap::iterator iterEquityModel =(*modelMapPtr)[EquityName];
			EquityModelPtr=(*iterEquityModel).Model();
			ARM_ModelNameMap::iterator iterDiscountingModel =(*modelMapPtr)[CurveName];
			DiscountingModelPtr=(*iterDiscountingModel).Model();
			zerocurvePtr=DiscountingModelPtr->GetZeroCurve();
		
			
		}
		else
		{
			zerocurvePtr=mod->GetZeroCurve();
			ARM_PricingModelPtr EquityModelPtr1(mod);
			EquityModelPtr=EquityModelPtr1;
		}

		/// for now the zero coupon curve sees just its date adjusted to the necessary date
		ARM_ZeroCurve*  zerocurve=&(*zerocurvePtr);
		ARM_ZeroCurve*  forward_zerocurve =dynamic_cast<ARM_ZeroCurve*>(zerocurve->Clone());
		forward_zerocurve->SetAsOfDate(ARM_Date::ARM_Date(evalDate));
		ARM_SABR_Eq* SABR_Model= dynamic_cast<ARM_SABR_Eq*>(&(*EquityModelPtr));


		size_t modelNb		= SABR_Model->GetModelNb();
		double new_f;		
		double new_Alpha;
		ARM_ModelParam* new_Alpha_param;
		ARM_SABR_Eq* SABR_Model_clone=static_cast<ARM_SABR_Eq *> (SABR_Model->Clone());
		ARM_GP_Vector* breakPointTimes = new ARM_GP_Vector(1,0.0);
		SABR_Model_clone->SetModelRank(0);
		SABR_Model_clone->SetModelNb(0);


		///  int NbIterations=(*(SABR_Model_clone->GetNumMethod())).
		ARM_TimeStepPerYearScheduler scheduler(MC_NB_INTER_STEPS);
		ARM_NormalCentredSamplerND sampler(&scheduler);

		int NbIterations=hermite_or_MC_size;
		ARM_MCMethod* mcMethod = new ARM_MCMethod(NbIterations,numRandGen,&sampler);

		SABR_Model_clone->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
		ARM_ModelParamsSABR_Eq* mdps;

		/// Construction d'un multiasset
		if (multiassetPtr)
		{
			multiassetPtr2=dynamic_cast<ARM_MultiAssetsModel*>(multiassetPtr->Clone());
			ARM_StringVector names (2);
			vector<ARM_PricingModelPtr> models (2);
			ARM_StringVectorVector depends(2);
			names[0] = CurveName; // plus généralement, le petit nom de la ccy
			names[1] = EquityName;
			models[0] = DiscountingModelPtr;
			models[1] = EquityModelPtr;
			ARM_StringVector sv0 (1,names[0]);
			ARM_StringVector sv1 (1,names[0]);
			depends[0] = sv0;
			depends[1] = sv1;
			ARM_ModelNameMap* modelMap  = new ARM_ModelNameMap (names, models, depends);
			multiassetPtr2->SetModelMap(modelMap);
			/// association de la methode numerique et du numeraire
			multiassetPtr2->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
			multiassetPtr2->SetNumeraire(multiassetPtr->GetNumeraire());
		}
		
		switch( distType )
		{
			case ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
				for( i=0; i<resultSize; ++i )
				{	
					new_f=		states->GetModelState(i,modelNb);
					new_Alpha=	states->GetModelState(i,modelNb+1);
					ARM_GP_Vector new_f_values(1,new_Alpha);
					ARM_GP_Vector new_Alpha_values(1,new_Alpha);
					new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
					SABR_Model_clone->GetModelParams()->SetModelParam(new_Alpha_param,0);

					(*result)[i] = mepi_VanillaOption_STOBS(evalDate,(*spot)[i],(*strike)[i],relative_maturity,forward_zerocurve, CurveName,
														borrowing_spread, Yearly_Fees, CashSpread,
														SABR_Model_clone,EquityName,
														Min_Exposure, Max_Exposure, Risk_Factor, 
														(*Initial_Protection)[i], Final_Protection, 
														PortMin, PortMax,
														(int) avg_period_number,
														(*AlreadyAsianReset)[i],
														(int) callOrPut,(int) periods_number,
														(int) port_discret_size,(int) hermite_or_MC_size);
					delete new_Alpha_param;
														
				}
					
				break;
			case ARM_CF_MepiDispatcher::K_SABR_DIST:
				for( i=0; i<resultSize; ++i )
				{
					new_f=		states->GetModelState(i,modelNb);
					new_Alpha=	states->GetModelState(i,modelNb+1);
					// test
			//		new_f=1.0;
			//		new_Alpha=0.05;
					// test
					ARM_GP_Vector new_Alpha_values(1,new_Alpha);
					new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
					// iniçtialization du spot par new_f
					mdps=dynamic_cast<ARM_ModelParamsSABR_Eq*>(SABR_Model_clone->GetModelParams());
					mdps->SetSpot(new_f);
					
					if (multiassetPtr)
					{
						(*result)[i] = mepi_VanillaOption_SABR(evalDate,(*spot)[i],(*strike)[i],ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
							borrowing_spread, Yearly_Fees, CashSpread,
							multiassetPtr2,EquityName,
							Min_Exposure, Max_Exposure, Risk_Factor, 
							(*Initial_Protection)[i], Final_Protection, 
							PortMin, PortMax,
							(int) avg_period_number,
							(*AlreadyAsianReset)[i],(int) periods_number);
					}
					else
					{
						(*result)[i] = mepi_VanillaOption_SABR(evalDate,(*spot)[i],(*strike)[i],ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
							borrowing_spread, Yearly_Fees, CashSpread,
							SABR_Model_clone,EquityName,
							Min_Exposure, Max_Exposure, Risk_Factor, 
							(*Initial_Protection)[i], Final_Protection, 
							PortMin, PortMax,
							(int) avg_period_number,
							(*AlreadyAsianReset)[i],(int) periods_number);
					}
					delete new_Alpha_param;
				}
					
				break;

			case ARM_CF_MepiDispatcher::K_SABR_MC_DIST:
				for( i=0; i<resultSize; ++i )
				{
					/// used in the MACRO ARM_RESULT
					
					ARM_GenSecurity * gensec = NULL;
					(*result)[i] = 0;
				}
				
				break;
			default:
				{
					delete forward_zerocurve;
					delete SABR_Model_clone;
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
				}
		}
		delete forward_zerocurve;
		delete SABR_Model_clone;
		return ARM_GramFctorArg(result );
	}
	else
	{
		double spot					= arg[0].GetDouble();
		double strike				= arg[1].GetDouble();
		double relative_maturity	= arg[2].GetDouble();
		double Risk_Factor			= arg[3].GetDouble();
		double Min_Exposure			= arg[4].GetDouble();
		double Max_Exposure			= arg[5].GetDouble();
		double Initial_Protection 	= arg[6].GetDouble();
		double Final_Protection		= arg[7].GetDouble();
		double borrowing_spread		= arg[8].GetDouble();
		double Yearly_Fees			= arg[9].GetDouble();
		double CashSpread			= arg[10].GetDouble();
		double periods_number		= arg[11].GetDouble();		/// in fact integer coded as a double
		double avg_period_number	= arg[12].GetDouble();		/// in fact integer coded as a double
		double AlreadyAsianReset	= arg[13].GetDouble();
		double callOrPut			= ARM_ArgConv_CallPut.GetNumber( arg[14].GetString() );
		double PortMin				= arg[15].GetDouble();
		double PortMax				= arg[16].GetDouble();
		double port_discret_size	= arg[17].GetDouble();		/// in fact integer coded as a double
		double hermite_or_MC_size	= arg[18].GetDouble();		/// in fact integer coded as a double
		
		string CurveName				= arg[19].GetString();		
		string EquityName				= arg[20].GetString();	
 
		double resultDouble = 0.0;

		ARM_MultiAssetsModel* multiassetPtr=dynamic_cast<ARM_MultiAssetsModel*>(mod);
		if (multiassetPtr)
		{
			ARM_ModelNameMap* modelMapPtr=multiassetPtr->GetModelMap();

			ARM_ModelNameMap::iterator iterEquityModel =(*modelMapPtr)[EquityName];
			EquityModelPtr=(*iterEquityModel).Model();

			ARM_ModelNameMap::iterator iterDiscountingModel =(*modelMapPtr)[CurveName];
			DiscountingModelPtr=(*iterDiscountingModel).Model();
			zerocurvePtr=DiscountingModelPtr->GetZeroCurve();
		
			
		}
		else
		{
			zerocurvePtr=mod->GetZeroCurve();
			ARM_PricingModelPtr EquityModelPtr1(mod);
			EquityModelPtr=EquityModelPtr1;
		}

		/// for now the zero coupon curve sees just its date adjusted to the necessary date
		ARM_ZeroCurve*  zerocurve=&(*zerocurvePtr);
		ARM_ZeroCurve*  forward_zerocurve =dynamic_cast<ARM_ZeroCurve*>(zerocurve->Clone());
		forward_zerocurve->SetAsOfDate(ARM_Date::ARM_Date(evalDate));
		
		ARM_SABR_Eq* SABR_Model= dynamic_cast<ARM_SABR_Eq*>(&(*EquityModelPtr));
		size_t modelNb		= SABR_Model->GetModelNb();
		double new_f;		
		double new_Alpha;
		ARM_ModelParam* new_Alpha_param;
		ARM_SABR_Eq* SABR_Model_clone=static_cast<ARM_SABR_Eq *> (SABR_Model->Clone());
		ARM_GP_Vector* breakPointTimes = new ARM_GP_Vector(1,0.0);	
		SABR_Model_clone->SetModelRank(0);
		SABR_Model_clone->SetModelNb(0);

		ARM_TimeStepPerYearScheduler scheduler(MC_NB_INTER_STEPS);
		ARM_NormalCentredSamplerND sampler(&scheduler);

		int NbIterations=hermite_or_MC_size;
		ARM_MCMethod* mcMethod = new ARM_MCMethod(NbIterations,numRandGen,&sampler);

		SABR_Model_clone->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
		ARM_ModelParamsSABR_Eq* mdps;
		/// Construction d'un multiasset
		if (multiassetPtr)
		{
			multiassetPtr2=dynamic_cast<ARM_MultiAssetsModel*>(multiassetPtr->Clone());
			ARM_StringVector names (2);
			vector<ARM_PricingModelPtr> models (2);
			ARM_StringVectorVector depends(2);
			names[0] = CurveName; // plus généralement, le petit nom de la ccy
			names[1] = EquityName;
			models[0] = DiscountingModelPtr;
			models[1] = EquityModelPtr;
			ARM_StringVector sv0 (1,"");
			ARM_StringVector sv1 (1,names[0]);
			depends[0] = sv0;
			depends[1] = sv1;
			ARM_ModelNameMap* modelMap  = new ARM_ModelNameMap (names, models, depends);
			multiassetPtr2->SetModelMap(modelMap);
			/// association de la methode numerique et du numeraire
			multiassetPtr2->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
			multiassetPtr2->SetNumeraire(multiassetPtr->GetNumeraire());
			
		}


		switch( distType )
		{
			case ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
				{
				new_f=		states->GetModelState(0,modelNb);
				new_Alpha=	states->GetModelState(0,modelNb+1);
				ARM_GP_Vector new_Alpha_values(1,new_Alpha);
				new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
				SABR_Model_clone->GetModelParams()->SetModelParam(new_Alpha_param,0);

				resultDouble = mepi_VanillaOption_STOBS(evalDate,spot,strike,relative_maturity,forward_zerocurve, CurveName,
														borrowing_spread, Yearly_Fees, CashSpread,
														SABR_Model_clone,EquityName,
														Min_Exposure, Max_Exposure, Risk_Factor, 
														Initial_Protection, Final_Protection, 
														PortMin, PortMax,
														(int) avg_period_number,
														AlreadyAsianReset,
														(int) callOrPut,(int) periods_number,
														(int) port_discret_size,(int) hermite_or_MC_size);
				delete new_Alpha_param;
														
				}
				break;
			case ARM_CF_MepiDispatcher::K_SABR_DIST:
				{	
				new_f=		states->GetModelState(i,modelNb);
				new_Alpha=	states->GetModelState(i,modelNb+1);
				ARM_GP_Vector new_Alpha_values(1,new_Alpha);
				new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
				// initialization du spot par new_f
				mdps=dynamic_cast<ARM_ModelParamsSABR_Eq*>(SABR_Model_clone->GetModelParams());
				mdps->SetSpot(new_f);
				if (multiassetPtr)
				{
					resultDouble = mepi_VanillaOption_SABR(evalDate,spot,strike,ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
						borrowing_spread, Yearly_Fees, CashSpread,
						multiassetPtr2,EquityName,
						Min_Exposure, Max_Exposure, Risk_Factor, 
						Initial_Protection, Final_Protection, 
						PortMin, PortMax,
						(int) avg_period_number,
						AlreadyAsianReset,(int) periods_number);
				}
				else
				{
					
					resultDouble = mepi_VanillaOption_SABR(evalDate,spot,strike,ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
						borrowing_spread, Yearly_Fees, CashSpread,
						SABR_Model_clone,EquityName,
						Min_Exposure, Max_Exposure, Risk_Factor, 
						Initial_Protection, Final_Protection, 
						PortMin, PortMax,
						(int) avg_period_number,
						AlreadyAsianReset,(int) periods_number);
				}

				delete new_Alpha_param;
														
				}
				break;


			default:
				{
					delete forward_zerocurve;
					delete SABR_Model_clone;
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
				}
		}
		delete forward_zerocurve;
		delete SABR_Model_clone;
		return ARM_GramFctorArg(resultDouble );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_CF_MepiCallFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_CF_MepiCallFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr( new ARM_GP_Vector( 0 )),ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
////////////////////// ARM_CF_MepiGreekFctor
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CF_MepiGreekFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the closed form function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_CF_MepiGreekFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 24, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_VECTOR_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[1], GFAT_VECTOR_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[2], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[3], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[4], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[5], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_VECTOR_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[7], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[8], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[9], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[10], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[11], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[12], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[13], GFAT_VECTOR_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[14], GFAT_STRING_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[15], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[16], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[17], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[18], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[19], GFAT_DOUBLE_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[20], GFAT_STRING_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[21], GFAT_STRING_TYPE, ARM_CF_MepiGreekFctor::itsFuncName );

	ARM_VectorPtr result;
	bool vectorialFormula = false;
	size_t vecSize = 0;

	for( size_t i=0; i<13; ++i )
	{
		if(	GFAT_VECTOR_TYPE == arg[i].GetType() )
		{
			vectorialFormula = true;
			vecSize = arg[i].GetVector()->size();
		}
	}
	
	// fix for the unability ofthe parser to see that the mot clé is stochastic : its depends on the states via the model
	vectorialFormula = true;
	ARM_PricingStates* pt=&(*states);
	vecSize = pt->size();
	// fix


	ARM_CountedPtr<ARM_ZeroCurve> zerocurvePtr;
	ARM_PricingModelPtr EquityModelPtr;
	
	ARM_PricingModelPtr  DiscountingModelPtr;
	ARM_MultiAssetsModel* multiassetPtr2;
	ARM_ModelParamsSABR_Eq* mdps;
	
	// MRGK5 Random Generator
	ARM_RandomGeneratorPtr  pBaseRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::MRGK5,
		ARM_RandGenFactoryImp::UnknownTransformAlgo ) );
	
	/// antithetic box muller
	ARM_RandomGeneratorPtr normRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::BoxMuller,
		pBaseRandomGen ) );
	
	/// antithetic variates!
	ARM_RandomGeneratorPtr numRandGen = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::AntitheticOne,
		normRandGen ) );		
	
	ARM_CF_MepiDispatcher::DistributionType distType = (ARM_CF_MepiDispatcher::DistributionType) ARM_ArgConv_cfMepiDistType.GetNumber( arg[22].GetString() );
	ARM_CF_MepiDispatcher::GreekType greekType		= (ARM_CF_MepiDispatcher::GreekType) ARM_ArgConv_cfGreekType.GetNumber( arg[23].GetString() );

	if ( vectorialFormula )
	{
		ARM_VectorPtr spot				= ARM_GramFctorArg_ConvertToVector( arg[0], vecSize );
		ARM_VectorPtr strike			= ARM_GramFctorArg_ConvertToVector( arg[1], vecSize );
		double relative_maturity		= arg[2].GetDouble();
		double Risk_Factor				= arg[3].GetDouble();
		double Min_Exposure				= arg[4].GetDouble();
		double Max_Exposure				= arg[5].GetDouble();
		ARM_VectorPtr Initial_Protection		= ARM_GramFctorArg_ConvertToVector( arg[6], vecSize );
		double Final_Protection			= arg[7].GetDouble();
		double borrowing_spread			= arg[8].GetDouble();
		double Yearly_Fees				= arg[9].GetDouble();
		double CashSpread				= arg[10].GetDouble();
		double periods_number			= arg[11].GetDouble();		/// in fact integer coded as a double
		double avg_period_number		= arg[12].GetDouble();		/// in fact integer coded as a double
		ARM_VectorPtr AlreadyAsianReset	= ARM_GramFctorArg_ConvertToVector( arg[13], vecSize );
		double callOrPut				= ARM_ArgConv_CallPut.GetNumber( arg[14].GetString() );
		double PortMin					= arg[15].GetDouble();
		double PortMax					= arg[16].GetDouble();
		double port_discret_size		= arg[17].GetDouble();		/// in fact integer coded as a double
		double hermite_or_MC_size		= arg[18].GetDouble();		/// in fact integer coded as a double
		double Shift_Size				= arg[19].GetDouble();		/// in fact integer coded as a double
		string CurveName				= arg[20].GetString();		
		string EquityName				= arg[21].GetString();	
		
		size_t resultSize		= spot->size();
		result = ARM_VectorPtr( new ARM_GP_Vector( resultSize ) );

		ARM_TimeStepPerYearScheduler scheduler(MC_NB_INTER_STEPS);
		ARM_NormalCentredSamplerND sampler(&scheduler);
		
		int NbIterations=hermite_or_MC_size;
		ARM_MCMethod* mcMethod = new ARM_MCMethod(NbIterations, numRandGen, &sampler);
		
		
		ARM_MultiAssetsModel* multiassetPtr=dynamic_cast<ARM_MultiAssetsModel*>(mod);
		if (multiassetPtr)
		{
			ARM_ModelNameMap* modelMapPtr=multiassetPtr->GetModelMap();
			
			ARM_ModelNameMap::iterator iterEquityModel =(*modelMapPtr)[EquityName];
			EquityModelPtr=(*iterEquityModel).Model();
			
			ARM_ModelNameMap::iterator iterDiscountingModel =(*modelMapPtr)[CurveName];
			DiscountingModelPtr=(*iterDiscountingModel).Model();
			zerocurvePtr=DiscountingModelPtr->GetZeroCurve();
			
			
		}
		else
		{
			zerocurvePtr=mod->GetZeroCurve();
			ARM_PricingModelPtr EquityModelPtr1(mod);
			EquityModelPtr=EquityModelPtr1;
		}
		/// for now the zero coupon curve sees just its date adjusted to the necessary date
		ARM_ZeroCurve*  zerocurve=&(*zerocurvePtr);
		ARM_ZeroCurve*  forward_zerocurve =dynamic_cast<ARM_ZeroCurve*>(zerocurve->Clone());
		forward_zerocurve->SetAsOfDate(ARM_Date::ARM_Date(evalDate));


		ARM_SABR_Eq* SABR_Model= dynamic_cast<ARM_SABR_Eq*>(&(*EquityModelPtr));
		size_t modelNb		= SABR_Model->GetModelNb();
		double new_f;		
		double new_Alpha;
		ARM_ModelParam* new_Alpha_param;
		ARM_SABR_Eq* SABR_Model_clone=static_cast<ARM_SABR_Eq *> (SABR_Model->Clone());
		ARM_GP_Vector* breakPointTimes = new ARM_GP_Vector(1,0.0);	
		SABR_Model_clone->SetModelRank(0);
		SABR_Model_clone->SetModelNb(0);

		/// ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		/// SABR_Model_clone->SetNumeraire(numeraire);

		SABR_Model_clone->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
		
		
		ARM_ModelParamsSABR_Eq* mdps;
		/// Construction d'un multiasset
		if (multiassetPtr)
		{
			multiassetPtr2=dynamic_cast<ARM_MultiAssetsModel*>(multiassetPtr->Clone());
			ARM_StringVector names (2);
			vector<ARM_PricingModelPtr> models (2);
			ARM_StringVectorVector depends(2);
			names[0] = CurveName; // plus généralement, le petit nom de la ccy
			names[1] = EquityName;
			models[0] = DiscountingModelPtr;
			models[1] = EquityModelPtr;
			ARM_StringVector sv0 (1,names[0]);
			ARM_StringVector sv1 (1,names[0]);
			depends[0] = sv0;
			depends[1] = sv1;
			ARM_ModelNameMap* modelMap  = new ARM_ModelNameMap (names, models, depends);
			multiassetPtr2->SetModelMap(modelMap);
			/// association de la methode numerique et du numeraire
			multiassetPtr2->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
			multiassetPtr2->SetNumeraire(multiassetPtr->GetNumeraire());
		}
		
		switch( distType )
		{
		case ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
			switch( greekType )
			{
			case ARM_CF_MepiDispatcher::K_Delta:
				for( i=0; i<resultSize; ++i )
				{
					new_f=		states->GetModelState(i,modelNb);
					new_Alpha=	states->GetModelState(i,modelNb+1);
					ARM_GP_Vector new_Alpha_values(1,new_Alpha);
					new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
					// initialization du spot par new_f
					mdps=dynamic_cast<ARM_ModelParamsSABR_Eq*>(SABR_Model_clone->GetModelParams());
					mdps->SetSpot(new_f);

					(*result)[i] = mepi_VanillaOption_STOBS_delta(
														evalDate,(*spot)[i],(*strike)[i],relative_maturity,forward_zerocurve, CurveName,
														borrowing_spread, Yearly_Fees, CashSpread,
														SABR_Model_clone,EquityName,
														Min_Exposure, Max_Exposure, Risk_Factor, 
														(*Initial_Protection)[i], Final_Protection, 
														PortMin, PortMax,
														(int) avg_period_number,
														(*AlreadyAsianReset)[i],
														(int) callOrPut,(int) periods_number,
														(int) port_discret_size,(int) hermite_or_MC_size);
					delete new_Alpha_param;
				}
				break;


			default:
				{
					delete forward_zerocurve;
					delete SABR_Model_clone;
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
				}
			}
			break;
		
	case ARM_CF_MepiDispatcher::K_SABR_DIST:
			switch( greekType )
			{
			case ARM_CF_MepiDispatcher::K_Delta:
				for( i=0; i<resultSize; ++i )
				{
					new_f=		states->GetModelState(i,modelNb);
					new_Alpha=	states->GetModelState(i,modelNb+1);
					ARM_GP_Vector new_Alpha_values(1,new_Alpha);
					new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
					// initialization du spot par new_f
					mdps=dynamic_cast<ARM_ModelParamsSABR_Eq*>(SABR_Model_clone->GetModelParams());
					mdps->SetSpot(new_f);
					
					if (multiassetPtr)
					{
						(*result)[i] = mepi_VanillaOption_SABR_delta(
							evalDate,(*spot)[i],(*strike)[i],ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
							borrowing_spread, Yearly_Fees, CashSpread,
							multiassetPtr2,EquityName,
							Min_Exposure, Max_Exposure, Risk_Factor, 
							(*Initial_Protection)[i], Final_Protection, 
							PortMin, PortMax,
							(int) avg_period_number,
							(*AlreadyAsianReset)[i],(int) periods_number,Shift_Size);
					}
					else
					{
						(*result)[i] = mepi_VanillaOption_SABR_delta(
							evalDate,(*spot)[i],(*strike)[i],ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
							borrowing_spread, Yearly_Fees, CashSpread,
							SABR_Model_clone,EquityName,
							Min_Exposure, Max_Exposure, Risk_Factor, 
							(*Initial_Protection)[i], Final_Protection, 
							PortMin, PortMax,
							(int) avg_period_number,
							(*AlreadyAsianReset)[i],(int) periods_number,Shift_Size);
					}
					delete new_Alpha_param;
					
				}
				break;


			default:
				{
					delete forward_zerocurve;
					delete SABR_Model_clone;
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
				}
			}
			break;
		default:
			{
				delete forward_zerocurve;
				delete SABR_Model_clone;
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
			}
		}
		delete forward_zerocurve;
		delete SABR_Model_clone;
		return ARM_GramFctorArg(result );
	}
	else
	{
		double spot					= arg[0].GetDouble();
		double strike				= arg[1].GetDouble();
		double relative_maturity	= arg[2].GetDouble();
		double Risk_Factor			= arg[3].GetDouble();
		double Min_Exposure			= arg[4].GetDouble();
		double Max_Exposure			= arg[5].GetDouble();
		double Initial_Protection 	= arg[6].GetDouble();
		double Final_Protection		= arg[7].GetDouble();
		double borrowing_spread		= arg[8].GetDouble();
		double Yearly_Fees			= arg[9].GetDouble();
		double CashSpread			= arg[10].GetDouble();
		double periods_number		= arg[11].GetDouble();		/// in fact integer coded as a double
		double avg_period_number	= arg[12].GetDouble();		/// in fact integer coded as a double
		double AlreadyAsianReset	= arg[13].GetDouble();
		double callOrPut			= ARM_ArgConv_CallPut.GetNumber( arg[14].GetString() );
		double PortMin				= arg[15].GetDouble();
		double PortMax				= arg[16].GetDouble();
		double port_discret_size	= arg[17].GetDouble();		/// in fact integer coded as a double
		double hermite_or_MC_size	= arg[18].GetDouble();		/// in fact integer coded as a double
		double Shift_Size			= arg[19].GetDouble();		/// in fact integer coded as a double
		
		string CurveName				= arg[20].GetString();		
		string EquityName				= arg[21].GetString();	
		
		double resultDouble = 0.0;

		ARM_TimeStepPerYearScheduler scheduler(MC_NB_INTER_STEPS);
		ARM_NormalCentredSamplerND sampler(&scheduler);

		int NbIterations=hermite_or_MC_size;
		ARM_MCMethod* mcMethod = new ARM_MCMethod(NbIterations,numRandGen,&sampler);
		
		ARM_MultiAssetsModel* multiassetPtr=dynamic_cast<ARM_MultiAssetsModel*>(mod);
		if (multiassetPtr)
		{
			ARM_ModelNameMap* modelMapPtr=multiassetPtr->GetModelMap();

			ARM_ModelNameMap::iterator iterEquityModel =(*modelMapPtr)[EquityName];
			EquityModelPtr=(*iterEquityModel).Model();

			ARM_ModelNameMap::iterator iterDiscountingModel =(*modelMapPtr)[CurveName];
			DiscountingModelPtr=(*iterDiscountingModel).Model();
			zerocurvePtr=DiscountingModelPtr->GetZeroCurve();
		
			
		}
		else
		{
			zerocurvePtr=mod->GetZeroCurve();
			ARM_PricingModelPtr EquityModelPtr1(mod);
			EquityModelPtr=EquityModelPtr1;
		}

		/// for now the zero coupon curve sees just its date adjusted to the necessary date
		ARM_ZeroCurve*  zerocurve=&(*zerocurvePtr);
		ARM_ZeroCurve*  forward_zerocurve =dynamic_cast<ARM_ZeroCurve*>(zerocurve->Clone());
		forward_zerocurve->SetAsOfDate(ARM_Date::ARM_Date(evalDate));
		

		ARM_SABR_Eq* SABR_Model= dynamic_cast<ARM_SABR_Eq*>(&(*EquityModelPtr));
		size_t modelNb		= SABR_Model->GetModelNb();
		double new_f;		
		double new_Alpha;
		ARM_ModelParam* new_Alpha_param;
		ARM_SABR_Eq* SABR_Model_clone=static_cast<ARM_SABR_Eq *> (SABR_Model->Clone());
		ARM_GP_Vector* breakPointTimes = new ARM_GP_Vector(1,0.0);
		SABR_Model_clone->SetModelRank(0);
		SABR_Model_clone->SetModelNb(0);

		/// ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		/// SABR_Model_clone->SetNumeraire(numeraire);

		/// Construction d'un multiasset
		if (multiassetPtr)
		{
			multiassetPtr2=dynamic_cast<ARM_MultiAssetsModel*>(multiassetPtr->Clone());
			ARM_StringVector names (2);
			vector<ARM_PricingModelPtr> models (2);
			ARM_StringVectorVector depends(2);
			names[0] = CurveName; // plus généralement, le petit nom de la ccy
			names[1] = EquityName;
			models[0] = DiscountingModelPtr;
			models[1] = EquityModelPtr;
			ARM_StringVector sv0 (1,names[0]);
			ARM_StringVector sv1 (1,names[0]);
			depends[0] = sv0;
			depends[1] = sv1;
			ARM_ModelNameMap* modelMap  = new ARM_ModelNameMap (names, models, depends);
			multiassetPtr2->SetModelMap(modelMap);
			/// association de la methode numerique et du numeraire
			multiassetPtr2->SetNumMethod(ARM_NumMethodPtr( mcMethod ));
			multiassetPtr2->SetNumeraire(multiassetPtr->GetNumeraire());
		}

		switch( distType )
		{
		case ARM_CF_MepiDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
			switch( greekType )
			{
			case ARM_CF_MepiDispatcher::K_Delta:
				{
						new_f=		states->GetModelState(i,modelNb);
					new_Alpha=	states->GetModelState(i,modelNb+1);
					ARM_GP_Vector new_Alpha_values(1,new_Alpha);
					new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
					// initialization du spot par new_f
					mdps=dynamic_cast<ARM_ModelParamsSABR_Eq*>(SABR_Model_clone->GetModelParams());
					mdps->SetSpot(new_f);

					resultDouble = mepi_VanillaOption_STOBS_delta(
														evalDate,spot,strike,relative_maturity,forward_zerocurve, CurveName,
														borrowing_spread, Yearly_Fees, CashSpread,
														SABR_Model_clone,EquityName,
														Min_Exposure, Max_Exposure, Risk_Factor, 
														Initial_Protection, Final_Protection, 
														PortMin, PortMax,
														(int) avg_period_number,
														AlreadyAsianReset,
														(int) callOrPut,(int) periods_number,
														(int) port_discret_size,(int) hermite_or_MC_size);
					delete new_Alpha_param;
				}
				break;

			default:
				{
					delete forward_zerocurve;
					delete SABR_Model_clone;
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
				}
			}
			break;
		case ARM_CF_MepiDispatcher::K_SABR_DIST:
			switch( greekType )
			{
			case ARM_CF_MepiDispatcher::K_Delta:
				{
						new_f=		states->GetModelState(i,modelNb);
					new_Alpha=	states->GetModelState(i,modelNb+1);
					ARM_GP_Vector new_Alpha_values(1,new_Alpha);
					new_Alpha_param= new ARM_CurveModelParam( ARM_ModelParamType::Alpha,&new_Alpha_values, breakPointTimes );
					// initialization du spot par new_f
					mdps=dynamic_cast<ARM_ModelParamsSABR_Eq*>(SABR_Model_clone->GetModelParams());
					mdps->SetSpot(new_f);

					resultDouble = mepi_VanillaOption_SABR_delta(
														evalDate,spot,strike,ConvertXLDateToJulian(relative_maturity),&(*zerocurvePtr), CurveName,
														borrowing_spread, Yearly_Fees, CashSpread,
														SABR_Model_clone,EquityName,
														Min_Exposure, Max_Exposure, Risk_Factor, 
														Initial_Protection, Final_Protection, 
														PortMin, PortMax,
														(int) avg_period_number,
														AlreadyAsianReset,
														(int) periods_number,Shift_Size);
					delete new_Alpha_param;
				}
				break;

			default:
				{
					delete forward_zerocurve;
					delete SABR_Model_clone;
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
				}
			}
			break;
	
		}
		delete forward_zerocurve;
		delete SABR_Model_clone;
		return ARM_GramFctorArg(resultDouble );
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_CF_MepiGreekFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_CF_MepiGreekFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 0 )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///	Class  : ARM_NumericalCallOnMepiDeltaFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : computes time lag given the input!
////////////////////////////////////////////////////
void ARM_NumericalCallOnMepiDeltaFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel *mod, double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
}

////////////////////////////////////////////////////
///	Class  : ARM_NumericalCallOnMepiDeltaFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr
///	Action : computes time lag given the input!
////////////////////////////////////////////////////
ARM_NodeInfo ARM_NumericalCallOnMepiDeltaFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector(1,1)),ARM_AdditionalTimeInfoPtr(NULL) );
}


ARM_GramFctorArg ARM_NumericalCallOnMepiDeltaFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_GramFctorArg(NULL);
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

