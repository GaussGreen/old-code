/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 08/15/2006
 *
 *  basic functions for the closed form framework 
 *
 *	\file fxconvertdata.cpp
 *
 *  \brief
 *
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date August 2006
 */
#include <glob/firsttoinc.h>
#include <glob/expt.h>   // for the exceptions
#include "gpbase/port.h"
#include "gpbase/utilityport.h"
#include "gpbase/gpmatrix.h"

#include "gpclosedforms/fxconvertdata.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/optimization1.h"

#include "inst/portfolio.h"
#include "mod/bsmodel.h"

#include "nag.h"
#include "nage04.h"

using namespace std;

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine:  Copy constructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_FXMktDataToTotemFormat::ARM_FXMktDataToTotemFormat(const ARM_FXMktDataToTotemFormat& rhs)
: 
  itsPortfolio(CreateClone(rhs.itsPortfolio)),
  itsIndexDeltaCalls(rhs.itsIndexDeltaCalls),
  itsIndexDeltaPuts(rhs.itsIndexDeltaPuts),
  itsATmVol(rhs.itsATmVol),
  itsUnkown(rhs.itsUnkown),
  itsLowBound(rhs.itsLowBound),
  itsUpBound(rhs.itsUpBound),
  itsFctObj(rhs.itsFctObj),
  itsMaturityIndex(rhs.itsMaturityIndex),
  itsMktModel(CreateClone(rhs.itsMktModel)),
  itsMaxIter(rhs.itsMaxIter),
  itsAlgoType(rhs.itsAlgoType),
  itsTolerance(rhs.itsTolerance)
{}

  ////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine:  destructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_FXMktDataToTotemFormat::~ARM_FXMktDataToTotemFormat() 
{
	delete itsPortfolio; 
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine:  Constructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_FXMktDataToTotemFormat::ARM_FXMktDataToTotemFormat(
	const ARM_StdPortfolio& portfolio,
	const ARM_GP_Vector& deltaCalls,
	const ARM_GP_Vector& deltaPuts,
	double atmVol,
	const ARM_BSModel& model,
	const ARM_GP_Vector& InitPoint,
	const ARM_GP_Vector& LowerBound,
	const ARM_GP_Vector& UpperBound,
	double maxIter,
	double algoType,
	double Tolerance)
: 
  itsIndexDeltaCalls(),
  itsIndexDeltaPuts(),
  itsATmVol(atmVol),
  itsUnkown(&InitPoint ? InitPoint :ARM_GP_Vector()),
  itsLowBound(&LowerBound ? LowerBound :ARM_GP_Vector()),
  itsUpBound(&UpperBound ? UpperBound :ARM_GP_Vector()),
  itsMaxIter(maxIter),
  itsAlgoType(algoType),
  itsTolerance(Tolerance),
  itsFctObj(0),
  itsMaturityIndex()
{
	itsMktModel = CreateClone(const_cast< ARM_BSModel*>(&model) );
	itsPortfolio = CreateClone(const_cast< ARM_StdPortfolio*>(&portfolio));
	ARM_FXVolCurve* vol = static_cast<ARM_FXVolCurve*>(itsMktModel->GetVolatility());
	ARM_Matrix RR = vol->GetRR();
	ARM_Matrix STR = vol->GetSTR();

	double asofdate = itsMktModel->GetStartDateJul();

	ARM_StdPortfolio* pf = dynamic_cast<ARM_StdPortfolio*>(itsPortfolio->GetAsset(0));
	ARM_Date expiry = pf ? pf->GetAsset(0)->GetMaturity() : itsPortfolio->GetAsset(0)->GetMaturity();

	double mat = (expiry.GetJulian()- asofdate)/ K_YEAR_LEN;

	itsMaturityIndex = vol->GetOptionsMatus().find(mat);

	const ARM_Vector* pivotVol = vol->GetPivotVols();
	(*pivotVol)[itsMaturityIndex] = atmVol;

	size_t size = 0;
	double x;

	ARM_GP_Vector unknown;
	ARM_GP_Vector lowbound;
	ARM_GP_Vector upbound;

	/// Risk Reversal
	if(!deltaPuts.empty())
	{
		if(deltaPuts.size() == 2)
		{
			itsIndexDeltaPuts.push_back(0);
			x = RR.Elt(itsMaturityIndex,0);
			unknown.push_back(x);
			lowbound.push_back(x - 0.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));

			itsIndexDeltaPuts.push_back(1);
			x = RR.Elt(itsMaturityIndex,1);
			unknown.push_back(x);
			lowbound.push_back(x - 0.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
		}
		else if(deltaPuts.size() == 1 && fabs(deltaPuts[0] + 10.00) < ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE){
			itsIndexDeltaPuts.push_back(0);
			x = RR.Elt(itsMaturityIndex,0);
			unknown.push_back(x);
			lowbound.push_back(x - 0.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
		}
		
		else if(deltaPuts.size() == 1 && fabs(deltaPuts[0] + 25.00) < ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE){
			itsIndexDeltaPuts.push_back(1);
			x = RR.Elt(itsMaturityIndex,1);
			unknown.push_back(x);
			lowbound.push_back(x - 0.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
		}

		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Delta Puts are not valid");
	}
	if(!deltaCalls.empty())
	{
		if(deltaCalls.size() == 2)
		{
			itsIndexDeltaCalls.push_back(0);
			x = STR.Elt(itsMaturityIndex,0);
			unknown.push_back(x);
			lowbound.push_back(x - 1.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
			
			itsIndexDeltaCalls.push_back(1);	
			x = STR.Elt(itsMaturityIndex,1);
			unknown.push_back(x);
			lowbound.push_back(x - 1.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
		}
		else if(deltaCalls.size() == 1 && fabs(deltaCalls[0] - 10.00) < ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE){
			itsIndexDeltaCalls.push_back(0);
			x = STR.Elt(itsMaturityIndex,0);
			unknown.push_back(x);
			lowbound.push_back(x - 1.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
		}
		
		else if(deltaCalls.size() == 1 && fabs(deltaCalls[0] - 25.00) < ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE){
			itsIndexDeltaCalls.push_back(1);
			x = STR.Elt(itsMaturityIndex,0);
			unknown.push_back(x);
			lowbound.push_back(x - 1.5*fabs(x));
			upbound.push_back(x + 0.5*fabs(x));
		}

		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Delta Calls are not valid");
	}
	if(itsUnkown.empty())
		itsUnkown = unknown;
	if(itsLowBound.empty())
		itsLowBound = lowbound;
	if(itsUpBound.empty())
		itsUpBound = upbound;

	Calibrate();
	UpdateParams(itsUnkown);
	size_t m = itsPortfolio->size();
	ARM_GP_Vector tmpf(m);
	WeightedSquareFunc( itsUnkown,tmpf);

}

////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine:  UpdateParams
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_FXMktDataToTotemFormat::UpdateParams(const ARM_GP_Vector& unknown) 
{
	ARM_FXVolCurve* vol = static_cast<ARM_FXVolCurve*>(itsMktModel->GetVolatility());
	ARM_Matrix RR = vol->GetRR();
	ARM_Matrix STR = vol->GetSTR();
	size_t j=0;
	if(!itsIndexDeltaPuts.empty())
		for(int i(0); i<itsIndexDeltaPuts.size(); ++i)
			RR.Elt(itsMaturityIndex,itsIndexDeltaPuts[i]) = unknown[j++];

	if(!itsIndexDeltaCalls.empty())
		for(int i(0); i<itsIndexDeltaCalls.size(); ++i)
			STR.Elt(itsMaturityIndex,itsIndexDeltaCalls[i]) = unknown[j++];

	vol->GenerateFXVols(STR,RR);
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine: WeightedSquareFunc
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_FXMktDataToTotemFormat::WeightedSquareFunc(const ARM_GP_Vector& unknown, ARM_GP_Vector& result)
{
	UpdateParams(unknown);
	double sumWeight = itsPortfolio->GetWeights()->sum();
	size_t m = itsPortfolio->size();

	for(int i(0); i<m; i++)
	{
        double precision = (*itsPortfolio->GetPrecision())[i];
		double weight = (*itsPortfolio->GetWeights())[i];
        double sqrtWeight = sqrt(weight/sumWeight);
        double a = sqrtWeight/precision;
		itsPortfolio->SetModel(itsMktModel);
		double fx = itsPortfolio->GetAsset(i)->ComputePrice();		
		double target = (*itsPortfolio->GetMktPrices())[i];
		fx-=target;
		fx*=a;

        result[i]=fx;
	}
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine: Derivatives
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_FXMktDataToTotemFormat::Derivatives(const ARM_GP_Vector& var,ARM_GP_Matrix& fjac)
{
	size_t m = itsPortfolio->size();
	ARM_GP_Vector result(m);
	ARM_GP_Vector resultDown(m);
	size_t n = var.size();
	ARM_GP_Vector unknown(var);
	for(int i(0); i<n; i++)
	{
		double x =  unknown[i];
		double epsilon = 0.001 * CC_Max<double>(x, 1.0 );
		double twiceepsilon = 2.0*epsilon;
		unknown[i]+=epsilon;
		WeightedSquareFunc(unknown, result);
		unknown[i]-=twiceepsilon;
		WeightedSquareFunc(unknown, resultDown);
		///reset unknown
		unknown[i]+=epsilon;
		result-=resultDown;

		result/=twiceepsilon;

		fjac.push_backColumn(result);
	}
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_FXMktDataToTotemFormat
///	Routine: Calibrate
///	Returns: void
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_FXMktDataToTotemFormat::Calibrate()
{
	class ObjFunction : public Optimization_ObjectiveFuntion
	{
	private:
		ARM_FXMktDataToTotemFormat itsMktDataToTotemFormat;
	  
	public:
		virtual void NAG_CALL operator() (
			Integer m,
			Integer n, 
			double x[],		/// input
			double f[],		/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac,
			Nag_Comm *comm )
		{
			ARM_GP_Vector var(n,x);
			ARM_GP_Vector tmpf(m);
			itsMktDataToTotemFormat.WeightedSquareFunc(var,tmpf);
			ARM_GP_Matrix tmpfjac(n,0);
			itsMktDataToTotemFormat.Derivatives(var,tmpfjac);
			for(int i(0); i<m; ++i){
				f[i] = tmpf[i];
				for(int j(0); j<n; ++j){
					fjac[m*i+j] = tmpfjac(i,j);
				}
			}
		}
		ObjFunction(const ARM_FXMktDataToTotemFormat& mktData)
		:
		Optimization_ObjectiveFuntion(),
		itsMktDataToTotemFormat(mktData)
        {}	
		~ObjFunction() {};
	};
	ObjFunction func(*this);

	int n = itsUnkown.size();
	size_t m = itsPortfolio->size();

	Optimization_Result_Set* result= OptimizeWithDerivatives(&ARM_GP_Vector(m,0.0),
		&ARM_GP_Vector(m,0.0),
		&ARM_GP_Vector(m,1.0),
		&func,
		&itsUnkown,
		&itsLowBound,
		&itsUpBound,
		itsAlgoType,
		false,
		//"C:\\tmp\\FxMktToTotem.txt",
		"",
		itsTolerance,  // Tolerance
		itsMaxIter);	// max_iter*/

	itsUnkown = *result->OptimalParamSet;

}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
