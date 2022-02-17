/*
 * Copyright (c) IXIS CIB February 2007 Paris
 *


/*! \file EqModelParams.cpp
 *
 *  \brief 
 *
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date February 2007
 *
 *	\class ARM_EqModelParams contains the parameters of HW model in the contex of Equity Model with Stochastic Volatility
 */


#include "gpmodels/EqModelParams.h"

// c++ standard headers
#include <map>
#include  <stdlib.h>

using std::pair;

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_ModelParams
///	Routine: ARM_EQHWSV_ModelParams 
///	Returns: void
///	Action : constructor
////////////////////////////////////////////////////

ARM_EQHWSV_ModelParams::ARM_EQHWSV_ModelParams( const ARM_ModelParamVector& params ): ARM_ModelParams(params ){

	const ARM_GP_Vector& schedAt	= GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam().GetCurve()->GetAbscisses();
 
	if( schedAt.size() != 1)
		ARMTHROW(ERR_INVALID_ARGUMENT, "Mean Reversion should be a scalar" );
	
	MergeSchedule();
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_ModelParams
///	Routine: MergeSchedule
///	Returns: void
///	Action : merge all model parameter schedules
////////////////////////////////////////////////////

void ARM_EQHWSV_ModelParams::MergeSchedule(){
	int i;

    const ARM_GP_Vector& schedAt	= GetModelParam( ARM_ModelParamType::Volatility			).ToCurveModelParam().GetCurve()->GetAbscisses();
    const ARM_GP_Vector& schedBt	= GetModelParam( ARM_ModelParamType::Correlation		).ToCurveModelParam().GetCurve()->GetAbscisses();
    const ARM_GP_Vector& schedCt	= GetModelParam( ARM_ModelParamType::CompoundVol		).ToCurveModelParam().GetCurve()->GetAbscisses();
 	const ARM_GP_Vector& schedDt	= GetModelParam( ARM_ModelParamType::VolOfVol			).ToCurveModelParam().GetCurve()->GetAbscisses();
 	const ARM_GP_Vector& schedEt	= GetModelParam( ARM_ModelParamType::VolMeanReversion	).ToCurveModelParam().GetCurve()->GetAbscisses();
	const ARM_GP_Vector& schedFt	= GetModelParam( ARM_ModelParamType::LongTermVol		).ToCurveModelParam().GetCurve()->GetAbscisses();
 
	map<int,int> tmpMap;

	for( i=0; i< schedAt.size(); i++)
		tmpMap.insert( pair<int,int> (schedAt[i],0) );

	for( i=0; i< schedBt.size(); i++)
		tmpMap.insert( pair<int,int> (schedBt[i],0) );

	for( i=0; i< schedCt.size(); i++)
		tmpMap.insert( pair<int,int> (schedCt[i],0) );

	for( i=0; i< schedDt.size(); i++)
		tmpMap.insert( pair<int,int> (schedDt[i],0) );

	for( i=0; i< schedEt.size(); i++)
		tmpMap.insert( pair<int,int> (schedEt[i],0) );

	for( i=0; i< schedFt.size(); i++)
		tmpMap.insert( pair<int,int> (schedFt[i],0) );

	std::map<int,int>::iterator it;
	for(it= tmpMap.begin(); it!=tmpMap.end(); it++)
		itsSchedule.push_back( it->first );

	if(itsSchedule[0] > K_NEW_DOUBLE_TOL)
        itsSchedule.insert(itsSchedule.begin(),0.0);
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_ModelParams
///	Routine: CptSchedule
///	Returns: ARM_GP_Vector
///	Action : return the pertinent date between t and T
////////////////////////////////////////////////////
const	ARM_GP_Vector & ARM_EQHWSV_ModelParams::CptSchedule( const double& t, const double & T) const{

	double tmin= t<T? t:T;
	double tmax= t<T? T:t;
	ARM_GP_Vector res(1,tmin);

	size_t nb = itsSchedule.size();
	for( int i = 0; i< nb; ++i){
		if ( itsSchedule[i]< tmax && itsSchedule[i]> tmin )
			res.push_back ( itsSchedule[i] );
	}
	return res;
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////

ARM_EQHWSV_ModelParams& ARM_EQHWSV_ModelParams::operator = (const ARM_EQHWSV_ModelParams & rhs){
	if(this != &rhs)
		ARM_ModelParams::operator=(rhs);
	itsSchedule = rhs.itsSchedule;

	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_NumMethods
///	Routine: ARM_EQHWSV_NumMethods
///	Returns: void
///	Action : constructor
////////////////////////////////////////////////////

ARM_EQHWSV_NumMethods::ARM_EQHWSV_NumMethods( const double & maxdecay, const double & imaxis, const int & intstep ):itsMaxDecay(maxdecay),itsImAxis(imaxis),itsIntStep(intstep){}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_NumMethods
///	Routine: ARM_EQHWSV_NumMethods
///	Returns: void
///	Action : recopy constructor
////////////////////////////////////////////////////

ARM_EQHWSV_NumMethods::ARM_EQHWSV_NumMethods( const ARM_EQHWSV_NumMethods & rhs ):itsMaxDecay(rhs.itsMaxDecay),itsImAxis(rhs.itsImAxis),itsIntStep(rhs.itsIntStep){}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_NumMethods
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////

ARM_EQHWSV_NumMethods& ARM_EQHWSV_NumMethods::operator = (const ARM_EQHWSV_NumMethods & rhs){
	if(this != &rhs)
		*this =  ARM_EQHWSV_NumMethods(rhs);
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV_NumMethods
///	Routine: toString
///	Returns: string
///	Action : provide the view of the class
////////////////////////////////////////////////////
string ARM_EQHWSV_NumMethods::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Numerical parameters involved in the Green function computation method\n";
    os << indent << "---------------------------------------------------------------------\n\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<CC_NS(std,right)<< "IntStep";
	os <<"  :  "<< itsIntStep<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<CC_NS(std,right)<< "ImAxis";
	os <<"  :  "<< itsImAxis<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<CC_NS(std,right)<< "MaxDecay";
	os <<"  :  "<< itsMaxDecay<<"\n";
    return os.str();
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

