
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Option PayOff Visitor													 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


#include "gpinflation\infoptionspreadvisitor.h"

#include "gpclosedforms/normal.h"
#include "gpinflation\infvisitor.h"

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridCap
///	Routine: ARM_InfHybridCap
///	Returns: void
///	Action : constructor

/////////////////////////////////////////////////////

ARM_InfHybridCap::ARM_InfHybridCap(	const ARM_MAP_Curve &	optCurve):ARM_InfHybridPayOff(){
	
	ARM_MAP_Curve::const_iterator it;

	for ( it = optCurve.begin(); it != optCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		std::pair<string,ARM_GP_CurvePtr> q ( it->first, CreateClonedPtr( &*it->second ) );
		itsOptCurve.insert(p);
		itsCpnCurve.insert(q);
		
		itsVolType.insert(std::pair<string,VolType> (it->first, CAP) );
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridCapValue
///	Routine: operator()
///	Returns: double
///	Action : compute the payoff

/////////////////////////////////////////////////////

double call( double fwd, double str, double var){

	double d1,d2;
	d1= log(fwd/str);
	d1/=var;

	d2 = d1-0.5*var;
	d1+= 0.5*var;

	return fwd * NormalCDF(d1) - str * NormalCDF(d2);
}

double ARM_InfHybridCapValue::operator() ( const  double &	lag, const  ARM_MAP_Double & index, const  ARM_MAP_Double & param   ){


	ARM_MAP_Double::const_iterator it=param.begin();
	string stg = it->first;
	double var = it->second;
	ARM_MAP_Double::const_iterator kt =index.find(it->first);
	double fwd =kt->second;


	ARM_MAP_Curve tmp =itsPayOff->GetOptCurve();
	double str = tmp["NO"]->Interpolate(lag);
	for (Iter_Curve jt = tmp.begin(); jt!= tmp.end(); jt++){
		if( jt->first!=stg && jt->first != "NO" ){
			kt =index.find(jt->first);
			str += tmp[jt->first]->Interpolate	(lag) * kt->second;
		}
	}

	fwd *= tmp[stg]->Interpolate(lag);
	str *= -1.0;

	if		( fwd >= 0.0 && str <= 0.0) 
		return fwd-str;
	else if ( fwd <= 0.0 && str <= 0.0 )
		return call( -fwd,-str,var)-(str - fwd);
	else if ( fwd >= 0.0 && str >= 0.0 )
		return call( fwd,str,var);
	else
		return 0.0;
}




/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridDigit
///	Routine: ARM_InfHybridDigit
///	Returns: void
///	Action : constructor

/////////////////////////////////////////////////////

ARM_InfHybridDigit::ARM_InfHybridDigit(	const ARM_MAP_Curve &	optCurve, const ARM_GP_CurvePtr & notional):ARM_InfHybridPayOff(){
	
	ARM_MAP_Curve::const_iterator it;

	for ( it = optCurve.begin(); it != optCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsOptCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first, DIG) );
	}
	itsCpnCurve.insert(std::pair<string,ARM_GP_CurvePtr>( "NO", notional ) );
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridDigit
///	Routine: operator()
///	Returns: double
///	Action : compute the payoff

/////////////////////////////////////////////////////

double Digit( double fwd, double str, double var){

	double d1;
	d1= log(fwd/str);
	d1/=var;
	d1-= 0.5*var;

	return NormalCDF(d1);
}


double ARM_InfHybridDigitValue::operator() ( const  double &	lag, const  ARM_MAP_Double & index, const  ARM_MAP_Double & param   ){


	ARM_MAP_Double::const_iterator it=param.begin();
	string stg = it->first;
	double var = it->second;
	ARM_MAP_Double::const_iterator kt =index.find(it->first);
	double fwd =kt->second;


	ARM_MAP_Curve tmp =itsPayOff->GetOptCurve();
	double str = tmp["NO"]->Interpolate(lag);
	for (Iter_Curve jt = tmp.begin(); jt!= tmp.end(); jt++){
		if( jt->first!=stg && jt->first != "NO" ){
			kt =index.find(jt->first);
			str += tmp[jt->first]->Interpolate	(lag) * kt->second;
		}
	}
	fwd *= tmp[stg]->Interpolate(lag);
	str *= -1.0;

	tmp =itsPayOff->GetCpnCurve();
	double notional = tmp["NO"]->Interpolate(lag);

	if		( fwd >= 0.0 && str <= 0.0) 
		return notional;
	else if ( fwd <= 0.0 && str <= 0.0 )
		return notional*(1.0-Digit( -fwd, -str, var) );
	else if ( fwd >= 0.0 && str >= 0.0 )
		return notional*Digit( fwd, str, var);
	else
		return 0.0;
}





CC_END_NAMESPACE()





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















