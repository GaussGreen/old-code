
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Corridor Visitor														 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: July, 5th 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


#include "gpinflation/infcorridorvisitor.h"
#include "gpinflation/infhybridpayoffvisitor.h"
#include "gpinflation/argconvdefault.h"

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////

///	Class  : ARM_InfCorridor
///	Routine: ARM_InfCorridor
///	Returns: void
///	Action : constructor

/////////////////////////////////////////////////////

ARM_InfCorridor::ARM_InfCorridor( const ARM_MAP_Curve &	cpnCurve, const ARM_MAP_Curve &	optCurve):ARM_InfHybridPayOff( cpnCurve ,optCurve){

	ARM_MAP_Curve::const_iterator it;

	itsCpnCurve.clear();
	for ( it = cpnCurve.begin(); it != cpnCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsCpnCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}

	itsOptCurve.clear();
	for ( it = optCurve.begin(); it != optCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsOptCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}

}


/////////////////////////////////////////////////////

///	Class  : ARM_InfCorridor
///	Routine: ARM_InfCorridor
///	Returns: void
///	Action : copy constructor

/////////////////////////////////////////////////////

ARM_InfCorridor::ARM_InfCorridor(	const ARM_InfCorridor &	rhs):ARM_InfHybridPayOff(rhs){

	ARM_MAP_Curve::const_iterator it;

	itsOptCurve.clear();
	for ( it = rhs.itsOptCurve.begin(); it != rhs.itsOptCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsOptCurve.insert(p);
	}
}



/////////////////////////////////////////////////////

///	Class  : ARM_InfCorridor
///	Routine: toString
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

string ARM_InfCorridor::toString(	const string& indent, const string& nextIndent) const{


	CC_Ostringstream	os;

	const int LAG = 13;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(58)<<"HYBRID CORRIDOR\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<<"LOWER BOUND DESCRIPTION\n";
	os	<<"================="<<"\n";

	ARM_GP_Vector	tmp;
	ARM_GP_CurvePtr	cur;
	Iter_Curve		it;

	it  = itsCpnCurve.begin();
	cur = it->second;
	tmp = cur->GetAbscisses();

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Lag";
	for ( it = itsCpnCurve.begin(); it != itsCpnCurve.end(); it++ ){

		tmp  = MergeCurve ( tmp, it->second );
		if ( it->first =="NO")
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Strike";
		else
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<it->first;
	}
	os	<< "\n";
	for ( int i = 0; i< tmp.size(); i++){

		double lag = tmp[i];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< lag;
		for ( it = itsCpnCurve.begin(); it != itsCpnCurve.end(); it++ )
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<it->second->Interpolate(lag);
	}

	os	<< "\n\n";	
	os	<<"UPPER BOUND DESCRIPTION\n";
	os	<<"================="<<"\n";

	it  = itsCpnCurve.begin();
	cur = it->second;
	tmp = cur->GetAbscisses();

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"LAG";
	for ( it = itsOptCurve.begin(); it != itsOptCurve.end(); it++ ){
		tmp  = MergeCurve ( tmp, it->second );
		if ( it->first =="NO")
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"strike";
		else
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<it->first;
	}
	os	<< "\n";
	for ( i = 0; i< tmp.size(); i++){

		double lag = tmp[i];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< lag;
		for ( it = itsOptCurve.begin(); it != itsOptCurve.end(); it++ )
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<it->second->Interpolate(lag);
	}	

	return os.str();

}
/////////////////////////////////////////////////////

///	Class  : ARM_InfCorridorValue
///	Routine: operator()
///	Returns: double
///	Action : compute the payoff

//			 We recall that	cpnCurve corresponds to the Lower Bound 
//			 and optCurve corresponds to the Upper Bound 
/////////////////////////////////////////////////////

double ARM_InfCorridorValue::operator()	( const  double &	lag, const  ARM_MAP_Double & index, const  ARM_MAP_Double & param   ){

	ARM_MAP_Double:: const_iterator jt;

	ARM_MAP_Curve tmp = itsPayOff->GetCpnCurve();
	double cpn =0.0;
	for (Iter_Curve it = tmp.begin(); it!= tmp.end(); it++){
		if ( it->first != "NO") {
			jt = index.find(it->first);
			cpn += tmp[it->first]->Interpolate	(lag) * jt->second;	
		}
	}

	double lower = cpn + tmp["NO"]->Interpolate	(lag) ;
	lower = lower>0.0?1.0:0.0;

	tmp = itsPayOff->GetCpnCurve();
	double upper = cpn + tmp["NO"]->Interpolate	(lag) ;
	upper = upper>0.0?1.0:0.0;

	return lower-upper;

}

CC_END_NAMESPACE()





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















