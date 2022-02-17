
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Hybrid PayOff Visitor													 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


#include "gpinflation\infhybridpayoffvisitor.h"
#include "gpinflation/argconvdefault.h"

#include "gpclosedforms/normal.h"



CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridPayOff
///	Routine: ARM_InfHybridPayOff
///	Returns: void
///	Action : constructor

/////////////////////////////////////////////////////

ARM_InfHybridPayOff::ARM_InfHybridPayOff(	const ARM_MAP_Curve &	cpnCurve,
											const ARM_MAP_Curve &	optCurve):ARM_InfPayOff(){
	
	ARM_MAP_Curve::const_iterator it;

	itsOptCurve.clear();
	for ( it = optCurve.begin(); it != optCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsOptCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}

	itsCpnCurve.clear();
	for ( it = cpnCurve.begin(); it != cpnCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second) );
		itsCpnCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridPayOff
///	Routine: ARM_InfHybridPayOff
///	Returns: void
///	Action : copy constructor

/////////////////////////////////////////////////////

ARM_InfHybridPayOff::ARM_InfHybridPayOff ( const ARM_InfHybridPayOff & rhs){

	
	ARM_MAP_Curve::const_iterator it;

	itsOptCurve.clear();
	for ( it = rhs.itsOptCurve.begin(); it != rhs.itsOptCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsOptCurve.insert(p);
	}

	itsCpnCurve.clear();
	for ( it = rhs.itsCpnCurve.begin(); it != rhs.itsCpnCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second) );
		itsCpnCurve.insert(p);
	}

	itsVolType.clear();
	itsVolType = rhs.itsVolType;
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridPayOff
///	Routine: ValidateSecurity
///	Returns: void
///	Action : validation of the underlying (from securty) with respect the payoff

/////////////////////////////////////////////////////
void ARM_InfHybridPayOff::ValidateSecurity ( const vector<string> & vec){

	ARM_MAP_Curve tmp = itsOptCurve;
	ARM_MAP_Curve ::iterator it;
	std::vector<std::string>::const_iterator is;

	for ( it = tmp.begin(); it != tmp.end(); it++){
		is = std::find( vec.begin(),vec.end(), it->first );
		if ( is == vec.end() )
			ARMTHROW(ERR_INVALID_ARGUMENT," payoff is incoherent with the security");	
	}

	tmp = itsCpnCurve;
	for ( it = tmp.begin(); it != tmp.end(); it++){
		is = std::find( vec.begin(),vec.end(), it->first );
		if ( is == vec.end() )
			ARMTHROW(ERR_INVALID_ARGUMENT," payoff is incoherent with the security");	
	}
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridPayOff
///	Routine: toString
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

ARM_GP_Vector MergeCurve ( ARM_GP_Vector abs, ARM_GP_CurvePtr cur){

	std::map<double,int>		tmpMap;
	std::map<double,int>::iterator	it;
	ARM_GP_Vector				tmp;

	for ( int i = 0; i< abs.size(); i++)
		tmpMap.insert(std::pair<double,int> ( abs[i],1) );

	tmp = cur->GetAbscisses();
	for ( i = 0; i< tmp.size(); i++)
		tmpMap.insert(std::pair<double,int> ( tmp[i],1) );

	tmp.clear();
	for( it=tmpMap.begin(); it!=tmpMap.end(); it++)
		tmp.push_back(it->first);

	return tmp;
}


string ARM_InfHybridPayOff::toString(	const string& indent, const string& nextIndent) const{


	CC_Ostringstream	os;

	const int LAG = 13;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(58)<<"HYBRID INFLATION PAYOFF\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<<"COUPON DESCRIPTION\n";
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
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<(int) it->second->Interpolate(lag);
	}

	os	<< "\n\n";	
	os	<<"OPTION DESCRIPTION\n";
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
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<(int) it->second->Interpolate(lag);
	}	

	return os.str();
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridPayOffValue
///	Routine: operator()
///	Returns: double
///	Action : compute the payoff

/////////////////////////////////////////////////////

double ARM_InfHybridPayOffValue::operator() ( const  double &	lag, const ARM_MAP_Double & index, const ARM_MAP_Double & param ){

	ARM_MAP_Double:: const_iterator jt;

	ARM_MAP_Curve tmp =itsPayOff->GetCpnCurve();
	double cpn = tmp["NO"]->Interpolate	(lag) ;

	for (Iter_Curve it = tmp.begin(); it!= tmp.end(); it++){
		if ( it->first != "NO") {
			jt = index.find(it->first);
			cpn += tmp[it->first]->Interpolate	(lag) * jt->second;	
		}
	}

	double opt = CptOpt (lag,  index);	
	opt = opt>0?1.0:0.0;
	return cpn*opt;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfHybridPayOffValue
///	Routine: CptOpt
///	Returns: double
///	Action : compute the argument of the option

/////////////////////////////////////////////////////

double ARM_InfHybridPayOffValue::CptOpt	(	const double & lag , const ARM_MAP_Double &  index){

	ARM_MAP_Double:: const_iterator jt;

	ARM_MAP_Curve tmp =itsPayOff->GetOptCurve();
	double opt = tmp["NO"]->Interpolate	(lag) ;

	for (Iter_Curve it = tmp.begin(); it!= tmp.end(); it++){
		if ( it->first != "NO") {
			jt = index.find(it->first);
			opt += tmp[it->first]->Interpolate	(lag) * jt->second;	
		}
	}	
	
	return opt;

}




CC_END_NAMESPACE()





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















