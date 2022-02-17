
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Double Digital PayOff Visitor											 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


#include "gpinflation/infdoubledigitalvisitor.h"
#include "gpinflation/infhybridpayoffvisitor.h"
#include "gpinflation/argconvdefault.h"

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////

///	Class  : ARM_InfDoubleDigit
///	Routine: ARM_InfDoubleDigit
///	Returns: void
///	Action : constructor

/////////////////////////////////////////////////////

ARM_InfDoubleDigit::ARM_InfDoubleDigit(	const ARM_MAP_Curve		&	cpnCurve,
										const ARM_MAP_Curve		&	m_optCurve,
										const ARM_MAP_Curve		&	s_optCurve):ARM_InfHybridPayOff(){

	ARM_MAP_Curve::const_iterator it;

	itsCpnCurve.clear();
	for ( it = cpnCurve.begin(); it != cpnCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second) );
		itsCpnCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}

	itsM_OptCurve.clear();
	for ( it = m_optCurve.begin(); it != m_optCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsM_OptCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}
	
	itsS_OptCurve.clear();
	for ( it = s_optCurve.begin(); it != s_optCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsS_OptCurve.insert(p);
		itsVolType.insert(std::pair<string,VolType> (it->first,DIG) );
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfDoubleDigit
///	Routine: ARM_InfDoubleDigit
///	Returns: void
///	Action : copyconstructor

/////////////////////////////////////////////////////

ARM_InfDoubleDigit::ARM_InfDoubleDigit(	const ARM_InfDoubleDigit &	rhs):ARM_InfHybridPayOff(rhs){

	ARM_MAP_Curve::const_iterator it;

	itsM_OptCurve.clear();
	for ( it = rhs.itsM_OptCurve.begin(); it != rhs.itsM_OptCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsM_OptCurve.insert(p);
	}

	itsS_OptCurve.clear();
	for ( it = rhs.itsS_OptCurve.begin(); it != rhs.itsS_OptCurve.end(); it++){
		std::pair<string,ARM_GP_CurvePtr> p ( it->first, CreateClonedPtr( &*it->second ) );
		itsS_OptCurve.insert(p);
	}
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfDoubleDigit
///	Routine: ValidateSecurity
///	Returns: void
///	Action : validation of the underlying (from securty) with respect the payoff

/////////////////////////////////////////////////////

void ARM_InfDoubleDigit::ValidateSecurity ( const vector<string> & vec){

	ARM_InfHybridPayOff::ValidateSecurity(vec);

	ARM_MAP_Curve tmp = itsM_OptCurve;
	ARM_MAP_Curve ::iterator it;
	std::vector<std::string>::const_iterator is;

	for ( it = tmp.begin(); it != tmp.end(); it++){
		is = std::find( vec.begin(),vec.end(), it->first );
		if ( is == vec.end() )
			ARMTHROW(ERR_INVALID_ARGUMENT," payoff is incoherent with the security");	
	}

	tmp = itsS_OptCurve;
	for ( it = tmp.begin(); it != tmp.end(); it++){
		is = std::find( vec.begin(),vec.end(), it->first );
		if ( is == vec.end() )
			ARMTHROW(ERR_INVALID_ARGUMENT," payoff is incoherent with the security");	
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfDoubleDigit
///	Routine: toString
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

string ARM_InfDoubleDigit::toString(	const string& indent, const string& nextIndent) const{


	CC_Ostringstream	os;

	const int LAG = 13;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(58)<<"DOUBLE DIGITAL INFLATION\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<<"COUPON DESCRIPTION\n";
	os	<<"======================"<<"\n";

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
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Spread";
		else
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<it->first;
	}
	os	<< "\n";
	for ( int i = 0; i< tmp.size(); i++){

		double lag = tmp[i];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< lag;
		for ( it = itsCpnCurve.begin(); it != itsCpnCurve.end(); it++ )
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< (int) it->second->Interpolate(lag);
	}

	os	<< "\n\n";	
	os	<<"MAIN DIGITAL DESCRIPTION\n";
	os	<<"====================="<<"\n";

	it  = itsM_OptCurve.begin();
	cur = it->second;
	tmp = cur->GetAbscisses();

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"LAG";
	for ( it = itsM_OptCurve.begin(); it != itsM_OptCurve.end(); it++ ){
		tmp  = MergeCurve ( tmp , it->second );
		if ( it->first =="NO")
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"strike";
		else
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<it->first;
	}
	os	<< "\n";
	for ( i = 0; i< tmp.size(); i++){

		double lag = tmp[i];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< lag;
		for ( it = itsM_OptCurve.begin(); it != itsM_OptCurve.end(); it++ )
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< it->second->Interpolate(lag);
	}	

	os	<< "\n\n";	
	os	<<"SUB DIGITAL DESCRIPTION\n";
	os	<<"====================="<<"\n";

	it  = itsS_OptCurve.begin();
	cur = it->second;
	tmp = cur->GetAbscisses();

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"LAG";
	for ( it = itsS_OptCurve.begin(); it != itsS_OptCurve.end(); it++ ){
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
		for ( it = itsS_OptCurve.begin(); it != itsS_OptCurve.end(); it++ )
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< it->second->Interpolate(lag);
	}	

	
	return os.str();
}



/////////////////////////////////////////////////////

///	Class  : ARM_InfDoubleDigitValue
///	Routine: operator()
///	Returns: double
///	Action : compute the payoff

/////////////////////////////////////////////////////

double ARM_InfDoubleDigitValue::operator()	( const  double &	lag, const  ARM_MAP_Double & index, const  ARM_MAP_Double & param   ){

	ARM_MAP_Double:: const_iterator jt;

	ARM_MAP_Curve tmp = itsPayOff->GetCpnCurve();
	double cpn = tmp["NO"]->Interpolate	(lag) ;
	for (Iter_Curve it = tmp.begin(); it!= tmp.end(); it++){
		if ( it->first != "NO") {
			jt = index.find(it->first);
			cpn += tmp[it->first]->Interpolate	(lag) * jt->second;	
		}
	}


	tmp =dynamic_cast<ARM_InfDoubleDigit*> ( &*itsPayOff )->GetM_OptCurve();
	double m_opt = tmp["NO"]->Interpolate	(lag) ;
	for (it = tmp.begin(); it!= tmp.end(); it++){
		if ( it->first != "NO") {
			jt = index.find(it->first);
			m_opt += tmp[it->first]->Interpolate(lag) * jt->second;	
		}
	}
	m_opt = m_opt>0?1.0:0.0;

	tmp =dynamic_cast<ARM_InfDoubleDigit*> ( &*itsPayOff )->GetS_OptCurve();
	double s_opt = tmp["NO"]->Interpolate	(lag) ;
	for (it = tmp.begin(); it!= tmp.end(); it++){
		if ( it->first != "NO") {
			jt = index.find(it->first);
			s_opt += tmp[it->first]->Interpolate(lag) * jt->second;	
		}
	}
	s_opt = s_opt>0?1.0:0.0;

	return cpn * m_opt * s_opt;

}

CC_END_NAMESPACE()





/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


















