
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Model Visitor															 *
 *																							 *
 *			This class builds a generic design pattern: visitor								 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 23rd 2007														 *	
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#include "gpinflation/infpayoffvisitor.h"
#include "gpinflation/infhybridpayoffvisitor.h"
#include "gpinflation/infoptionspreadvisitor.h"

#include "gpinflation/infmodelvisitor.h"
#include "gpinflation/infhkvisitor.h"
#include "gpinflation/infbilogvisitor.h"

#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal.h"


CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////

///	Class  : ARM_InfNoModel
///	Routine: toString
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

string ARM_InfNoModel::toString(	const string& indent, const string& nextIndent) const{

	CC_Ostringstream	os;

	const int LAG = 13;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(58)<<"HYBRIID INFLATION MODEL\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<<"DESCRIPTION\n";
	os	<<"================="<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Name";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< itsName;
	os <<"\n";

	return os.str();
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIntModel
///	Routine: toString
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

string ARM_InfIntModel::toString(	const string& indent, const string& nextIndent) const{

	CC_Ostringstream	os;

	const int LAG = 13;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(58)<<"HYBRID INFLATION MODEL\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<<"DESCRIPTION\n";
	os	<<"==========="<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Name";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< itsName;
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Int Disc";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< itsDis;
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Int Domain";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< itsDom;
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Epsilon";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< itsEps;
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Center";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<< itsCen;
	os <<"\n";


	return os.str();
}



/////////////////////////////////////////////////////

///	Class  : ARM_InfIntModel
///	Routine: ARM_InfIntModel
///	Returns: void
///	Action : Constructor

/////////////////////////////////////////////////////

ARM_InfIntModel::ARM_InfIntModel(	const string & name, 
									const double & dis, 
									const double & dom, 
									const double & eps,
									const double & cen):ARM_InfModel(name),
														itsDis(dis), 
														itsDom(dom), 
														itsEps(eps),
														itsCen(cen){


	GaussLegendre_Coefficients c( (int) dis, cen-dom, cen+dom);
	itsPosit= ARM_GP_VectorPtr ( new ARM_GP_Vector (  (int) dis ) );
	itsWeigt= ARM_GP_VectorPtr ( new ARM_GP_Vector (  (int) dis ) );

	for ( int i = 0; i< (int) dis; i++){
		itsPosit->Elt(i) =	c.get_point(i);
		itsWeigt->Elt(i) =	c.get_weight(i)*exp(-0.5*itsPosit->Elt(i)*itsPosit->Elt(i))/sqrt(2.0*PI);
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_InfIntModel
///	Routine: ARM_InfIntModel
///	Returns: void
///	Action : copy Constructor

/////////////////////////////////////////////////////
ARM_InfIntModel::ARM_InfIntModel( const ARM_InfIntModel & rhs):ARM_InfModel(rhs),itsDis(rhs.itsDis),itsDom(rhs.itsDom),itsEps(rhs.itsEps),itsCen(rhs.itsCen){
	itsPosit = CreateClonedPtr( &*rhs.itsPosit) ;
	itsWeigt = CreateClonedPtr( &*rhs.itsWeigt) ;
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfNoModelValue
///	Routine: operator()
///	Returns: double
///	Action : evaluation of the payoff

/////////////////////////////////////////////////////

double ARM_InfNoModelValue::operator() (	const ARM_InfPayOffValuePtr		& payOff, 	
											const double					& res,
											const ARM_MAP_Double			& fwd,
											const ARM_MAP_PairDb			& cor,											
											const ARM_MAP_VolPar			& vol ){
	ARM_MAP_Double par;
	return (*payOff)(res, fwd, par );	
}

/////////////////////////////////////////////////////

///	Class  : ARM_InfModelVisitor
///	Routine: Visit()
///	Returns: void
///	Action : attribuate the appropriate model

/////////////////////////////////////////////////////

void ARM_InfModelVisitor::Visit( const ARM_InfNumBiLog & payOff){
		ARM_InfNumBiLogValue* tmp = new ARM_InfNumBiLogValue(payOff);
		itsInfModelValue	= ARM_InfModelValuePtr ( tmp );		
}

void ARM_InfModelVisitor::Visit( const ARM_InfBiLog & payOff){
		ARM_InfBiLogValue* tmp = new ARM_InfBiLogValue(payOff);
		itsInfModelValue	= ARM_InfModelValuePtr ( tmp );		
}


void ARM_InfModelVisitor::Visit( const ARM_InfHK & payOff){
		ARM_InfHKValue* tmp = new ARM_InfHKValue(payOff);
		itsInfModelValue	= ARM_InfModelValuePtr ( tmp );		
} 

void ARM_InfModelVisitor::Visit( const ARM_InfNoModel & payOff){
		ARM_InfNoModelValue* tmp = new ARM_InfNoModelValue(payOff);
		itsInfModelValue	= ARM_InfModelValuePtr ( tmp );		
} 


CC_END_NAMESPACE()






/*--------------------------------------------------------------------------*/
/*---- End of file ----*/



