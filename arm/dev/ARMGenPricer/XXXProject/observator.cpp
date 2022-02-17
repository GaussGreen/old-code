/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file observator.cpp
 *  \brief file for uniformized mkt
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date Agust 2006
 */

#include <gpinfra/pricingmodel.h>
#include <gpinfra\curvemodelparam.h>
#include <gpmodels\ModelParams_EqFxBase.h>
#include <gpmodels\Mixture_FX.h>
#include <gpbase/countedptr.h>
#include <gpbase/curvematrix.h>
#include <gpinflation/infcurv.h>
#include <gpinflation/infbsmodel.h>

#include <mod/bssmiled.h>
#include <crv/correlmanager.h>
#include <crv/zerocurv.h>
#include <crv/volcube.h>
#include <crv/volint.h>
#include <crv/zerointspreaded.h>
#include <inst/forex.h>
#include "xxxproject/observator.h"
#include "xxxproject/argconvdefault.h"
#include "xxxproject/mktdatas.h"

#include  <stdlib.h>

string ConvertDoubleToStringMatu( double d ){
	if( d==1.0/365)
		return "1D";
	else if( d == 2.0/365)
		return "2D";
	else if( d == 1.0/12)
		return "1M";
	else
		return YearTermToStringMatu(d);
}

CC_BEGIN_NAMESPACE( ARM )


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*					Static: ARM Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

vector<ARM_Subject* >	ARM_Subject::CreateSubject( ){

	vector<ARM_Subject* > tmpSub(25);

	tmpSub[0]		=	new ARM_YC_Subject;

	tmpSub[1]		=	new ARM_CAP_Subject< 0 >;	//	ATM
	tmpSub[2]		=	new ARM_CAP_Subject< 1 >;	//	RHO
	tmpSub[3]		=	new ARM_CAP_Subject< 2 >;	//	NU
	tmpSub[4]		=	new ARM_CAP_Subject< 3 >;	//	BETA

	tmpSub[5]		=	new ARM_OSW_Subject< 0 >;	//	ATM
	tmpSub[6]		=	new ARM_OSW_Subject< 1 >;	//	RHO
	tmpSub[7]		=	new ARM_OSW_Subject< 2 >;	//	NU
	tmpSub[8]		=	new ARM_OSW_Subject< 3 >;	//	BETA

	tmpSub[9]		=	new ARM_BS_Subject;
	tmpSub[10]		=	new	ARM_FX_Subject;

	tmpSub[11]		=	new	ARM_VFX_Subject< 4 >;	//	PIV
	tmpSub[12]		=	new	ARM_VFX_Subject< 5 >;	//	RR
	tmpSub[13]		=	new	ARM_VFX_Subject< 6 >;	//	STR

	tmpSub[14]		=	new	ARM_MIX_Subject< 7 >;	//	VOL
	tmpSub[15]		=	new	ARM_MIX_Subject< 8 >;	//	SMILE
	tmpSub[16]		=	new	ARM_MIX_Subject< 9 >;	//	SHIFT
	tmpSub[17]		=	new	ARM_MIX_Subject< 10>;	//	Q

	tmpSub[18]		=	new	ARM_COR_Subject;
	tmpSub[19]		=	new	ARM_DAT_Subject;

	tmpSub[20]		=	new	ARM_SO_Subject< 0 >;	// ATM
	tmpSub[21]		=	new	ARM_SO_Subject< 11 >;	// ADJ

	tmpSub[22]		=	new	ARM_INF_Subject;	

	tmpSub[23]		=	new	ARM_IBS_Subject<12>;	//	CPI
	tmpSub[24]		=	new	ARM_IBS_Subject<13>;	//	YOY

	return tmpSub;
}

vector<ARM_Subject* >  ARM_Subject::AttributeSubject( const string & str, vector<ARM_Subject* > * sub, ARM_Object * obj){


	vector<ARM_Subject* > tmpSub;

	if ( str == "DATE" ){
		tmpSub.resize(1);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[19];
	}

	else if ( str == "YC" ){
		tmpSub.resize(2);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[0];
		tmpSub[1]	=	(*sub)[19];
	}

	else if ( str == "CAPMOD"){
		if ( dynamic_cast<const ARM_BSSmiledModel*> (obj) ){
			tmpSub.resize(6);
			itsCcyDim	=	 1;
			tmpSub[0]	=	(*sub)[0];
			tmpSub[1]	=	(*sub)[1];
			tmpSub[2]	=	(*sub)[2];
			tmpSub[3]	=	(*sub)[3];
			tmpSub[4]	=   (*sub)[4];
			tmpSub[5]	=	(*sub)[19];
		}
		else if (  dynamic_cast<const ARM_BSModel*> (obj)  ){
			tmpSub.resize(3);
			itsCcyDim	=	 1;
			tmpSub[0]	=	(*sub)[0];
			tmpSub[1]	=	(*sub)[1];
			tmpSub[2]	=	(*sub)[19];
		}
	}

	else if ( str == "OSWMOD"){
		if ( dynamic_cast<const ARM_BSSmiledModel*> (obj) ){
			tmpSub.resize(6);
			itsCcyDim	=	 1;
			tmpSub[0]	=	(*sub)[0];
			tmpSub[1]	=	(*sub)[5];
			tmpSub[2]	=	(*sub)[6];
			tmpSub[3]	=	(*sub)[7];
			tmpSub[4]	=   (*sub)[8];
			tmpSub[5]	=	(*sub)[19];
		}
		else if (  dynamic_cast<const ARM_BSModel*> (obj)  ){
			tmpSub.resize(3);
			itsCcyDim	=	 1;
			tmpSub[0]	=	(*sub)[0];
			tmpSub[1]	=	(*sub)[5];
			tmpSub[2]	=	(*sub)[19];
		}
	}

	else if ( str == "YC_BASIS" ){
		ARM_BasisCurve * tmp= dynamic_cast<ARM_BasisCurve*> (obj);
		if ( tmp && tmp->GetBasisCurve()  ) {
			tmpSub.resize(3);
			itsCcyDim	=	 1;
			tmpSub[0]	=	(*sub)[0];
			tmpSub[1]	=	(*sub)[19];
			tmpSub[2]	=	(*sub)[9];
		}
	else{
		tmpSub.resize(2);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[0];
		tmpSub[1]	=	(*sub)[19];
		}
	}

	else if ( str == "FOREX" ){
		tmpSub.resize(1);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[10];
	}
	
	else if ( str == "FXMOD" ){
		if (  dynamic_cast<const ARM_BSModel*> (obj)  ){
			tmpSub.resize(7);
			itsCcyDim	=	 2;
			tmpSub[0]	=	 (*sub)[0];
			tmpSub[1]	=	 (*sub)[9];
			tmpSub[2]	=	 (*sub)[10];
			tmpSub[3]	=	 (*sub)[11];
			tmpSub[4]	=	 (*sub)[12];
			tmpSub[5]	=	 (*sub)[13];
			tmpSub[6]	=	 (*sub)[19];
		}
		else if (  dynamic_cast<const ARM_MixtureModel_Fx*> (obj)  ){
			tmpSub.resize(8);
			itsCcyDim	=	 2;
			tmpSub[0]	=	 (*sub)[0];
			tmpSub[1]	=	 (*sub)[9];
			tmpSub[2]	=	 (*sub)[10];
			tmpSub[3]	=	 (*sub)[14];
			tmpSub[4]	=	 (*sub)[15];
			tmpSub[5]	=	 (*sub)[16];
			tmpSub[6]	=	 (*sub)[17];
			tmpSub[7]	=	 (*sub)[19];
		}
	}

	else if ( str == "CORREL" ){
		tmpSub.resize(1);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[18];
	}

	else if ( str == "SOMOD" ){
		tmpSub.resize(8);
		itsCcyDim	=	 1;
			tmpSub[0]	=	(*sub)[0];
			tmpSub[1]	=	(*sub)[1];
			tmpSub[2]	=	(*sub)[2];
			tmpSub[3]	=	(*sub)[3];
			tmpSub[4]	=   (*sub)[4];
			tmpSub[5]	=   (*sub)[20];
			tmpSub[6]	=   (*sub)[21];
			tmpSub[7]	=	(*sub)[19];
	}

	else if ( str == "INF" ){
		tmpSub.resize(2);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[19];
		tmpSub[1]	=   (*sub)[22];
	}
	
	else if ( str == "INFMOD" ){
		tmpSub.resize(8);
		itsCcyDim	=	 1;
		tmpSub[0]	=	(*sub)[0];
		tmpSub[1]	=	(*sub)[1];
		tmpSub[2]	=	(*sub)[5];
		tmpSub[3]	=	(*sub)[18];
		tmpSub[4]	=	(*sub)[19];
		tmpSub[5]	=   (*sub)[22];
		tmpSub[6]	=	(*sub)[23];
		tmpSub[7]	=	(*sub)[24];
	}

	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "This market data is not already integrate to xxx project." );

	return tmpSub;
}



			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*					ARM Observator						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

int		ARM_Subject::itsCcyDim;

ARM_Observator::ARM_Observator( ARM_Object*	obj, V_Sub	sub, string & str) :itsObject(obj)  {

	It_Sub it;
	itsSubject.clear();
	itsCurrency.clear();

	for( it = sub.begin(); it != sub.end() ; it++)
		itsSubject.push_back( *it );

	for( it = sub.begin(); it != sub.end() ; it++)
		(*it)->Attach(this,str);

	itsCurrency = GetCurrency();
}

map<string, int>	ARM_Observator::GetCurrency( ){

		
	V_Str						tmp;
	M_Str						tmpCcy;
	M_Str::iterator				it;

	if ( itsCurrency.size() != 0)
		return itsCurrency;
	else {
		tmpCcy.clear();
		for(It_Sub  its	 = itsSubject.begin(); its!= itsSubject.end() ; its++){
	
			tmp.clear();
			tmp = (*its) -> GetCurrency( this->GetObject() );

			if (tmp.size() ==0 )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT, 	"the input market data did not have the same currency" );
			else {
				for( int i=0; i< tmp.size(); i++){
					pair<string, int>	p( tmp[i] , 1 );
					tmpCcy.insert( p );
				}
			}
		}
	}
	return tmpCcy;
}

void	ARM_Observator::SetDomCcy( string domCcy){

	if ( itsCurrency.find(domCcy) != itsCurrency.end() )
		itsCurrency.find(domCcy)->second = 2;
}

ARM_Object*		ARM_Observator::ReduceObject	( ){

	for(int i = 0; i < itsSubject.size(); i++){
		if( dynamic_cast<ARM_COR_Subject*> ( itsSubject[i] ) && dynamic_cast<ARM_VolCurve*> ( itsObject) )
			return ARM_COR_Subject::Convert( dynamic_cast<ARM_VolCurve*> ( itsObject) );
	}
	return itsObject;
}


ARM_Observator::~ARM_Observator(){

	for( It_Sub it= itsSubject.begin(); it!= itsSubject.end() ; it++)
		(*it)->Detach(this);
	itsSubject.clear();
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						ARM Subject						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


vector<string>	ARM_Subject::GetCurrency ( ARM_Object* obj )	const	{ 

	Ct_Mod itm = itsModel.find(obj);

	if (itm == itsModel.end())
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the model.")

	M_Mkt	tmp = itm->second;

	V_Str	tmpCcy;

	tmpCcy.clear();
	for ( It_Mkt it = tmp.begin(); it != tmp.end() ; it ++ ){
		tmpCcy.push_back( it->first );
	}

	return tmpCcy;
}

ARM_Observator* ARM_Subject::GetModel ( const string & ccy ) const{

	Ct_Inv ito = itsModelInv.find( ccy );
	if ( ito != itsModelInv.end() ){
		return (ito -> second )[0];
	}
	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the model "<<ccy<<".")
	return NULL;
	
}

void ARM_Subject::Reset(){	
	
	It_Inv	itv;
	It_Mkt	itm;
	
	for( itm = itsCurrentMkt.begin(); itm != itsCurrentMkt.end();	itm++)
		if( itm->second ) { delete( itm->second); itm->second=NULL; }
	itsCurrentMkt.clear(); 

	for( itm = itsInitialMkt.begin(); itm != itsInitialMkt.end();	itm++)
		if( itm->second ) { delete( itm->second); itm->second=NULL; }
	itsInitialMkt.clear();

	for( itv = itsModelInv.begin(); itv != itsModelInv.end(); itv ++){ 
		for(int i=0; i< (itv-> second).size(); i++)
			Detach( (itv -> second )[i] );
	}

}

void ARM_Subject::Equalize( const string&	ccy	){	

	It_Mkt	  itm	=	itsCurrentMkt.find(ccy);
	if (itm == itsCurrentMkt.end())
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the ccy " << ccy << ".")
	itm ->	second	=	(itsInitialMkt.find(ccy)->second)->Clone();
	Notify( ccy );
} 

void ARM_Subject::Map( const string & ccy, const ARM_Object* obj){ 
	
	It_Mkt	itm	=	itsCurrentMkt.find(ccy);

	if ( itm == itsCurrentMkt.end() ){
		pair<string, ARM_Object*>	p( ccy , obj->Clone( ) );
		itsInitialMkt.insert( p	);
		}
	pair<string, ARM_Object*>	p( ccy , obj->Clone( ) );
	itsCurrentMkt.insert( p	);
}


ARM_Subject::ARM_Subject( const ARM_Subject & sub ){

	string		ccy;
	It_Mkt		itm;
	M_Mkt		tmpMkt;
	M_Mkt		sCurrentMkt = sub.itsCurrentMkt;

	ARM_Object*	tmp;
	M_Mkt		tmpMod;

	this->Reset();

	tmpMkt		=	sub.itsInitialMkt;
	for( itm	=	tmpMkt.begin(); itm != tmpMkt.end(); itm ++){
		ccy		=	itm->first;
		tmp		=	itm->second->Clone();

		this	->	itsInitialMkt.insert( pair<string, ARM_Object*>	( ccy , tmp ) );
	}

	tmpMkt		=	sub.itsCurrentMkt;
	for( itm	=	tmpMkt.begin(); itm != tmpMkt.end(); itm ++){
		ccy		=	itm->first;
		tmp		=	itm->second->Clone();

		this		->	itsCurrentMkt.insert( pair<string, ARM_Object*>	( ccy , tmp ) );
	}

	for( Ct_Mod it	=	sub.itsModel.begin()	; it  != sub.itsModel.end() ; it++ ){
		tmpMod.clear();
		tmpMkt		=	it->second;
		for ( itm	=	tmpMkt.begin()	; itm != tmpMkt.end()	; itm ++)
			tmpMod.insert( pair<string, ARM_Object*> (itm->first, itm -> second -> Clone() ) );

		this		->	itsModel.insert( pair<ARM_Object*, M_Mkt> (it->first -> Clone(), tmpMod ) );
	}

	itsCcyDim	= sub.itsCcyDim;

}

ARM_Object*	 ARM_Subject::Get	( const string & ccy) const {

 	Ct_Mkt	itm	=	itsCurrentMkt.find(ccy);

	if (itm == itsCurrentMkt.end())
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the ccy " << ccy << ".")

	return  itm	->	second;
}


void	ARM_Subject::Set	( const string & ccy, const ARM_Object* obj){

	It_Mkt	  itm	=	itsCurrentMkt.find(ccy); 

	if (itm == itsCurrentMkt.end())
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the ccy " << ccy << ".")

	itm ->	second	=	obj -> Clone() ;
	Notify( ccy );
}

void ARM_Subject::Notify( const string & ccy ){

	int				Ind;
	V_Obs			tmpV;
	M_Mkt			tmpM;
	It_Inv			it;
	It_Mkt			itm;

	ARM_Object*		tmp		=	Get(ccy);
	it = itsModelInv.find( ccy );
	if ( it != itsModelInv.end() ){
		tmpV = it->second;
		for ( int i=0 ; i< tmpV.size() ; i++){
			Ind		  =	tmpV[i]->GetCurrency().find(ccy)->second;
			From( tmpV[i]->GetObject() , tmp , Ind);
		}
	}
	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the ccy " << ccy << ".")
}

void ARM_Subject::Attach( ARM_Observator*  obs, string & str ){	
	
	string			ccy;
	M_Mkt			M_Obs;
	V_Obs			tmp;
	It_Inv			it;

	BuildModel(obs,str);
	ARM_Object*		obj		=	obs -> GetObject();
	ARM_Object*		tmpObj;

	M_Obs.clear();

	for ( int i=0; i< itsCcyDim ; i++){
		
		tmpObj	=	To( obj, i+1);
		ccy		=	GetCcy( tmpObj );
		pair<string, ARM_Object*>	p( ccy , tmpObj );
		M_Obs.insert( p );

		Map( ccy , tmpObj );

		it = itsModelInv.find ( ccy );
		if ( it != itsModelInv.end() )
			it->second.push_back( obs );
		else{
			tmp.clear();
			tmp.push_back( obs );
			pair<string, V_Obs>	t( ccy , tmp );	
			itsModelInv.insert (t);
		}
	}
	pair<ARM_Object*, M_Mkt>	q( obj, M_Obs );
	itsModel.insert( q ) ;
}

void ARM_Subject::Detach( ARM_Observator*  obs){	

	It_Obs  tmp;
	ARM_Object*		obj = obs->GetObject();

	It_Mod			ito = itsModel.find( obj );

	if (ito == itsModel.end())
		ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the observator.")

	( ito->second ).clear();
	itsModel.erase( ito );

	for( It_Inv it= itsModelInv.begin(); it !=itsModelInv.end(); it++){
		tmp = find( ( it-> second ).begin(), ( it-> second ).end(), obs);
		if(tmp == ( it-> second ).end())
		{
			ARMTHROW(ERR_INVALID_ARGUMENT, "Cannot find the subject.")
		}
		( it->second ).erase( tmp );
		if ( ( it->second) .size ==0) 
			itsModelInv.erase( it );
	}

}

V_Str	ARM_Subject::GetCurrency( )		const{
	V_Str tmp;

	for( Ct_Mkt it= itsCurrentMkt.begin(); it!= itsCurrentMkt.end(); it++)
		tmp.push_back(it->first);
	return tmp;
}


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM 1D Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


string		ARM_1D_Subject::GetCcy	( ARM_Object* obj)	const{
	
	return dynamic_cast<ARM_ZeroCurve*> ( To( obj, 1 ) ) -> GetCurrencyUnit() -> GetCcyName();
}

void	ARM_1D_Subject::Bump	( ARM_Object* curve, const ARM_Matrix  & shift, const int & isRelative){	
	
	ARM_CRV_TERMS		psMatu;
	ARM_ZeroCurve*		tmp		=	dynamic_cast< ARM_ZeroCurve* > ( curve );
	int					dim		=	tmp->GetMktData() -> itsMktValue -> GetSize();
	ARM_Vector*			epsilon	=	new ARM_Vector(dim,0.0 );
	double				rate;
	
	for( int i=0; i<dim;i++) {
		sprintf(psMatu[i], "%s", (const char*)( tmp->GetMktData()->itsMktTerms[i]));
		rate = tmp->GetMktData()->itsMktValue->Elt(i);

		if ( isRelative )
			epsilon->Elt(i)= rate* shift.Elt(i,0);
		else
			epsilon->Elt(i)= shift.Elt(i,0);
	}
	tmp->Copy( tmp->GenerateShiftCurve(psMatu,epsilon) );	
}

string		ARM_1D_Subject::toString		( const	string & ccy, const string  & indent, const string & nextIndent) const {
	CC_Ostringstream	os;
	int					dim;
	double				val;
	string				buc;
	ARM_ZeroCurve*		tmp = NULL;

	Ct_Mkt it = itsCurrentMkt.find(ccy);

	if ( it!=itsCurrentMkt.end() ){

		os << "\n"				<< indent	<<"\n";
		os << GetSubjectInfo()	<<" :"		<<"\n";

		tmp = dynamic_cast< ARM_ZeroCurve*> (it->second);
		dim = tmp->GetMktData()->itsMktValue->size(); 

		for ( int i=0 ; i< dim; i++){
			val = tmp->GetMktData()->itsMktValue->Elt(i);
			buc = (string)( tmp->GetMktData()->itsMktTerms[i]);
			os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(9)<<  buc;
			os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(9)<<CC_NS(std,setprecision)(3)<<val<<"\n";
		}
	}
	os << "\n";
	return os.str();
}
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM YC Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


ARM_Object*		ARM_YC_Subject::To	( ARM_Object * obj, int dim)	const{

	ARM_ZeroCurve*		tmpObj	=	NULL;

	switch ( dim ){

	case 1:
		if	( dynamic_cast< ARM_ZeroCurve* >		( obj ) ){

			if	( dynamic_cast< ARM_BasisCurve* >	( obj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( obj ) -> GetZeroCurve();
			else 
				tmpObj	=	dynamic_cast< ARM_ZeroCurve* >	( obj );
		}
	
		else if ( dynamic_cast< ARM_BSSmiledModel* >		( obj ) ){
			tmpObj	=	dynamic_cast<  ARM_BSSmiledModel* >	( obj ) -> GetZeroCurve();
			if ( dynamic_cast<  ARM_BasisCurve*  >	( tmpObj ) )
					tmpObj	=	dynamic_cast<  ARM_BasisCurve*  >	( tmpObj )	-> GetZeroCurve();
		}

		else if ( dynamic_cast< ARM_InfBSModel *>	(obj ) ) 
			tmpObj		= dynamic_cast< ARM_InfBSModel *> ( obj )->GetIRModel()-> GetYCModel() -> GetZeroCurve();

		else if ( dynamic_cast< ARM_BSModel*  >				( obj ) ){
			tmpObj	=	dynamic_cast<  ARM_BSModel*		>	( obj ) -> GetDividend();						// for fx
			if ( !tmpObj)
				tmpObj	=	dynamic_cast<  ARM_BSModel*   >		( obj ) -> GetYCModel() -> GetZeroCurve();  // for osw
			if ( dynamic_cast<  ARM_BasisCurve*  >	( tmpObj ) )
				tmpObj	=	dynamic_cast<  ARM_BasisCurve*  >	( tmpObj )	-> GetZeroCurve();
		}

		else if ( dynamic_cast< ARM_MixtureModel_Fx*	  >		( obj ) ){
			ARM_ModelParams*	tmp			=	dynamic_cast< ARM_MixtureModel_Fx* >		( obj ) ->GetModelParams();
			tmpObj	= &*dynamic_cast< ARM_ModelParams_Fx* >	( tmp	) -> GetForCurve();
			if	( dynamic_cast< ARM_BasisCurve* >	( tmpObj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( tmpObj ) -> GetZeroCurve();
		}

		break;

	case 2:
		if	( dynamic_cast< ARM_ZeroCurve* >		( obj ) ){
			if	( dynamic_cast< ARM_BasisCurve* >	( obj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( obj ) -> GetZeroCurve();
			else 
				tmpObj	=	dynamic_cast< ARM_ZeroCurve* >	( obj );
		}
		else if ( dynamic_cast< ARM_BSModel* >			( obj ) ){
			if( !dynamic_cast< ARM_BSSmiledModel* >			( obj ) ){
				tmpObj	=	dynamic_cast<  ARM_BSModel*   >		( obj ) -> GetYCModel() -> GetZeroCurve();
				if ( dynamic_cast<  ARM_BasisCurve*  >	( tmpObj ) ) 
					tmpObj	=	dynamic_cast<  ARM_BasisCurve*  >	( tmpObj )	-> GetZeroCurve();
			}
		}
		else if ( dynamic_cast< ARM_MixtureModel_Fx*	>		( obj ) ){
			ARM_ModelParams*	tmp			=	dynamic_cast< ARM_MixtureModel_Fx* >		( obj ) ->GetModelParams();	
			tmpObj	= &*dynamic_cast< ARM_ModelParams_Fx* >	( tmp	) -> GetDomCurve();
			if	( dynamic_cast< ARM_BasisCurve* >	( tmpObj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( tmpObj ) -> GetZeroCurve();
		}

		break;
	}
	return tmpObj;
}


// This function recalculates the curve (a zero curve or a basis curve)

void GenerateZeroCurve(ARM_Object* fromObj, ARM_Object* toObj)
{
	ARM_ZeroCurve* fromZc = dynamic_cast<ARM_ZeroCurve*>(fromObj);
	ARM_ZeroCurve* toZc = dynamic_cast<ARM_ZeroCurve*>(toObj);

	if (!fromZc)
		ARMTHROW(ERR_INVALID_ARGUMENT, "GenerateZeroCurve: From Object is not a zero curve." );

	if (!toZc)
		ARMTHROW(ERR_INVALID_ARGUMENT, "GenerateZeroCurve: To Object is not a zero curve.");

	ARM_BasisCurve* fromBc;
	if (fromBc = dynamic_cast<ARM_BasisCurve*>(fromZc))
		fromBc->GenerateZeroCurve(toZc);
	else
		fromZc->Copy(toZc);
}

void		ARM_YC_Subject::From			( ARM_Object * obs,	ARM_Object * obj, int  dim	){

	ARM_ZeroCurve*		cur		=	NULL;
	ARM_ModelParams*	tmp		=	NULL;

	if ( dynamic_cast< ARM_ZeroCurve *> ( obj ) ){
		ARM_ZeroCurve* cur = dynamic_cast< ARM_ZeroCurve *> ( obj );

		switch ( dim ){
		case 1:

			if( dynamic_cast< ARM_ZeroCurve *> ( obs ) )
				GenerateZeroCurve(obs,cur);
			
			else if ( dynamic_cast< ARM_BSSmiledModel *> ( obs )){
				ARM_BSSmiledModel* tmp	=	dynamic_cast< ARM_BSSmiledModel *>	( obs );
				tmp		->	SetZeroCurve( cur );
			}

			else if ( dynamic_cast< ARM_InfBSModel *>	(obs) ){
				dynamic_cast< ARM_InfBSModel *>	(obs) -> SetZeroCurve( cur );
				dynamic_cast< ARM_InfBSModel *>	(obs) -> GetIRModel ( )->SetZeroCurve( cur );
			}

			else if ( dynamic_cast< ARM_BSModel *> ( obs )){
				ARM_BSModel*	tmp			=	dynamic_cast< ARM_BSModel*>	( obs );

				if ( dynamic_cast< ARM_FXVolCurve*> ( tmp -> GetVolatility() ) ){
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*> ( tmp -> GetVolatility() );
					ARM_ZeroCurve* tmpDomCurve =	dynamic_cast< ARM_ZeroCurve*> ( tmp -> GetYCModel() -> GetZeroCurve()->Clone() );
					ARM_ZeroCurve*	tmpForCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp -> GetDividend()->Clone() );	
					GenerateZeroCurve(tmpForCurve,cur);
					tmp			->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
	
				else{
					ARM_ZeroCurve* tmpObj = dynamic_cast< ARM_ZeroCurve* >(tmp -> GetDividend());
					if (!tmpObj)
						tmpObj = dynamic_cast< ARM_ZeroCurve* >(tmp-> GetYCModel()->GetZeroCurve());
					ARM_ZeroCurve*	tmpForCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmpObj->Clone() );
					ARM_VolCurve*	volCurve	=	dynamic_cast< ARM_VolCurve*> ( tmp -> GetVolatility()->Clone() );
					GenerateZeroCurve(tmpForCurve,cur);
					
					if( dynamic_cast< ARM_BSModel *> ( obs )->GetSabrModel() ){
						int				modelType		= dynamic_cast< ARM_BSModel *> ( obs ) -> GetModelType();	
						int				spreadVolType	= dynamic_cast< ARM_BSModel *> ( obs ) -> GetSpreadVolType();

						ARM_VolCurve*			spreadLockCurve = dynamic_cast< ARM_BSModel *> ( obs ) -> GetSpreadLock() ;
						if(	spreadLockCurve )	spreadLockCurve = dynamic_cast< ARM_VolCurve*> ( spreadLockCurve -> Clone() );

						ARM_VolCurve*			capVolCurve		= dynamic_cast< ARM_BSModel *> ( obs ) -> GetCvxAdjVolatility() ;
						if( capVolCurve )		capVolCurve		= dynamic_cast< ARM_VolCurve*> ( capVolCurve -> Clone() );
						
						ARM_VolCurve*			correlationCurve= dynamic_cast< ARM_BSModel *> ( obs ) -> GetCorrelationMatrix();
						if( correlationCurve )	correlationCurve= dynamic_cast< ARM_VolCurve*> ( correlationCurve -> Clone() );

						ARM_VolCurve*			spreadVolCurve	= dynamic_cast< ARM_BSModel *> ( obs ) -> GetSpreadVolCurve();
						if( spreadVolCurve )	spreadVolCurve	= dynamic_cast< ARM_VolCurve*> ( spreadVolCurve -> Clone() );

						ARM_VolCurve*			cashVolCurve	= dynamic_cast< ARM_BSModel *> ( obs ) -> GetCashVolMatrix();
						if( cashVolCurve )		cashVolCurve	= dynamic_cast< ARM_VolCurve*> ( cashVolCurve -> Clone() );

						tmp	-> Set(	tmpForCurve, 		tmpForCurve,	spreadLockCurve, 	capVolCurve ,	volCurve, 		
									correlationCurve, 	cashVolCurve, 	spreadVolCurve,		modelType,		spreadVolType);

						ARM_BSSmiledModel* tmp	=	dynamic_cast< ARM_BSModel *> ( obs )->GetSabrModel();
						tmp		->	SetZeroCurve( cur );
						dynamic_cast< ARM_BSModel *> ( obs )->SetSabrModel( new ARM_BSSmiledModel( *tmp ) );
					}
					else
						tmp	->  Set(0, tmpForCurve, tmpForCurve, volCurve, tmp->GetType());
				}
			}
			
			else if ( dynamic_cast< ARM_MixtureModel_Fx*	>	( obs ) ){

				ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( obs ) -> GetModelParams();
				ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast<ARM_ModelParamsMixture_Fx*>( tmp );
				ARM_ZeroCurve*				tmpForCurve	=	dynamic_cast< ARM_ZeroCurve *> ( &*tmpFx -> GetForCurve() );
				GenerateZeroCurve(tmpForCurve,cur);
			}

			else
				throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET, "the object argument is not a ARM_ZeroCurve" );
			break;

		case 2:

			if ( dynamic_cast< ARM_BSModel *> ( obs )){
				ARM_BSModel*	tmp			=	dynamic_cast< ARM_BSModel*>	( obs );
				if ( dynamic_cast< ARM_BasisCurve*> ( tmp -> GetDividend() ) ) {
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*>( tmp -> GetVolatility() -> Clone() );
					ARM_BasisCurve* tmpDomCurve	=	dynamic_cast< ARM_BasisCurve*> ( tmp -> GetYCModel() -> GetZeroCurve()->Clone() );
					tmpDomCurve	->	GenerateZeroCurve( cur );
					ARM_BasisCurve*	tmpForCurve	=	dynamic_cast< ARM_BasisCurve*> ( tmp -> GetDividend() -> Clone() );	
					tmp		->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
				if ( dynamic_cast< ARM_FXVolCurve*> ( tmp -> GetVolatility() ) ){
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*>( tmp -> GetVolatility() -> Clone() );
					ARM_ZeroCurve* tmpDomCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp -> GetYCModel() -> GetZeroCurve()->Clone() );
					GenerateZeroCurve(tmpDomCurve,cur);
					ARM_ZeroCurve*	tmpForCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp -> GetDividend() -> Clone() );	
					tmp		->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
			}

			else if ( dynamic_cast< ARM_MixtureModel_Fx*	>	( obs ) ){
				ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( obs ) -> GetModelParams();
				ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast<ARM_ModelParamsMixture_Fx*>( tmp );
				ARM_ZeroCurve*				tmpDomCurve	=	dynamic_cast< ARM_ZeroCurve *> ( &*tmpFx -> GetDomCurve() );
			
				GenerateZeroCurve(tmpDomCurve,cur);
				tmpFx		->	SetDomCurve ( CreateClonedPtr<ARM_ZeroCurve>(tmpDomCurve) );
			}

			break;

		default:
			throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET, "the object argument is not a ARM_ZeroCurve" );
			break;
		}
	}
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM BS Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

ARM_Object*		ARM_BS_Subject::To	( ARM_Object * obj, int dim)	const{

	ARM_ZeroCurve*		tmpObj	=	NULL;
	ARM_ModelParams*	tmp		=	NULL;

	switch( dim ){
	case 1:

		if	( dynamic_cast< ARM_ZeroCurve* >				( obj ) ){

			if		( dynamic_cast< ARM_BasisCurve* >		( obj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( obj )->GetBasisCurve();
			
			else
				tmpObj = dynamic_cast< ARM_ZeroCurve* >(obj);
		}

		else if ( dynamic_cast< ARM_BSModel*  >				( obj	) ){
			if( ! dynamic_cast< ARM_BSSmiledModel *> ( obj ) ){

				tmpObj	=	dynamic_cast<  ARM_BSModel*		>	( obj	) -> GetDividend();
				if (dynamic_cast<  ARM_BasisCurve*  >	( tmpObj) )
					tmpObj	=	dynamic_cast<  ARM_BasisCurve*  >	( tmpObj) -> GetBasisCurve();
				if ( dynamic_cast< ARM_BasisCurve* >		( tmpObj ) )
					tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( tmpObj )->GetBasisCurve();
			}
		}

		else if ( dynamic_cast< ARM_MixtureModel_Fx*	>	( obj	) ){
			ARM_ModelParams*	tmp			=	dynamic_cast< ARM_MixtureModel_Fx* >		( obj ) ->GetModelParams();	
			tmpObj	= &*dynamic_cast< ARM_ModelParams_Fx* >	( tmp	) -> GetForCurve();
			if	( dynamic_cast< ARM_BasisCurve* >		( tmpObj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( tmpObj )->GetBasisCurve();
		}
		break;

	case 2:
		if	( dynamic_cast< ARM_ZeroCurve* >				( obj ) )
			if		( dynamic_cast< ARM_BasisCurve* >		( obj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( obj )->GetBasisCurve();
			else
				tmpObj = dynamic_cast< ARM_ZeroCurve* >(obj);
		else if ( dynamic_cast< ARM_BSModel* >				( obj	) ){
			if( ! dynamic_cast< ARM_BSSmiledModel *> ( obj ) ){
				tmpObj	=	dynamic_cast<  ARM_BSModel*   >		( obj	) -> GetYCModel() -> GetZeroCurve();
				if ( dynamic_cast<  ARM_BasisCurve*  >	( tmpObj) ) 
					tmpObj	=	dynamic_cast<  ARM_BasisCurve*  >	( tmpObj) -> GetBasisCurve();
				if	( dynamic_cast< ARM_BasisCurve* >		( tmpObj ) )
					tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( tmpObj )->GetBasisCurve();
			}
		}

		else if ( dynamic_cast< ARM_MixtureModel_Fx*	>	( obj	) ){
			ARM_ModelParams*	tmp			=	dynamic_cast< ARM_MixtureModel_Fx* >		( obj ) ->GetModelParams();	
			tmpObj	= &*dynamic_cast< ARM_ModelParams_Fx* >	( tmp	) -> GetDomCurve();
			if	( dynamic_cast< ARM_BasisCurve* >		( tmpObj ) )
				tmpObj	=	dynamic_cast< ARM_BasisCurve* > ( tmpObj )->GetBasisCurve();
		}
		break;
	}
	return tmpObj;
}


// This function recalculates the curve (a zero curve or a basis curve)

void GenerateBasisCurve(ARM_Object* fromObj, ARM_Object* toObj)
{
	ARM_ZeroCurve* fromZc = dynamic_cast<ARM_ZeroCurve*>(fromObj);
	ARM_ZeroCurve* toZc = dynamic_cast<ARM_ZeroCurve*>(toObj);

	if (!fromZc)
		ARMTHROW(ERR_INVALID_ARGUMENT, "GenerateZeroCurve: From Object is not a zero curve." );

	if (!toZc)
		ARMTHROW(ERR_INVALID_ARGUMENT, "GenerateZeroCurve: To Object is not a zero curve.");

	ARM_BasisCurve* fromBc;
	if (fromBc = dynamic_cast<ARM_BasisCurve*>(fromZc))
		fromBc->GenerateBasisCurve(toZc);
	else
		fromZc->Copy(toZc);
}

void		ARM_BS_Subject::From			( ARM_Object * obs,	ARM_Object * obj,  int dim	){

	ARM_ZeroCurve*		cur		=	NULL;
	ARM_BasisCurve*		tmpObj	=	NULL;
	ARM_ModelParams*	tmp		=	NULL;

	if ( dynamic_cast< ARM_ZeroCurve *> ( obj ) ){
		ARM_ZeroCurve* cur = dynamic_cast< ARM_ZeroCurve *> ( obj );

		switch( dim){
		case 1:
			if (  dynamic_cast< ARM_ZeroCurve *> ( obs )){
				GenerateBasisCurve(obs,cur);
			}
			
			else if ( dynamic_cast< ARM_BSModel *> ( obs )){
				ARM_BSModel*	tmp			=	dynamic_cast< ARM_BSModel*>	( obs );
				if (dynamic_cast< ARM_BasisCurve*> (tmp-> GetDividend())){
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*> ( tmp-> GetVolatility()->Clone() );
					ARM_BasisCurve* tmpDomCurve =	dynamic_cast< ARM_BasisCurve*> ( tmp-> GetYCModel() -> GetZeroCurve()->Clone() );
					ARM_BasisCurve*	tmpForCurve	=	dynamic_cast< ARM_BasisCurve*> ( tmp-> GetDividend()->Clone() );	
					tmpForCurve	->	GenerateBasisCurve( cur );
					tmp			->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
				if (dynamic_cast< ARM_BasisCurve*> (tmp-> GetDividend()))
				{
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*> ( tmp-> GetVolatility()->Clone() );
					ARM_ZeroCurve* tmpDomCurve =	dynamic_cast< ARM_ZeroCurve*> ( tmp-> GetYCModel() -> GetZeroCurve()->Clone() );
					ARM_ZeroCurve*	tmpForCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp-> GetDividend()->Clone() );	
					GenerateBasisCurve(tmpForCurve,cur);
					tmp			->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
			}

			else if ( dynamic_cast< ARM_MixtureModel_Fx*	>	( obs ) ){
				ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( obs ) -> GetModelParams();
				ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast<ARM_ModelParamsMixture_Fx*>( tmp );
				ARM_ZeroCurve*				tmpForCurve	=	dynamic_cast< ARM_ZeroCurve *> ( &*tmpFx -> GetForCurve() );		
			
				GenerateBasisCurve(tmpForCurve,cur);
			}

			else
				throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET, "the object argument is not a ARM_ZeroCurve" );
			break;

		case 2:
			
			if ( dynamic_cast< ARM_BSModel *> ( obs )){
				ARM_BSModel*	tmp			=	dynamic_cast< ARM_BSModel*>	( obs );
				if (dynamic_cast< ARM_BasisCurve*> (tmp-> GetDividend())){
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*>( tmp -> GetVolatility() -> Clone() );
					ARM_BasisCurve* tmpDomCurve	=	dynamic_cast< ARM_BasisCurve*> ( tmp -> GetYCModel() -> GetZeroCurve()->Clone() );
					tmpDomCurve	->	GenerateBasisCurve( cur );
					ARM_BasisCurve*	tmpForCurve	=	dynamic_cast< ARM_BasisCurve*> ( tmp -> GetDividend() -> Clone() );	
					tmp		->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
				if (dynamic_cast< ARM_BasisCurve*> (tmp-> GetDividend()))
				{
					ARM_BSModel*	tmp			=	dynamic_cast< ARM_BSModel*>	( obs );
					ARM_FXVolCurve*	tmpVol		=	dynamic_cast< ARM_FXVolCurve*>( tmp -> GetVolatility() -> Clone() );
					ARM_ZeroCurve* tmpDomCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp -> GetYCModel() -> GetZeroCurve()->Clone() );
					GenerateBasisCurve(tmpDomCurve,cur);
					ARM_ZeroCurve*	tmpForCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp -> GetDividend() -> Clone() );	
					tmp		->  Set( tmp->GetSpot(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );
				}
			}

			else if ( dynamic_cast< ARM_MixtureModel_Fx*	>	( obs ) ){
				ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( obs ) -> GetModelParams();
				ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast<ARM_ModelParamsMixture_Fx*>( tmp );
				ARM_ZeroCurve*				tmpDomCurve	=	dynamic_cast< ARM_ZeroCurve *> ( &*tmpFx -> GetDomCurve() );
			
				GenerateBasisCurve(tmpDomCurve,cur);
				tmpFx		->	SetDomCurve ( CreateClonedPtr<ARM_ZeroCurve>(tmpDomCurve) );
				}

			break;

		default:
			ARMTHROW(ERR_INVALID_ARGUMENT, "the object argument is not a ARM_ZeroCurve" );
			break;

		}
	}
}


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						ARM INF Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


ARM_Object*		ARM_INF_Subject::To	( ARM_Object * obj, int dim)	const{

	ARM_ZeroCurve*	tmp = NULL;
	ARM_InfCurv*	cur = NULL;
	vector<string>  str;

	ARM_CRV_TERMS Terms;


	if(  dynamic_cast< ARM_ZeroCurve* >	( obj ) )
		tmp = dynamic_cast< ARM_ZeroCurve* >	( obj );


	else if ( dynamic_cast< ARM_InfBSModel*> ( obj ) )
		cur = dynamic_cast< ARM_InfBSModel*> ( obj ) -> GetInfFwdCurv();

	else if ( dynamic_cast< ARM_InfCurv* >	( obj )  )
		cur = dynamic_cast< ARM_InfCurv* >	( obj );

	if( cur) {
		str = cur->GetMktTerms();

		for( int k = 0 ; k< cur->GetMktValues().size(); k++){
			double	buc = (double) atoi(str[k].c_str());
			if( buc>2000000)
				sprintf(Terms[k], "%s", ARM_Date(buc).GetStrDate() );
			else
				sprintf(Terms[k], "%s", str[k].c_str() );
		}

		tmp =  new ARM_ZeroCurve();
		tmp -> SetAsOfDate( cur->GetAsOf() );
		tmp -> SetCurrencyUnit( new ARM_Currency( cur->GetCalendar().c_str() ));
		tmp -> SetMktData(	Terms,	new ARM_Vector(cur->GetMktValues() ), 0, 0,	0,	0 );
	}

	return tmp;

}

void	ARM_INF_Subject::From( ARM_Object * obs,	ARM_Object * obj,  int dim	){

	ARM_ZeroCurve*		tmp		=	NULL;

	if ( dynamic_cast< ARM_ZeroCurve *> ( obj ) ){
		tmp = dynamic_cast< ARM_ZeroCurve *> ( obj );
		
		
		if (  dynamic_cast< ARM_InfCurv *> ( obs ) ){
			ARM_InfCurv*		cur = dynamic_cast< ARM_InfCurv *> ( obs );
			ARM_Vector*			vec = tmp->GetMktData()->itsMktValue;
			vector<double>		mkt(vec->size());
			for ( int i =0 ; i< mkt.size(); i++) mkt[i] = vec->Elt(i);
			cur->SetCurve( cur -> GetMktTerms(), mkt);
		}
		else if( dynamic_cast< ARM_InfBSModel *> ( obs ) ){
			ARM_InfBSModel *	mod =  dynamic_cast< ARM_InfBSModel *> ( obs );
			ARM_InfCurv*		cur = mod->GetInfFwdCurv();
			ARM_Vector*			vec = tmp->GetMktData()->itsMktValue;
			vector<double>		mkt(vec->size());
			for ( int i =0 ; i< mkt.size(); i++) mkt[i] = vec->Elt(i);
			cur->SetCurve( cur -> GetMktTerms(), mkt);
			mod->SetInfFwdCurv(dynamic_cast< ARM_InfCurv* > (cur->Clone()) );
		}
	}
	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "the object argument is not a ARM_InfZeroCurve" );

}

void	ARM_INF_Subject::Bump( ARM_Object* curve, const ARM_Matrix  & shift, const int & isRelative){	
	
	ARM_ZeroCurve*		tmp		=	dynamic_cast< ARM_ZeroCurve* > ( curve );
	ARM_Vector*			rate	=	tmp -> GetMktData() -> itsMktValue;
	int					dim		=	rate-> size();

	for( int i=0; i<dim; i++) {
		if(rate->Elt(i) >50)
			shift.Elt(i,0) *=100;
	}

	ARM_CRV_TERMS		mktTerms;
	ARM_Vector*			mktValues =	new ARM_Vector(dim);
	
	for( i=0; i<dim;i++) {
		sprintf(mktTerms[i], "%s", (const char*)( tmp->GetMktData()->itsMktTerms[i]));
		mktValues ->Elt(i) = tmp->GetMktData()->itsMktValue->Elt(i);

		if ( isRelative )
			mktValues ->Elt(i) += mktValues ->Elt(i) * shift.Elt(i,0);
		else
			mktValues ->Elt(i) += shift.Elt(i,0);
	}

	tmp -> SetMktData( mktTerms, mktValues,0,0,0,0);           
}
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						ARM FX Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

 
ARM_Object*		ARM_FX_Subject::To	( ARM_Object * obj, int)	const{

	ARM_Currency*			tmpMainCcy	=	NULL;   
	ARM_Currency*			tmpMoneyCcy	=	NULL;
	ARM_ModelParams*		tmpFx		=	NULL;
	ARM_ZeroCurve*			tmpBs		=	NULL;

	double				tmpSpot;				

	ARM_Forex*		tmpObj	=	NULL;

	if		( dynamic_cast< ARM_Forex* >	( obj ) )
		tmpObj	=	dynamic_cast<  ARM_Forex* > ( obj );

	else if ( dynamic_cast< ARM_BSModel*>	( obj ) ){

		tmpSpot		=	dynamic_cast<  ARM_BSModel* > ( obj )-> GetSpot() ;
		tmpMainCcy	=	dynamic_cast<  ARM_BSModel* > ( obj )-> GetDividend() -> GetCurrencyUnit();
		tmpMoneyCcy =	dynamic_cast<  ARM_BSModel* > ( obj )-> GetYCModel()  -> GetZeroCurve() -> GetCurrencyUnit();

		tmpObj		=	new ARM_Forex(tmpMainCcy, tmpMoneyCcy, tmpSpot);
	}
	
	else if ( dynamic_cast< ARM_MixtureModel_Fx*	  >		( obj ) ){
		ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( obj ) -> GetModelParams();
		ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast<ARM_ModelParamsMixture_Fx*>( tmp );
		ARM_Currency*				tmpMoneyCcy	=	dynamic_cast< ARM_ZeroCurve *> ( &*tmpFx -> GetDomCurve() )-> GetCurrencyUnit();
		ARM_Currency*				tmpMainCcy	=	dynamic_cast< ARM_ZeroCurve *> ( &*tmpFx -> GetForCurve() )-> GetCurrencyUnit();
		
		tmpSpot		=	tmpFx->GetSpot();
		tmpObj		=	new ARM_Forex(tmpMainCcy, tmpMoneyCcy, tmpSpot);		
		}

	return tmpObj;
}

void		ARM_FX_Subject::From	( ARM_Object * obs,	ARM_Object * obj, int	){

	ARM_Forex*				cur		=	NULL;
	ARM_Forex*				tmpObj	=	NULL;

	if ( dynamic_cast< ARM_Forex *> ( obj ) ){
		ARM_Forex* cur = dynamic_cast< ARM_Forex *> ( obj );

		if ( dynamic_cast< ARM_Forex *> ( obs )){
			ARM_Forex* tmp	=	dynamic_cast< ARM_Forex *>	( obs );
			(*tmp) = (*cur);
		}
		
		else if ( dynamic_cast< ARM_BSModel*> ( obs )){
			ARM_BSModel*   tmp			=	dynamic_cast< ARM_BSModel*>	( obs );
			ARM_FXVolCurve*  tmpVol		=	dynamic_cast< ARM_FXVolCurve*> (tmp->GetVolatility() );
			ARM_ZeroCurve* tmpDomCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp-> GetYCModel() -> GetZeroCurve()->Clone() );
			ARM_ZeroCurve* tmpForCurve	=	dynamic_cast< ARM_ZeroCurve*> ( tmp-> GetDividend()->Clone() );
			tmpVol	->	BumpFxVol(0.0, 0, 0, cur -> GetMarketPrice(),tmpDomCurve,tmpForCurve );
			tmpVol  =   dynamic_cast< ARM_FXVolCurve*> ( tmpVol->Clone());
			tmp		->  Set( cur -> GetMarketPrice(), tmpForCurve, tmpDomCurve, tmpVol , tmp->GetType() );	
		}

		else if ( dynamic_cast< ARM_MixtureModel_Fx*	  >		( obs ) ){
			ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( obs ) -> GetModelParams();
			ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast<ARM_ModelParamsMixture_Fx*>( tmp );

			tmpFx 	->	SetSpot( cur->  GetMarketPrice() );
		}

		else
			ARMTHROW(ERR_CONDITION_NOT_MEET, "the object argument is not a ARM_ZeroCurve" );
	}

}

string		ARM_FX_Subject::GetCcy	( ARM_Object* obj )	const{
	string tmp;
	tmp = dynamic_cast<ARM_Forex*> ( To( obj ) ) -> GetMainCurrency() -> 	GetCcyName();
	tmp +="/";
	tmp +=dynamic_cast<ARM_Forex*> ( To( obj ) ) -> GetMoneyCurrency() -> 	GetCcyName();
	return tmp;
}

void	ARM_FX_Subject::Bump( ARM_Object* curve,  const ARM_Matrix & shift, const int & isRelative){
	
	double		val	= shift.Elt(0,0);
	ARM_Forex * tmp = dynamic_cast<ARM_Forex *> ( curve); 
	if ( isRelative ) 
		tmp->SetMarketPrice( ((ARM_Security *) curve )->GetMarketPrice() * (1+val) );
	else
		tmp->SetMarketPrice( ((ARM_Security *) curve )->GetMarketPrice() + val );
}

string		ARM_FX_Subject::toString		( const	string & ccy , const string  & indent, const string & nextIndent) const{
	CC_Ostringstream	os;
	Ct_Mkt it = itsCurrentMkt.find(ccy);

	if ( it!=itsCurrentMkt.end() ){

		os << "\n"				<< indent	<<"\n";
		os << GetSubjectInfo()	<<" :"		<<"\n";

		double tmp = dynamic_cast< ARM_Forex*> (it->second)->GetMarketPrice();

		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(9)<<"Spot";
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(9)<< CC_NS(std,setprecision)(3)<<tmp;
		os<<"\n";
	}
	os << "\n";
	return os.str();
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						ARM DAT Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
 
ARM_Object*		ARM_DAT_Subject::To	( ARM_Object * obj, int)	const{

	ARM_Date*	tmpObj = NULL;

	if		( dynamic_cast< ARM_Date* >	( obj ) )
		tmpObj	=	dynamic_cast<  ARM_Date* > ( obj );

	else if		( dynamic_cast< ARM_ZeroCurve* >	( obj ) )
		tmpObj	=	new  ARM_Date( dynamic_cast< ARM_ZeroCurve* >( obj ) ->GetAsOfDate( ) );

	else if		( dynamic_cast< ARM_VolCurve* >	( obj ) )
		tmpObj	=	new  ARM_Date( dynamic_cast< ARM_VolCurve* >( obj ) ->GetAsOfDate( ) );  

	else if		( dynamic_cast< ARM_InfCurv* >	( obj ) )
		tmpObj	=	new  ARM_Date( dynamic_cast< ARM_InfCurv* >( obj ) ->GetAsOf( ) );

	else if ( dynamic_cast< ARM_Model* >	( obj ) )
		tmpObj	=	new  ARM_Date( dynamic_cast< ARM_Model* >( obj ) ->GetAsOfDate( ) );

	else if ( dynamic_cast< ARM_PricingModel* >	( obj ) )
		tmpObj	=	new  ARM_Date( dynamic_cast< ARM_PricingModel* >( obj ) ->GetAsOfDate( ) );

	return	tmpObj;
}

void		ARM_DAT_Subject::From	( ARM_Object * obs,	ARM_Object * obj, int	){

	ARM_Date*			tmpObj	=	NULL;
	int					Index	=	1;

	if ( dynamic_cast< ARM_Date *> ( obj ) ){
		tmpObj = dynamic_cast< ARM_Date *> ( obj );
		if ( dynamic_cast< ARM_Date *> (obs)  )
			obs = tmpObj;

		else{
			if( dynamic_cast<ARM_ZeroCurve* > ( ARM_INF_Subject::To(obs,1)  ) ){
				ARM_ZeroCurve*		cur		=	dynamic_cast<ARM_ZeroCurve* > ( ARM_INF_Subject::To(obs,1)->Clone() );
				cur	->SetAsOfDate(*tmpObj);
				if( cur->GetMktData() ) {
					ARM_INF_Subject::Bump(cur, ARM_Matrix(cur->GetMktData() -> itsMktValue -> GetSize(),1,0.0), true);
					ARM_INF_Subject::From(obs,cur->Clone(),1);
				}
			}
			for( int i = 0; i<2;i++){
				if( dynamic_cast<ARM_ZeroCurve* > ( ARM_YC_Subject::To(obs,i+1)  ) ){
					ARM_ZeroCurve*		cur		=	dynamic_cast<ARM_ZeroCurve* > ( ARM_YC_Subject::To(obs,i+1)->Clone() );
					cur	->SetAsOfDate(*tmpObj);
					if( cur->GetMktData()) {
						ARM_YC_Subject::Bump(cur, ARM_Matrix(cur->GetMktData() -> itsMktValue -> GetSize(),1,0.0), true);
						ARM_YC_Subject::From(obs,cur->Clone(),i+1);
					}
				}
				if( dynamic_cast<ARM_ZeroCurve* > ( ARM_BS_Subject::To(obs,i+1)  ) ){
					ARM_ZeroCurve*		cur		=	dynamic_cast<ARM_ZeroCurve* > ( ARM_BS_Subject::To(obs,i+1)->Clone() );
					cur	->SetAsOfDate(*tmpObj);
					if( cur->GetMktData()) {
						ARM_YC_Subject::Bump(cur, ARM_Matrix(cur->GetMktData() -> itsMktValue -> GetSize(),1,0.0), true);
						ARM_BS_Subject::From(obs,cur->Clone(),i+1);
					}
				}
			}
		}
	}
	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "The method cannot integrate AsofDate" );
}




string		ARM_DAT_Subject::GetCcy	( ARM_Object* obj )	const{
	return string("DATE");
}

void	ARM_DAT_Subject::Bump( ARM_Object* curve,  const ARM_Matrix & shift, const int &){
	dynamic_cast< ARM_Date *> (curve)->SetJulian( shift.Elt(0,0) ); 
}

string		ARM_DAT_Subject::toString		( const	string & ccy , const string  & indent, const string & nextIndent) const{
	CC_Ostringstream	os;
	char d[20];

	Ct_Mkt it=itsCurrentMkt.begin(); 
	ARM_Date* tmp = dynamic_cast<ARM_Date*> ( it->second );
	tmp->JulianToStrDate(d);

	os<< "\n"				<< indent	<<"\n";
	os<< GetSubjectInfo()	<<" : ";
	os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(9)<<(string) d;
	os<< "\n\n";

	return os.str();
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM 2D Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template<int n>
V_Str	ARM_2D_Subject<n>::GetLabel		( const string & ccy, const string & lab ) const{
		
	V_Str			tmp;
	ARM_VolCurve*	mod			= GetVol ( dynamic_cast<ARM_Object*> ( GetModel( ccy )->GetObject() ) );
	int				NbPlots;			

	if		( lab == "X"){
		NbPlots	= mod->GetVolatilities()->GetNumLines();
		for (int i = 0; i < NbPlots ; i++)
			tmp.push_back(	mod->itsYearTermsX[i] );

	}
	
	else if ( lab == "Y" ){
		NbPlots	= mod->GetVolatilities()->GetNumCols();
		for (int i = 0; i < NbPlots ; i++)
			tmp.push_back(	mod->itsYearTermsY[i] );
		
	}

	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "string must be X or Y" );
	return tmp;
}

template<int n>
ARM_Object*		ARM_2D_Subject<n>::To	( ARM_Object * obj, int)	const{

	ARM_VolCurve*		tmpObj	=	NULL;

	if		( dynamic_cast< ARM_VolCurve* >	( obj ) )
		tmpObj	=	dynamic_cast< ARM_VolCurve* > ( obj );

	else if ( dynamic_cast	< ARM_Object * >	( obj ) ){
			tmpObj		=	GetVol ( dynamic_cast	< ARM_Object * > ( obj ) );
	}
	
	return tmpObj;
}

template<int n>
void	ARM_2D_Subject<n>::From	( ARM_Object * obs,	ARM_Object * obj, int	){

	ARM_VolCurve*		vol		= NULL;
	ARM_Object*			tmp		= NULL;
	ARM_PricingModel*	tmpFx	= NULL;

	if ( dynamic_cast< ARM_VolCurve *> ( obj ) ){
		vol = dynamic_cast< ARM_VolCurve *> ( obj );

		if ( dynamic_cast< ARM_VolCurve *> ( obs) )	
			obs -> Copy( vol );

		else if ( dynamic_cast< ARM_Object *> ( obs) ){
			tmp = dynamic_cast< ARM_Object *> ( obs );
			SetVol( tmp, vol );
		}
		
		else
			ARMTHROW(ERR_INVALID_ARGUMENT,"the object argument is not a ARM_ZeroCurve" );
		}
}

template<int n>
string		ARM_2D_Subject<n>::GetCcy	( ARM_Object* obj )	const{

		return dynamic_cast<ARM_VolCurve*> ( To( obj ) ) -> GetCurrency() -> GetCcyName();
}

template< int n>
ARM_VolCurve*	ARM_2D_Subject<n>::GetVol	( ARM_Object *	mod		) const{

	ARM_VolCurve*		tmp		= NULL;
	ARM_BSSmiledModel *	tmpMod	= NULL;
	ARM_BSModel *		tmpBs	= NULL;
	bool				tmpSabr	= false;

	if( dynamic_cast< ARM_VolCurve *>	(mod ) )
		tmp = dynamic_cast< ARM_VolCurve *>	(mod );
	
	else if(  dynamic_cast< ARM_BSModel *>	(mod ) ){
		tmpBs	=	dynamic_cast< ARM_BSModel *>	(mod );
		tmp		=	tmpBs -> GetVolatility();
		if ( dynamic_cast<ARM_VolCube*> (tmp) ){

			if( dynamic_cast<ARM_VolCube*> (tmp)->GetATMVol() )
				tmp =	dynamic_cast<ARM_VolCube*> (tmp)->GetATMVol();

			else if ( dynamic_cast< ARM_InfBSModel *>	(mod ) ) {
				tmp = dynamic_cast< ARM_InfBSModel *>	(mod )->GetIRModel()-> GetVolatility();
				tmp = dynamic_cast<ARM_VolCube*> (tmp)->GetATMVol();
			}

			else if ( dynamic_cast< ARM_BSSmiledModel *>( tmpBs->GetSabrModel()) ){
				tmpMod	=	tmpBs->GetSabrModel();
				tmpSabr	=	true;
			}
		}
		else if(  dynamic_cast< ARM_BSSmiledModel *>	(mod ) ){
			tmpMod	=	dynamic_cast< ARM_BSSmiledModel *>	(mod );
			tmpSabr	=	true;
		}
		if( tmpSabr ){
			switch(n){
				case	ATM:
					tmp = tmpMod -> GetVolatility();
					break;
				case	RHO:
					tmp = tmpMod -> GetRho();
					break;
				case	NU:
					tmp = tmpMod -> GetNu();
					break;
				case	BETA:
					tmp = tmpMod -> GetBeta();
					break;
				default:
					ARMTHROW(ERR_INVALID_ARGUMENT, "The SABR volatility surfaces are missing." );
					break;
			}
		}
	}
	return  tmp;
}


template< int n> 
ARM_VolCurve*	ARM_CAP_Subject<n>::GetVol	( ARM_Object * mod	)	const{

	ARM_VolCurve*		tmp		= NULL;
	if ( dynamic_cast< ARM_InfBSModel *>	(mod ) ){
		tmp = dynamic_cast< ARM_InfBSModel *>	(mod )->GetIRModel()-> GetVolatility();
		tmp = dynamic_cast<ARM_VolCube*> (tmp)->GetATMVol();
	}
	else
		tmp = ARM_2D_Subject<n>::GetVol(mod);

	return tmp;
}

template< int n> 
ARM_VolCurve*	ARM_OSW_Subject<n>::GetVol	( ARM_Object * mod	)	const{

	ARM_VolCurve*		tmp		= NULL;
	if ( dynamic_cast< ARM_InfBSModel *>	(mod ) ){
		tmp = dynamic_cast< ARM_InfBSModel *>	(mod )->GetIRSwoptVolCurve();
	}
	else
		tmp = ARM_2D_Subject<n>::GetVol(mod);

	return tmp;
}

template< int n>
void	ARM_2D_Subject<n>::SetVol	(ARM_Object *	mod, ARM_VolCurve *	vol	){

	ARM_VolCurve*		tmp		=	NULL;
	ARM_BSSmiledModel *	tmpMod	=	NULL;
	ARM_BSModel *		tmpBs	=	NULL;
	
	bool				tmpSabr = false;

	if( dynamic_cast<ARM_VolCurve *> (mod) ){
		mod->Copy(vol);
	}
	else if(  dynamic_cast< ARM_BSModel *>	( mod ) ){
		tmpBs	=	dynamic_cast< ARM_BSModel *>	( mod );
		tmp		=	tmpBs -> GetVolatility();

		if ( dynamic_cast< ARM_VolCube *>	( tmp ) ){
			tmpMod	= tmpBs->GetSabrModel();
			if( tmpMod )
				tmpSabr	= true;
			else
				tmp = vol;
		}
		
		else if(  dynamic_cast< ARM_BSSmiledModel *>	(mod ) ){
			tmpMod = dynamic_cast< ARM_BSSmiledModel *>	(mod );
			tmpSabr= true;
		}
		else if ( dynamic_cast< ARM_BSModel* >  (mod) ){
			tmp->Copy(vol);
		}
		if( tmpSabr ){
			switch(n){
				case	ATM:
					tmpMod -> Set_XXX_Volatility( vol );
					break;

				case	RHO:
					tmpMod -> Set_XXX_Rho( vol );
					break;

				case	NU:
					tmpMod -> Set_XXX_Nu( vol );
					break;

				case	BETA:
					tmpMod -> Set_XXX_Beta( vol );
					break;

				default:
					ARMTHROW(ERR_INVALID_ARGUMENT, "The SABR volatility surfaces are missing." );
					break;
			}
		}
	}
}

template< int n> 
void	ARM_CAP_Subject<n>::SetVol(	ARM_Object *	mod, ARM_VolCurve *	vol) {

	ARM_VolCurve*		tmp		= NULL;
	if ( dynamic_cast< ARM_InfBSModel *>	(mod ) ){

		ARM_InfBSModel*	infMod  = dynamic_cast< ARM_InfBSModel *>	(mod );
		ARM_BSModel*	irMod	= infMod->GetIRModel();
		tmp = irMod -> GetVolatility();

		dynamic_cast<ARM_VolCube*> (tmp)->SetATMVol(vol);
		ARM_ZeroCurve* zc= dynamic_cast<ARM_ZeroCurve*> (irMod-> GetYCModel() -> GetZeroCurve()->Clone());
		irMod= new   ARM_BSModel( zc, tmp, K_INDEX_RATE);
		infMod->SetIRModel(irMod);
	}
	else
		ARM_2D_Subject<n>::SetVol(mod ,vol );
}

template< int n> 
void	ARM_OSW_Subject<n>::SetVol(	ARM_Object*	mod, ARM_VolCurve *	vol) {

	ARM_VolCurve*		tmp		= NULL;
	if ( dynamic_cast< ARM_InfBSModel *>	(mod ) ){
		tmp = dynamic_cast< ARM_InfBSModel *>	(mod )->GetIRSwoptVolCurve();
		tmp = dynamic_cast<ARM_VolCube*> (vol->Clone() );
		dynamic_cast< ARM_InfBSModel *>	(mod )->SetIRSwoptVolCurve( tmp );
	}
	else
		ARM_2D_Subject<n>::SetVol(mod ,vol );
}


template< int n>
string		ARM_2D_Subject<n>::toString		( const	string & ccy, const string  & indent, const string & nextIndent) const {
	CC_Ostringstream	os;
	int					dimX;
	int					dimY;
	V_Str				labX;
	V_Str				labY;

	ARM_Matrix*			tmp = NULL;

	Ct_Mkt it = itsCurrentMkt.find(ccy);

	if ( it!=itsCurrentMkt.end() ){

		os << "\n"				<< indent	<<"\n";
		os << GetSubjectInfo()	<<" :"		<<"\n";

		tmp = dynamic_cast< ARM_VolCurve*> (To(it->second) )->GetVolatilities();
		labX= GetLabel(ccy,"X");
		labY= GetLabel(ccy,"Y");
		dimX = labX.size(); 
		dimY = labY.size(); 
		
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(12)<<"";
		for ( int j=0 ; j< dimY; j++)
			os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(12)<<labY[j];

		os << "\n";
		for ( int i =0; i<dimX; i++){
			os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(12)<<labX[i];
			for ( j=0 ; j< dimY; j++)
				os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(12)<< CC_NS(std,setprecision)(3)<<tmp->Elt(i,j);
			os<<"\n";
		}
	}
	os << "\n";
	return os.str();
} 


			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM VFX Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


template< int n> 
string		ARM_VFX_Subject<n>::GetCcy	( ARM_Object* obj	)	const{

	string tmp;
	tmp= dynamic_cast<ARM_FXVolCurve*> ( To( obj ) ) -> GetForCcy() -> GetCcyName();
	tmp += "/";
	tmp += dynamic_cast<ARM_FXVolCurve*> ( To( obj ) ) -> GetDomCcy() -> GetCcyName();
	return tmp;
}


template< int n> 
V_Str	ARM_VFX_Subject<n>::GetLabel		( const string & ccy, const string & lab ) const{
	
	int				NbPlots;
	V_Str			tmp;
	ARM_FXVolCurve*	mod	 = dynamic_cast<ARM_FXVolCurve*> (GetVol ( dynamic_cast<ARM_Object*> ( GetModel( ccy )->GetObject() ) ) );
	

	if		( lab == "X"){
		NbPlots = 	mod -> GetOptionsMatus().GetSize();
		for (int i = 0; i < NbPlots ; i++)
			tmp.push_back(	mod->itsYearTermsX[i] );
	}

	if	( lab == "Y" ){
		switch(n){
			case	PIV:
				tmp.push_back( mod->itsYearTermsY[0] );
				break;
			case	RR:{
				for (int i = 0; i < 2 ; i++)
					tmp.push_back(	mod->itsYearTermsY[1+i] );
				break;
				}
			case	STR:{
				for (int i = 0; i < 2 ; i++)
					tmp.push_back(	mod->itsYearTermsY[3+i] );
				break;
			}
		}	
	}

	return tmp;	
}	
	
template< int n > 
ARM_VolCurve*	ARM_VFX_Subject<n>::GetVol	( ARM_Object *	mod		) const{

	ARM_Matrix*			Mtmp	= NULL;
	ARM_Vector*			Vtmp	= NULL;

	ARM_FXVolCurve*		tmpVol	= NULL;

	tmpVol = dynamic_cast<ARM_FXVolCurve*> (dynamic_cast<ARM_BSModel*> (mod) -> GetVolatility() );

	switch( n ){

		case	PIV:
			Vtmp = const_cast<ARM_Vector* > (tmpVol -> GetPivotVols());
			Mtmp= new ARM_Matrix (*Vtmp);
			break;

		case	RR:
			Mtmp = new ARM_Matrix (tmpVol -> GetRR());
			break;

		case	STR:
			Mtmp = new ARM_Matrix (tmpVol -> GetSTR());
			break;
		}
	
	tmpVol->SetVolatilities( Mtmp);
	return tmpVol;

}

template< int n > 
void	ARM_VFX_Subject<n>::SetVol	(ARM_Object *	mod, ARM_VolCurve *	vol	) {

	ARM_Matrix			Mtmp;
	ARM_Vector*			Vtmp	= NULL;


	ARM_VolCurve*		tmp		= NULL;
	ARM_FXVolCurve*		tmpVol	= NULL;

	tmp = dynamic_cast<ARM_BSModel* > (mod) -> GetVolatility();
	if ( !dynamic_cast<ARM_FXVolCurve*> (tmp) )
		throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "This model do not contain a FX Vol Curve." );
	

	if ( !dynamic_cast<ARM_FXVolCurve*> (vol) )
		throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "This is not a FX Vol Curve." );
	else
		tmpVol = dynamic_cast<ARM_FXVolCurve*> (vol);

	switch(n){

		case	PIV:
			Vtmp = const_cast<ARM_Vector* > (tmpVol -> GetPivotVols());
			dynamic_cast<ARM_FXVolCurve*>(tmp) -> BumpFxVol( ARM_Vector(Vtmp) );
			break;

		case	RR:
			Mtmp = tmpVol -> GetRR();
			dynamic_cast<ARM_FXVolCurve*>(tmp) -> FXBumpRRorSTR( Mtmp, K_YES);
			break;

		case	STR:
			Mtmp = tmpVol -> GetSTR();
			dynamic_cast<ARM_FXVolCurve*>(tmp) ->FXBumpRRorSTR(	Mtmp, K_NO);
			break;
		}
	dynamic_cast< ARM_BSModel* > (mod) -> Set_XXX_Volatility( tmp );
	
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*				Template ARM MIX Subject				*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
 
template< int n > 
string		ARM_MIX_Subject<n>::GetCcy	( ARM_Object* obj	)	const{

	ARM_FXVolCurve*			tmp	 = dynamic_cast<ARM_FXVolCurve*> ( obj ); 
	string					tmpForCcy (tmp->GetForCcy()->GetCcyName() ) ;
	string					tmpDomCcy (tmp->GetDomCcy()->GetCcyName() ) ;
	return tmpForCcy +"/"+ tmpDomCcy;
}

	
template< int n > 
ARM_VolCurve*	ARM_MIX_Subject<n>::GetVol	( ARM_Object *	mod	) const{

	ARM_CurveModelParam*		tmpVol		=	NULL;
	ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( mod ) -> GetModelParams();
	ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast< ARM_ModelParamsMixture_Fx*>( tmp );


	switch(n){

		case	VOL:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::Volatility) );
			break;

		case	SMILE:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::Smile) );
			break;
		
		case	SHIFT:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::Shift) );
			break;

		case	Q:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::QParameter) );
			break;
	}	
	
	ARM_FXVolCurve*	vol =new ARM_FXVolCurve;

	ARM_GP_Vector	gpExp	= tmpVol->GetCurve()->GetAbscisses();
	ARM_GP_Vector	gpVal	= tmpVol->GetCurve()->GetOrdinates();
	ARM_Matrix*		valTerm	= new ARM_Matrix( gpVal.size(),1 );

	sprintf(vol->itsYearTermsY[0], "%d", 0 );
	
	switch(n){
		case	VOL:
			sprintf(vol->itsYearTermsY[0], "%s", "VOL" );
			break;

		case	SMILE:
			sprintf(vol->itsYearTermsY[0], "%s", "SMILE" );
			break;
		
		case	SHIFT:
			sprintf(vol->itsYearTermsY[0], "%s", "SHIFT" );		
			break;

		case	Q:
			sprintf(vol->itsYearTermsY[0], "%s", "Q" );		
			break;
	}


	for ( int i= 0; i< gpVal.size(); i++){
		sprintf(vol->itsYearTermsX[i], "%d", (int) gpExp.Elt(i) );	
		valTerm->Elt(i,0) = gpVal.Elt(i);
	}
	vol->SetVolatilities(valTerm);

	vol->SetDomCcy( dynamic_cast<ARM_ZeroCurve *> (&*tmpFx->GetDomCurve())-> GetCurrencyUnit() -> GetCcyName());
	vol->SetForCcy( dynamic_cast<ARM_ZeroCurve *> (&*tmpFx->GetForCurve())-> GetCurrencyUnit() -> GetCcyName());
	
	return  vol;
}

template< int n> 
void	ARM_MIX_Subject<n>::SetVol	( ARM_Object *	mod, ARM_VolCurve *	vol	){

	ARM_Matrix*		mat		= vol->GetVolatilities();
	ARM_GP_Vector	gpMat(mat->GetNumLines());

	for ( int i= 0; i< mat->GetNumLines(); i++)
		gpMat.Elt(i)= mat->Elt(i,0);

	ARM_CurveModelParam*		tmpVol		=	NULL;
	ARM_ModelParams*			tmp			=	dynamic_cast< ARM_MixtureModel_Fx*> ( mod ) -> GetModelParams();
	ARM_ModelParamsMixture_Fx*	tmpFx		=	dynamic_cast< ARM_ModelParamsMixture_Fx*>( tmp );
	
	switch( n ){

		case	VOL:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::Volatility) );
			break;
		case	SMILE:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::Smile) );
			break;
		case	SHIFT:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::Shift) );
			break;
		case	Q:
			tmpVol		=	dynamic_cast< ARM_CurveModelParam *> ( &tmpFx->GetModelParam(ARM_ModelParamType::QParameter) );
			break;
	}	

	ARM_Curve* cur = tmpVol->GetCurve();
	cur->SetOrdinates( gpMat);
	tmpVol->SetCurve( (ARM_Curve*) cur->Clone());
	tmpFx->SetModelParam( (ARM_CurveModelParam*) tmpVol ->Clone());
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						ARM IBS Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template< int n> 
ARM_VolCurve*	ARM_IBS_Subject<n>::GetVol	( ARM_Object *	obj	)	const{

	ARM_VolCurve* mod = NULL;

	if( dynamic_cast< ARM_VolCurve*> ( obj )  )
		mod = dynamic_cast< ARM_VolCurve*> ( obj );

	else if( dynamic_cast< ARM_InfBSModel*> ( obj )  ){

		ARM_VolCube*	vol		= dynamic_cast< ARM_VolCube *>(dynamic_cast< ARM_InfBSModel*> ( obj )->GetVolatility() );
		vector< ARM_VolCurve* >* vCurve =  vol->GetVols();
		
		ARM_Vector*		tenor	= vol->GetUnderLyings();
		ARM_Vector*		strike	= dynamic_cast<ARM_VolLInterpol*> (  dynamic_cast< ARM_VolCurve *>( (*vCurve)[0] )  )->GetStrikes();
		ARM_Matrix*		values	= new ARM_Matrix(tenor->size(), strike->size() );

		ARM_Matrix*		tmp		= NULL;

		switch (n) {

			case	CPI:{
				for( int i=0; i<tenor->size(); i++){
					tmp = (*vCurve)[i]->GetVolatilities();
					for( int j= 0; j< strike->size();j++)
						values->Elt(i,j)= tmp->Elt(0,j);
				}
			}
			break;
			case	YOY:{
				for( int i=0; i<tenor->size(); i++){
					tmp = (*vCurve)[i]->GetVolatilities();
					if (tmp->GetNumLines() == 1){
						for( int j= 0; j< strike->size();j++)
							values->Elt(i,j)= tmp->Elt(0,j);
					}
					else {
						for( int k=0; k<tmp->GetNumLines();k++){
							for( int j= 0; j< strike->size();j++)
								values->Elt(k+i,j)= tmp->Elt(k,j);
						}
						break;
					}
				}
			}
			break;
		}

		
		mod =  new  ARM_VolCurve(vol->GetAsOfDate(), new ARM_Currency (vol -> GetCurrency() -> GetCcyName() )  );
		for(int i =0; i< strike->size(); i++)
			sprintf(mod->itsYearTermsY[i], "%.2f", strike->Elt(i) );
		for( i =0; i< tenor->size(); i++)
			sprintf(mod->itsYearTermsX[i], "%s", ConvertDoubleToStringMatu(tenor->Elt(i)).c_str() );
		mod->SetVolatilities(values);
	}

	return mod;
}

template< int n> 
void	ARM_IBS_Subject<n>::SetVol	( ARM_Object * mod , ARM_VolCurve* vol	){

	ARM_VolCube*			volCube	= dynamic_cast< ARM_VolCube *>(dynamic_cast< ARM_InfBSModel*> ( mod )->GetVolatility() );
	vector<ARM_VolCurve* >* vCurve =  volCube->GetVols();
	int						dim		= volCube->GetUnderLyings()->size();
	int						nbStrike= dynamic_cast<ARM_VolLInterpol*> (  dynamic_cast< ARM_VolCurve *>( (*vCurve)[0] )  )->GetStrikes()->size();

	ARM_Matrix*				values	= vol -> GetVolatilities();
	ARM_Matrix*				tmp		= NULL;

	switch (n) {

		case	CPI:{
			for( int i=0; i< dim; i++){
				tmp = (*vCurve)[i]->GetVolatilities();
				for( int j= 0; j< nbStrike; j++)
					tmp->Elt(0,j) = values->Elt(i,j);
			}
		}
		break;
		case	YOY:{
			for( int i=0; i< dim; i++){
				tmp = (*vCurve)[i]->GetVolatilities();
				if (tmp->GetNumLines() == 1){
					for( int j= 0; j< nbStrike; j++)
						tmp->Elt(0,j)=values->Elt(i,j);
				}
				else {
					for( int k=0; k<tmp->GetNumLines();k++){
						for( int j= 0; j< nbStrike; j++)
							tmp->Elt(k,j)=values->Elt(k+i,j);
					}
					break;
				}
			}
		}
		break;
	}

}



			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*					ARM CORREL Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


string		ARM_COR_Subject::GetCcy	( ARM_Object* obj	)	const{
		return dynamic_cast<ARM_VolCurve*> ( obj )->GetIndexName();
}


void	ARM_COR_Subject::BuildModel( ARM_Observator * obs, const string & str ) {

	ARM_VolCurve*			vol		= NULL;;
	ARM_CurveCorrelMatrix*	correl	= NULL;
	ARM_GP_Vector			curve;
	ARM_Matrix*				valTerm;
	vector<string>			lab;
	int						pos;
	int						dim;
	int						Index	= 0;
	string					tmp		= str;


	if	( dynamic_cast< ARM_GP_Matrix* > ( obs->GetObject() ) ){
		pos = tmp.find(",");
		while( pos >= 0 && pos <= tmp.size() ){
			lab.push_back( tmp.substr(0,pos) );
			tmp = tmp.substr(pos+1,tmp.size() );
			pos = tmp.find(",");
		}
		if ( tmp.size() !=0)	lab.push_back( tmp );
		correl = new ARM_CurveCorrelMatrix(*dynamic_cast< ARM_GP_Matrix* > ( obs->GetObject() ),lab);
	}
	else if	( dynamic_cast< ARM_VolCurve* > ( obs->GetObject() ) ){
		vol = CreateClone<ARM_VolCurve>(dynamic_cast< ARM_VolCurve* > ( obs->GetObject() ));

		if( string(vol->itsYearTermsX[0]).size()==0 && string(vol->itsYearTermsY[0]).size()==0  ){

			ARM_Vector* expery = vol->GetExpiryTerms();
			ARM_Vector* strike = vol->GetStrikes();

			for (int i = 0; i<expery->size(); i++)
				sprintf(vol->itsYearTermsX[i], "%s", ConvertDoubleToStringMatu ( expery->Elt(i) ).c_str() );
			for (	 i = 0; i<strike->size(); i++)
				sprintf(vol->itsYearTermsY[i],"%s",  ConvertDoubleToStringMatu ( strike->Elt(i) ).c_str());

		}
		obs->SetObject(vol);
	}

	if	( dynamic_cast< ARM_CurveCorrelMatrix* > ( obs->GetObject() ) ){
		correl = dynamic_cast< ARM_CurveCorrelMatrix* > ( obs->GetObject() );
		
		if (!vol)	{
			vol = new ARM_VolCurve;

			Index	= 0;
			dim		= correl -> rows();
			lab		= correl -> GetLabel();

			Index	= 0;
			dim		= correl -> rows();
			lab		= correl -> GetLabel();

			curve = correl->GetCorrelCurves()[1].itsCurve.GetAbscisses();
			for( int k = 0 ; k< curve.size(); k++)
				sprintf(vol->itsYearTermsX[k], "%d", (int) curve.Elt(k) );

			valTerm = new ARM_Matrix(curve.size(), (int)  dim*(dim-1)/2, 0.0);
			tmp.resize(0);
			for( int i =0 ; i< dim; i++ ){
				if( tmp.size()==0 )	tmp = lab[i];
				else				tmp = tmp +"," + lab[i];
			
				for( int j =i+1 ; j< dim; j++ ){
					sprintf(vol->itsYearTermsY[Index], "%s", (correl->GetCorrelCurves()[i*dim+j].itsLabel).c_str() );
					curve = correl->GetCorrelCurves()[ i*dim+j ].itsCurve.GetOrdinates();			
					for( int k = 0 ; k< curve.size(); k++)
						valTerm->Elt(k,Index) = curve.Elt(k);
					Index++;
				}
			}
			vol->SetVolatilities(valTerm);
			vol->SetIndexName( tmp.c_str() );
			obs->SetObject(vol);
		}
	}
}

ARM_VolCurve*	ARM_COR_Subject::GetVol	( ARM_Object *	obj	)	const{
	ARM_VolCurve* vol = NULL;

	if ( dynamic_cast< ARM_VolCurve*> ( obj) )
		vol = dynamic_cast< ARM_VolCurve*> (obj);

	else if ( dynamic_cast< ARM_InfBSModel*> ( obj) ){
		string			key;

		ARM_InfBSModel*		mod		= dynamic_cast< ARM_InfBSModel*> ( obj );
		ARM_CorrelManager* correl	=  mod->GetCorrelManager();
		ARM_CorrelManager::AllMarketCorrelsMap	correlMap = correl->GetMktData();

		ARM_CorrelManager::AllMarketCorrelsMap::const_iterator iter;
		ARM_CorrelManager::intraMarketCorrelsMap::const_iterator it;

		for( iter=correlMap.begin(); iter != correlMap.end(); iter++ ){
			key = iter->first;
			if ( key == "INF_YOY/IR_FWD"){
				for( it = iter->second.begin();  it != iter->second.end() ; it++ ){
					vol = const_cast<ARM_CorrelMatrix*>(&it->second);
					vol = dynamic_cast<ARM_CorrelMatrix* > (vol->Clone() );
					vol ->SetIndexName(key);
					ARM_Vector* expery = vol->GetExpiryTerms();
					ARM_Vector* strike = vol->GetStrikes();


					for (int i = 0; i<expery->size(); i++)
						sprintf(vol->itsYearTermsX[i], "%s", ConvertDoubleToStringMatu ( expery->Elt(i) ).c_str() );
					for (	 i = 0; i<strike->size(); i++)
						sprintf(vol->itsYearTermsY[i],"%s",  ConvertDoubleToStringMatu ( strike->Elt(i) ).c_str());
				}
			}	
		}
	}
	return vol;
}

void	ARM_COR_Subject::SetVol	( ARM_Object *	obj, ARM_VolCurve* vol	){
	ARM_VolCurve* tmp = NULL;

	if(  dynamic_cast<ARM_CorrelMatrix*>(vol) ){
		if ( dynamic_cast< ARM_InfBSModel*> ( obj) ){

			ARM_InfBSModel*		mod		= dynamic_cast< ARM_InfBSModel*> ( obj );
			ARM_CorrelManager* correl	=  mod->GetCorrelManager();
			ARM_CorrelManager::AllMarketCorrelsMap	correlMap = correl->GetMktData();

			ARM_CorrelManager::AllMarketCorrelsMap::iterator iter;
			ARM_CorrelManager::intraMarketCorrelsMap::iterator it;

			for( iter=correlMap.begin(); iter != correlMap.end(); iter++ ){
				if ( iter->first == "INF_YOY/IR_FWD"){
					for( it = iter -> second.begin();  it != iter->second.end() ; it++ ){
						if( it->first.size()!=0){
							ARM_CorrelMatrix* matrix = correl->ComputeCorrelData( "INF_YOY/IR_FWD",  it->first ); 
							matrix ->Copy( dynamic_cast<ARM_CorrelMatrix*>(vol) );
						}
					}	
				}	
			}
		}	
	}
	else
		ARMTHROW(ERR_INVALID_ARGUMENT, "this should be a ARM_CorrelMatrix" ); 
}


void	ARM_COR_Subject::Bump( ARM_Object* vol , const ARM_Matrix & shift, const int & isRelative){
	
	ARM_Matrix*		tmp	= dynamic_cast<ARM_VolCurve * > ( vol )->GetVolatilities();
	ARM_Matrix*		val = new ARM_Matrix ( shift.GetNumLines(), shift.GetNumCols()  );

	for ( int i=0; i<shift.GetNumLines(); i++){
		for( int j=0; j< shift.GetNumCols(); j++){
			if ( isRelative ) 
				val->Elt(i,j) = tmp->Elt(i,j)*( 1 + shift.Elt(i,j) );
			else
				val->Elt(i,j) = tmp->Elt(i,j) + shift.Elt(i,j);
		}
	}
	dynamic_cast<ARM_VolCurve * > ( vol )->SetVolatilities( val );
}


ARM_Object*		ARM_COR_Subject::Convert ( ARM_VolCurve* vol) {
	
	int				dim;
	ARM_Matrix*		tmp	= dynamic_cast<ARM_VolCurve * > ( vol )->GetVolatilities();

	int nbCols = tmp->GetNumCols();
	double nb = 0.5*(1+sqrt(1+8*nbCols));
	if( nb-(int) nb <1e-8 )
		dim = (int) nb;
	else 
		ARMTHROW(ERR_INVALID_ARGUMENT, "it can't be reduced to a correlation matrix" ); 

	ARM_GP_Matrix*	correl= new ARM_GP_Matrix(dim,dim,0.0);

	int Index=0;
	for( int i=0; i< dim ; i++){
		for( int j= i+1; j< dim; j++){
			correl->Elt(i,j) = tmp->Elt(0,Index);
			correl->Elt(j,i) = tmp->Elt(0,Index);
			Index++;
		}
		correl->Elt(i,i) = 1.0;
	}
	return correl;
}





			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*					ARM SO Subject					*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

template< int n>
ARM_VolCurve*	ARM_SO_Subject<n>::GetVol	( ARM_Object *	mod		) const{

	ARM_VolCurve*		tmp		= NULL;
	ARM_BSModel *		tmpBs	= NULL;

	if( dynamic_cast< ARM_VolCurve *>	(mod ) )
		tmp = dynamic_cast< ARM_VolCurve *>	(mod );
	
	else if(  dynamic_cast< ARM_BSModel *>	(mod ) ){
		tmpBs	=	dynamic_cast< ARM_BSModel *>	(mod );
		tmp = tmpBs -> GetVolatility();
		if( dynamic_cast< ARM_VolCube *>( tmp ) ){
			switch( n ){
				case 	ATM:
					tmp = dynamic_cast< ARM_VolCube *> ( tmp )-> GetATMVol( );
					break;
				case	ADJ:
					tmp = tmpBs -> GetCvxAdjVolatility( );
					break;
			}
		}
	}
	return  tmp;
}

template< int n>
void	ARM_SO_Subject<n>::SetVol	(ARM_Object *	mod, ARM_VolCurve *	vol	){

	ARM_VolCurve*		tmp		=	NULL;
	ARM_BSModel *		tmpBs	=	NULL;

	if( dynamic_cast<ARM_VolCurve *> (mod) ){
		mod->Copy(vol);
	}
	else if(  dynamic_cast< ARM_BSModel *>	( mod ) ){
		tmpBs	=	dynamic_cast< ARM_BSModel *>	( mod );
		tmp		=	tmpBs -> GetVolatility();

		if ( dynamic_cast< ARM_VolCube *>	( tmp ) ){
			switch( n ){
				case	ATM:
					dynamic_cast< ARM_VolCube *> ( tmp )-> SetATMVol( vol);
					tmpBs -> SetVolatility( new ARM_VolCurve(*vol) );
					break;
				case	ADJ:
					tmpBs -> SetCvxAdjVolatility( vol );
					break;
			}
		}
		else 
			tmpBs	-> Set_XXX_Volatility( vol );
	}
}


CC_END_NAMESPACE( )
