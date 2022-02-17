
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Hybride Inflation Leg														 *
 *																							 *
 *			This class builds a corridor leg from swap legs	and inherits from Inf Swap Leg	 *									 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 1st 2007														 *																											 *
 *			Copyright (c) Natixis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
#ifndef _INGPINFLATION_HYBRIDINFLEG_CPP
#define _INGPINFLATION_HYBRIDINFLEG_CPP

#pragma warning(disable : 4786)
#include <gpbase\assignop.h>

#include "gpinfra/argconvdefault.h"

#include "gpbase/typedef.h"
#include "gpbase/datestrip.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/cloneutilityfunc.h"

#include "gpinflation/hybridinfleg.h"
#include "gpinflation/infleg.h"			
#include "gpinflation/infdata.h"

#include "gpinflation/infbssmiledmodel.h"
#include "gpinflation/argconvdefault.h"

#include "mod/bssmiled.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: ARM_HybridInfIrLeg
///	Returns: Void
///	Action : Construction

////////////////////////////////////////////////////

ARM_HybridInfIrLeg::ARM_HybridInfIrLeg(	map<string,ARM_Date>	&	mDate, 
										map<string,string>		&	mString,
										map<string,int>			&	mInt){

	typedef std::pair<IndexType,ARM_InfIrIndex> indexPair;
	indexPair tmp;

	if( mString [ "mainIndex"] != "N"){
		tmp =	indexPair( MAI,	ARM_InfIrIndex(mString [ "mainIndex"], mDate, mString, mInt ) );
		itsIdx.insert(tmp);
	}

	if( mString [ "subIndex" ]!= "N"){
		tmp =	indexPair( SUB,	ARM_InfIrIndex(mString [ "subIndex" ], mDate, mString, mInt ) );
		itsIdx.insert(tmp);
	}

	if( mString [ "supIndex" ]!= "N"){
		tmp =	indexPair( SUP,	ARM_InfIrIndex(mString [ "supIndex" ], mDate, mString, mInt ) );
		itsIdx.insert(tmp);
	}

	itsAsOfDate		= mDate		[	"asOfDate"			];
	itsStartDate	= mDate		[	"startDate"			];
	itsEndDate		= mDate		[	"endDate"			];
	itsCurrency		= mString	[	"currency"			];
	isInit			= false;

	BuildSchedule( );
}

////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: BuildSchedule
///	Returns: Void
///	Action : Construction

////////////////////////////////////////////////////


void ARM_HybridInfIrLeg::BuildSchedule(){
	
	if ( itsIdx[MAI].itsType == NO){
		ARMTHROW(ERR_INVALID_ARGUMENT," The Main Index should be initialized");	
	}

	idxIter it=itsIdx.find(SUP);
	if ( it!=itsIdx.end() ) {
		if (  itsIdx[SUP].itsType != NO && itsIdx[SUB].itsType== NO )
			ARMTHROW(ERR_INVALID_ARGUMENT," The Sub Index should be initialized before the Sup Index");	
	}

	itsNbFlows		= itsIdx[MAI].itsIndex ->GetFlowStartDates()->GetSize();

	itsBegDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsEndDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsResDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		 
	itsPayDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	
	itsIntTerms		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsIntrDays		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsDiscFact		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	
	itsIntTerms     = ARM_InfIrIndex::ConvertToARM_GP_Vector (  itsIdx[MAI].itsIndex -> GetInterestTerms() );
	itsIntrDays     = ARM_InfIrIndex::ConvertToARM_GP_Vector (  itsIdx[MAI].itsIndex ->GetInterestDays() );

	for ( it= itsIdx.begin(); it != itsIdx.end(); it++){

		if ( it->first == MAI){
			itsBegDates		= ARM_InfIrIndex::ConvertToARM_GP_Vector ( it->second.itsIndex->GetFlowStartDates() );
			itsEndDates		= ARM_InfIrIndex::ConvertToARM_GP_Vector ( it->second.itsIndex->GetFlowEndDates()	);
			itsPayDates		= ARM_InfIrIndex::ConvertToARM_GP_Vector ( it->second.itsIndex->GetPaymentDates()	);
			itsResDates		= ARM_InfIrIndex::ConvertToARM_GP_Vector ( it->second.itsIndex->GetResetDates()		);
		}

		if ( it->second.itsType == IN ){
			itsNumDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsDemDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			isHybrid		= true;
		}

		if ( it->second.itsType != NO){
			itsRes[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsTen[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsBeg[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsEnd[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsPay[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsFwd[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsAdj[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
			itsVol[it->first] = ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );

			itsRes[it->first] =	it->second.GetResLags();
			itsTen[it->first] =	it->second.GetTenLags();
			itsBeg[it->first] =	it->second.GetBegDates();
			itsEnd[it->first] =	it->second.GetEndDates();

			for( int i = 0; i<itsNbFlows; i++)
				itsPay[it->first]->Elt(i) = itsPayDates->Elt(i);

			if ( it->second.itsType == IN ){
				itsNumDates= it->second.GetNumDates();
				itsDemDates= it->second.GetDemDates();
				isHybrid = true;
			}
		}
	}

	char	d[20];

	ARM_GramFctorArg*	tmpFunctor = new ARM_GramFctorArg;
	ARM_StringVector*	vec =  new ARM_StringVector( itsNbFlows );
	for ( int i = 0; i<itsNbFlows; i++){
		( (ARM_Date) itsBegDates->Elt(i)).JulianToSpeDateDay(d);		
		(*vec)[i] = (string) d;
	}
	tmpFunctor->SetStringVector( ARM_StringVectorPtr( vec ) );
	itsFunctor.InsertData( "STARTDATES", *tmpFunctor	);


 	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}
	tmpFunctor	=	new ARM_GramFctorArg;
	vec			=	new ARM_StringVector( itsNbFlows );
	for ( i = 0; i<itsNbFlows; i++){
		( (ARM_Date) itsEndDates->Elt(i)).JulianToSpeDateDay(d);		
		(*vec)[i] = (string) d;
	}
	tmpFunctor->SetStringVector( ARM_StringVectorPtr( vec ) );
	itsFunctor.InsertData( "ENDDATES", *tmpFunctor	);

	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}
	tmpFunctor	=	new ARM_GramFctorArg;
	vec			=	new ARM_StringVector( itsNbFlows );
	for ( i = 0; i<itsNbFlows; i++){
		( (ARM_Date) itsResDates->Elt(i)).JulianToSpeDateDay(d);		
		(*vec)[i] = (string) d;
	}
	tmpFunctor->SetStringVector( ARM_StringVectorPtr( vec ) );
	itsFunctor.InsertData( "RESDATES", *tmpFunctor	);

	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}
	tmpFunctor	=	new ARM_GramFctorArg;
	vec			=	new ARM_StringVector( itsNbFlows );
	for ( i = 0; i<itsNbFlows; i++){
		( (ARM_Date) itsPayDates->Elt(i)).JulianToSpeDateDay(d);		
		(*vec)[i] = (string) d;
	}
	tmpFunctor->SetStringVector( ARM_StringVectorPtr( vec ) );
	itsFunctor.InsertData( "PAYDATES", *tmpFunctor	);

}


////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: ARM_HybridInfIrLeg
///	Returns: Void
///	Action : Copy construction

////////////////////////////////////////////////////


ARM_HybridInfIrLeg::ARM_HybridInfIrLeg( const ARM_HybridInfIrLeg& 	rhs	){

	itsNbFlows			= rhs.itsNbFlows;
	itsCurrency			= rhs.itsCurrency;
	isHybrid			= rhs.isHybrid;
	isInit				= rhs.isInit;

	itsAsOfDate			= rhs.itsAsOfDate;
	itsStartDate		= rhs.itsStartDate;
	itsEndDate			= rhs.itsEndDate;

	itsCorObj			= rhs.itsCorObj;
	itsIdx				= rhs.itsIdx;

	itsBegDates			= CreateClonedPtr ( &*	rhs.itsBegDates			);	
	itsEndDates			= CreateClonedPtr ( &*	rhs.itsEndDates			);
	itsResDates			= CreateClonedPtr ( &*	rhs.itsResDates			);		 
	itsPayDates			= CreateClonedPtr ( &*	rhs.itsPayDates			);
	
	itsIntTerms			= CreateClonedPtr ( &*	rhs.itsIntTerms			);
	itsDiscFact			= CreateClonedPtr ( &*	rhs.itsDiscFact			);

	itsNumDates			= CreateClonedPtr ( &*	rhs.itsNumDates			);		
	itsDemDates			= CreateClonedPtr ( &*	rhs.itsDemDates			);	

	ARM_GP_MapPtr::const_iterator	vec;
	ARM_GP_CorPtr::const_iterator	cor;
	ARM_GP_MapIdx::const_iterator	idx;
	

	for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){
		vec = rhs.itsRes.find(it->first);
		itsRes[it->first] = CreateClonedPtr ( &*	vec->second	);

		vec = rhs.itsTen.find(it->first);
		itsTen[it->first] = CreateClonedPtr ( &*	vec->second	);

		vec = rhs.itsBeg.find(it->first);
		itsBeg[it->first] = CreateClonedPtr ( &*	vec->second	);

		vec = rhs.itsEnd.find(it->first);
		itsEnd[it->first] = CreateClonedPtr ( &*	vec->second	);

		vec = rhs.itsPay.find(it->first);
		itsPay[it->first] = CreateClonedPtr ( &*	vec->second	);

		vec = rhs.itsFwd.find(it->first);
		itsFwd[it->first] = CreateClonedPtr ( &*	vec->second	);

		vec = rhs.itsAdj.find(it->first);
		itsAdj[it->first] = CreateClonedPtr ( &*	vec->second	);
		
		vec = rhs.itsVol.find(it->first);
		itsVol[it->first] = CreateClonedPtr ( &*	vec->second	);
	}

	itsCor = rhs.itsCor; 
	for ( corIter ip=itsCor.begin(); ip!=itsCor.end(); ip++){
		ARM_GP_CorPtr::const_iterator	ipTmp=(rhs.itsCor).find(ip->first);
		itsCor[ip->first]= CreateClonedPtr ( &*ipTmp->second );
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: Validate
///	Returns: Void
///	Action : Validation of the mkt data manager.

/////////////////////////////////////////////////////


void	ARM_HybridInfIrLeg::Validate( ARM_MarketData_ManagerRep* mkt ){

	string	tmp;
	bool	isYc = false;
	ARM_MarketData_ManagerRepPtr MktManager =	CreateClonedPtr(&*mkt);
	ARM_MarketData_ManagerImp*   MktMap		=	MktManager->GetMktDataManager();
	ARM_MarketData_ManagerImp::iterator itm;

	for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){

		if ( it->second.itsType != NO) {
			tmp = it->second.GetName( );

			if( it->second.itsType==IR){
				isYc = true;
				ARM_BSSmiledModel* cap = dynamic_cast<ARM_BSSmiledModel*>( MktMap->GetData("CAPMOD_" + tmp) );
				if(!cap)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : BS Smiled Model for key=" + tmp + " is expected in the Market Data Manager");

				ARM_ZeroLInterpol* yc = dynamic_cast<ARM_ZeroLInterpol*>( MktMap->GetData("YC_" + tmp) );
				if(!yc)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Zero Curve for key=" + tmp + " is expected in the Market Data Manager");
			}
		
			else if ( it->second.itsType==IN){
				ARM_InfBSSmiledModel* cap = dynamic_cast<ARM_InfBSSmiledModel*>( MktMap->GetData("CAPMOD_" + tmp) );
				if(!cap)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Inf BS Smiled Model for key=" + tmp + " is expected in the Market Data Manager");

				ARM_InfCurv* yc = dynamic_cast<ARM_InfCurv*>( MktMap->GetData("INF_" + tmp) );
				if(!yc)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Inf Curve for key=" + tmp + " is expected in the Market Data Manager");
			}
		}

		for ( idxIter is= itsIdx.begin(); is != itsIdx.end(); is++){
			if ( is->second.itsType != NO && it->first != is->first ){
				string str = is->second.GetName( );
				ARM_VolCurve* cor = dynamic_cast<ARM_VolCurve*>( MktMap->GetData("CORREL_"+tmp +"/"+ str) );
				if ( !cor ){
					cor = dynamic_cast<ARM_VolCurve*>( MktMap->GetData("CORREL_"+str+"/"+tmp) );
					if ( !cor ){
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correl Curve for key="+ str+"/"+tmp+ " is expected in the Market Data Manager");
					}
				}
			}
		}
	}
}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: Init
///	Returns: Void
///	Action : build latente quantities like Discount Factor at payment dates

/////////////////////////////////////////////////////


void ARM_HybridInfIrLeg::Init( ARM_MarketData_ManagerRep* mkt){

	Validate(mkt);
	double tmp, tmpStr;
	ARM_GP_VecPtr tmpDates = itsPayDates; 

	ARM_MarketData_ManagerImp*  mktMap	=	mkt->GetMktDataManager();
 	ARM_ZeroLInterpol*			dfCurve	=	dynamic_cast<ARM_ZeroLInterpol*> (mktMap->GetData("YC_"+itsCurrency) );
	
	double initDate = itsAsOfDate.GetJulian();
	for ( int i =0; i< itsNbFlows; i++){
		tmp = ( itsPayDates->Elt(i) - initDate )/ K_YEAR_LEN;
		itsDiscFact->Elt(i) =dfCurve->nvDiscountPrice( tmp );
	}

	for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){
		string				str	= it->second.GetName( );
		ARM_BSSmiledModel*	mod	=	dynamic_cast<ARM_BSSmiledModel*> (mktMap->GetData("CAPMOD_"+str ) );
	
		it->second.SetMod(mod);

		ARM_GP_VecPtr beg = itsBeg[it->first];
		ARM_GP_VecPtr end = itsEnd[it->first];
		ARM_GP_VecPtr pay = itsPay[it->first];
		ARM_GP_VecPtr res = itsRes[it->first];
		ARM_GP_VecPtr ten = itsTen[it->first];

		for ( int i =0; i<itsNbFlows; i++){	
			itsFwd[it->first]->Elt(i) = it->second.CptFwd(beg->Elt(i), end->Elt(i), pay->Elt(i));
			tmp = itsFwd[it->first]->Elt(i);
			tmpStr = pow(1.0+tmp/100,ten->Elt(i));

			itsAdj[it->first]->Elt(i) = it->second.CptAdj(beg->Elt(i), end->Elt(i), pay->Elt(i), tmp);
			itsVol[it->first]->Elt(i) = it->second.CptVol(res->Elt(i), ten->Elt(i), 1.0+tmp/100.0, tmpStr);
		}

		for ( idxIter is= itsIdx.begin(); is != itsIdx.end(); is++){
			if ( is->second.itsType != NO && it->first != is->first ){
				string	stg			= is->second.GetName( );
				ARM_VolLInterpol* mod	= dynamic_cast<ARM_VolLInterpol*> (mktMap->GetData("CORREL_"+stg +"/"+ str) );
				if ( mod )
					itsCorObj.SetCorrel( stg, str, mod);
			}
		}
	}

	itsCor.clear();
	for ( it= itsIdx.begin(); it != itsIdx.end(); it++){
		for ( idxIter is= itsIdx.begin(); is != itsIdx.end(); is++){

			string	str			= it->second.GetName( );
			string	stg			= is->second.GetName( );
				
			ARM_GP_VecPtr	res	= itsRes[it->first];
			ARM_GP_VecPtr	ten	= itsTen[is->first];
			ARM_VolCurvePtr	vol	= itsCorObj.GetCorrel(str,stg);
			if ( &*vol ) {
				ARM_GP_VecPtr vec	= ARM_GP_VecPtr (new ARM_GP_Vector);
				for ( int i =0; i<itsNbFlows; i++){	
					tmp = res->Elt(i)<1.0/K_YEAR_LEN?1.0/K_YEAR_LEN:res->Elt(i);
					vec->push_back(vol->ComputeVolatility( tmp,ten->Elt(i)) );
				}
				pair<IndexType,IndexType> p(it->first,is->first);
				itsCor[p] = vec;
			}
		}
	}
	isInit = true;
}

/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: View
///	Returns: Void
///	Action : General View

/////////////////////////////////////////////////////

void	ARM_HybridInfIrLeg::View( char* id , FILE* ficOut ){

    FILE* fOut;
    char fOutName[200];	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	CC_Ostringstream	os;
	os << ViewGlob();

	fprintf(fOut, "%s" ,(os.str()).c_str() );
	if ( ficOut == NULL )
		fclose(fOut);
}

/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: ViewGlob
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

string	ARM_HybridInfIrLeg::ViewGlob( ){

	CC_Ostringstream	os;
	Ostream<>			ss;
	Ostream<string,int> is;
	string				str;
	char				d[20];
	double				tmp;
	vector<string>		vc(4);
	int					nb = 4;


	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(55)<<"HYBRIID INFLATION LEG"<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n";	
	os	<<"INDEX DESCRIPTION\n";
	os	<<"================="<<"\n";

	os	<< ss.ToStream("Index Currency",	itsIdx[MAI].GetIndex()->GetCurrencyUnit()->GetCcyName() );

	for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){
		if ( it->second.itsType != NO){
			str = ARM_ArgConvReverse_SubordIndex.GetString(it->first);
			os	<< ss.ToStream( str + " Index",	it->second.GetIndexName() );
		}
	}

	os	<< "\n";
	os	<<"SCHEDULE DESCRIPTION\n";
	os	<<"===================="<<"\n";

	itsAsOfDate.JulianToStrDate(d);
	os	<< ss.ToStream("As of Date",(string) d);
	itsStartDate.JulianToStrDate(d);
	os	<< ss.ToStream("Start Date",(string) d);
	itsEndDate.JulianToStrDate(d);
	os	<< ss.ToStream("End Date",	(string) d);
	os	<< "\n";
	
	os	<< ss.ToStream("Reset Calendar",	itsIdx[MAI].itsIndex->GetResetCalName() );
	os	<< ss.ToStream("Reset Frequency",	ARM_ParamView::GetMappingName(	S_FREQUENCY,	itsIdx[MAI].GetIndex()->GetResetFrequency() ));
	os	<< ss.ToStream("Reset Timing",		ARM_ParamView::GetMappingName(	S_TIMING_MOD,	itsIdx[MAI].GetIndex()->GetResetTiming() ));    
	os	<< is.ToStream("Reset Gap",			itsIdx[MAI].GetIndex()->GetResetGap() );

	for ( it= itsIdx.begin(); it != itsIdx.end(); it++){
		if ( it->second.itsType != IN){
			os	<< ss.ToStream("Num Gap",	ARM_ParamView::GetMappingName(	S_INDEX_TYPES, it->second.GetIndex()->GetIndexType() ));
			os	<< ss.ToStream("Dem Gap",	ARM_ParamView::GetMappingName(	S_INDEX_TYPES, it->second.GetIndex()->GetIndexType() ));
			break;
		}
	}

	os	<< ss.ToStream("Pay Calendar",		itsIdx[MAI].itsIndex->GetPayCalName() );
	os	<< ss.ToStream("Pay Frequency",		ARM_ParamView::GetMappingName( S_FREQUENCY,		itsIdx[MAI].GetIndex()->GetPayFrequency() ));
	os	<< ss.ToStream("Reset Timing",		ARM_ParamView::GetMappingName( S_TIMING_MOD,	itsIdx[MAI].GetIndex()->GetPayTiming() ));    
	os	<< is.ToStream("Reset Gap",			itsIdx[MAI].GetIndex()->GetPayGap() );
  	os	<< "\n";
	os	<< ss.ToStream("Forward Rule",		ARM_ParamView::GetMappingName( S_FORWARD_RULES,	itsIdx[MAI].GetIndex()->GetFwdRule() ));
	os	<< ss.ToStream("Int Rule",			ARM_ParamView::GetMappingName( S_INTEREST_RULES,itsIdx[MAI].GetIndex()->GetIntRule() )); 
	os	<< ss.ToStream("Stub Rule",			ARM_ParamView::GetMappingName( S_STUB_RULES,	itsIdx[MAI].itsIndex->GetStubMeth() )); 
	os	<< ss.ToStream("Adj First Rule",	ARM_ParamView::GetMappingName( S_INTEREST_RULES,itsIdx[MAI].itsIndex->GetAdjStartDateFlag() )); 
	os	<< ss.ToStream("Day Count",			ARM_ParamView::GetMappingName( S_DAYCOUNT,		itsIdx[MAI].GetIndex()->GetDayCount() ));

	os	<< "\n";	
	os	<<"SCHEDULE MANAGER\n";
	os	<<"================"<<"\n";
	

	os	<< "\n" <<"General Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Start Dates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"End Dates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Res Dates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Pay Dates";
	
	if ( isHybrid ){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Num Dates";
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Dem Dates";
		vc.resize(6);
		nb += 2;
		}

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Int Terms";
	nb++;


	for ( it= itsIdx.begin(); it != itsIdx.end(); it++){
		if ( it->second.itsType != NO){
			str = ARM_ArgConvReverse_SubordIndex.GetString(it->first);
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<str+" Terms";
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<str+" Tenors";
			nb += 2;
		}
	}


	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(LAG*nb+2)<<" \n";

	os <<"\n";

	for ( int i = 0; i < itsNbFlows; ++i )
	{
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;

		( (ARM_Date) itsBegDates->Elt(i)).JulianToStrDate(d);		vc[0] = (string) d;
		( (ARM_Date) itsEndDates->Elt(i)).JulianToStrDate(d);		vc[1] = (string) d;
		( (ARM_Date) itsResDates->Elt(i)).JulianToStrDate(d);		vc[2] = (string) d;
		( (ARM_Date) itsPayDates->Elt(i)).JulianToStrDate(d);		vc[3] = (string) d;
		
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[3];
		
		if ( isHybrid ){
			( (ARM_Date) itsNumDates->Elt(i)).JulianToStrDate(d);		vc[3] = (string) d;
			( (ARM_Date) itsDemDates->Elt(i)).JulianToStrDate(d);		vc[4] = (string) d;
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[3];
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[4];
		}
		
		tmp =	itsIntTerms ->Elt(i);
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<tmp;

		for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){
			if ( it->second.itsType != NO){
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<itsRes[it->first]->Elt(i);
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<itsTen[it->first]->Elt(i);
			}
		}
		os <<"\n";
	}

	
	if ( isInit ){
	  	os	<< "\n";	
		os	<<"SCHEDULE INIT\n";
		os	<<"============="<<"\n";
		os  << ViewInit	( );
	}

	return os.str();  
}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfIrLeg
///	Routine: ViewInit
///	Returns: string
///	Action : return the initialized (by mkt) quantities like the discount factor

/////////////////////////////////////////////////////


string	ARM_HybridInfIrLeg::ViewInit	(  ){
   
	CC_Ostringstream	os;
	double				tmp;

	int	nb	= itsIdx.size();
	nb		= 3*nb + 1 +int( nb*(nb-1)/2 );

	os	<< "\n" <<"Initialized Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Disc Factors";


	for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){
		if ( it->second.itsType != NO){
			string str = ARM_ArgConvReverse_SubordIndex.GetString(it->first);
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<str +" Fwd";
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<str +" Adj";
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<str +" Vol";
		}
	}
	for ( it= itsIdx.begin(); it != itsIdx.end(); it++){
		for ( idxIter is= itsIdx.begin(); is != itsIdx.end(); is++){
			pair<IndexType,IndexType> p(it->first,is->first);
			corIter iter = itsCor.find( p );
			if ( iter !=itsCor.end() ){	
				string str = itsIdx[it->first].GetName();
				string stg = itsIdx[is->first].GetName();
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" Cor "+ str +"/"+stg;
			}
		}
	}

	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(nb*LAG+2)<<" \n";
	os	<< "\n"; 
	for ( int i = 0; i < itsNbFlows; ++i ){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;

		tmp =	itsDiscFact ->Elt(i);
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmp;	

		for ( idxIter it= itsIdx.begin(); it != itsIdx.end(); it++){
			if ( it->second.itsType != NO){
				tmp = itsFwd[it->first]->Elt(i);
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmp;
				tmp = itsAdj[it->first]->Elt(i);
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmp;
				tmp = itsVol[it->first]->Elt(i);
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmp;
			}
		}

		for ( it= itsIdx.begin(); it != itsIdx.end(); it++){
			for ( idxIter is= itsIdx.begin(); is != itsIdx.end(); is++){
				pair<IndexType,IndexType> p(it->first,is->first);
				corIter iter = itsCor.find( p );
				if ( iter !=itsCor.end() ){	
					tmp = itsCor[p]->Elt(i);
					os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmp;
				}
			}
		}
		os <<"\n";
	}
	return os.str();
}

CC_END_NAMESPACE()	
#endif
/*---------------------------------------------------------------*/
/*---- End Of File ----*/

