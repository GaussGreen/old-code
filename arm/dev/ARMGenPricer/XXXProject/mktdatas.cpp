/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file mktdata.cpp
 *  \brief file for the mkt datas container
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */


#include <gpbase/removeidentifiedwarning.h>

#include "xxxproject/mktdatas.h"
#include "xxxproject/observator.h"
#include "xxxproject/argconvdefault.h"

#include <gpbase/cloneutilityfunc.h>

#include <mod/bssmiled.h>
#include <crv/zerocurv.h>
#include <crv/zerointspreaded.h>
#include <inst/forex.h>
#include <utility>


CC_BEGIN_NAMESPACE( ARM )

ARM_MktData::ARM_MktData(	const ARM_Date & asOfDate, const V_Str	& key,	const V_Obj	& mkt ){

	InitMktData(key,mkt);
	Init(string("DATE_DATE"), new ARM_Date(asOfDate) );
}


void ARM_MktData::InitMktData(const V_Str	& key,	const V_Obj	& mkt)
{
	itsSub = ARM_Subject::CreateSubject();

	if (key.size() != mkt.size() )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						 "Keys and MktDatas should have the same size." );	
	itsObs.clear();

	for (int i = 0; i <key.size(); i++)
		Init(key[i], CreateClone(mkt[i]));
}


ARM_MktData* ARM_MktData::CreateCopy	( ) const
{
	// We try to reproduce the behaviour of the constructor
	M_Obs::const_iterator it;

	V_Str	key;
	V_Obj	mkt;

	ARM_MktData* copy = new ARM_MktData();

	for (it = itsObs.begin(); it != itsObs.end(); ++it)
	{
		key.push_back(it->first);
		mkt.push_back(CreateClone(it->second->GetObject()));
	}

	copy->InitMktData(key,mkt);

	return copy;
}

ARM_MktData::ARM_MktData(	const ARM_MktData & rhs ){

	itsSub = rhs.itsSub;
	itsObs = rhs.itsObs;
}


string	ARM_MktData::GetInstrument(const string & ins ){

	string	tmp;

	int pos		=	ins.rfind("_");
	if (pos<0	||	pos > ins.size())
		throw	Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "instrument name not available" );
	else 
		tmp	=	ins.substr(0,pos);

	return tmp;

}


string	ARM_MktData::GetDomCcy (	const string & dom){

	string	tmp;
	int		pos;

	pos			=	dom.rfind("_");
	if (pos<0	||	pos > dom.size())
		throw	Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "instrument name not available" );
	else
		tmp	= dom.substr(pos+1,dom.size() );

	pos			= dom.rfind("/");
	if (pos<0	||	pos > dom.size())
		tmp	= string("");
	else
		tmp	= dom.substr(pos+1,dom.size() );

	return tmp;

}

string	ARM_MktData::GetLabel (	const string & dom){

	string	tmp;
	int		pos;

	pos			=	dom.rfind("_");
	if (pos<0	||	pos > dom.size())
		throw	Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "instrument name not available" );
	else
		tmp	= dom.substr(pos+1,dom.size() );
	
	pos			=	tmp.find("(");
	if (pos>=0	&&	pos <= tmp.size())
		tmp	= tmp.substr(pos+1,tmp.size() );

	pos			=	tmp.find(")");
	if (pos>=0	&&	pos <= tmp.size())
		tmp	= tmp.substr(0,pos );

	return tmp;

}

void ARM_MktData::Reset(){

	It_Sub	its;

	for( It_Mbs ito = itsObs.begin(); ito != itsObs.end(); ito ++){
		for( its = itsSub.begin(); its != itsSub.begin(); its ++){
			(*its)->Detach( ito->second );
		}
	}

	for( its = itsSub.begin(); its != itsSub.begin(); its ++)
		if (*its){ delete( *its ); *its=NULL; }
}	

void ARM_MktData::Init(const string & ins, ARM_Object* obj ){

	string					str = ARM_MktData::GetInstrument(ins);
	string					domCcy = ARM_MktData::GetDomCcy(ins);
	string					label = ARM_MktData::GetLabel(ins);

	V_Sub					sub = ARM_Subject::AttributeSubject	( str, &itsSub, obj);
	ARM_Observator*			obs	= new ARM_Observator(obj,sub,label);

	if (!domCcy.empty())
		obs->SetDomCcy(domCcy);
	pair<string, ARM_Observator*> p( ins,  obs);
	itsObs.insert(p);
}



string  ARM_MktData::toString	(	const string & indent, 	const string & nextIndent) const{

	CC_Ostringstream	os;
	V_Str				ccy;

	os << "\n"<< indent <<"\n";
	os << "\t\t\t" <<"======================================================================"	<<"\n\n";
	os << "\t\t\t" <<"\t\t\t" <<"MARKET DATA MANAGER"<<"\n\n";
	os << "\t\t\t" <<"======================================================================"	<<"\n";
	os << "\n\n";

	os << itsSub[19]->toString( string(""), indent,   nextIndent) ;
	ccy = itsSub[0]->GetCurrency( );

	for ( int i =0; i<ccy.size(); i++){
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(70)<< "CURRENCY : "<<ccy[i]<<"\n";
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(70)<< "===========";
		os<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(ccy[i].size())<<"="<<"\n";

		os << itsSub[0]->toString( ccy[i], indent,   nextIndent) ;
		os << itsSub[9]->toString( ccy[i], indent,   nextIndent) ;
		for ( int j =1; j<9; j++)
			os << itsSub[j]->toString( ccy[i], indent,   nextIndent) ;

		os << itsSub[20]->toString( ccy[i], indent,   nextIndent) ;
		os << itsSub[21]->toString( ccy[i], indent,   nextIndent) ;
		os<<"\n";	
		}

	ccy = itsSub[10]->GetCurrency( );
	for ( i =0; i<ccy.size(); i++){
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(70)<< "CHANGE : "<<ccy[i]<<"\n";
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(70)<< "=========";
		os<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(ccy[i].size())<<"="<<"\n";
		
		os << itsSub[10]->toString( ccy[i], indent,   nextIndent) ;
		for ( int j =11; j<18; j++)
			os << itsSub[j]->toString( ccy[i], indent,   nextIndent) ;

		os<<"\n";
	}

	ccy = itsSub[18]->GetCurrency( );
	for ( i =0; i<ccy.size(); i++){
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(70)<< "CORRELATION : "<<ccy[i]<<"\n";
		os<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(70)<< "==============";
		os<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(ccy[i].size())<<"="<<"\n";
		
		os << itsSub[18]->toString( ccy[i], indent,   nextIndent) ;
		os<<"\n";
	}

	return os.str();
}

CC_END_NAMESPACE( )







