/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: InfHybModel.cpp,v $
 * Revision 1.1  2004/04/19 12:15:13  ebenhamou
 * Initial revision
 *
 *
 *
 */

/*! \file InfHybModel.cpp
 *
 *  \brief inflation hybrid model
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date July 2004
 */

/*
#pragma warning(disable : 4786) 

#include "gpinflation/infhybmodel.h"
#include "mod/bssmiled.h"
#include "crv/zeroint.h"
#include "crv/volint.h"
#include <gpinflation/infcurv.h>
#include <gpinflation/infbsmodel.h>
#include "gpbase/cloneutilityfunc.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_InfHybridModel::ARM_InfHybridModel( const map< string,ARM_Object*> & arg):ARM_Model(){

	map< string, ARM_Object*>::const_iterator it;
	for( it = arg.begin(); it!= arg.end(); it++){
		pair<string, ARM_ObjectPtr> p( it->first, ARM_ObjectPtr(it->second->Clone()) );
		itsMktMap.insert(p);
	}
	CheckMkt();
}


ARM_InfHybridModel::ARM_InfHybridModel( const ARM_InfHybridModel & rhs ){

	map< string, ARM_ObjectPtr>::const_iterator it;

	for( it = rhs.itsMktMap.begin(); it!= rhs.itsMktMap.end(); it++){
		pair<string, ARM_ObjectPtr> p( it->first, CreateClonedPtr(&*it->second) );
		itsMktMap.insert(p);
	}
}


ARM_ObjectPtr ARM_InfHybridModel::Get( const string & key) const{

	map< string, ARM_ObjectPtr>::const_iterator it;
	it= itsMktMap.find(key);
	if ( it == itsMktMap.end() ){
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The market data" + key + "doesn't belong to the Mkt set");
	}
  	else
		return ARM_ObjectPtr( (it->second)->Clone() );
}

string ExtractInstrument( const string& key){

	int pos = key.find("_");
	if (pos<0	||	pos > key.size() ){
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : key is incorrect" );
	}
	else
		return ( key.substr(1,pos) );
}

string ExtractLabel( const string& key){

	int pos = key.find("_");

	if (pos<0	||	pos > key.size() ){
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : key is incorrect" );
	}
	else
		return ( key.substr(pos+1,key.size() ) );
}

void ARM_InfHybridModel::CheckMkt(){

	map< string, ARM_ObjectPtr>::iterator it;
	map< string, ARM_ObjectPtr>::iterator tmpIrgIt;
	map< string, ARM_ObjectPtr>::iterator tmpOswIt;

	vector<string> vec;
	vector<string>::iterator itv;
	for( it = itsMktMap.begin(); it!= itsMktMap.end(); it++){

		if ( ExtractInstrument(it->first) =="YC") {
			string lab = ExtractLabel(it->first);
			itv =find( vec.begin(), vec.end(), lab);
			if ( itv ==vec.end() )
				vec.push_back(lab);
			if ( !dynamic_cast< ARM_ZeroCurve*> ( &*it->second )  ){
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Yield Curve YC_" + lab + " is not a ARM_ZeroCurve");
			}
			else{
				string ccy = dynamic_cast<ARM_ZeroCurve*> ( &*it->second ) -> GetCurrencyUnit() -> GetCcyName();
				if ( ccy!= lab) 
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Label incoherent with the ral currency :" + ccy );


				tmpIrgIt = itsMktMap.find("CAPMOD"+lab);
				tmpOswIt = itsMktMap.find("OSWMOD"+lab);

				if ( tmpIrgIt == itsMktMap.end() && tmpOswIt == itsMktMap.end()){
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : We need a lol market data for :" + ccy );
				}
				else if ( tmpIrgIt != itsMktMap.end() ){
					if ( !dynamic_cast< ARM_BSSmiledModel* > ( &*tmpIrgIt->second )  )	
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BS Model " + lab + " is not a BSSmiledModel");
				}
				else if ( tmpOswIt != itsMktMap.end() ){
					if ( !dynamic_cast< ARM_BSSmiledModel* > ( &*tmpOswIt->second )  )	
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BS Model " + lab + " is not a BSSmiledModel");
				}
			}
		}

		if ( ExtractInstrument(it->first) =="INF") {
			string lab = ExtractLabel(it->first);
			if ( itv ==vec.end() )
				vec.push_back(lab);
			if ( !dynamic_cast< ARM_InfCurv *> ( &*it->second )  ){
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Yield Curve INF_" + lab + " is not a ARM_ZeroCurve");
			}
			else{
				string ccy = dynamic_cast<ARM_InfCurv *> ( &*it->second ) -> GetInfIdxName();
				if ( ccy!= lab) 
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : label incoherent with the ral currency :" + ccy );
			
				tmpIrgIt = itsMktMap.find("CAPMOD"+lab);
				tmpOswIt = itsMktMap.find("OSWMOD"+lab);

				if ( tmpIrgIt == itsMktMap.end() && tmpOswIt == itsMktMap.end() ){
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : We need a vol market data for :" + ccy );
				}
				else if ( tmpIrgIt != itsMktMap.end() ){
					if ( !dynamic_cast< ARM_InfBSModel * > ( &*tmpIrgIt->second )  )	
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BS Model " + lab + " is not a InfBsModel");
				}
				else if ( tmpOswIt != itsMktMap.end() ){
					if ( !dynamic_cast< ARM_InfBSModel * > ( &*tmpOswIt->second )  )	
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": BS Model " + lab + " is not a InfBsModel");
				}
			}
		}
	}
	for ( int i=0; i< vec.size(); i++){
		for ( int j=i+1; j< vec.size(); j++){
			tmpIrgIt = itsMktMap.find("Correl_" +vec[i]+"/"+vec[j] );
			tmpOswIt = itsMktMap.find("Correl_" +vec[j]+"/"+vec[i] );
			if ( tmpIrgIt == itsMktMap.end() && tmpOswIt == itsMktMap.end() ){
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : We need a correlation market data for :" + vec[i] +" and " + vec[j] );
			}
			else if ( tmpIrgIt != itsMktMap.end() ){
				if (  !dynamic_cast< ARM_VolCurve * > ( &*tmpIrgIt->second ) ) 
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": The correlation "+ vec[i] +" and " + vec[j] + " is not a vol curve");
			}
			else if ( tmpIrgIt != itsMktMap.end() )
				if (  !dynamic_cast< ARM_VolCurve * > ( &*tmpOswIt->second ) ) {
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": The correlation "+ vec[i] +" and " + vec[j] + " is not a vol curve");
			}
		}
	}
}



CC_END_NAMESPACE()
*/
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

