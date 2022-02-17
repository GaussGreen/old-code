/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file hedge.cpp
 *  \brief file for pricer builder
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */
#include <pricer\ipricer.h>


#include <gpinfra/gramfunctorargdict.h>
#include <gpinfra/gramfunctorarg.h>
#include "xxxproject/hedge.h"
#include "xxxproject/scenario.h"
#include "xxxproject/mktdatas.h"
#include "xxxproject/pricer.h"
#include "xxxproject/pricerbuilder.h"

#include <mod/bssmiled.h>
#include <crv/zerocurv.h>

CC_BEGIN_NAMESPACE( ARM )

ARM_Hedge::ARM_Hedge(const ARM_Hedge& rhs):itsInitPrice(rhs.itsInitPrice),itsFunctor(rhs.itsFunctor){

	itsPricer		=	rhs.itsPricer;
	itsMktData		=	rhs.itsMktData;
	itsFunctor		=	rhs.itsFunctor;
}


void ARM_Hedge::Init(	ARM_Object* theSecurity, 
						ARM_MktData*  theMktData){

	if (theSecurity)
		itsPricer		=	ARM_PricerBuilder::BuildPricer	(theSecurity);
	else
		itsPricer		=	NULL;
	itsMktData		=	theMktData;
	if (itsPricer)
	{
		itsPricer		->	SetMkt( itsMktData );
		itsInitPrice	=	itsPricer->Price();
	}
}



void ARM_Hedge::ComputeHedge( ARM_Scenario* theScenario){

	theScenario -> InitMktData (  itsMktData	);
	if (  dynamic_cast< ARM_ND_Scenario* > ( theScenario ) )
		Apply_ND_Scenario( theScenario, theScenario->GetNbShift(), true	);
	else if ( theScenario -> GetDim() == 0)
		Apply_0D_Scenario( theScenario, true	);
	else if	( theScenario -> GetDim() == 1)
		Apply_1D_Scenario( theScenario, theScenario->GetNbShift(), true );
	else if	( theScenario -> GetDim() == 2)
		Apply_2D_Scenario( theScenario, theScenario->GetNbShift(), true	);
}

void ARM_Hedge::ApplyScenario( ARM_Scenario* theScenario, int nbShift){

	theScenario -> InitMktData (  itsMktData	);
	if (  dynamic_cast< ARM_ND_Scenario* > ( theScenario ) )
		Apply_ND_Scenario( theScenario, nbShift, false	);
	else if ( theScenario -> GetDim() == 0)
		Apply_0D_Scenario( theScenario, false	);
	else if	( theScenario -> GetDim() == 1)
		Apply_1D_Scenario( theScenario, nbShift, false );
	else if	( theScenario -> GetDim() == 2)
		Apply_2D_Scenario( theScenario, nbShift, false	);
}

void ARM_Hedge::Apply_0D_Scenario( ARM_Scenario* theScenario, bool computeHedge ){

	ARM_GramFctorArg			tmpFunctor;

	ARM_GP_VectorPtr			Container	(new ARM_GP_Vector		( 1, 0.0	)	);
	ARM_StringVectorPtr			expLabel	(new ARM_StringVector	( 1 )	);	

	(*expLabel)[0]=	theScenario -> GetPlotsLabel()[0][0];
	theScenario -> InitPosition();

	theScenario	->	ShiftNextPlot();
	theScenario ->	ShiftMktData( itsMktData );

	if (computeHedge)
	{
		(*Container).Elt(0)	=	itsPricer->Price()- itsInitPrice;

		theScenario ->	Finalize( itsMktData );
		
		tmpFunctor.SetVector(Container);
		itsFunctor.InsertData( theScenario->GetScenarioKey(), tmpFunctor	);

		tmpFunctor.SetStringVector(expLabel);
		itsFunctor.InsertData( string("EXPIRY"), tmpFunctor	);
	}
}


void ARM_Hedge::Apply_1D_Scenario( ARM_Scenario* theScenario, int nbShift, bool computeHedge){

	int							Ind;
	double						tmpPrice;
	double						tmpPriceRef = itsInitPrice;

	ARM_GramFctorArg			tmpFunctor;

	vector<ARM_StringVector>	tmpPlotsPos		=	theScenario -> GetPlotsPosition();
	vector<ARM_StringVector>	tmpPlotsLab		=	theScenario -> GetPlotsLabel();
	map<string, ScenarioType>	tmpPriceCritera	=	theScenario -> GetPriceCritera();

	ARM_GP_VectorPtr			Container	(new ARM_GP_Vector		( tmpPlotsPos[0].size(), 0.0	)	);
	ARM_StringVectorPtr			expLabel	(new ARM_StringVector	( tmpPlotsPos[0].size()		)	);	

	Ind =0;
	for(int  i=0; i<tmpPlotsLab[0].size();	i++){
		if ( find(	tmpPlotsPos[0].begin(),	tmpPlotsPos[0].end(), tmpPlotsLab[0][i] ) != tmpPlotsPos[0].end() ){
			(*expLabel)[Ind]= tmpPlotsLab[0][i];
			Ind++;
		}
	}

	int		row = 0;
	theScenario -> InitPosition();

	for( i=0; (i<theScenario->GetNbShift()) && (i < nbShift);i++){

		theScenario	->	ShiftNextPlot();
		theScenario				->	ShiftMktData( itsMktData );

		if (computeHedge)
		{
			switch ( tmpPriceCritera[tmpPlotsLab[0][i]] ){

				case REF:{
					tmpPriceRef				=	itsPricer->Price();
					break;
					}

				case VALUE:{
					tmpPrice				=	itsPricer->Price();
					(*Container).Elt(row)	=	tmpPrice - tmpPriceRef;
					row++;
					break;
				}

				case REF_VALUE:{
					tmpPrice				=	itsPricer->Price();
					(*Container).Elt(row)	=	tmpPrice - tmpPriceRef;
					tmpPriceRef				=	tmpPrice;
					row++;
					break;
				}
			}
		}

		if (computeHedge || ((i < nbShift-1) && (i < theScenario->GetNbShift()-1)))
			theScenario				->	Finalize( itsMktData );
	}

	if (computeHedge)
	{
		theScenario -> InitPosition();

		theScenario -> OrderData( Container);
		tmpFunctor.SetVector(Container);
		itsFunctor.InsertData( theScenario->GetScenarioKey(), tmpFunctor	);

		tmpFunctor.SetStringVector(expLabel);
		itsFunctor.InsertData( string("EXPIRY"), tmpFunctor	);
	}

}



void ARM_Hedge::Apply_2D_Scenario( ARM_Scenario* theScenario, int nbShift, bool computeHedge ){

	int							Ind;
	vector<It_Str >				its			( 2 );
	ARM_IntVector				coordEff	( 2 );
	ARM_IntVector				pos			( 2 );
	ARM_IntVector				coord		( 2, 0 );	
	ARM_IntVector				nbCoord		( 2, 0 );
	vector<ARM_StringVector>	tmpPlotsPos	( 2 );
	vector<ARM_StringVector>	tmpPlotsLab	( 2 );
	map<string, ScenarioType>	tmpPriceCritera	=	theScenario -> GetPriceCritera();

	ARM_IntVector				id			=	theScenario -> GetPermOrder();	

	for(int  j=0 ;	j<2 ;	j++){
		tmpPlotsPos[j]	=	theScenario -> GetPlotsPosition()[j];
		tmpPlotsLab[j]	=	theScenario -> GetPlotsLabel()[j];
		for( int i=0; i<tmpPlotsLab[j].size();	i++){
			if ( find(	tmpPlotsPos[j].begin(),	tmpPlotsPos[j].end(), tmpPlotsLab[j][i] )	!=	tmpPlotsPos[j].end()	)	
				nbCoord[j]++;
		}
	}

	ARM_GP_MatrixPtr			Container	( new ARM_GP_Matrix		( nbCoord[0], nbCoord[1], 0.0) );
	ARM_StringVector*			expLabel	= new ARM_StringVector	( nbCoord[0], ""	);
	ARM_StringVector*			tenLabel	= new ARM_StringVector	( nbCoord[1], ""	);
	
	double						tmpPrice;
	double						tmpPriceRef = itsInitPrice;

	ARM_GramFctorArg			tmpFunctor;

	theScenario -> InitPosition();

	Ind =0;
	for(int i=0; i< tmpPlotsLab[0].size();	i++){
		if ( find(	tmpPlotsPos[0].begin(),	tmpPlotsPos[0].end(),tmpPlotsLab[0][i]	)	!=	tmpPlotsPos[0].end()){
			(*expLabel)[Ind]= tmpPlotsLab[0][i];
			Ind++;
		}
	}

	Ind =0;
	for( i=0; i< tmpPlotsLab[1].size();	i++){
		if ( find(	tmpPlotsPos[1].begin(),	tmpPlotsPos[1].end(),tmpPlotsLab[1][i]	)	!=	tmpPlotsPos[1].end()){
			(*tenLabel)[Ind]= tmpPlotsLab[1][i];
			Ind++;
		}
	}
	
	for( i=0; (i<theScenario->GetNbShift()) && (i < nbShift);i++){

		theScenario	->	ShiftNextPlot();
		pos			=	theScenario->GetPosition();

		theScenario				->	ShiftMktData( itsMktData );

		if (computeHedge)
		{
			switch ( tmpPriceCritera[ tmpPlotsLab[0][pos[0] ] + tmpPlotsLab[1][pos[1]] ] ){
				case REF:{
					tmpPriceRef				=	itsPricer->Price();
					break;
				}

				case VALUE:{
					tmpPrice				=	itsPricer->Price();
					coordEff[0]				=	coord[1]*id[0] + coord[0]*id[1];
					coordEff[1]				=	coord[0]*id[0] + coord[1]*id[1];
					Container->Elt(coordEff[0], coordEff[1])	=	tmpPrice - tmpPriceRef;
					if ( coord[0]<nbCoord[1]* id[0] + nbCoord[0] * id[1] -1 )
						coord[0]++ ;
					else { 
						coord[1]++; 
						coord[0] = 0;
						}
					break;
				}

				case REF_VALUE:{
					tmpPrice				=	itsPricer->Price();
					coordEff[0]				=	coord[1]*id[0] + coord[0]*id[1];
					coordEff[1]				=	coord[0]*id[0] + coord[1]*id[1];
					Container->Elt(coordEff[0], coordEff[1])	=	tmpPrice - tmpPriceRef;
					if ( coord[0]<nbCoord[1]* id[0] + nbCoord[0] * id[1] -1 )
						coord[0]++ ;
					else { 
						coord[1]++; 
						coord[0] = 0;
						}
					tmpPriceRef				=	tmpPrice;
					break;
				}
			}
		}

		if (computeHedge || ((i < nbShift-1) && (i < theScenario->GetNbShift()-1)))
			theScenario				->	Finalize( itsMktData );
	}

	if (computeHedge)
	{
  		theScenario -> InitPosition();

		theScenario -> OrderData(Container);
		tmpFunctor.SetMatrix( Container );
		itsFunctor.InsertData( theScenario->GetScenarioKey(), tmpFunctor	);

		tmpFunctor.SetStringVector( ARM_StringVectorPtr( expLabel ) );
		itsFunctor.InsertData( string("EXPIRY"), tmpFunctor	);

		tmpFunctor.SetStringVectorTrans(ARM_StringVectorPtr( tenLabel ) );
		itsFunctor.InsertData( string("TENOR"), tmpFunctor	);
	}
}


void ARM_Hedge::Apply_ND_Scenario( ARM_Scenario* theScenario, int nbShift, bool computeHedge){

	string						key;
	ARM_IntVector				tmpPos;
	double						tmp;
	double						tmpPriceRef;
	if (computeHedge)
		tmpPriceRef =	itsPricer->Price();
	int							dim				=	theScenario -> GetDim();
	vector<ARM_StringVector>	tmpPlotsLab		=	theScenario -> GetPlotsLabel();
	map<string, ScenarioType>	tmpPriceCritera	=	theScenario -> GetPriceCritera();

	ARM_GP_VectorPtr			container	( new ARM_GP_Vector );
	vector<ARM_StringVector*>	expLabel(dim);
	ARM_GramFctorArg			tmpFunctor;

	for(int j=0; j<dim; j++)		expLabel[j]= new ARM_StringVector;
	
	theScenario -> InitPosition();
	for(int i=0; (i<theScenario->GetNbShift()) && (i<nbShift);i++){

		theScenario	->	ShiftNextPlot();
		key.resize(0);
		tmpPos	=	theScenario->GetPosition();
		for(int j=0; j<dim; j++){
			if ( key.size()!=0 )	key = key + string(":");
			key = key + tmpPlotsLab[j][tmpPos[j] ];
		}

		theScenario	->	ShiftMktData( itsMktData );

		if (computeHedge)
		{
			switch ( tmpPriceCritera[ key ] ){
				
				case REF:{
					tmpPriceRef	= itsPricer->Price();
					theScenario	->	Finalize( itsMktData );
					break;
				}
				case VALUE:{
					container	->	push_back(itsPricer->Price()-tmpPriceRef );
					for(int j=0; j<dim; j++)
						expLabel[j]->push_back(tmpPlotsLab[j][tmpPos[j] ] );
						
					break;
				}
				case REF_VALUE:{
					tmp			=	itsPricer->Price();
					container	->	push_back(tmp-tmpPriceRef );
					tmpPriceRef = tmp;
					for(int j=0; j<dim; j++)
						expLabel[j]->push_back(tmpPlotsLab[j][tmpPos[j] ] );
					break;
				}
			}
		}

		if (computeHedge || ((i < nbShift-1) && (i < theScenario->GetNbShift()-1)))
			theScenario	->	Finalize( itsMktData );
	}
	theScenario -> InitPosition();

	theScenario -> OrderData(container);
	tmpFunctor.SetVector( container );
	key = theScenario->GetScenarioKey();
	itsFunctor.InsertData( key, tmpFunctor	);

	for( j=0; j<dim; j++){
		tmpFunctor.SetStringVector( ARM_StringVectorPtr( expLabel[j] ) );
		key = string("EXPIRY_")+ dynamic_cast<ARM_ND_Scenario*> (theScenario)->GetSubScenario()[j];
		itsFunctor.InsertData( key, tmpFunctor	);
	}
}


string ARM_Hedge::toString(	const string& indent, const string& nextIndent) const
{
	return itsFunctor.toString();
}


CC_END_NAMESPACE()













