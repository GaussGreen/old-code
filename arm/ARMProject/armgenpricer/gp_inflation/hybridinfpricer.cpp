
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Hybrid Inf Ir Pricer														 *
 *																							 *
 *			This class builds a hybrid inf ir leg from swap legs and  Inf Swap Leg			 *									 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: April, 2nd 2007														 *																											 *
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#include "gpinflation/hybridinfpricer.h"
#include "gpinflation/infpayoffvisitor.h"
#include "gpinflation/infhybridpayoffvisitor.h"
#include "gpinflation/infoptionspreadvisitor.h"
#include "gpinflation/infdoubledigitalvisitor.h"

#include <gpinflation/infbilogvisitor.h>

#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"

#pragma warning(disable : 4786) 

#ifndef _INGPINFLATION_HYBRIDINFLEGPRICER_CPP
#define _INGPINFLATION_HYBRIDINFLEGPRICER_CPP

CC_BEGIN_NAMESPACE( ARM )

 

/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: ARM_HybridInfLegPricer
///	Returns: void
///	Action : constructor

/////////////////////////////////////////////////////

 
ARM_HybridInfLegPricer::ARM_HybridInfLegPricer( ARM_HybridInfIrLeg*		security, 
												ARM_InfHybridPayOff*	payOff, 
												ARM_InfModel*			model ){

	itsPayOff		= CreateClonedPtr( payOff	);
	itsSecurity		= CreateClonedPtr( security	);
	itsModel		= CreateClonedPtr( model	);

	vector<string> vec;
	ARM_GP_MapIdx idx = itsSecurity->GetIdx();
	for (idxIter  iter= idx.begin(); iter!= idx.end(); iter++)
		vec.push_back(iter->second.GetName() );
	vec.push_back("NO");
		
	itsPayOff->ValidateSecurity ( vec );
	isComputed = false;
}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: ARM_HybridInfLegPricer
///	Returns: void
///	Action : copy constructor

/////////////////////////////////////////////////////

ARM_HybridInfLegPricer::ARM_HybridInfLegPricer(	const ARM_HybridInfLegPricer & rhs):isComputed(false){
	itsPayOff		= CreateClonedPtr( &*rhs.itsPayOff		);
	itsSecurity		= CreateClonedPtr( &*rhs.itsSecurity	);
	itsModel		= CreateClonedPtr( &*rhs.itsModel		);

}

/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: toString
///	Returns: string
///	Action : Viewer

/////////////////////////////////////////////////////

string ARM_HybridInfLegPricer::toString(	const string& indent, const string& nextIndent) const{

	CC_Ostringstream	os;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(57)<<"HYBRIID INFLATION PRICER\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<< "SECURITY DESCRIPTION\n";
	os	<< "===================="<<"\n";

	os	<< itsSecurity->ViewGlob();

	os	<< "\n\n";	
	os	<<"PAYOFF DESCRIPTION\n";
	os	<<"=================="<<"\n";

	os	<< itsPayOff->toString(	indent,  nextIndent);


	os	<< "\n\n";	
	os	<<"MODEL DESCRIPTION\n";
	os	<<"================="<<"\n";

	os	<< itsModel->toString(	indent,  nextIndent);

	return os.str();

}

/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: Map
///	Returns: void
///	Action : map the models with payoff

/////////////////////////////////////////////////////

void ARM_HybridInfLegPricer::Map(ARM_InfPayOffValuePtr& pValue , ARM_InfModelValuePtr& mValue){

	ARM_InfModelVisitor modelVisitor;
	if ( !&*itsModel) itsModel = ARM_CountedPtr<ARM_InfModel>(new ARM_InfNoModel);
	itsModel->Accept(modelVisitor);
	mValue = ARM_InfModelValuePtr (modelVisitor.GetInfModelValue() );

	if (	!dynamic_cast<ARM_InfBiLog*> (&*itsModel ) && dynamic_cast<ARM_InfHybridCap*>	(&*itsPayOff) ){
		ARM_MAP_VolType tmp = itsPayOff->GetVolType();
		itsPayOff=ARM_CountedPtr<ARM_InfHybridPayOff> ( new ARM_InfHybridPayOff(itsPayOff->GetOptCurve(),itsPayOff->GetOptCurve() ) );
		itsPayOff->SetVolType(tmp);
	}

	if (	!dynamic_cast<ARM_InfBiLog*> (&*itsModel ) && dynamic_cast<ARM_InfHybridDigit*> (&*itsPayOff) ){

		ARM_MAP_VolType tmp = itsPayOff->GetVolType();
		itsPayOff=ARM_CountedPtr<ARM_InfHybridPayOff> ( new ARM_InfHybridPayOff(itsPayOff->GetCpnCurve(),itsPayOff->GetOptCurve() ) );
		itsPayOff->SetVolType(tmp);
	}


	ARM_InfPayOffVisitor payOffVisitor;
	itsPayOff->Accept(payOffVisitor);
	pValue = ARM_InfPayOffValuePtr (payOffVisitor.GetInfPayOffValue() );

}

/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: Cpt
///	Returns: void
///	Action : return price and its ventilation

/////////////////////////////////////////////////////

void ARM_HybridInfLegPricer::Compute(){

	if ( !isComputed ) CptFlows( );
	isComputed = true;
}


void ARM_HybridInfLegPricer::CptFlows(  ){

	ARM_InfPayOffValuePtr pValue;
	ARM_InfModelValuePtr  mValue;	

	Map(pValue, mValue);

	ARM_GP_MapIdx	tmpIdx  = itsSecurity->GetIdx	( );

	ARM_GP_VecPtr	intTerm = itsSecurity->GetIntTerm( );
	ARM_GP_VecPtr	disFact = itsSecurity->GetDisFact( );
	ARM_GP_VecPtr	resDate = itsSecurity->GetDisFact( );

	ARM_MAP_Double	fwdArg;
	ARM_MAP_PairDb	corArg;
	ARM_MAP_VolPar	volArg;

	double				tmp;
	double				dat;
	double				sum = 0.0;
	int					dim = intTerm->size();

	ARM_GP_VectorPtr	rat = ARM_GP_VectorPtr(new ARM_GP_Vector);
	ARM_GP_VectorPtr	dif = ARM_GP_VectorPtr(new ARM_GP_Vector);
	ARM_GP_VectorPtr	val = ARM_GP_VectorPtr(new ARM_GP_Vector);
	ARM_GP_VectorPtr	cor = ARM_GP_VectorPtr(new ARM_GP_Vector);

	for( int i	= 0; i< dim ; i++){
		fwdArg	=	CptFwd(i);
		volArg	=	CptVol(i);
		corArg	=	CptCor(i);

		dat		=	intTerm->Elt(i);
		rat		->	push_back( dat );
		dif		->	push_back( disFact->Elt(i) );

		tmp		=	(*mValue)(pValue, dat, fwdArg, corArg, volArg)/RateBase;
		val		->	push_back( tmp );
		tmp		*=  dat*disFact->Elt(i);

		sum		+=tmp;
	}

	ARM_GramFctorArg*	tmpFunctor = new ARM_GramFctorArg;

	tmpFunctor->SetDouble( sum );
	itsFunctor.InsertData( "PRICE", *tmpFunctor	);

	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}
	tmpFunctor = new ARM_GramFctorArg;
	tmpFunctor->SetVector( rat );
	itsFunctor.InsertData( "RATIO", *tmpFunctor	);

	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}
	tmpFunctor = new ARM_GramFctorArg;
	tmpFunctor->SetVector( dif );
	itsFunctor.InsertData( "DISCOUNT_FACTOR", *tmpFunctor	);

	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}
	tmpFunctor = new ARM_GramFctorArg;
	tmpFunctor->SetVector( val );
	itsFunctor.InsertData( "FLOWS", *tmpFunctor	);

	if ( tmpFunctor)	{ delete tmpFunctor;	tmpFunctor	=NULL;}

}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: CptFwd
///	Returns: ARM_MAP_Double
///	Action : return the Fwd for each index

/////////////////////////////////////////////////////


ARM_MAP_Double ARM_HybridInfLegPricer::CptFwd( const int & i ){

	ARM_GP_MapIdx	mapIdx  = itsSecurity->GetIdx	( );

	ARM_GP_MapPtr	mapBeg	= itsSecurity->GetBeg( );
	ARM_GP_MapPtr	mapEnd	= itsSecurity->GetEnd( );
	ARM_GP_MapPtr	mapPay	= itsSecurity->GetPay( );
	
	ARM_MAP_Double	fwdArg;

	string			str;
	ARM_InfIrIndex	idx;
	double			beg, end, pay, fwd;

	for( idxIter iter = mapIdx.begin(); iter!= mapIdx.end(); iter++){
		idx	=	iter->second;
		str	=	idx.GetName();
		beg	=	mapBeg[iter->first]->Elt(i);
		end	=	mapEnd[iter->first]->Elt(i);
		pay	=	mapPay[iter->first]->Elt(i);
		fwd	=	idx.CptFwd(beg,end,pay);
	
		if ( idx.itsType == IN )
			fwd  = fwd/100+1.0;
		fwdArg.insert( pair<string,double> (str,fwd) );
	}
	return fwdArg;
}




/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: CptVol
///	Returns: ARM_MAP_VolPar
///	Action : return the ARM_VolParam for each index

/////////////////////////////////////////////////////


ARM_MAP_VolPar ARM_HybridInfLegPricer::CptVol( const int & i ){


	ARM_GP_MapIdx	mapIdx	= itsSecurity->GetIdx	( );

	ARM_GP_MapPtr	mapRes	= itsSecurity->GetRes( );
	ARM_GP_MapPtr	mapTen	= itsSecurity->GetTen( );
	
	ARM_MAP_VolPar	volArg;

	string			str;
	ARM_InfIrIndex	idx;
	double			res, ten;

	for( idxIter iter = mapIdx.begin(); iter!= mapIdx.end(); iter++){
		idx	=	iter->second;
			
		str	=	idx.GetName();
		res	=	mapRes[iter->first]->Elt(i);
		ten	=	mapTen[iter->first]->Elt(i);

		ARM_VolParam vol( idx, res, ten );
		volArg.insert( pair<string,ARM_VolParam> (str,vol) );
	}

	return volArg;
}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: CptCor
///	Returns: ARM_MAP_PairDb
///	Action : return the correlation for each pair of indexes 

/////////////////////////////////////////////////////


ARM_MAP_PairDb ARM_HybridInfLegPricer::CptCor( const int & i ){


	ARM_GP_MapIdx	mapIdx  = itsSecurity->GetIdx	( );

	ARM_GP_MapPtr	mapRes	= itsSecurity->GetRes( );
	ARM_GP_MapPtr	mapTen	= itsSecurity->GetTen( );
	
	ARM_InfIrCorrel	corObj	= itsSecurity->GetCorObj( );	

	ARM_MAP_PairDb	corArg;

	string			str1, str2;
	ARM_InfIrIndex	idx1, idx2;
	double			res,  ten, cor;

	for( idxIter iter = mapIdx.begin(); iter!= mapIdx.end(); iter++){
		for( idxIter jter = mapIdx.begin(); jter!= mapIdx.end(); jter++){
			idx1	=	iter->second;
			idx2	=	jter->second;
				
			str1	=	idx1.GetName();
			str2	=	idx2.GetName();
			res		=	mapRes[iter->first]->Elt(i);
			ten		=	mapTen[jter->first]->Elt(i);

			pair<string,string> p(str1,str2);

			if ( idx1.itsType == IR && idx2.itsType == IN){
				res  =	mapRes[jter->first]->Elt(i);
				ten  =	mapTen[iter->first]->Elt(i);
				str1 =	idx2.GetName();
				str2 =	idx1.GetName();
			}

			if( str1 == str2)
				cor = 1.0;
			else
				cor  =	corObj.GetCorrel( str1, str2, res, ten);

			corArg.insert( pair<pair<string,string>, double>	(p,cor) );
		}
	}
	return corArg;
}


/////////////////////////////////////////////////////

///	Class  : ARM_HybridInfLegPricer
///	Routine: ImplicitCorrel
///	Returns: double
///	Action : return the equivalent correlation for the bilog model 

/////////////////////////////////////////////////////

class CptPrice: public ARM_GP::UnaryFunc<double,double>{

public:
	CptPrice(	 const double					& target,
				 ARM_InfPayOffValuePtr			& pValue,
				 ARM_InfModelValuePtr			& mValue,
				 const double					& lag,
				 const ARM_MAP_Double			& fwdArg,
				 const ARM_MAP_PairDb			& corArg,
				 const ARM_MAP_VolPar			& volArg):	itsTarget(target),
													itsPValue(pValue),
													itsMValue(mValue),
													itsLag(lag),
													itsFwdArg(fwdArg),
													itsCorArg(corArg),
													itsVolArg(volArg){}
	virtual double operator() (double x) const;


private:
	double					itsTarget;
	ARM_InfPayOffValuePtr	itsPValue;
	ARM_InfModelValuePtr	itsMValue;
	double					itsLag;
	ARM_MAP_Double			itsFwdArg;
	ARM_MAP_PairDb			itsCorArg;
	ARM_MAP_VolPar			itsVolArg;
};

double CptPrice::operator() (double x) const{ 

	ARM_MAP_PairDb cor = itsCorArg;
	ARM_MAP_PairDb::iterator it;
	for ( it = cor.begin(); it!=cor.end(); it++){
		if( it->second!=1){
			it->second = x;
		}
	}
	double tmp = (*itsMValue)(itsPValue,  itsLag, itsFwdArg, cor, itsVolArg);
	return tmp - itsTarget;		
}



CC_END_NAMESPACE()
#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/



