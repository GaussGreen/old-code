
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Corridor Leg															 *
 *																							 *
 *			This class builds a corridor leg from swap legs	and inherits from Inf Swap Leg	 *									 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: October, 25th 2006													 *																											 *
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

#ifndef _INGPINFLATION_INFCORRIDORLEG_H
#define _INGPINFLATION_INFCORRIDORLEG_H

#define LAG			13

#define RateBase	CC_NS( ARM_Constants, rateBase )
#define VolBase		CC_NS( ARM_Constants, volBase )


#include <gpbase/gplinalgtypedef.h>
#include <gpbase/gpvector.h>
#include <gpbase/rootobject.h>
#include <gpbase/typedef.h>
#include <inst/swapleg.h>
#include <gpbase/curvetypedef.h>
#include <gpbase/stringconvert.h>
#include <map>

class ARM_SwapLeg;
class ARM_Date;
class ARM_IrIndex;
class ARM_Currency;

CC_BEGIN_NAMESPACE( ARM )

class ARM_InfLeg;

class ARM_InfCorridorLeg: public ARM_Security {

public:
	ARM_InfCorridorLeg( ){ }
	ARM_InfCorridorLeg(		map<string,ARM_Date>	&	mDate, 
							map<string,string>		&	mString,
							map<string,ARM_Curve*>	&	mCurve,
							map<string,int>			&	mInt,
							map<string,double>		&	mDouble);



	ARM_InfCorridorLeg						( const ARM_InfCorridorLeg &	);
	virtual ARM_Object*	Clone				( ) const =0;
	virtual ~ARM_InfCorridorLeg				( ) { };

	virtual void			View			( char* id = NULL, FILE* ficOut = NULL);

	virtual	void			BuildSchedule	( );

	void					InitSchedule	(	ARM_Curve	*	notional,
												ARM_Curve	*	irLeverage,	
												ARM_Curve	*	infLeverage,	
												ARM_Curve	*	constant,						
												ARM_Curve	*	multipleUp,						
												ARM_Curve	*	multipleDown,						
												ARM_Curve	*	rangeUp,						
												ARM_Curve	*	rangeDown	);


	virtual void			CptExpFwdRates	( );
	virtual double			ComputePrice	( int mode=0 );

	ARM_CountedPtr<ARM_InfLeg>		GetInfLeg		( ) const {	return itsInfLeg;	}	
	ARM_CountedPtr<ARM_SwapLeg>		GetIrLeg		( )	const {	return itsIrLeg;	}
	
protected:

	virtual void			InitModel			( )=0;
	virtual void			CptExpFwdCorrel		( )=0;
	virtual	void			CptExpFwdInfVol		( )=0;
	virtual	void			CptExpFwdIrVol		( )=0;

	virtual string			ViewModelFeatures	( )=0;
	virtual string			ViewPricerFeatures	( )=0;
	virtual	void			StripCompute		( const int & )=0;

private:
	string					ViewGeneralFeatures		( );
	string					ViewInflationFeatures	( );
	string					ViewRateFeatures		( );
	string					ViewDealFeatures		( );


	template< typename A = string, typename B = string>
	struct Ostream{
		string	ToStream( const A & str1, const B & str2)	{
			CC_Ostringstream	os;
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
			os <<CC_NS(std,right)<<str1;
			os <<"  :  "<< str2<<"\n";
			return os.str();	
		};
	};

public:

	ARM_GP_VectorPtr  itsFlowStartDates;		
	ARM_GP_VectorPtr  itsFlowEndDates;		
	ARM_GP_VectorPtr  itsPaymentDates;		 
	ARM_GP_VectorPtr  itsInterestDays;		
	ARM_GP_VectorPtr  itsInterestTerms;		
	ARM_GP_VectorPtr  itsDiscFactor;			

	ARM_GP_VectorPtr  itsFwdRateStartDates;	
	ARM_GP_VectorPtr  itsFwdRateEndDates	;	
	ARM_GP_VectorPtr  itsIrResetDates;		 
	ARM_GP_VectorPtr  itsIrResetTerms;
	ARM_GP_VectorPtr  itsIrInterTerms;	
	ARM_GP_VectorPtr  itsFwdIrRates;			
	ARM_GP_VectorPtr  itsIrAdjConv;				

	ARM_GP_VectorPtr  itsNumResetDates;		
	ARM_GP_VectorPtr  itsDemResetDates;		
	ARM_GP_VectorPtr  itsInfResetDates;		 
	ARM_GP_VectorPtr  itsInfResetTerms;
	ARM_GP_VectorPtr  itsInfInterTerms;	
	ARM_GP_VectorPtr  itsNumCPIRates;		
	ARM_GP_VectorPtr  itsDemCPIRates;	
	ARM_GP_VectorPtr  itsNumPublishDates;	
	ARM_GP_VectorPtr  itsDemPublishDates;
	ARM_GP_VectorPtr  itsFwdInfRates;			
	ARM_GP_VectorPtr  itsInfAdjConv;				

	ARM_GP_VectorPtr  itsIrLeverage;
	ARM_GP_VectorPtr  itsInfLeverage;	
	ARM_GP_VectorPtr  itsConstant;			
	ARM_GP_VectorPtr  itsUpMultiple;			
	ARM_GP_VectorPtr  itsDownMultiple;		
	ARM_GP_VectorPtr  itsRangeUp;			
	ARM_GP_VectorPtr  itsRangeDown;			
	ARM_GP_VectorPtr  itsNotional;			

	ARM_GP_VectorPtr  itsUpPrice;
	ARM_GP_VectorPtr  itsDownPrice;

protected:

	ARM_CountedPtr<ARM_InfLeg>		itsInfLeg;
	ARM_CountedPtr<ARM_SwapLeg>		itsIrLeg;

	ARM_Date		itsAsOfDate;
	ARM_Date		itsStartDate;
	ARM_Date		itsEndDate;

	int				itsNbFlows;
	int				itsRecOrPay;
	int				isIrCritera;
	int				isModulable;
	bool			isComputed;

	double			itsEpsilon;
	double			itsNbGaussLeg;
	double			itsPtGaussLeg;
};

typedef ARM_CountedPtr<ARM_VolCurve> ARM_VolCurvePtr;

class ARM_Copula{
	
public:
	ARM_Copula( ){ }
	ARM_Copula( ARM_VolCurvePtr, ARM_VolCurvePtr, map<string,double> &);
	ARM_Copula( const ARM_Copula & );
	virtual ARM_Copula*	Clone	( ) const	{ return new ARM_Copula (*this);	}
	virtual ARM_Copula*	Clone	( )			{ return new ARM_Copula (*this);	}
	ASSIGN_OPERATOR	( ARM_Copula );
	virtual ~ARM_Copula( ){ };

	void	SetMarginal	( double (*f1)(ARM_GP_Vector, ARM_VolCurvePtr, double), double (*f2)(ARM_GP_Vector, ARM_VolCurvePtr, double) ) ;
	void	SetPayoff	( const ARM_GP_Vector &);
	void	SetVolParam ( const ARM_GP_Vector &, const ARM_GP_Vector & );
	void	InitGaussLeg( );	
	virtual double	Compute	( ){ return 0.0; }

public:

	double itsStrike1;
	double itsStrike2;
	
protected:
	double itsFwd1;
	double itsTerms1;
	double itsLev1;

	double itsFwd2;
	double itsTerms2;
	double itsLev2;

	double itsCorrel;
	double itsStrike;
	double itsMultiple;
	double itsConstant;

	double itsEpsilon;
	double itsNbGaussLeg;
	double itsPtGaussLeg;
	
	ARM_GP_Vector	itsParam1;
	ARM_VolCurvePtr	itsVol1;
	double			(*CptVol1)(ARM_GP_Vector, ARM_VolCurvePtr, double);

	ARM_GP_Vector	itsParam2;
	ARM_VolCurvePtr	itsVol2;
	double			(*CptVol2)(ARM_GP_Vector, ARM_VolCurvePtr, double);

	ARM_GP_Vector	itsPosit;
	ARM_GP_Vector	itsWeigth;

};


class ARM_NaiveCopula: public ARM_Copula{

public:
	ARM_NaiveCopula( ):ARM_Copula( ){ }
	ARM_NaiveCopula( ARM_VolCurvePtr vol1, ARM_VolCurvePtr vol2, map<string,double> & param):ARM_Copula(vol1, vol2, param){ }
	ARM_NaiveCopula( const ARM_NaiveCopula & cop):ARM_Copula(cop) { }
	virtual ARM_Copula*	Clone	( ) const	{ return new ARM_NaiveCopula (*this);	}
	virtual ARM_Copula*	Clone	( )			{ return new ARM_NaiveCopula (*this);	}
	ASSIGN_OPERATOR				( ARM_NaiveCopula );

	double	Compute		( ) { return -999;}
	static string	ViewModelFeatures(){ return "NAIVE COPULA"; }
};


class ARM_GaussCopula: public ARM_Copula{

public:
	ARM_GaussCopula( ):ARM_Copula( ){ }
	ARM_GaussCopula( ARM_VolCurvePtr vol1, ARM_VolCurvePtr vol2, map<string,double> & param):ARM_Copula(vol1, vol2, param){ }
	ARM_GaussCopula( const ARM_GaussCopula & cop):ARM_Copula(cop) { }
	virtual ARM_Copula*	Clone	( ) const	{ return new ARM_GaussCopula (*this);	}
	virtual ARM_Copula*	Clone	( )			{ return new ARM_GaussCopula (*this);	}
	ASSIGN_OPERATOR				( ARM_GaussCopula );

	double Compute		( );
	static string	ViewModelFeatures(){ return "GAUSS COPULA"; }

public:
	static int nbStrike;
};

template< class T=ARM_Copula > 
class ARM_InfCorridorLegWithModel : public ARM_InfCorridorLeg{

public:
	ARM_InfCorridorLegWithModel(	map<string,ARM_Date>	&	mDate, 
									map<string,string>		&	mString,
									map<string,ARM_Curve*>	&	mCurve,
									map<string,int>			&	mInt,
									map<string,double>		&	mDouble);

	ARM_InfCorridorLegWithModel	( const ARM_InfCorridorLegWithModel & );
	virtual ~ARM_InfCorridorLegWithModel( ){ };
	virtual ARM_Object*	Clone	( ) const		{ return new ARM_InfCorridorLegWithModel<T> (*this);	}
	virtual ARM_Object*	Clone	( )				{ return new ARM_InfCorridorLegWithModel<T> (*this);	}
	ASSIGN_OPERATOR				( ARM_InfCorridorLegWithModel );

	virtual	void	BuildSchedule	( );
	virtual void	InitModel		( );
	virtual void	InitInfVolParam	( );
	virtual void	InitIrVolParam	( );

	virtual	void	CptExpFwdInfVol	( );
	virtual	void	CptExpFwdIrVol  ( );
	virtual void	CptExpFwdCorrel	( );

protected:

	static double	GetInfVol		( ARM_GP_Vector, ARM_VolCurvePtr, double );
	static double	GetIrVol		( ARM_GP_Vector, ARM_VolCurvePtr, double );
		
	virtual	void	StripCompute		( const int & );
	virtual string	ViewModelFeatures	( ){ return T::ViewModelFeatures();	}
	virtual string	ViewPricerFeatures	( );

	virtual ARM_GP_Vector	GetInfParam( const int & ); 
	virtual ARM_GP_Vector	GetIrParam( const int & );
public:

	ARM_GP_VectorPtr  itsFwdInfVolUp;			
	ARM_GP_VectorPtr  itsFwdInfVolDown;		
	ARM_GP_VectorPtr  itsFwdIrVolUp;				
	ARM_GP_VectorPtr  itsFwdIrVolDown;
	ARM_GP_VectorPtr  itsFwdCorrel;


private:
	ARM_GP_VectorPtr	itsInfTenor;
	ARM_GP_VectorPtr	itsInfRenor;
	ARM_GP_VectorPtr	itsInfExper;
	ARM_GP_VectorPtr	itsInfInter;
	ARM_GP_VectorPtr	itsIrTenor;
	ARM_GP_VectorPtr	itsIrExper;
	ARM_GP_VectorPtr	itsIrForwd;

	ARM_VolCurvePtr		itsInfVol;
	ARM_VolCurvePtr		itsIrVol;

	ARM_CountedPtr<T>	itsModel;
};






		/********************************************************************************/
		/*																				*/
		/*								SOURCE TEMPLATE									*/
		/*																				*/
		/********************************************************************************/



//======>	ARM_InfCorridorLegWithModel
template< class T >
ARM_InfCorridorLegWithModel<T>::ARM_InfCorridorLegWithModel(	map<string,ARM_Date>	&	mDate, 
																map<string,string>		&	mString,
																map<string,ARM_Curve*>	&	mCurve,
																map<string,int>			&	mInt,
																map<string,double>		&	mDouble):	
																ARM_InfCorridorLeg(	mDate,	mString, mCurve,	mInt, mDouble	)
																{	BuildSchedule( );	}

//======>	ARM_InfCorridorLegWithModel

template< class T >
ARM_InfCorridorLegWithModel<T>::ARM_InfCorridorLegWithModel( const ARM_InfCorridorLegWithModel & leg):ARM_InfCorridorLeg( leg ){

	itsFwdInfVolUp		= CreateClonedPtr ( &*	leg.itsFwdInfVolUp		);
	itsFwdInfVolDown	= CreateClonedPtr ( &*	leg.itsFwdInfVolDown	);
	itsFwdIrVolUp		= CreateClonedPtr ( &*	leg.itsFwdIrVolUp		);
	itsFwdIrVolDown		= CreateClonedPtr ( &*	leg.itsFwdIrVolDown		);
	itsFwdCorrel		= CreateClonedPtr ( &*	leg.itsFwdCorrel		);

	itsInfTenor			= CreateClonedPtr ( &*	leg.itsInfTenor			);
	itsInfRenor			= CreateClonedPtr ( &*	leg.itsInfRenor			);
	itsInfExper			= CreateClonedPtr ( &*	leg.itsInfExper			);
	itsInfInter			= CreateClonedPtr ( &*	leg.itsInfInter			);
	itsIrTenor			= CreateClonedPtr ( &*	leg.itsIrTenor			);
	itsIrExper			= CreateClonedPtr ( &*	leg.itsIrExper			);
	itsIrForwd 			= CreateClonedPtr ( &*	leg.itsIrForwd			);

	itsInfVol			= CreateClonedPtr ( &*	leg.itsInfVol			);
	itsIrVol			= CreateClonedPtr ( &*	leg.itsIrVol			);
	itsModel			= CreateClonedPtr ( &*  leg.itsModel			);
}

//======>	BuildSchedule

template< class T >
void	ARM_InfCorridorLegWithModel<T>::BuildSchedule	( ){

	itsFwdInfVolUp	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			
	itsFwdInfVolDown= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsFwdIrVolUp	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsFwdIrVolDown	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsFwdCorrel	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );

	itsInfTenor		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsInfRenor		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsInfExper		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsInfInter		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsIrTenor		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsIrExper		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsIrForwd		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
}

//======>	InitModel

template< class T >
void ARM_InfCorridorLegWithModel<T>::InitModel(){

	if ( !dynamic_cast<ARM_InfBSModel*	>(GetModel()) )
		ARMTHROW(ERR_INVALID_ARGUMENT,"Model is of bad type... Its should contain a Swap Leg and a Inflation leg "); 
	if ( !dynamic_cast<ARM_IRIndex*>(itsInfLeg->GetIRIndex()) )
		ARMTHROW(ERR_INVALID_ARGUMENT,"Model is of bad type... Its should contain a Swap Leg and a Inflation leg "); 
	if ( !dynamic_cast<ARM_IRIndex*>(itsIrLeg->GetIRIndex()) )
		ARMTHROW(ERR_INVALID_ARGUMENT,"Model is of bad type... Its should contain a Swap Leg and a Inflation leg "); 

	ARM_InfBSModel*	pricingMod	= dynamic_cast<ARM_InfBSModel*	>(GetModel());
	ARM_IRIndex*	tmpIdx		= dynamic_cast<ARM_IRIndex*>(itsInfLeg->GetIRIndex() );
	double			InitDate	= GetModel()->GetStartDate().GetJulian();
	
	for ( int i =0; i< itsNbFlows; i++){
		itsDemPublishDates->Elt(i) = pricingMod -> GetModelDateWPublishLag	( ARM_Date( itsDemResetDates->Elt(i) ),	(ARM_InfIdx*)	tmpIdx	).GetJulian();
		itsNumPublishDates->Elt(i) = pricingMod -> GetModelDateWPublishLag	( ARM_Date( itsNumResetDates->Elt(i) ),	(ARM_InfIdx*)	tmpIdx	).GetJulian();
		itsInfResetTerms  ->Elt(i) = ( itsNumPublishDates ->Elt(i) -  InitDate )/K_YEAR_LEN;
		itsIrResetTerms	  ->Elt(i) = ( itsIrLeg ->GetResetDates() ->Elt(i) -  InitDate )/K_YEAR_LEN;

		if ( fabs(itsInfResetTerms ->Elt(i) - itsIrResetTerms ->Elt(i) ) >0.1) 
			ARMTHROW(ERR_INVALID_ARGUMENT,"The Ir and Inf Schedule should be synchronized. The difference is more than 1 month"); 
	}

	InitInfVolParam();
	InitIrVolParam();

	map<string, double> modelParam;
	modelParam.insert(pair<string, double> ("epsilon",		itsEpsilon		) );
	modelParam.insert(pair<string, double> ("nbGaussLeg",	itsNbGaussLeg	) );
	modelParam.insert(pair<string, double> ("ptGaussLeg",	itsPtGaussLeg	) );
	
	if ( isIrCritera ) {
		itsModel = ARM_CountedPtr<T> ( new T(itsIrVol,itsInfVol, modelParam) );
		itsModel -> SetMarginal( GetIrVol, GetInfVol);
	}
	else{
		itsModel = ARM_CountedPtr<T> (  new T(itsInfVol,itsIrVol, modelParam) );
		itsModel -> SetMarginal( GetInfVol, GetIrVol);
	}
}

//======>	CptExpFwdCorrel

template< class T >
void	ARM_InfCorridorLegWithModel<T>::CptExpFwdCorrel(){

	double			tenor;
	double			TiForMaturity;
	ARM_Date		tmpDate;

	ARM_InfBSModel*	pricingMod		= dynamic_cast<ARM_InfBSModel*	>(GetModel());
	ARM_IRIndex*	infIdx			= dynamic_cast<ARM_IRIndex*>(itsInfLeg->GetIRIndex() );
	ARM_IRIndex*	irIdx			= dynamic_cast<ARM_IRIndex*>(itsIrLeg->GetIRIndex() );
		
	double			dayCount		= pricingMod -> GetInfFwdCurv() -> GetMonthlyInterpType();
	ARM_Date		lastKnownDate	= pricingMod -> GetVolatility() -> GetLastKnownDate();
	ARM_Date		modelAsOf		= GetModel() -> GetStartDate();

	for( int i=0 ; i<itsNbFlows; i++){
		tmpDate				= pricingMod->GetModelDateWPublishLag		( ARM_Date( itsNumResetDates->Elt(i) ),	(ARM_InfIdx*)	infIdx	);
		TiForMaturity		= CountYearsWithoutException				( dayCount, 	modelAsOf,		tmpDate );
		tenor				= irIdx -> GetYearTerm();

		if( TiForMaturity <= 0.0 )	{
			tenor				+=	TiForMaturity;
			TiForMaturity		=	StringMaturityToYearTerm( "1d" );
			}

		itsFwdCorrel->Elt(i) = pricingMod->GetInfIRCorrel( TiForMaturity, tenor, (ARM_InfIdx*) infIdx, irIdx, "INF_YoY/IR_FWD" );
	}
}

//======>	InitIrVolParam

template< class T >
void	ARM_InfCorridorLegWithModel<T>::InitIrVolParam	( ){

	ARM_InfBSModel*	pricingMod = dynamic_cast<ARM_InfBSModel*	>(GetModel());
	itsIrVol	=   CreateClonedPtr ( &* pricingMod->GetIRModel()->GetVolatility()  );

	for ( int i= 0; i< itsNbFlows; i++){
		itsIrTenor -> Elt(i) = itsIrLeg -> GetIRIndex() -> GetYearTerm();
		itsIrExper -> Elt(i) = itsIrResetTerms->Elt(i);
		itsIrForwd -> Elt(i) = itsFwdIrRates->Elt(i);	
	}
}

//======>	InitInfVolParam

template< class T >
void	ARM_InfCorridorLegWithModel<T>::InitInfVolParam	( ){
	
	double			renormFactor;
	double 			timeToStart;	
	double			timeToExpiry;
	double			TiVolLookUp;
	double			TjVolLookUp;	
	double 			tenor;		
	double			TjForMaturity;
	double			TiForMaturity;

	ARM_Date		numDate;
	ARM_Date		demDate;


	ARM_InfBSModel*	pricingMod		= dynamic_cast<ARM_InfBSModel*	>(GetModel());
	ARM_IRIndex*	tmpIdx			= dynamic_cast<ARM_IRIndex*>(itsInfLeg->GetIRIndex() );
	
	int				dayCount		= pricingMod -> GetInfFwdCurv() -> GetMonthlyInterpType();
	ARM_Date		lastKnownDate	= pricingMod -> GetVolatility() -> GetLastKnownDate();
	ARM_Date		modelAsOf		= GetModel() -> GetStartDate();
	
	itsInfVol = CreateClonedPtr ( &* pricingMod -> GetVolatility() );
				
	for ( int i=0; i< itsNbFlows; i++){
		renormFactor	= 1.0;
		numDate			= ARM_Date( itsNumResetDates->Elt(i) );
		demDate			= ARM_Date( itsDemResetDates->Elt(i) );
		timeToStart		= pricingMod -> GetModelTimeWPublishLag	( demDate,	(ARM_InfIdx*)	tmpIdx	);
		timeToExpiry	= pricingMod -> GetModelTimeWPublishLag	( numDate,  (ARM_InfIdx*)	tmpIdx	);	
		TjForMaturity	= CountYearsWithoutException			( dayCount,	modelAsOf,	itsDemPublishDates->Elt(i) );
		TiForMaturity	= CountYearsWithoutException			( dayCount,	modelAsOf,	itsNumPublishDates->Elt(i) );
		TiVolLookUp		= CountYearsWithoutException			( dayCount,	lastKnownDate,	numDate ); 
		TjVolLookUp		= CountYearsWithoutException			( dayCount,	lastKnownDate,	demDate ); 

		tenor			= TiVolLookUp	- TjVolLookUp;
		if( timeToStart < 0.0 )	
			renormFactor				=	itsDemCPIRates -> Elt(i)/ pricingMod -> GetCPIIndexValue();
			
		if (TiForMaturity > 0.0	){
			if( TjForMaturity	<= 0.0 ){
				tenor			+= TjForMaturity;
				TjVolLookUp		=  StringMaturityToYearTerm( "1d" );
				}

			if( TjVolLookUp		<= 0.0 )
				TjVolLookUp		=	StringMaturityToYearTerm( "1d" );

		}
		else
			TjVolLookUp = - 1.0;

		itsInfTenor -> Elt(i) = tenor;
		itsInfRenor -> Elt(i) = renormFactor;
		if ( TjVolLookUp<0)
			itsInfExper -> Elt(i) =	StringMaturityToYearTerm( "1d" );
		else
			itsInfExper -> Elt(i) = TjVolLookUp;
		itsInfInter -> Elt(i) = itsInfResetTerms -> Elt(i);
	}
}

//======>	GetInfParam

template< class T >
ARM_GP_Vector	ARM_InfCorridorLegWithModel<T>::GetInfParam( const int & i) {
	
	ARM_GP_Vector vec(4);
	vec[0] = itsInfTenor->Elt(i);
	vec[1] = itsInfRenor->Elt(i);
	vec[2] = itsInfExper->Elt(i);
	vec[3] = itsInfInter->Elt(i);

	return vec;
}

//======>	GetInfVol

template< class T >  
double ARM_InfCorridorLegWithModel<T>::GetInfVol( ARM_GP_Vector input, ARM_VolCurvePtr infVol, double strike){

	double tmp;
	double vol;

	double tenor = input[0];
	double renor = input[1];
	double exper = input[2];
	double inter = input[3];

	if ( tenor > 0.0 && inter > 0.0){
		tmp	= pow( renor * strike/VolBase, 1.0/tenor )-1.0;
		vol	= infVol->ComputeVolatility( exper, tmp, tenor )/VolBase;
		vol = vol  *  sqrt( tenor )/ sqrt(inter );  
	}
	else
		vol	= K_NEW_DOUBLE_TOL;
	return vol;
}

//======>	GetIrParam

template< class T >
ARM_GP_Vector	ARM_InfCorridorLegWithModel<T>::GetIrParam( const int & i) {
	ARM_GP_Vector vec(3);
	vec[0] = itsIrTenor->Elt(i);
	vec[1] = itsIrExper->Elt(i);
	vec[2] = itsIrForwd->Elt(i);

	return vec;
}

//======>	GetIrVol

template< class T >  
double ARM_InfCorridorLegWithModel<T>::GetIrVol( ARM_GP_Vector input, ARM_VolCurvePtr irVol, double strike){

	double tenor = input[0];
	double exper = input[1];
	double forwd = input[2];

	return irVol->ComputeVolatility( exper, strike - forwd, tenor )/VolBase;
}

//======>	CptExpFwdInfVol

template< class T > 
void	ARM_InfCorridorLegWithModel<T>::CptExpFwdInfVol(){

	double			strike;

	for( int i=0 ; i<itsNbFlows; i++){

		if ( isIrCritera ){
			if ( itsDownMultiple->Elt(i) != 0.0  )											
				strike	=	RateBase - ( itsFwdIrRates->Elt(i)- itsRangeDown->Elt(i) )/ itsDownMultiple->Elt(i) ;
			else
				strike	=	K_NEW_DOUBLE_TOL;
		}
		else
			strike	=	RateBase + itsDownMultiple->Elt(i)*itsFwdIrRates->Elt(i)+itsRangeDown->Elt(i);
		itsFwdInfVolDown->Elt(i) = GetInfVol(GetInfParam(i), itsInfVol, strike);


		if ( isIrCritera ){
			if ( itsUpMultiple->Elt(i) != 0.0  )											
				strike	=	RateBase - ( itsFwdIrRates->Elt(i)- itsRangeUp->Elt(i) )/ itsUpMultiple->Elt(i) ;
			else
				strike	=	K_NEW_DOUBLE_TOL;
		}
		else
			strike	=	RateBase + itsUpMultiple->Elt(i)*itsFwdIrRates->Elt(i)+itsRangeUp->Elt(i);
		itsFwdInfVolUp->Elt(i) = GetInfVol(GetInfParam(i), itsInfVol, strike);
	}
}

//======>	CptExpFwdIrVol

template< class T >
void	ARM_InfCorridorLegWithModel<T>::CptExpFwdIrVol(){

	double			strike;

	for( int i=0 ; i<itsNbFlows; i++){

		strike	=	itsFwdInfRates->Elt(i) * itsUpMultiple -> Elt(i) + itsRangeUp->Elt(i);
		if ( !isIrCritera ){
			if ( itsUpMultiple -> Elt(i)!= 0.0 )
				strike	=	( itsFwdInfRates->Elt(i) -itsRangeUp->Elt(i) ) / itsUpMultiple -> Elt(i);
			else
				strike	= K_NEW_DOUBLE_TOL;
		}
		itsFwdIrVolUp ->Elt(i)	= GetIrVol(GetIrParam(i),itsIrVol, strike);

		strike	=	itsFwdInfRates->Elt(i) * itsDownMultiple -> Elt(i) + itsRangeDown->Elt(i);
		if ( !isIrCritera ){
			if ( itsDownMultiple -> Elt(i)!= 0.0 )
				strike	=	( itsFwdInfRates->Elt(i) -itsRangeDown->Elt(i) ) / itsDownMultiple -> Elt(i); 
			else
				strike	= K_NEW_DOUBLE_TOL;
		}
		itsFwdIrVolDown ->Elt(i) = GetIrVol(GetIrParam(i), itsIrVol, strike);
	} 
}

//======>	StripCompute

template< class T > 
void	ARM_InfCorridorLegWithModel<T>::StripCompute( const int & k ){
	
	ARM_GP_Vector input(10);

	input[6] = itsFwdCorrel		->	Elt(k);
	input[9] = itsConstant		->	Elt(k);

	double	initDate	= GetModel()->GetStartDate().GetJulian();
	double	refDate		= itsPaymentDates->Elt(k);

	if ( refDate-initDate<0){
		itsUpPrice		->	Elt(k)	= 0.0;
		itsDownPrice	->	Elt(k)	= 0.0;
	}
	else{
		if( isIrCritera){

			itsModel->SetVolParam( GetIrParam(k), GetInfParam(k) );

			input[0] = itsFwdIrRates	->	Elt(k);
			input[1] = itsIrResetTerms	->	Elt(k);
	
			input[2] = itsFwdInfRates	->	Elt(k) + RateBase;
			input[3] = itsInfResetTerms	->	Elt(k);

			input[7] = itsIrLeverage	->	Elt(k);
			input[8] = itsInfLeverage	->	Elt(k);

			input[4] = itsRangeUp		->	Elt(k) - itsUpMultiple	->	Elt(k) * RateBase;
			input[5] = itsUpMultiple	->	Elt(k);

			itsModel	->SetPayoff( input );
			itsUpPrice	->	Elt(k)	=	itsModel->Compute( );

			input[4] = itsRangeDown		->	Elt(k) - itsDownMultiple->	Elt(k) * RateBase;
			input[5] = itsDownMultiple	->	Elt(k);

			itsModel	->SetPayoff( input );
			itsDownPrice->	Elt(k)	=	itsModel->Compute( );
	
		}
		else{

			itsModel->SetVolParam( GetInfParam(k), GetIrParam(k) );

			input[0] = itsFwdInfRates	->	Elt(k) + RateBase;
			input[1] = itsInfResetTerms	->	Elt(k);
	
			input[2] = itsFwdIrRates	->	Elt(k);
			input[3] = itsIrResetTerms	->	Elt(k);

			input[7] = itsInfLeverage	->	Elt(k);
			input[8] = itsIrLeverage	->	Elt(k);

			input[4] = itsRangeUp		->	Elt(k) + RateBase;
			input[5] = itsUpMultiple	->	Elt(k) ;

			itsModel	->SetPayoff( input );
			itsUpPrice	->	Elt(k)	=	itsModel->Compute( );

			input[4] = itsRangeDown		->	Elt(k) + RateBase;
			input[5] = itsDownMultiple	->	Elt(k) ;

			itsModel	->SetPayoff( input );
			itsDownPrice->	Elt(k)	=	itsModel->Compute( );
		}
	}
}


//======>	ViewPricerFeatures


template< class T >
string	ARM_InfCorridorLegWithModel<T>::ViewPricerFeatures ( ){

	CC_Ostringstream	os;
	vector<string>		vc;
	vector<double>		vd(7);

	os	<< "\n" <<"Pricer Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"UpInfVol";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DownInfVol";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"UpIrVol";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DownIrVol";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"CorrelIr/Inf";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"UpPrice";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DownPrice";	
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed)<< CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(7*LAG+2)<<" \n";
	
	for (int i = 0; i < itsNbFlows; ++i )
	{
		vd[0]=	itsFwdInfVolUp	->Elt(i) * sqrt( fabs( itsInfResetTerms->Elt(i) ) ) * VolBase;
		vd[1]=	itsFwdInfVolDown->Elt(i) * sqrt( fabs( itsInfResetTerms->Elt(i) ) ) * VolBase;
		vd[2]=	itsFwdIrVolUp	->Elt(i) * VolBase;
		vd[3]=	itsFwdIrVolDown ->Elt(i) * VolBase;
		vd[4]=	itsFwdCorrel	->Elt(i);
		vd[5]=	itsUpPrice		->Elt(i);
		vd[6]=	itsDownPrice	->Elt(i);
		
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[3];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[4];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[5];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[6];
		os <<"\n";
	}
	return os.str( );
}

CC_END_NAMESPACE()

#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/























