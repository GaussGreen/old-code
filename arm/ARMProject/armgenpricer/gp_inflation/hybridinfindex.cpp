
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


 /*

	Ceci n'est pas un index 
	Magritte

*/

#pragma warning(disable : 4786)

#include "gpinfra/argconvdefault.h"

#include "gpbase/typedef.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/cloneutilityfunc.h"

#include "gpinflation/hybridinfindex.h"

//#include "gpinflation/infleg.h"			
#include "gpinflation/infComposedLeg.h"
#include "gpinflation/infPayOff.h"


#include "gpinflation/infdata.h"
#include "gpinflation/infbsmodel.h"
#include "gpinflation/infbssmiledmodel.h"

#include "gpclosedforms/smile_sabr.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal.h"

#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/brent.h"
#include "gpnumlib/numfunction.h"

#include "mod/bssmiled.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: ConvertToARM_GP_Vector
///	Returns: ARM_GP_VectorPtr
///	Action : Convert a ARM_Vector* into a ARM_GP_VectorPtr

////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_InfIrIndex::ConvertToARM_GP_Vector ( const ARM_Vector * vec){
	
	ARM_GP_VectorPtr tmp( new ARM_GP_Vector);

	double dim = vec->size();
	if( dim >0){
		tmp->resize(dim);
		for ( int i = 0; i< dim; i++)
			tmp->Elt(i) = vec->Elt(i);
	}
	else
		tmp = ARM_GP_VectorPtr(NULL);
	return tmp;
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: ARM_InfIrIndex
///	Returns: Void
///	Action : Construction

////////////////////////////////////////////////////
ARM_InfIrIndex::ARM_InfIrIndex( const string			&	indexName,
								map<string,ARM_Date>	&	mDate, 
								map<string,string>		&	mString,
								map<string,int>			&	mInt){

	itsAsOfDate = mDate["asOfDate"];

	if( indexName =="N"){
		itsIndex = ARM_CountedPtr<ARM_SwapLeg>	(NULL);
		itsType  = NO;
	}
	else if ( indexName =="EMU" || indexName =="IFRF"){

		int				tmpRec			= 0;
		char *			tmpPayCal		= const_cast< char *> (mString["payCal"].c_str() );
		char *			tmpInfResetCal	= "INF";
		ARM_CountedPtr<ARM_Currency>	ccy(new ARM_Currency	( mString["currency"].c_str() ));
		double			tmpSpread		= 0.0;
		double			tmpMultiple		= 1.0;
		double			tmpComultiple	= 1.0;
		int				tmpDayCount		= ARM_ArgConv_LgNameDayCount.GetNumber("30/360"); //PFFFFFF ...
		int				tmpAdjFirstRule	= ARM_ArgConv_IntRule.GetNumber("UNADJ");
		int				resetFreq		= mInt["resetFreq"] ;
		int				payFreq			= mInt["payFreq"] ;


		ARM_GP_Vector nominals(1,1.);
		bool monthBegin(false) ;
		Period numGap(mInt["resetNumGap"],K_DAILY) ;
		Period demGap(mInt["resetDemGap"],K_DAILY) ;
		ARM_CountedPtr <YOYPayOff> payoff( new  YOYPayOff());
		int resetRollConv = mInt["fwdRule"], payRollConv = mInt["fwdRule"]  ;//PFFFF
		Period tenor(1, K_ANNUAL) ;
		
		ARM_CountedPtr<ARM_InfIdx> cpi(new ARM_InfIdx(indexName));
		
		ARM_CountedPtr<InfPerformance> cpiPerf( new InfPerformance(cpi,tenor));
		


		ARM_GP_T_Vector<ARM_CountedPtr < CashFlow > >
		vCahFlows = InfFloatingCouponVector(nominals,payoff, cpiPerf, mDate["startDate"],
											mDate["endDate"], resetFreq,
											payFreq, tmpDayCount,
											mInt["fwdRule"], mInt["stubRule"],
											mInt["intRule"], resetRollConv,
											payRollConv, tmpPayCal, monthBegin,
											demGap,numGap);

		//TODO : faudrait trouver mieux (faire un wrapper)
		if (payFreq == resetFreq)
			itsIndex = ARM_CountedPtr<InfLeg>(new InfLeg(vCahFlows, mInt["infInterType"], mInt["firstReset"]));	
		else
			itsIndex = ARM_CountedPtr<InfComposedLeg>(new InfComposedLeg(vCahFlows, mInt["infInterType"], mInt["firstReset"]));	
		
		itsType = IN;
	}
	else
	{

		ARM_Currency*	tmpCcy			= new ARM_Currency	( mString["currency"].c_str() );
		ARM_INDEX_TYPE  tmpIndex		= (ARM_INDEX_TYPE)	ARM_ArgConv_IndexType.GetNumber(indexName.c_str());		
		int				tmpRec			= 0;
		double			tmpSpread		= 0.0;
		char *			tmpResetCal		= const_cast< char *> (mString["resetCal"].c_str() );
		char *			tmpPayCal		= const_cast< char *> (mString["payCal"	 ].c_str() );
		char *			tmpRefDate		= const_cast< char *> (mString["refDate" ].c_str() );;

		ARM_IRIndex* tmpIrIndex	=	new ARM_IRIndex(	
									mInt	[	"dayCount"			], 
									mInt	[	"resetFreq"			], 
									mInt	[	"payFreq"			], 
									mInt	[	"maturity"			], 
									mInt	[	"compMeth"			], 
									mInt	[	"fwdRule"			], 
									mInt	[	"intRule"			], 
									mInt	[	"resetTiming"		],  
									mInt	[	"resetGap"			],
									mInt	[	"payTiming"			],
									mInt	[	"payGap"			], 
									tmpCcy,
									tmpIndex,
									mInt	[	"decompFreq"		]);	
		
		itsIndex = ARM_CountedPtr<ARM_SwapLeg> (	new ARM_SwapLeg(	
									mDate	[	"startDate"			],			
									mDate	[	"endDate"			],	
									tmpIrIndex,
									tmpRec,
									tmpSpread ,
									mInt	[	"stubRule"			],
									mInt	[	"decompFreq"		],  
									tmpCcy,
									mInt	[	"dayCount"		],
									10000,					
									tmpResetCal,
									tmpPayCal,
									mInt	[	"decompPriceFlag"	],
									mInt	[	"finNotioType"		],
									tmpRefDate,
									mInt	[	"adjFirstRule"		]) );
		int toto = mInt	[	"adjFirstRule"		];



		itsType = IR;
		if ( tmpIrIndex){ delete tmpIrIndex;	tmpIrIndex=NULL; }
		if ( tmpCcy    ){ delete tmpCcy;		tmpCcy=NULL;	 }
	}
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: ARM_InfIrIndex
///	Returns: Void
///	Action : Copy Construction

////////////////////////////////////////////////////
ARM_InfIrIndex::ARM_InfIrIndex( const ARM_InfIrIndex & rhs ):itsType( rhs.itsType){ 
	itsIndex			= rhs.itsIndex;
	itsMod				= rhs.itsMod;
	itsAsOfDate			= rhs.itsAsOfDate;
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetNumDates
///	Returns: ARM_GP_VectorPtr
///	Action : return Num ( Pub ) Date if index is inf

////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfIrIndex::GetNumDates( ) const { 

	ARM_CountedPtr<InfLegBase> tmpLeg = itsIndex ;
	if (tmpLeg.IsNull()) ARMTHROW(ERR_INVALID_ARGUMENT," ARM_InfIrIndex::GetDemDates( ) You cannot use this method on IR leg");
	ARM_GP_VectorPtr tmpDates ( new ARM_GP_Vector(*tmpLeg->GetNumJulianDates() )) ;
	return tmpDates;
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetDemDates
///	Returns: ARM_GP_VectorPtr
///	Action : return Dem Date if index is inf

////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_InfIrIndex::GetDemDates( ) const { 

	ARM_CountedPtr<InfLegBase> tmpLeg = itsIndex ;
	if (tmpLeg.IsNull()) ARMTHROW(ERR_INVALID_ARGUMENT," ARM_InfIrIndex::GetDemDates( ) You cannot use this method on IR leg");

	ARM_GP_VectorPtr tmpDates (new ARM_GP_Vector(*tmpLeg->GetDemJulianDates() )) ;
	return tmpDates;

}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetIndex
///	Returns: ARM_IRIndex
///	Action : return the Type ofInfex

////////////////////////////////////////////////////

ARM_CountedPtr<ARM_IRIndex>		ARM_InfIrIndex::GetIndex( ) const{

	if( dynamic_cast<InfLegBase*> (&*itsIndex) )
		return dynamic_cast<InfLegBase*> (&*itsIndex)->GetIndex();
	else
		return ARM_CountedPtr<ARM_IRIndex> ( CreateClone(itsIndex->GetIRIndex()) );
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetBegDates
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the Reset Dates 

////////////////////////////////////////////////////

ARM_GP_VectorPtr  ARM_InfIrIndex::GetBegDates( ) const {

	ARM_GP_VectorPtr tmp;
	double initDate = itsAsOfDate.GetJulian();
	if  ( itsType == IR)	tmp = ConvertToARM_GP_Vector (itsIndex->GetFlowStartDates() );
	else					
		tmp = GetDemDates( );

	return tmp;
}


////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetEndDates
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the Reset Dates 

////////////////////////////////////////////////////

ARM_GP_VectorPtr  ARM_InfIrIndex::GetEndDates( ) const {

	ARM_GP_VectorPtr tmp;
	double initDate = itsAsOfDate.GetJulian();
	if  ( itsType == IR)	tmp = ConvertToARM_GP_Vector (itsIndex->GetFlowEndDates() );
	else					
		tmp = GetNumDates( );

	return tmp;
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetResDates
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the Reset Dates 

////////////////////////////////////////////////////

ARM_GP_VectorPtr  ARM_InfIrIndex::GetResLags( ) const {

	ARM_GP_Vector resetTimes ;
	ARM_CountedPtr<InfLegBase> infLeg = itsIndex ;
	ARM_CountedPtr<ARM_SwapLeg> irLeg = itsIndex;
	int dc ;
	double time;
	if ( !infLeg.IsNull())
	{
		ARM_GP_VectorPtr demJulianDates = GetDemDates();
		dc = infLeg->GetIndex()->GetDayCount();
		string	lag	= InfData::GetPublishLag( infLeg->GetIndexName().c_str() );
		ARM_Date demDate ;
		for ( int i =0; i< demJulianDates->size(); ++i)
		{	
			demDate = ARM_Date(demJulianDates->Elt(i) );
			time = CountYearsWithoutException(dc, itsAsOfDate, demDate.AddPeriod(lag)) ;
			resetTimes.push_back(time);
		}
//ARM_GP_VectorPtr toto = ARM_GP_VectorPtr(new ARM_GP_Vector( *GetDemDates()));


	}
	else if  ( !irLeg.IsNull())
	{
		ARM_Vector* resetDates = itsIndex->GetResetDates();
		dc = irLeg->GetIRIndex()->GetDayCount();
		for ( int i = 0; i < resetDates->size() ; ++i)
		{	
			time = CountYearsWithoutException(dc, itsAsOfDate, resetDates->Elt(i) ) ;
			resetTimes.push_back(time);
		}
	}
	else ARMTHROW(ERR_INVALID_ARGUMENT," ARM_InfIrIndex::itsIndex must be an InfComposedLeg or ARM_SwapLeg");
	
	ARM_GP_VectorPtr resetTimesPtr ( new ARM_GP_Vector(resetTimes)) ;

	//	if ( !infLeg.IsNull())
	//ARM_GP_VectorPtr toto = GetDemDates();

	return resetTimesPtr;
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetTenDates
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the Tenor Dates 

////////////////////////////////////////////////////

ARM_GP_VectorPtr  ARM_InfIrIndex::GetTenLags( ) const{

	ARM_GP_Vector tenorTimes ;
	ARM_CountedPtr<InfLegBase> infLeg = itsIndex ;
	ARM_CountedPtr<ARM_SwapLeg> irLeg = itsIndex;

	if ( !infLeg.IsNull())
	{
		ARM_GP_VectorPtr numJulianDates =	GetNumDates( );
		ARM_GP_VectorPtr demJulianDates =	GetDemDates( );
	
		const double	dayCount= infLeg->GetIndex()->GetDayCount();
		const string	pubLag	= InfData::GetPublishLag( infLeg->GetIndex()->GetIndexName().c_str() );

		double time;
		ARM_Date demDate, numDate ;
		for ( int i =0; i< numJulianDates->size(); ++i)
		{	
			demDate = ARM_Date(demJulianDates->Elt(i) );
			numDate = ARM_Date(numJulianDates->Elt(i) );
			demDate.AddPeriod(pubLag);
			numDate.AddPeriod(pubLag);
			if( demDate<itsAsOfDate) demDate=itsAsOfDate;
			time = CountYearsWithoutException(dayCount, demDate, numDate) ;
			tenorTimes.push_back(time);
		}
	}
	else if  ( !irLeg.IsNull())
	{
		int dim = itsIndex->GetFlowStartDates()->size();
		tenorTimes = ARM_GP_Vector(dim, itsIndex->GetIRIndex() -> GetYearTerm()) ;
	}
	else ARMTHROW(ERR_INVALID_ARGUMENT," ARM_InfIrIndex::itsIndex must be an InfComposedLeg or ARM_SwapLeg");

	ARM_GP_VectorPtr tenorTimesPtr ( new ARM_GP_Vector(tenorTimes)) ;
	return tenorTimesPtr;

}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetIndexName
///	Returns: string
///	Action : return index name

////////////////////////////////////////////////////
string ARM_InfIrIndex::GetIndexName( ) const { 

	if ( dynamic_cast<InfLegBase*> ( &*itsIndex ) )
		return dynamic_cast<InfLegBase*> ( &*itsIndex )->GetIndexName();
	else 
		return ARM_ParamView::GetMappingName(	S_INDEX_TYPES, itsIndex->GetIRIndex()->GetIndexType() );
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: GetName
///	Returns: string
///	Action : return name

////////////////////////////////////////////////////
string ARM_InfIrIndex::GetName( ) const { 

	if ( dynamic_cast<InfLegBase*> ( &*itsIndex ) )
		return dynamic_cast<InfLegBase*> ( &*itsIndex )->GetIndexName();
	else 
		return itsIndex->GetCurrencyUnit() -> GetCcyName();
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptQua
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the quantile 

////////////////////////////////////////////////////

void	ARM_InfIrIndex::SmoothQua( const ARM_GP_Vector& abs, ARM_GP_Vector& ord){

/*	const double	tol		= 1.2;
	double			slope	= (ord[1]-ord[0])/(abs[1]-abs[0]);
	int				nb		= abs.size();
	double			b_Abs, b_Ord;


	ARM_GP_Vector	tmpAbs;
	ARM_GP_Vector	tmpOrd;
	
	tmpAbs.push_back(abs[0]);
	tmpOrd.push_back(ord[0]);


	b_Abs	=	 abs[0];
	b_Ord	=	 ord[0];

	for ( int i = 1; i<nb; ++i){
		tmp = b_Ord+slope*(abs[i]-b_Abs);

		if ( vec[i]/tmp < tol && vec[i]/tmp>1/tol ){
			tmpAbs.push_back(abs[i]);
			tmpOrd.push_back(ord[i]);			
		}
	
	slope= (ord[i]-ord[i-1])/(abs[i]-abs[i-1]);
	}*/
}


////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptQua
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the quantile 

////////////////////////////////////////////////////

double	ARM_InfIrIndex::CptQua(	const double & exp,
								const double & ten,
								const double & fwd, 
								const double & gue,
								const double & str) const{

	if (itsType ==IR)
		return CptIrQuantile( exp, ten, fwd, str );
	else
		return CptInQuantile( exp, ten, fwd, gue, str );
}



////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptQuantile
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the quantile 

////////////////////////////////////////////////////

double	ARM_InfIrIndex::CptIrQuantile(	const double & exp,
										const double & ten,
										const double & fwd, 
										const double & str) const{

	double tmpExp = exp;
	double tmpTen = ten;
	double tmpFwd = fwd;

	double sigma	= itsMod->GetVolatility()	->ComputeVolatility(tmpExp,  tmpTen);
	double beta		= itsMod->GetBeta()			->ComputeVolatility(tmpExp,  tmpTen);
	double rho		= itsMod->GetRho()			->ComputeVolatility(tmpExp,  tmpTen);
	double nu		= itsMod->GetNu()			->ComputeVolatility(tmpExp,  tmpTen);
	

	if( itsType==IR) tmpFwd = fwd/RateBase;

    double alpha	= itsMod->ComputeAlpha(tmpFwd, tmpFwd, tmpExp, sigma/RateBase, 
                                rho, nu, beta, itsMod->GetWeight(), K_SABR_IMPLNVOL);

/*
double SABR_smile::inverse_distribution(double f,double proba,double tex, double alpha, double beta, double rho, 
										double nu,int flag,int nbsteps,
										double alpha_exp,double alpha_tanh,double kb_tanh)
*/
	double result= SABR_smile::inverse_distribution( tmpFwd, str, tmpExp, alpha, beta, rho, nu, K_SABR_IMPLNVOL,200,0.01,1.5,0.025);
	result *= RateBase;

	return result;
}




////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptInQuantile
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the quantile 

////////////////////////////////////////////////////


class CptDist: public ARM_GP::UnaryFunc<double,double>{

public:
	CptDist	(	ARM_InfIrIndex		  idx,
				const double		& res,
				const double		& ten,
				const double		& fwd, 
				const double		& pro): itsRes(res), itsTen(ten), itsFwd(fwd), itsPro(pro){	
				itsEpsilon	= 1e-6;	
				itsIdx		= idx;
				}

	virtual double operator() (double x) const;

private:
	virtual double CptCall( const double & fwd, const double & str, const double & var) const;

private:
	double			itsRes;
	double			itsTen; 
	double			itsFwd;
	double			itsPro; 

	ARM_InfIrIndex		itsIdx;
	double				itsEpsilon;
};

double CptDist::CptCall( const double & fwd, const double & str, const double & var)const{
	double d1,d2;

	d1 = fwd;
	d1 /=str;
	d1 = log(d1);
	d1 /=var;
	
	d2 = d1 - 0.5*var;
	d1 +=0.5*var;

	return fwd * NormalCDF(d1) - str * NormalCDF(d2);
}

double CptDist::operator() (double x) const{ 
	
	double var	= itsIdx.CptVol( itsRes, itsTen, itsFwd, x);
	var			*=sqrt(fabs(itsRes));
	double tmp	= CptCall( itsFwd, x, var);

	var	= itsIdx.CptVol( itsRes, itsTen, itsFwd, x+itsEpsilon);
	var *=sqrt(fabs(itsRes));
	tmp	-= CptCall( itsFwd, x+itsEpsilon, var);
	tmp /= itsEpsilon;

	tmp = 1-tmp;
	tmp -= itsPro;
	return tmp;		
}


double	ARM_InfIrIndex::CptInQuantile(	const double & exp,
										const double & ten,
										const double & fwd, 
										const double & gue,
										const double & str) const{

	double result, result1;

	CptDist func( *this, exp, ten, fwd, str );
	UnaryFuncWithNumDerivative<double> funcWithDeriv(func);
	T_NewtonRaphsonSolverNoThrow<UnaryFuncWithNumDerivative<double> > solver(funcWithDeriv, 1e-6, 1e-6 ,1e-6);
	solver.setInitialGuess(gue);
	result= solver.Solve();
/*

	double tmpExp = exp+ten ;
	double tmpTen = ten;
	double tmpFwd = fwd;

	double sigma	= itsMod->GetVolatility()	->ComputeVolatility(tmpExp,  tmpTen);
	double beta		= itsMod->GetBeta()			->ComputeVolatility(tmpExp,  tmpTen);
	double rho		= itsMod->GetRho()			->ComputeVolatility(tmpExp,  tmpTen);
	double nu		= itsMod->GetNu()			->ComputeVolatility(tmpExp,  tmpTen);
	

    double alpha	= itsMod->ComputeAlpha(tmpFwd, tmpFwd, tmpExp, sigma/RateBase, 
                                rho, nu, beta, itsMod->GetWeight(), K_SABR_IMPLNVOL);

	int		nbsteps		= 200;
	double	alpha_exp	= 0.01;
	double	alpha_tanh	= 1.5;
	double	kb_tanh		= 0.025;

	result1= SABR_smile::inverse_distribution( tmpFwd, str, tmpExp, alpha, beta, rho, nu, K_SABR_IMPLNVOL,nbsteps,alpha_exp,alpha_tanh,kb_tanh);
*/
	return result ;
}

////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptAdj
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the convexity adjustment ( private ) 

////////////////////////////////////////////////////


double ARM_InfIrIndex::CptAdj(	const double & res,
								const double & mat,
								const double & pay,
								const double & fwd) const{

	ARM_CountedPtr<InfLegBase> infLeg = itsIndex ;
	ARM_CountedPtr<ARM_SwapLeg> irLeg = itsIndex;
	if ( !infLeg.IsNull())
	{
		ARM_GP_Vector rat;
		ARM_CountedPtr<ARM_InfIdx>	idx(infLeg->GetIndex()->GetIndex()) ;
		ARM_CountedPtr<ARM_InfBSSmiledModel> mod = itsMod;
		if (mod.IsNull()) ARMTHROW(ERR_INVALID_ARGUMENT,"ARM_InfIrIndex::CptAdj : ARM_InfBSSmiledModel has to be used");
		rat = mod->CptFwdCPIRatio(	mat, //numdate
									res, //demdate
									pay,
									infLeg->GetInterpType(), 
									infLeg->GetFirstReset(), 
									idx.operator->());

		return  (1+rat[0])*rat[2]/rat[1];
	}
	else if  ( !irLeg.IsNull())
	{
		return 	 fwd/itsMod->ExpectedFwdYieldNonAdj( res , mat , pay );
	}
	else
		ARMTHROW(ERR_INVALID_ARGUMENT," This Index is not available");
	
}



////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptFwd
///	Returns: ARM_GP_VectorPtr
///	Action : Compute the Fwd ( public ) 

////////////////////////////////////////////////////

double ARM_InfIrIndex::CptFwd(	const double & res,
								const double & mat,
								const double & pay) const{ 

	ARM_CountedPtr<InfLegBase> infLeg = itsIndex ;
	ARM_CountedPtr<ARM_SwapLeg> irLeg = itsIndex;
	if ( !infLeg.IsNull())
	{
		ARM_GP_Vector rat;
		ARM_CountedPtr<ARM_InfIdx>	idx(infLeg->GetIndex()->GetIndex()) ;
		ARM_CountedPtr<ARM_InfBSSmiledModel> mod = itsMod;
		if (mod.IsNull()) ARMTHROW(ERR_INVALID_ARGUMENT,"ARM_InfIrIndex::CptFwd : ARM_InfBSSmiledModel has to be used");

		rat = mod->CptFwdCPIRatio(	mat,	//numer 
									res,	//denom
									pay,	//pay
									infLeg->GetInterpType(), 
									infLeg->GetFirstReset(), 
									idx.operator->());

		return  rat[0]*RateBase;
	}
	else if  ( !irLeg.IsNull())
	{
		return  itsMod->ExpectedFwdYield( res , mat , pay );
	}
	else ARMTHROW(ERR_INVALID_ARGUMENT," ARM_InfIrIndex::itsIndex must be an InfComposedLeg or ARM_SwapLeg");
}


////////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: CptVol
///	Returns: double
///	Action : compute the ATM Vol

////////////////////////////////////////////////////


double	ARM_InfIrIndex::CptVol(	const double & exp, 
								const double & ten,  
								const double & fwd, 
								const double & str) const{

	double vol= K_NEW_DOUBLE_TOL;

	ARM_BSSmiledModel * mod = dynamic_cast< ARM_BSSmiledModel*> ( &*itsMod ); 
	double tmp = exp<0?1.0/K_YEAR_LEN:exp;

	if ( itsType == IR && str<= 0.0)
		vol = mod ->  ComputeVol( tmp, ten, fwd, fwd);
	else if ( itsType == IN && str<= -RateBase)
		vol = mod ->  ComputeVol( tmp, ten, fwd, fwd);
	else
		vol = mod ->  ComputeVol( tmp, ten, fwd, str);
	
	if ( itsType == IN ) {
		if (exp< 0)
			vol *= sqrt( fabs(ten) )/sqrt(fabs(exp) );
		else
			vol *= sqrt( fabs(tmp+ten) )/sqrt(fabs(exp) );
	}
	return vol/VolBase;
}

///////////////////////////////////////////////////

///	Class  : ARM_InfIrIndex
///	Routine: SetMod
///	Returns: void
///	Action : set the model

////////////////////////////////////////////////////
void	ARM_InfIrIndex::SetMod( ARM_BSSmiledModel* mod){
	itsMod = ARM_CountedPtr<ARM_BSSmiledModel> ( mod );
}

///////////////////////////////////////////////////

///	Class  : ARM_InfIrCorrel
///	Routine: ARM_InfIrCorrel
///	Returns: void
///	Action : copy consructor


////////////////////////////////////////////////////
ARM_InfIrCorrel::ARM_InfIrCorrel( const ARM_InfIrCorrel & rhs ){

	map< pair<string,string>,  ARM_VolCurvePtr>::const_iterator it;
	for ( it =rhs.itsCorrel.begin() ; it!=rhs.itsCorrel.end();it++)
		itsCorrel[it->first] = CreateClonedPtr(&*it->second );
}

///////////////////////////////////////////////////

///	Class  : ARM_InfIrCorrel
///	Routine: SetCorrel
///	Returns: void
///	Action : fill the double map of correlation

////////////////////////////////////////////////////
void	ARM_InfIrCorrel::SetCorrel	( const string & idx1, const string & idx2, ARM_VolCurve* vol){

	pair<string,string> q(idx1,idx2);
	pair<pair<string,string>, ARM_VolCurvePtr>  p(q, ARM_VolCurvePtr(vol));
	itsCorrel.insert(p);
}

///////////////////////////////////////////////////

///	Class  : ARM_InfIrCorrel
///	Routine: GetCorrel
///	Returns: double
///	Action : return the value of the correlation for index of same expery and the same tenor

////////////////////////////////////////////////////

double ARM_InfIrCorrel::GetCorrel	( const string & idx1, const string & idx2, const double & exp, const double & ten) const {
	
	ARM_VolCurvePtr vol = GetCorrel( idx1, idx2);

	double	tmp = exp<1.0/K_YEAR_LEN?1.0/K_YEAR_LEN:exp;

	if ( &*vol)
		return vol->ComputeVolatility( tmp, ten)/ VolBase;

	vol = GetCorrel( idx2, idx1);
	if ( &*vol)
		return vol->ComputeVolatility( tmp, ten)/ VolBase;
	else
		ARMTHROW(ERR_INVALID_ARGUMENT," bad argument in ARM_InfIrCorrel::GetCorrel");
	

	return -99.9;
}

///////////////////////////////////////////////////

///	Class  : ARM_InfIrCorrel
///	Routine: GetCorrel
///	Returns: ARM_VolCurvePtr
///	Action : return the corresponding vol curve

////////////////////////////////////////////////////

ARM_VolCurvePtr	ARM_InfIrCorrel::GetCorrel	( const string &idx1, const string &idx2) const { 

	ARM_VolCurvePtr vol=ARM_VolCurvePtr(NULL);

	pair<string,string> p ( idx1, idx2);
	map< pair<string,string>,  ARM_VolCurvePtr>::const_iterator it = itsCorrel.find(p);
	if ( it != itsCorrel.end() )
		vol= CreateClonedPtr(&*it->second);

	return vol;
}


CC_END_NAMESPACE()	
/*---------------------------------------------------------------*/
/*---- End Of File ----*/

