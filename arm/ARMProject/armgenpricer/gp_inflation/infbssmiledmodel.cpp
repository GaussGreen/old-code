/*----------------------------------------------------------------------------*/
 
/*! \file infbssmiledmodel.cpp
 * Copyright (c) NATIXIS CM April 2007 Paris
 *
 *  \brief Sabr model for the inflation volatility 
 *  + vol Atm for the interest rate
 *  + corelation between 2 atm volatilities ( Ir and Inf) 
 *
 *	\author  Mathieu Bernardo
 *	\version 1.0
 *	\date April 2007
 */

/*----------------------------------------------------------------------------*/
#include <gpbase\assignop.h>

#include "gpinflation/infbssmiledmodel.h"

/// gpbase
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/stringconvert.h"

#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpmatrix.h"

/// kernel
#include <mod/bsconvadjust.h>

/// gpinflation
#include "gpinflation/infdata.h"
#include "gpinflation/infcapfloor.h"
#include "gpinflation/infswaption.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
//#include "gpclosedforms\sabrbdiff1.h"

/// STL
#include <functional>
#include <crv/correlmanager.h>
#include <crv/volcube.h>
//////






#include "gpinflation/infbsmodel.h"

/// gpbase
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpmatrixlinalg.h"
#include <gpbase/cloneutilityfunc.h>

/// gpinflation
#include "gpinflation/infcurv.h"
#include "gpinflation/infcapfloor.h"
#include "gpinflation/infswaption.h"
#include "gpinflation/infswopvol.h"
#include "gpinflation/infdata.h"
#include "gpinflation/sparsevolcube.h"
#include "gpinflation/infidx.h"
#include "gpinflation/assetinfo.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/spreadoption_lognormal_interface.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_shiftedlognormal_interface.h"
#include "gpclosedforms/vanilla_shifted_lognormal.h"
#include "gpclosedforms/vanille_normal_interface.h"
#include "gpclosedforms/vanilla_normal.h"


/// kernel
#include <mod/bsconvadjust.h>
#include <crv/correlmanager.h>
#include <crv/volcube.h>


//////////

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ARM_InfBSSmiledModel
///	Returns: void
///	Action : constructor
////////////////////////////////////////////////////

ARM_InfBSSmiledModel::ARM_InfBSSmiledModel(	ARM_Date&			asOfDate,
											ARM_ZeroCurve*		discountCurve, 
											ARM_InfCurv*		infFwdCurve,
											ARM_VolCurve*		irCapVolCurve,
											ARM_VolCurve*		correlVolCurve,
											ARM_VolCurve*		correlAdjVolCurve,
											ARM_VolCurve*		infCapVolCurve,
											ARM_VolCurve*		nuCurve,
											ARM_VolCurve*		rhoCurve,
											ARM_VolCurve*		betaCurve,
											int sabrMode,
											int alphaOrSigmaInput):
					InfOptionModel		(infFwdCurve) ,
					ARM_BSSmiledModel	(	asOfDate, 0.0,	
											discountCurve,
											discountCurve,
											infCapVolCurve, 
											K_PRICE, 
											rhoCurve, 
											nuCurve,
											sabrMode,
											betaCurve),
					itsCorrelVolCurve	( correlVolCurve? (ARM_VolCurve*) correlVolCurve->Clone() : NULL  ), 
					itsIrCapVolCurve	( irCapVolCurve?  (ARM_VolCurve*) irCapVolCurve->Clone()  : NULL  ){
						if ( !correlAdjVolCurve )
							itsCorrelAdjVolCurve = GenerateAdjCorrelMatrix( );	
						else
							itsCorrelAdjVolCurve = (ARM_VolCurve*) correlAdjVolCurve->Clone();
					}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ARM_InfBSSmiledModel
///	Returns: void
///	Action : Copy constructor
////////////////////////////////////////////////////

ARM_InfBSSmiledModel::ARM_InfBSSmiledModel(const ARM_InfBSSmiledModel& rhs): 
									InfOptionModel(rhs.GetInfFwdCurv() ),
									ARM_BSSmiledModel(rhs){
	itsCorrelVolCurve	= CreateClone( rhs.itsCorrelVolCurve );
	itsIrCapVolCurve	= CreateClone( rhs.itsIrCapVolCurve  );	
	itsCorrelAdjVolCurve= CreateClone( rhs.itsCorrelAdjVolCurve  );	
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ~ARM_InfBSSmiledModel
///	Returns: void
///	Action : destructor
////////////////////////////////////////////////////

ARM_InfBSSmiledModel::~ARM_InfBSSmiledModel(){
		if ( itsCorrelVolCurve )	{ delete itsCorrelVolCurve;		itsCorrelVolCurve = NULL; }
		if ( itsIrCapVolCurve )		{ delete itsIrCapVolCurve;		itsIrCapVolCurve  = NULL; }
		if ( itsCorrelAdjVolCurve ) { delete itsCorrelAdjVolCurve;  itsCorrelAdjVolCurve=NULL;}
	}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ARM_InfBSSmiledModel
///	Returns: FwdCPIRatio
///	Action : return yoy adjusted
////////////////////////////////////////////////////

ARM_GP_Vector ARM_InfBSSmiledModel::FwdCPIRatio(	const ARM_Date& numDate,
													const ARM_Date& demDate,
													const ARM_Date& paymentDate,
													double multiple,
													double spread,
													long dailyInterpType,
													double denomFixing,
													ARM_InfIdx* infIdx)
{
	return CptFwdCPIRatio(numDate.GetJulian(), demDate.GetJulian(), paymentDate.GetJulian(), dailyInterpType, denomFixing, infIdx) ;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: GenerateAdjCorrelMatrix
///	Returns: ARM_VolCurve*
///	Action : return the correlation matrix used in the computation 
///			 of the convexity adjustement ( return - matrix) due to 
///			 the implemantation of FwdCPIRatioAdjust
////////////////////////////////////////////////////

ARM_VolCurve*	ARM_InfBSSmiledModel::GenerateAdjCorrelMatrix( ){

	const int		nb			= 40;
	const double	D_1			= 1.0/K_YEAR_LEN;

	const double	RateBase	= CC_NS( ARM_Constants, rateBase );
	const double	VolBase		= CC_NS( ARM_Constants, volBase);

	ARM_Date		asOfDate	= GetStartDate();
	double			asOfTime	= asOfDate.GetJulian();
	double tmp, tmpYtY, tmpDem, tmpNum;
	double volYtY, volDem, volNum;

	ARM_Vector* X_label = new ARM_Vector(nb);		// input the label of the correl matrix: [abs] denom pub date & [ord] reset date or Eurib 1Y 
	ARM_Vector* Y_label = new ARM_Vector(nb);
	ARM_Matrix* correl  = new ARM_Matrix (nb,nb);
	ARM_GP_Vector volBond(nb);
	
	ARM_BSModel* mod= new ARM_BSModel(GetYCModel()->GetZeroCurve(), itsIrCapVolCurve );

	for( int i = 0; i < nb; i++){
		correl->Elt(0,i) = itsCorrelVolCurve->ComputeVolatility(D_1,i+1);// correspond to the YtY with denom equal 0 

		tmp			= mod->ExpectedFwdYield( asOfTime + (i+1)*K_YEAR_LEN, asOfTime + (i+2)*K_YEAR_LEN, asOfTime + (i+2)*K_YEAR_LEN )/RateBase;
		tmp			= tmp/(1+tmp);
		volBond[i]	= tmp * itsIrCapVolCurve->ComputeVolatility( i+1, 1) / VolBase;
		X_label->Elt(i)= i>0?i:D_1;
		Y_label->Elt(i)= (i+1);
	}

	if ( mod ) { delete mod; mod = NULL; }

	for ( i = 1; i< nb; i++){
		volYtY	= GetVolatility()->ComputeVolatility( i	 , 1	) / VolBase ;	
		volDem	= GetVolatility()->ComputeVolatility( D_1, i	) / VolBase ;
		volNum	= GetVolatility()->ComputeVolatility( D_1, i+1	) / VolBase ;

		for( int j=0; j< nb; j++){
			tmpYtY = 0;
			tmpDem = 0;
			tmpNum = 0;
			for( int k=0; k<j; k++){
				tmpYtY += volBond[k] * itsCorrelVolCurve->ComputeVolatility(i,k+1);
				tmpDem += volBond[k] * correl->Elt( i-1, k);
				tmpNum += volBond[k] * correl->Elt( i,   k);
			}

			tmpYtY += volBond[j] * itsCorrelVolCurve->ComputeVolatility(i,j+1);
			tmpDem += volBond[j] * correl->Elt( i-1, j);

			tmpYtY *= (i+1.0)*volYtY;
			tmpNum *= (i+1.0)*volNum;
			tmpDem *= i*volDem;

			tmp		= tmpYtY+tmpDem-tmpNum;
			tmp		/= (i+1.0) * volNum * volBond[j];

			if ( tmp > 99.9 ) 
				ARM_THROW( ERR_INVALID_ARGUMENT, "Pb on computation of correlation matrix for adjustement of convexity ");

			correl->Elt( i, j) = tmp;
		}
	}

	return new 	ARM_VolLInterpol(asOfDate, X_label, Y_label, correl);
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ARM_InfBSSmiledModel
///	Returns: CptFwdCPIRatio
///	Action : return the inflation performance reajusted 
////////////////////////////////////////////////////

ARM_GP_Vector ARM_InfBSSmiledModel::CptFwdCPIRatio(	const double&	numJulian,
													const double&	demJulian,
													const double&	payJulian,
													long			dailyInterpType,
													double			denomFixing,
													ARM_InfIdx*		infIdx){


	const double	rateBase= CC_NS( ARM_Constants, rateBase );
	const double	volBase	= CC_NS( ARM_Constants, volBase);
	const double	dayCount= infIdx->GetDayCount();
	const string	pubLag	= InfData::GetPublishLag( infIdx->GetIndexName().c_str() );


	ARM_Date numDate	(numJulian);
	ARM_Date demDate	(demJulian);
	ARM_Date payDate	(payJulian);

	double	 D_1(1.0/K_YEAR_LEN);
	double	 adj(1.0);
	
	ARM_Date asOfDate	= GetStartDate();
	ARM_Date demPub		= demDate;
	ARM_Date numPub		= numDate;
	demPub.AddPeriod(pubLag);
	if( demDate<asOfDate )  demPub = asOfDate;
	numPub.AddPeriod(pubLag);

	double demLag	= CountYearsWithoutException( dayCount, asOfDate, demPub ); 
	demLag	= demLag>0? demLag: D_1;
	double numLag	= CountYearsWithoutException( dayCount, asOfDate, numPub ); 	
	double payLag	= CountYearsWithoutException( dayCount, asOfDate, payDate ); 

	if( payLag < numLag )
			ARM_THROW( ERR_INVALID_ARGUMENT, "The payment date should be above the Num Publishing date ");

	ARM_GP_Vector result(3,0.0);

	if( demDate < asOfDate )
		result[2] = denomFixing == GETDEFAULTVALUE? FwdCPI( demDate, dailyInterpType ) : denomFixing;
	
	else {


		double volTj		= GetVolatility()->ComputeVolatility( D_1,			demLag )/volBase;
		double volTi		= GetVolatility()->ComputeVolatility( D_1,			numLag )/volBase;
		double volYtY		= GetVolatility()->ComputeVolatility( demLag, numLag-demLag)/volBase;

		double RhoInfIR			= itsCorrelAdjVolCurve	-> ComputeVolatility(demLag, demLag)/volBase;
		double RhoInfIRPaymenti	= itsCorrelAdjVolCurve	-> ComputeVolatility(numLag, numLag)/volBase;
		double RhoInfIRPaymentj	= itsCorrelAdjVolCurve	-> ComputeVolatility(demLag, numLag)/volBase;

		double forward = 0.0, forwardPayment = 0.0;
		double volIR = 0.0, volIRPayment = 0.0;

		ARM_BSModel* IRModel= new ARM_BSModel(GetYCModel()->GetZeroCurve(), itsIrCapVolCurve );
		if( IRModel )
		{
			forward			= IRModel->ExpectedFwdYield( demPub, numPub, numPub )/rateBase;
			forwardPayment	= IRModel->ExpectedFwdYield( numPub, payDate, payDate )/ rateBase;
			volIR			= IRModel->GetVolatility()->ComputeVolatility( demLag, numLag-demLag)/volBase;
			volIRPayment	= IRModel->GetVolatility()->ComputeVolatility( numLag , payLag-numLag ) /volBase;
		}
		
		ARM_GP_Vector Input(14);
		if(fabs(numPub.GetJulian()-demPub.GetJulian()-K_YEAR_LEN) < 7.0)		{			
			Input[0] = demLag;
			Input[1] = numLag;
			Input[2] = numLag-demLag;
			Input[3] = volTj;
			Input[4] = volTi;
			Input[5] = volYtY;
			Input[6] = forward;
			Input[7] = volIR;
			Input[8] = RhoInfIR;
			
			Input[9]  = payLag-numLag;
			Input[10] = forwardPayment;
			Input[11] = volIRPayment;
			Input[12] = RhoInfIRPaymentj;
			Input[13] = RhoInfIRPaymenti;			
		}
		else
		{

			double RhoInfInf = 1.0;
			RhoInfInf = itsCorrelAdjVolCurve -> ComputeVolatility(demLag, demLag);
			
			volYtY	= numLag*volTi*volTi-2*RhoInfInf*sqrt(numLag*demLag)*volTi*volTj+demLag*volTj*volTj;
		
			double tenor		= numLag-demLag;
			double tenorPayment = payLag-numLag;

			if(IRModel)
			{
				double volBond	= 0.0;
				double YF_TERM  = 1.0;
				int k=0;
				double T_kmoins1 = demLag;
				ARM_Date date_kmoins1 = demPub;
				double T_k = T_kmoins1;
				ARM_Date date_k = date_kmoins1;
				date_k.AddYears(YF_TERM);
				while (date_k.GetJulian() <numPub.GetJulian()+7.0)
				{			
					double delta_k = YF_TERM;
					double forward_kmoins1	= IRModel->ExpectedFwdYield( date_kmoins1, date_k, date_k )/ CC_NS( ARM_Constants, rateBase );
					double yf_kmoins1		= (date_kmoins1.GetJulian()	- asOfDate.GetJulian())/K_YEAR_LEN;
					double volIR_kmoins1	= IRModel->GetVolatility()	->ComputeVolatility( yf_kmoins1 , delta_k)/ CC_NS( ARM_Constants, volBase );
					double RhoInfIR_kmoins1 = itsCorrelAdjVolCurve		-> ComputeVolatility(demLag, yf_kmoins1 );

					volBond+= delta_k*forward_kmoins1/(1.0+delta_k*forward_kmoins1)*volIR_kmoins1*RhoInfIR_kmoins1;

					date_kmoins1 = date_k;
					k+=1;
					date_k.AddYears(YF_TERM);
				}

				if( ( numPub.GetJulian()-date_k.GetJulian() ) > 7.0 )
				{
					date_kmoins1 = date_k;
					k+=1;
					date_k=numPub;
					double delta_k = (date_k.GetJulian()-date_kmoins1.GetJulian())/K_YEAR_LEN;
					double forward_kmoins1	= IRModel->ExpectedFwdYield( date_kmoins1, date_k, date_k )/ CC_NS( ARM_Constants, rateBase );
					double yf_kmoins1		= (date_kmoins1.GetJulian() -  asOfDate.GetJulian())/K_YEAR_LEN;
					double volIR_kmoins1	= IRModel->GetVolatility()  -> ComputeVolatility( yf_kmoins1 , delta_k)/ CC_NS( ARM_Constants, volBase );
					double RhoInfIR_kmoins1 = itsCorrelAdjVolCurve		-> ComputeVolatility(demLag, yf_kmoins1 );

					volBond += delta_k*forward_kmoins1/(1.0+delta_k*forward_kmoins1)*volIR_kmoins1*RhoInfIR_kmoins1;
				}			
				RhoInfIR	= 1.0;
				volIR = ( 1.0 + tenor * forward )* volBond /(tenor * forward );
			}
		
			Input[0] = demLag;
			Input[1] = numLag;
			Input[2] = tenor;
			Input[3] = volTj;
			Input[4] = volTi;			
			Input[5] = sqrt(volYtY/tenor);

			Input[6] = forward;			
			Input[7] = volIR;
			Input[8] = RhoInfIR;
			
			/// part for the payment lag
			Input[9]  = tenorPayment;
			Input[10] = forwardPayment;
			Input[11] = volIRPayment;
			Input[12] = RhoInfIRPaymentj;
			Input[13] = RhoInfIRPaymenti;	
		}
		if ( IRModel ) { delete IRModel; IRModel =NULL ;}
		ARM_Vector tmpInput = To_ARM_Vector(Input);
		adj = GetConvAdjustManager()->FwdCPIRatioAdjust(this, &tmpInput);
	
		result[2] = FwdCPI( demDate, dailyInterpType );

	}
	result[1] = FwdCPI( numDate, dailyInterpType );
	result[0] = adj * result[1]/result[2] - 1.0;
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ComputeVolAndStrikeForCap
///	Returns: void
///	Action : Computes the total vol, the pricing strike and the tenor
////////////////////////////////////////////////////
void ARM_InfBSSmiledModel::ComputeVolAndStrikeForCap( 
	double			CPIForward,
	double			strike,
	int				callput,
	ARM_InfIdx*		infIdx,			
	StoreInfoObj&	storeInfo,
	const ARM_Date& numDate,	
	const ARM_Date& denomDate,	
	int				optionType,				
	double&			totalVol,
	double&			pricingStrike,
	double&			tenor ){

	double demLag, numLag, payLag, tenLag ;

	fillTimesForCalculation(	denomDate,  
								numDate, 
								GetModelDateWPublishLag(numDate, infIdx), 
								infIdx,
								demLag,
								numLag,
								payLag, 	
								tenLag);  // check if numdate< paydate
	
	/// test past cash flows or negative tenor
	if( tenLag< 0.0 )
	{
		double data[2] = { -1.0, -1.0 };
		storeInfo.Store( data );
		tenor = -1;
		return ;
	}
	strike /= CC_NS( ARM_Constants, rateBase ) ;
	pricingStrike = pow( 1.0 + strike, tenLag );

	if( pricingStrike<0){
		double data[2] = { -1.0, -1.0 };
		storeInfo.Store( data );
		return ;
	}

	double vol	= ComputeVol( demLag, tenLag, CPIForward, 1.0 + strike )/CC_NS( ARM_Constants, volBase );

	/// in order to track pricing information
	double data[2];
	data[0] = vol*CC_NS( ARM_Constants, volBase );
	data[1] = pricingStrike;
	storeInfo.Store( data );
	totalVol	= vol * sqrt( numLag );
}

void ARM_InfBSSmiledModel::ComputeVolAndStrikeForCap( 
										double CPIForward,
										double strike,
										int callput,
										ARM_InfIdx* infIdx,			
										StoreVolInfoWithNewCF* storeInfo,
										const ARM_Date& numDate,	
										const ARM_Date& denomDate,	
										double& totalVol,
										double& pricingStrike,
										double& tenor)
{
		
	double demLag, numLag, payLag, tenLag ;

	fillTimesForCalculation(	denomDate,  
								numDate, 
								GetModelDateWPublishLag(numDate, &*infIdx), 
								&*infIdx,
								demLag,
								numLag,
								payLag, 	
								tenLag);  // check if numdate< paydate
	
	/// test past cash flows or negative tenor
	if( tenLag< 0.0 )
	{
		double data[2] = { -1.0, -1.0 };
		storeInfo->Store( data );
		tenor = -1;
		return ;
	}
	strike /= CC_NS( ARM_Constants, rateBase ) ;
	pricingStrike = pow( 1.0 + strike, tenLag );

	if( pricingStrike<0){
		double data[2] = { -1.0, -1.0 };
		storeInfo->Store( data );
		return ;
	}

	double vol	= ComputeVol( demLag, tenLag, CPIForward, 1.0 + strike )/CC_NS( ARM_Constants, volBase );

	/// in order to track pricing information
	double data[2];
	data[0] = vol*CC_NS( ARM_Constants, volBase );
	data[1] = pricingStrike;
	storeInfo->Store( data );
	totalVol	= vol * sqrt( numLag );
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: SingleAssetOptionPrice
///	Returns: double
///	Action : compute the oplet price
////////////////////////////////////////////////////
double ARM_InfBSSmiledModel::SingleAssetOptionPrice(	double			CPIForward,
														double			strike,
														int				callput,
														double			discounting, 
														ARM_InfIdx*		infIdx,
														ARM_Object*		optionContext,
														StoreInfoObj&	storeInfo)
{	
	double totalVol = 0., pricingStrike = 0., tenor = 0.;
     
	ARM_InfCapFloorContext* infCapFloorContext = dynamic_cast<ARM_InfCapFloorContext*>(optionContext);

	if(infCapFloorContext)
	{
	
		ComputeVolAndStrikeForCap(		CPIForward,	
										strike,	
										callput,
										infIdx,	
										storeInfo,	
										infCapFloorContext->GetNumDate(), 
										infCapFloorContext->GetDenomDate(),
										callput,
										totalVol,
										pricingStrike, 
										tenor );
	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"unknown type of option... price only inflation caps or floors");

	/// test for past deal
	if( tenor < 0 )
		return 0.0;

	/// test for zero strike
	if( pricingStrike == 0 )
		pricingStrike = K_NEW_DOUBLE_TOL;

	/// test for negative strike that have no sense
	if( pricingStrike < 0 )
	{
		/// if the strike is negative... try to return the intrinsic value
		/// if positive!
		double intrinsic = (CPIForward-pricingStrike) * callput;
		if( intrinsic >= 0 )
			return intrinsic * CC_NS( ARM_Constants, rateBase )*discounting;

		/// otherwise return an exception
		else
		{
			char msg[255];
			sprintf( msg, "%s: strike is equal to %f should be at least positive with the model assumption Probably you want to use a swap!.. Please advise",
				ARM_USERNAME.c_str(), pricingStrike );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		}
	}
	
	/// compute the BS value
	double result;
	result = BlackSholes_Formula( CPIForward, totalVol, discounting, 1.0+strike / CC_NS( ARM_Constants, rateBase ) , callput );
	


	return CC_NS( ARM_Constants, rateBase ) * result;	
}

double ARM_InfBSSmiledModel::SingleAssetOptionPrice(	double CPIForward,
														double strike,
														int callput,
														double discounting, 
														ARM_InfIdx* infIdx,
														InfCapFloorContext* infCapFloorContext,
														StoreVolInfoWithNewCF* storeInfo)
{	
	double totalVol = 0., pricingStrike = 0., tenor = 0.;
     
		ComputeVolAndStrikeForCap(		CPIForward,	
										strike,	
										callput,
										infIdx,	
										storeInfo,	
										infCapFloorContext->GetNumDate(), 
										infCapFloorContext->GetDenomDate(),
										totalVol,
										pricingStrike, 
										tenor);
	
	/// test for past deal
	if( tenor < 0 )
		return 0.0;

	/// test for zero strike
	if( pricingStrike == 0 )
		pricingStrike = K_NEW_DOUBLE_TOL;

	/// test for negative strike that have no sense
	if( pricingStrike < 0 )
	{
		/// if the strike is negative... try to return the intrinsic value
		/// if positive!
		double intrinsic = (CPIForward-pricingStrike) * callput;
		if( intrinsic >= 0 )
			return intrinsic * CC_NS( ARM_Constants, rateBase )*discounting;

		/// otherwise return an exception
		else
		{
			char msg[255];
			sprintf( msg, "%s: strike is equal to %f should be at least positive with the model assumption Probably you want to use a swap!.. Please advise",
				ARM_USERNAME.c_str(), pricingStrike );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		}
	}
	
	/// compute the BS value
	double result;
	result = BlackSholes_Formula( CPIForward, totalVol, discounting, 1.0+strike / CC_NS( ARM_Constants, rateBase ), callput );
	


	return CC_NS( ARM_Constants, rateBase ) * result;	

}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: fillTimesForCalculation
///	Returns: void
///	Action : Computes some relevant dates
////////////////////////////////////////////////////

void ARM_InfBSSmiledModel::fillTimesForCalculation(	const ARM_Date& denomDate,  
													const ARM_Date& numDate, 
													const ARM_Date& paymentDate, 
													ARM_InfIdx* infIdx,
													double& demLag,
													double& numLag,
													double& payLag, 	
													double& tenLag)
{
	int dayCount					= infIdx->GetDayCount();
	ARM_Date lastKnownDate			= GetVolatility()->GetLastKnownDate();
	ARM_Date modelAsOfDate			= GetStartDate();
	ARM_Date denomDatewPublishLag	= GetModelDateWPublishLag( denomDate, infIdx );
	ARM_Date numDatewPublishLag		= GetModelDateWPublishLag( numDate, infIdx );
	if(paymentDate.GetJulian()  < numDatewPublishLag.GetJulian() )
			ARM_THROW( ERR_INVALID_ARGUMENT, "The payment date should be above the Num Publishing date ");
	demLag							= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );
	numLag							= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
	payLag							= CountYearsWithoutException( dayCount, modelAsOfDate, paymentDate );
	tenLag							= numLag-demLag;
	if ( denomDatewPublishLag < modelAsOfDate ){
		demLag	= 1.0/K_YEAR_LEN;
		tenLag	= numLag;
	}

};



////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: ComputeVol
///	Returns: double
///	Action : compute the implied lognormal vol
////////////////////////////////////////////////////
double ARM_InfBSSmiledModel::ComputeVol(	double	maturity, 
											double	tenor,
											double	fwd, 
											double	strike, 
											int		mode )
{
	double pricingStrike = pow(strike,tenor);

	if ( maturity<0 ){
		tenor	= tenor+maturity>0?tenor+maturity:1.0/K_YEAR_LEN;
		maturity= 1.0/K_YEAR_LEN;
	}
	return ARM_BSSmiledModel::ComputeVol(maturity, tenor, fwd, pricingStrike, mode ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: DiscountedCPI
///	Returns: double
///	Action : not implemented
////////////////////////////////////////////////////
double ARM_InfBSSmiledModel::DiscountedCPI(	const ARM_Date&		resetDate, 
											const ARM_Date&		paymentDate, 
											long				dailyInterpType,
											const string&		CPILag,
											const string&		DCFLag,
											ARM_InfIdx*			infIdx ){
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : DiscountedCPI not implemented" );
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: GetModelTimeWPublishLag
///	Returns: double
///	Action : not implemented
////////////////////////////////////////////////////

double			ARM_InfBSSmiledModel::GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx){
	int dayCount				= GetInfFwdCurv()->GetMonthlyInterpType();
	ARM_Date tmpDateWPublishLag	= GetModelDateWPublishLag( date, infIdx );
	return CountYearsWithoutException( dayCount, GetStartDate(), tmpDateWPublishLag	); 
}
////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: GetModelDateWPublishLag
///	Returns: ARM_Date
///	Action : not implemented
////////////////////////////////////////////////////

ARM_Date		ARM_InfBSSmiledModel::GetModelDateWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx){
	string indexName	= infIdx->GetIndexName();
	string publishLag	= InfData::GetPublishLag( indexName.c_str() );
	ARM_Date tmpDate( date );
	tmpDate.AddPeriod( publishLag );
	return tmpDate;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSSmiledModel
///	Routine: View
///	Returns: void
///	Action : return the viewer
////////////////////////////////////////////////////

void ARM_InfBSSmiledModel::View(char* id , FILE* ficOut ){

	 FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;


	int nb;
	CC_Ostringstream	os;
	const int LAG = 13;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(58)<<"INFLATION BSSMILED MODEL\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";	
	os	<<"Characteristics\n";
	os	<<"==============="<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<"As Of date:";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<GetStartDate();
	os	<< "\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<"Currency:";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
	os	<< "\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<"Rate Index:";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<GetZeroCurve()->GetExternalIndex();
	os	<< "\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<"Inflation Index:";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(20)<<GetInfFwdCurv()->GetIndexName();
	os	<< "\n";


	os	<< "\n\n";	
	os	<<"Zero Curve\n";
	os	<<"==========="<<"\n\n\n";

  	os	<<"\n==> Discount Factor Curve"<<"\n\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Dates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"NbDays";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Rates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DFs";	
	os	<<"\n";	

	nb = GetZeroCurve()->GetDateTerms()->size();
	for( int i = 0; i<nb ; i++){
		int			tmp		= (int) GetZeroCurve()->GetDateTerms()->Elt(i);
		ARM_Date	tmpDate = GetStartDate();
		tmpDate .AddDays(tmp);
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmpDate;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<tmp;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetZeroCurve()->GetZeroRates()->Elt(i);
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetZeroCurve()->GetDiscountFactors()->Elt(i);	
		os	<<"\n";	
	}
        

	os	<< "\n\n";	
	os	<<"Inflation Curve\n";
	os	<<"==============="<<"\n\n\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Dates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"CPIs";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"ZCs";	
	os	<<"\n";


	nb = GetInfFwdCurv()->GetExpiryTermsVec().size();
	for( i = 0; i < nb; ++i )	{
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<ARM_Date( GetInfFwdCurv()->GetExpiryTermsVec()[i] );
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetInfFwdCurv()->GetCPIValuesVec()[i];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetInfFwdCurv()->GetZCValuesVec()[i];
		os	<<"\n";
	}

	os	<< "\n\n";	
	os	<<"Sabr Inflation\n";
	os	<<"=============="<<"\n\n\n";

   	os	<<"==> Vol ATM"<<"\n\n";

	int xNb = GetVolatility()->GetExpiryTerms()->size();
	int yNb = GetVolatility()->GetStrikes()->size();

    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( int j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetVolatility()->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetVolatility()->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetVolatility()->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}
	os	<< "\n";
	
   	os	<<"==> Rho"<<"\n\n";
	xNb = GetRho()->GetExpiryTerms()->size();
	yNb = GetRho()->GetStrikes()->size();
    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetRho()->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetRho()->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetRho()->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}
	os	<< "\n";

   	os	<<"==> Nu"<<"\n\n";
	xNb = GetNu()->GetExpiryTerms()->size();
	yNb = GetNu()->GetStrikes()->size();
    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetNu()->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetNu()->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetNu()->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}
	os	<< "\n";

   	os	<<"==> Beta"<<"\n\n";
	xNb = GetBeta()->GetExpiryTerms()->size();
	yNb = GetBeta()->GetStrikes()->size();
    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetBeta()->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetBeta()->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<GetBeta()->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}

	os	<< "\n\n";	
	os	<<"Vol IR Cap Atm\n";
	os	<<"=============="<<"\n\n\n";

  	os	<<"==> Vol"<<"\n\n";
	xNb = itsIrCapVolCurve->GetExpiryTerms()->size();
	yNb = itsIrCapVolCurve->GetStrikes()->size();
    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsIrCapVolCurve->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsIrCapVolCurve->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsIrCapVolCurve->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}

	os	<< "\n\n";	
	os	<<"Correl INF/IR\n";
	os	<<"============="<<"\n\n\n";

  	os	<<"==> denom Pub Date: YtY / ResetDate: Libor"<<"\n\n";
	xNb = itsCorrelVolCurve->GetExpiryTerms()->size();
	yNb = itsCorrelVolCurve->GetStrikes()->size();
    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsCorrelVolCurve->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsCorrelVolCurve->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsCorrelVolCurve->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}
	os	<<"\n";

  	os	<<"==>Pub Date: Cpi / ResetDate: Bond "<<"\n\n";
	xNb = itsCorrelAdjVolCurve->GetExpiryTerms()->size();
	yNb = itsCorrelAdjVolCurve->GetStrikes()->size();
    os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"";
	for ( j =0; j< yNb; ++j)
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsCorrelAdjVolCurve->GetStrikes()->Elt(j);
	os	<<"\n";	
	for( i= 0; i< xNb; ++i){
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsCorrelAdjVolCurve->GetExpiryTerms()->Elt(i);
		for( j = 0; j< yNb; j++)
			os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<itsCorrelAdjVolCurve->GetVolatilities()->Elt(i,j);
		os	<<"\n";
	}

    fprintf(fOut, (os.str()).c_str()  );

	if ( ficOut == NULL )
		fclose(fOut);

}


CC_END_NAMESPACE()


/*---------------------------------------------------------------*/
/*---- End Of File ----*/

