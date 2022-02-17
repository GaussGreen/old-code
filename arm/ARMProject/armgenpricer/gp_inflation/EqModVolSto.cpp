/*----------------------------------------------------------------------------*/
 
/*! \file EqModVolSto.h
 * Copyright (c) CDC IXIS CM February 2007 Paris
 *
 *  \HW model on IRS Part... Heston model on each cpi.
 *
 *
 *	\author  Mathieu Bernardo
 *	\version 1.0
 *	\date February 2007
 */

/*----------------------------------------------------------------------------*/


#include "gpmodels\HWSV.h"
#include "gpinflation\EqModVolSto.h"

CC_BEGIN_NAMESPACE(ARM)
/*
ARM_ModifiedRiccati::ARM_ModifiedRiccati(	double eta,		double kappa_x,		double epsilon_x,	double vol_x,
											double kappa_y,	double epsilon_y,	double dt ){

	itsA  = epsilon_y*epsilon_y/2;
	itsB1 = - kappa_y;
	itsB2 =	- epsilon_x * epsilon_y;
	itsC  = 1.0;
	itsD  = 1.0;

	itsdt = dt;
	itsReV= -epsilon_x*epsilon_x/2;
	itsImV=	kappa_x;

	itsReG= -vol_x*epsilon_x*epsilon_x/2;
	itsImG= -eta;
};


*/


////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: ARM_EQHWSV
///	Returns: void
///	Action : Constructor
////////////////////////////////////////////////////

ARM_EQHWSV::ARM_EQHWSV( ARM_ZeroCurve*						zcCurve,	
					    const ARM_EQHWSV_ModelParams*		modelParams, 
					    const ARM_EQHWSV_NumMethods*		numMethods): 
						ARM_HWSV1F( ARM_ZeroCurvePtr(&*dynamic_cast< ARM_ZeroCurve*>(zcCurve->Clone() ) ), modelParams, numMethods->itsMaxDecay){
	itsModelParams		= ARM_EQHWSV_ModelParamsPtr (&* dynamic_cast<ARM_EQHWSV_ModelParams		*>	( modelParams	->Clone() ) );
	itsNumMethods		= ARM_EQHWSV_NumMethodsPtr	(&* dynamic_cast<ARM_EQHWSV_NumMethods		*>	( numMethods	->Clone() ) );

}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: ARM_EQHWSV
///	Returns: void
///	Action : Recopy Constructor
////////////////////////////////////////////////////

ARM_EQHWSV::ARM_EQHWSV( const  ARM_EQHWSV & rhs ):ARM_HWSV1F(rhs){

	itsModelParams	= rhs.itsModelParams;
	itsNumMethods	= rhs.itsNumMethods;
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: ValidateModelParams
///	Returns: boolean
///	Action : Check if all the model paramters are provided
////////////////////////////////////////////////////
bool ARM_EQHWSV::ValidateModelParams(const ARM_ModelParams& params) const
{
    if(!params.DoesModelParamExist(ARM_ModelParamType::MeanReversion)	||
	   !params.DoesModelParamExist(ARM_ModelParamType::Volatility)		||
	   !params.DoesModelParamExist(ARM_ModelParamType::CompoundVol)		||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolOfVol)		||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolMeanReversion)||
	   !params.DoesModelParamExist(ARM_ModelParamType::Correlation))
    {
       ARMTHROW(ERR_CONDITION_NOT_MEET, "At least 1 Model Param is not of a good type!");
    }
	return true;
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: toString
///	Returns: string
///	Action : provide the view of the class
////////////////////////////////////////////////////

string ARM_EQHWSV::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Equity Hull White Stochastic Volatility Model - EQHWSM \n";
    os << indent << "------------------------------------------------------\n\n";


	os << indent << "==> Model Parameters \n";
    os << itsModelParams->toString(indent);
    os << "\n\n";

	os << indent << "==> Numeric Methods \n";
    os << itsNumMethods->toString(indent);
    os << "\n\n";

    return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : compute a swaption price
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EQHWSV::VanillaSwaption(	const string& curveName,
											double evalTime,
											double swapResetTime,
											const ARM_GP_Vector& fixNotional,
											const ARM_GP_Vector& floatNotional,
											double floatStartTime,
											double floatEndTime,
											const ARM_GP_Vector& floatResetTimes,
											const ARM_GP_Vector& floatStartTimes,
											const ARM_GP_Vector& floatEndTimes,
											const ARM_GP_Vector& floatIntTerms,
											const ARM_GP_Vector& fixPayTimes,
											const ARM_GP_Vector& fixPayPeriods,
											const ARM_GP_Matrix& strikesPerState,
									        int callPut,
											const ARM_PricingStatesPtr& states,
											bool isConstantNotional,
											bool isConstantSpread,
											bool isConstantStrike) const{

	ReduceSwaptionSet(	evalTime,
						floatEndTime,
						fixPayTimes[fixPayTimes.size()-1 ],
						fixNotional,
						floatNotional,
						strikesPerState,
						isConstantNotional,
						isConstantSpread,
						isConstantStrike);
	

	double strike		= strikesPerState(0,0);
	double notional		= fixNotional[0];

	double freezedPrice	= CptOswDriftFreezedNonVolSto(	evalTime,
														swapResetTime,
														floatStartTime, 
														floatEndTime,
														strike,
														notional,
														callPut,
														fixPayTimes,
														fixPayPeriods);

	double exactPrice	= CptOswCorrelNull			(	evalTime,
														swapResetTime,
														floatStartTime, 
														floatEndTime,
														strike,
														notional,
														callPut,
														fixPayTimes,
														fixPayPeriods);

	ARM_VectorPtr values(&* new  ARM_GP_T_Vector<double>(1,freezedPrice - exactPrice));
	return values;

}
////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: ReduceSwaptionSet
///	Returns: void
///	Action : exclusion of the none standard swaption
////////////////////////////////////////////////////												
void ARM_EQHWSV::ReduceSwaptionSet	(	const double & evalTime,
									 	const double & floatEndTime,
										const double & fixPayTime,
										const ARM_GP_Vector& fixNotional,
										const ARM_GP_Vector& floatNotional,
										const ARM_GP_Matrix& strikesPerState,
										bool isConstantNotional,
										bool isConstantSpread,
										bool isConstantStrike ) const{
	if ( !isConstantNotional )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Notional should be constant for calibration" );
	if( fixNotional[0] != floatNotional[0] )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The notional of the 2 legs should be identical" );
	if ( !isConstantSpread )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : spread should be constant for calibration" );
	if ( !isConstantStrike )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : strike should be constant for calibration" );
	if(evalTime > K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : EQHWSV model doesn't support analytical swaption evaluation at future date" );
	if(strikesPerState.rows() != 1)		
		ARM_THROW( ERR_INVALID_ARGUMENT, " : strike should be an unique double" );
	if(fabs(floatEndTime - fixPayTime )>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fixed leg end dates are not matching in HWSV1F swaption pricing" );

	double strike = strikesPerState(0,0);
	for(int i=1;i<strikesPerState.cols(); i++){
		if( strikesPerState(0,i) != strike )
			ARM_THROW( ERR_INVALID_ARGUMENT, " : strike should be constant for pricing of standard swaption" );

	}
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptFreezedQuantities
///	Returns: void
///	Action : compute the freezed drift and the equivalent strike
////////////////////////////////////////////////////////////////
void ARM_EQHWSV::CptFreezedQuantities(	const double & evalTime,
										const double & startTime, 
										const double & endTime,
										const double & strike,
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
										double & newStrike,
										double & drift) const{	
/*	newStrike	= 0.0;
	drift		= 0.0;

	double				flow;
	double				lambda		= GetMrs();
	double				expT0		= exp(-lambda *(startTime - evalTime)/K_YEAR_LEN);
	ARM_ZeroCurvePtr	zcCurve		= GetZeroCurve();
	ARM_GP_Vector		coeff		= fixPayPeriods;
	int					nbPeriods	= fixPayTimes.size();

	coeff *= strike;
	coeff[nbPeriods-1] +=1;

	for(int i = 0; i<nbPeriods; i++)	{
		flow		= coeff[i] * zcCurve->DiscountPrice( (fixPayTimes[i]-evalTime)/K_YEAR_LEN);
		newStrike	+= flow;
		drift		+= flow * (exp (-lambda*(fixPayTimes[i]-evalTime)/K_YEAR_LEN) - expT0);
		}

	drift			=	drift/( lambda * newStrike );
	newStrike		=	zcCurve->DiscountPrice((startTime-evalTime)/K_YEAR_LEN) / newStrike;*/
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptBasketPrice
///	Returns: double
///	Action : Compute the price of the basket price
////////////////////////////////////////////////////////////////

double ARM_EQHWSV::CptBasketPrice(	const ARM_GP_Vector	& polCoef, 
									const ARM_GP_Vector	& expCoef, 
									const double		& x){	
	double value=0;
	for( int i	= 0; i<expCoef.size(); i++)	
		value	+= polCoef[i]*exp( expCoef[i]*x );
	return value;
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: Cpt_ExpPol_Coef
///	Returns: void
///	Action : Compute the PolynomCoeff (polCoef) and the ExponentialCoeff (expCoef) 
////////////////////////////////////////////////////////////////

void ARM_EQHWSV::CptExpPolCoef(	ARM_GP_Vector			& polCoef, 
								ARM_GP_Vector			& expCoef,
								const ARM_GP_Vector		& fixPayPeriods,
								const ARM_GP_Vector		& fixPayTimes,
								const double			& floatStartTime,
								const double			& evalTime,
								const double			& swapResetTime,
								const double			& strike,
								const double			& variance) const {	

/*	polCoef.clear();
	expCoef.clear();

	double				tmp;

	ARM_ZeroCurvePtr	zcCurve		= GetZeroCurve();
	double				DF			= zcCurve->DiscountPrice(( floatStartTime-evalTime)/K_YEAR_LEN);

	ARM_GP_Vector		coeff		= fixPayPeriods;
	int					nbPeriods	= fixPayTimes.size();
	coeff				*=	strike;
	coeff[nbPeriods-1]	+=	1;

	double				kappa		= GetMrs();
	double				ExpT0		= exp(-kappa*(floatStartTime-evalTime)/K_YEAR_LEN);
	double				ExpTi;

	for( int i	= 0; i<nbPeriods; i++)	{

		ExpTi	=	exp(-kappa*(fixPayTimes[i]-evalTime)/K_YEAR_LEN);
		tmp		=	(ExpT0-ExpTi)/kappa;
		tmp		*=	tmp;
		tmp		*=	variance;
		expCoef.push_back(sqrt(tmp) );
	
		tmp		=	exp(-0.5*tmp);
		tmp		*=	coeff[i] * zcCurve->DiscountPrice( (fixPayTimes[i]-evalTime) /K_YEAR_LEN)/DF;
		polCoef.push_back(tmp);
	}*/
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptFrontiere
///	Returns: double
///	Action : Compute the frontiere
////////////////////////////////////////////////////////////////
/*
class BasketPriceToInverse : public ARM_GP::UnaryFunc<double,double>{

public: 
	BasketPriceToInverse(	const ARM_GP_Vector& polCoef, 	const ARM_GP_Vector& expCoef):itsPolCoef(polCoef),itsExpCoef(expCoef){ };
	virtual double operator() (double z) const{ return ARM_EQHWSV::CptBasketPrice(itsPolCoef, itsExpCoef, z);	}

private:
	ARM_GP_Vector	itsPolCoef;
	ARM_GP_Vector	itsExpCoef;
};

double ARM_EQHWSV::CptFrontiere(	const ARM_GP_Vector & polCoef, 
									const ARM_GP_Vector & expCoef,
									double x ) const {
	double			root = 0.0;

	BasketPriceToInverse func( polCoef, expCoef );

	UnaryFuncWithNumDerivative<double> funcWithDeriv(func);

	T_NewtonRaphsonSolver<UnaryFuncWithNumDerivative<double> > solver(funcWithDeriv,x,DEF_PRECISION,DEF_PRECISION);

	solver.setInitialGuess(0.0);

	root= solver.Solve();
	
	return root;
}

*/
////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: GenerateCirModel 
///	Returns: ARM_CIRModelParams
///	Action : Generate the CIR model corresponding to 
	
///			Y_t		= p_t²*exp(2k.t)*V_t
///			dY_t	= H_t.( L_t-Y_t) +Q_t.Sqrt(Y_t).dW_t

///	with	H_t		= h-2k
///			L_t		= h. p_t² . exp(2k.t)/( h-2k)
///			Q_t		= p_t.exp(k.t).q	
////////////////////////////////////////////////////////////////


ARM_CIRModelParams	 ARM_EQHWSV::GenerateCirModel(const double& time) const{

/*	double				tp;
	double				tv;
	double				k			= GetMrs();
	double				h			= GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetValue(0.0);

	ARM_GP_VectorPtr	tmp( new ARM_GP_Vector (itsModelParams->GetSchedule()) );
	ARM_GP_VectorPtr	sch( new ARM_GP_Vector () );
	ARM_GP_Vector*		tmpMrs		= new ARM_GP_Vector();
	ARM_GP_Vector*		tmpVol		= new ARM_GP_Vector();
	ARM_GP_Vector*		tmpLtv		= new ARM_GP_Vector();

	for(int i = 0; i< tmp->size();i++){
		if( tmp->Elt(i) > time )
			sch->push_back(tmp->Elt(i)-time);
	}
	sch->insert(sch->begin(),0.0);


	ARM_CurvePtr vol ( GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve() );
	ARM_CurvePtr mrs ( GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam().GetCurve() );
	ARM_CurvePtr vov ( GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve() );

	tmp->clear();
	for( i = 0; i < sch->size(); i++){
		tp				=	sch->Elt(i)+time;
		tv				=	exp(k*tp/K_YEAR_LEN) * vol->Interpolate(tp);
		tmpVol			->	push_back( tv* vov->Interpolate(tp) );
		tmpMrs			->	push_back( h-2.0*k );
		tmpLtv			->	push_back( tv*tv*h/(h-2.0*k) );
	}

	ARM_ModelParamVector modelParams;
	modelParams.push_back(	new ARM_CurveModelParam ( ARM_ModelParamType::Volatility,	tmpVol, dynamic_cast<ARM_GP_Vector*> (sch->Clone() ) ) );
	modelParams.push_back(  new ARM_CurveModelParam ( ARM_ModelParamType::MeanReversion,tmpMrs, dynamic_cast<ARM_GP_Vector*> (sch->Clone() ) ) );
	modelParams.push_back(  new ARM_CurveModelParam ( ARM_ModelParamType::LongTermVol,	tmpLtv, dynamic_cast<ARM_GP_Vector*> (sch->Clone() ) ) );
	*/
	return ARM_CIRModelParams ();
}



double ARM_EQHWSV::CptLocalVariance( const double& Tstart, const double& Tend) const {

/*	ARM_CurveModelParam	vol			= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam();
	double				lambda		= GetMrs();

	double				timeStart,	timeEnd;
	double				expStart,	expEnd;
	double				variance,	flow;

	ARM_GP_Vector		tmp = ( (ARM_EQHWSV_ModelParams* ) GetModelParams() )->GetSchedule();
	ARM_GP_Vector		schedule;

	for( int i = 0 ; i<tmp.size(); i++){
		if( Tstart< tmp[i] && tmp[i]<Tend )
			schedule.push_back(tmp[i]);
		else if( tmp[i]>Tend)
			break;
	}
	if( schedule[schedule.size()-1]!= Tend)
		schedule.push_back(Tend);

	expStart	= 1;
	timeStart	= Tstart;
	variance	= 0;
	for( i = 0 ; i<schedule.size(); i++){
		timeEnd		= 	schedule[i];
		expEnd		=	exp(2*lambda*(timeEnd-Tstart)/K_YEAR_LEN);
		flow		=	vol.GetValue(timeStart);
		flow		*=	flow;
		flow		*=	expEnd - expStart ;
		variance	+=	flow;
		timeStart	=	timeEnd;
		expStart	=	expEnd;
	}*/
	return 1.0;// variance/(2*lambda);
}


///////////////////////////////////////////////////////////////////////////////////
/// Struct : VarianceDistribution
///	Routine: VarianceDistribution
///	Returns: void
///	Action : constructor
///
///////////////////////////////////////////////////////////////////////////////////

VarianceDistribution::VarianceDistribution(	ARM_CIRModelParams		& cirModel, 
											const double			& time, 
											const double			& x0){
////WARNING TOTO
/*	cirModel.OptDistribSet(time, x0,40);
//	DistribSet	tmp = cirModel.GetDistribSet(time, x0,40);

	itsWeight	=	tmp.itsWeight;
	itsPoint	=	tmp.itsPoint;
	for( int i	= 0; i< itsPoint.size(); i++)
		itsDensity.push_back(cirModel.CptDensity_Bond( time,x0,itsPoint[i]) );*/
}

VarianceDistribution::VarianceDistribution(	const ARM_GP_Vector	& weight, 
											const ARM_GP_Vector	& point, 
											const ARM_GP_Vector	& density):itsWeight(weight),itsPoint(point),itsDensity(density){
	
	if( itsWeight.size() != itsPoint.size() ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Nb point != Nb weight" );

	if( itsDensity.size() != itsPoint.size() ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Nb densiy != Nb weight" );
}


///////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptVarDist
///	Returns: double
///	Action : Compute the distribution of Int( p_s²*exp(2k.s)*V_s, s= tStart..tEnd, V_tStart= V0)
///
///	Remark : if the vol of vol is equal to zero v_s remains concentrate around 1
///////////////////////////////////////////////////////////////////////////////////
/*
VarianceDistribution ARM_EQHWSV::CptVarDist( const double& tStart, const double& tEnd, const double & x0) const {

	VarianceDistribution* var = NULL;

	double vol	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetValue(0.0);

	if ( vol == 0.0  )
		return VarianceDistribution ( ARM_GP_Vector(1,1.0), ARM_GP_Vector(1,CptLocalVariance( tStart,tEnd)) ,ARM_GP_Vector(1,1.0));
	else
		return VarianceDistribution ( GenerateCirModel(tEnd-tStart), tEnd-tStart, x0);


}
*/
///////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptOswCorrelNull
///	Returns: double
///	Action : Compute the distribution price of a swaption for a HW Vol sto in a zero correlation case
///
///////////////////////////////////////////////////////////////////////////////////

double ARM_EQHWSV::CptOswCorrelNull		(	const	double			& evalTime,
											const	double			& swapResetTime,
											const	double			& floatStartTime, 
											const	double			& floatEndTime,
											const	double			& strike,
											const	double			& notional,
											const	int				& callPut,
											const	ARM_GP_Vector	& fixPayTimes,
											const	ARM_GP_Vector	& fixPayPeriods) const{	
/*	double			price;
	double			sum; 
	double			target;
	double			tmpStrike;
	double			tmpVol;
	double			tmpFwd;

	ARM_GP_Vector	polCoef;
	ARM_GP_Vector	expCoef;

	double  DF		=	GetZeroCurve()->DiscountPrice( (floatStartTime-evalTime) /K_YEAR_LEN);
	double	kappa	=	GetMrs();
	double	initVar	=	GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetValue(evalTime);
	initVar			*=	exp(kappa*evalTime);
	initVar			*=	initVar;
	VarianceDistribution varDistrib = CptVarDist( evalTime, swapResetTime,initVar);

	ARM_GP_Vector weight	= varDistrib.itsWeight;
	ARM_GP_Vector point		= varDistrib.itsPoint;
	ARM_GP_Vector density	= varDistrib.itsDensity;

	price = 0.0;
	for( int i =0 ; i < point.size(); i++){
		
		CptExpPolCoef(	polCoef, 
						expCoef,
						fixPayPeriods,
						fixPayTimes,
						floatStartTime,
						evalTime,
						swapResetTime,
						strike,
						point[i] );	

		target	= CptFrontiere(	polCoef, expCoef,	1.0 );

		sum = 0.0;
		for ( int j= 0; j<fixPayTimes.size(); j++){
			tmpStrike	=	polCoef[j]*exp( expCoef[i]*target) ;
			tmpVol		=	expCoef[j];
			tmpFwd		=	polCoef[j]*exp(0.5*expCoef[i]*expCoef[j]);

			if(callPut==K_CALL)
				sum +=  BlackSholes_Formula(tmpFwd,tmpVol,DF,tmpStrike,K_PUT);
			else
				sum +=  BlackSholes_Formula(tmpFwd,tmpVol,DF,tmpStrike,K_CALL);
		}
		price += weight[i] * density[i] *sum;
	}
*/
	return 1.0;//notional*price;
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptOswDriftFreezedNonVolSto
///	Returns: double
///	Action : Compute the swaption price in the Hull white with a const variance process and the freezed drift1
////////////////////////////////////////////////////////////////

double ARM_EQHWSV::CptOswDriftFreezedNonVolSto(	const	double			& evalTime,
												const	double			& swapResetTime,
												const	double			& floatStartTime, 
												const	double			& floatEndTime,
												const	double			& strike,
												const	double			& notional,
												const	int				& callPut,
												const	ARM_GP_Vector	& fixPayTimes,
												const	ARM_GP_Vector	& fixPayPeriods) const{	
/*	double newStrike;
	double drift;
	CptFreezedQuantities(	evalTime, floatStartTime, floatEndTime,	strike, fixPayTimes, fixPayPeriods,	newStrike,	drift);

	double stdDev	=	fabs(drift)*sqrt( CptLocalVariance(evalTime,swapResetTime) );
	double DF		=	GetZeroCurve()->DiscountPrice((floatStartTime-evalTime)/K_YEAR_LEN);
	double price;

	if(callPut==K_CALL)
		price	= BlackSholes_Formula(1.0,stdDev,DF,newStrike,K_PUT);
	else
		price	= BlackSholes_Formula(1.0,stdDev,DF,newStrike,K_CALL);*/

	return 1.0;// notional*price/newStrike;
}





CC_END_NAMESPACE()


/*---------------------------------------------------------------*/
/*---- End Of File ----*/

