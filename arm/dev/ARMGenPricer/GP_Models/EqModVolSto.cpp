/*----------------------------------------------------------------------------*/
 
/*! \file EqModVolxto.h
 * Copyright (c) CDC IxIx CM February 2007 Parix
 *
 *  \HW model on IRx Part... Hexton model on each cpi.
 *
 *
 *	\author  Mathieu Bernardo
 *	\verxion 1.0
 *	\date February 2007
 */

/*----------------------------------------------------------------------------*/

#include "gpbase/functor.h"


#include "gpinfra\pricingmodelir.h"
#include "gpinfra\curvemodelparam.h"
#include "gpinfra\argconvdefault.h"
#include "gpinfra\zccrvfunctor.h"

#include "gpmodels\EqModVolsto.h"
#include "gpmodels\ModelParamsHW1F.h"
#include "gpmodels\HW.h"
#include "gpmodels\HWsV1F.h"

#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/ExpRiccati.h"

#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"

CC_BEGIN_NAMESPACE(ARM)

const int X_VARIABLE	= 0; // for xt
const int V_VARIABLE	= 1; // for Vt
const int PHI_VARIABLE	= 2; // for phit


const double DEF_PRECISION = 1e-8;

double	ARM_EQHWSV::itsDilatation;

ARM_ModifiedRiccati::ARM_ModifiedRiccati(	double eta,		double kappa_x,		double epsilon_x,	double vol_x,
											double kappa_y,	double epsilon_y,	double correl,		double dt ){

	itsA  = epsilon_y*epsilon_y/2;
	itsB1 = - kappa_y;
	itsB2 =	- epsilon_x * epsilon_y * sqrt( 1- correl* correl) ;
	itsC1 = 0.0;
	itsC2 = 1.0;
	itsD  = 1.0;

	itsdt = dt;
	itsReV= -epsilon_x*epsilon_x/2;
	itsImV=	kappa_x;

	itsReG= -vol_x*epsilon_x*epsilon_x/2;
	itsImG= -eta;
}

////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: ARM_EQHWSV
///	Returnx: void
///	Action : constructor
////////////////////////////////////////////////////

ARM_EQHWSV::ARM_EQHWSV( ARM_ZeroCurve*						zcCurve,	
					    const ARM_EQHWSV_ModelParams*		modelParams, 
					    const ARM_EQHWSV_NumMethods*		numMethods,
						const double&						dilatation): 
ARM_HWSV1F( ARM_ZeroCurvePtr(&*dynamic_cast< ARM_ZeroCurve*>(zcCurve->Clone() ) ), modelParams, numMethods->itsMaxDecay){
	itsModelParams		= ARM_EQHWSV_ModelParamsPtr (&* dynamic_cast<ARM_EQHWSV_ModelParams		*>	( modelParams	->Clone() ) );
	itsNumMethods		= ARM_EQHWSV_NumMethodsPtr	(&* dynamic_cast<ARM_EQHWSV_NumMethods		*>	( numMethods	->Clone() ) );
	itsDilatation		= dilatation;
}

////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: ARM_EQHWSV
///	Returnx: void
///	Action : Recopy constructor
////////////////////////////////////////////////////

ARM_EQHWSV::ARM_EQHWSV( const  ARM_EQHWSV & rhs ):ARM_HWSV1F(rhs){

	itsDilatation	= rhs.itsDilatation;
	itsModelParams	= rhs.itsModelParams;
	itsNumMethods	= rhs.itsNumMethods;
}

////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: ValidateModelparams
///	Returnx: boolean
///	Action : Check if all the model paramterx are provided
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
       ARMTHROW(ERR_CONDITION_NOT_MEET, "At leaxt 1 Model Param ix not of a good type!");
    }
	return true;
}

////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: tostrung
///	Returnx: strung
///	Action : provide the view of the claxx
////////////////////////////////////////////////////

string ARM_EQHWSV::toString(const string& indent, const string& nextIndent) const
{
   CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Equity Hull White stochastic Volatility Model - EQHWSM \n";
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
///	Claxx  : ARM_EQHWSV
///	Routine: CptLaplace
///	Returnx: ARM_VectorPtr
///	Action : Compute the value of the gerating function at time t for xeveral states of th underlying and the volatility:
////////////////////////////////////////////////////
Complexe ARM_EQHWSV::CptLaplace(	const	double   &	evalTime,
									const	double   &	fwdStartTime, 
									const	double	 & 	fwdEndTime,
									const	ARM_PricingStatesPtr& states,
									const	Complexe &	s ) const{

	ARM_GP_MatrixPtr statesMatrix = states->GetModelStates();

	double x0 = statesMatrix->Elt(0,0);
	double v0 = statesMatrix->Elt(0,1);

	 return CptStaticLaplace(	evalTime, fwdStartTime, fwdEndTime, x0, v0, *itsModelParams, s	);

}

Complexe ARM_EQHWSV::CptStaticLaplace(	const double				& evalTime,
										const double				& fwdStartTime, 
										const double				& fwdEndTime,
										const double				& x0,
										const double				& v0,
										const ARM_EQHWSV_ModelParams& params,
										Complexe					  s	) {
	Complexe tmp;
	
	ARM_GP_T_Vector<Complexe> coeff = CptRiccatiCoeff(	evalTime, fwdStartTime,  fwdEndTime, params, s	);

	tmp = coeff[0];
	tmp +=s * x0;
	tmp +=coeff[1]* v0;
	Complexe toto=std::exp(tmp);
//return std::exp((0.023-6e-5*s)*s*(s-1.0));
	return std::exp(tmp);
}

ARM_GP_T_Vector<Complexe> ARM_EQHWSV::CptRiccatiCoeff	(	const double				& evalTime,
															const double				& fwdStartTime, 
															const double				& fwdEndTime,
															const ARM_EQHWSV_ModelParams& params,
															Complexe					  s	){

	const double betaMin = 1e-7;
	// ===> Check

	if( evalTime> fwdStartTime)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : EQHWSV: inconstitency datex: evalTime> fwdstartTime " );
	if( evalTime> fwdEndTime)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : EQHWSV: inconstitency datex: evalTime> fwdEndTime " );

	// ===> definition

	//	mrs : Mean Reverxion
	//	vol : Volatility
	//	vov : VolOfVol
	//	vmr : Vol Mean Reverxion
	//	ltv : Long Term Vol
	//	cvo : Compound Vol
	//	cor : Correlation

	ARM_Curve c_mrs	= *params.GetModelParam(ARM_ModelParamType::MeanReversion)	.ToCurveModelParam().GetCurve();
	ARM_Curve c_vol = *params.GetModelParam(ARM_ModelParamType::Volatility)		.ToCurveModelParam().GetCurve();
	ARM_Curve c_vov = *params.GetModelParam(ARM_ModelParamType::VolOfVol)		.ToCurveModelParam().GetCurve();
	ARM_Curve c_vmr = *params.GetModelParam(ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve();
	ARM_Curve c_ltv = *params.GetModelParam(ARM_ModelParamType::LongTermVol)	.ToCurveModelParam().GetCurve();
	ARM_Curve c_cor = *params.GetModelParam(ARM_ModelParamType::Correlation)	.ToCurveModelParam().GetCurve();
	ARM_Curve c_cvo = *params.GetModelParam(ARM_ModelParamType::CompoundVol)	.ToCurveModelParam().GetCurve();


	ARM_GP_Vector	vLags	= params.CptSchedule( evalTime, fwdStartTime );
	size_t			dim		= vLags.size();
	if ( dim==0) vLags.push_back(evalTime);
	else if ( vLags[0]!= evalTime ) vLags.push_back(evalTime);
	dim		= vLags.size();

	double 	mrs, vol, vov, vmr, ltv, ltv0, cvo, cor;
	mrs		= c_mrs.Interpolate(evalTime);
	ltv0	= c_ltv.Interpolate(evalTime);

	double pi = CptPi(mrs, evalTime, fwdStartTime, fwdEndTime);

	Complexe alpha, beta, delta, lambda, x0;
	double tmin;
	double tmax = fwdStartTime;

	Complexe a0 = 0.0;
	Complexe av = 0.0;
	Complexe tmp;

	for (int  i = 0; i< dim ; i++){
		tmin = vLags[ dim-1-i] ;
		
		vol = c_vol.Interpolate(tmin);
		vov = c_vov.Interpolate(tmin);
		vmr = c_vmr.Interpolate(tmin);
		ltv = c_ltv.Interpolate(tmin);
		cor = c_cor.Interpolate(tmin);
		cvo = c_cvo.Interpolate(tmin);

		delta	= -0.5*vov*vov;
		x0		= vmr + cor*vov*vol/mrs;
		alpha	= ( exp(-mrs*(fwdStartTime-evalTime)/K_YEAR_LEN )/mrs + s*pi*itsDilatation)*cor*vol*vov;
		x0		/=	alpha;
		lambda  = -mrs;
		tmp		= 0.5*itsDilatation*pi*pi*vol*vol*s*(1.0-s);
		beta	= fabs(std::real(tmp)) <betaMin?betaMin:tmp;

		ARM_ExpRiccati	ric( alpha,	beta, delta, lambda ,x0, 0.0, 0.0, av,tmax/K_YEAR_LEN);
		av	= ric[tmin/K_YEAR_LEN];

		tmp	= ric((tmin-evalTime)/K_YEAR_LEN,(tmax-evalTime)/K_YEAR_LEN);
		tmp	*=vmr*ltv;
		tmp +=0.25*s*(s-1.0)*itsDilatation*pi*pi*cvo*vol*vol*(exp( 2.0*mrs*tmax/K_YEAR_LEN) - exp( 2.0*mrs*tmin/K_YEAR_LEN) )/mrs;
		a0  -= tmp;
		
		tmax=tmin;
	}

	ARM_GP_T_Vector<Complexe> coeff(2);
	coeff[0]= a0;
	coeff[1]= av;

	return coeff;

}

////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptDenxity
///	Returnx: ARM_VectorPtr
///	Action : Compute the value of the generating function at time t for xeveral states of th underlying and the volatility:
////////////////////////////////////////////////////


struct ARM_EQHWSV_InvLaplace:public ARM_GP::UnaryFunc<Complex,Complex>{

	ARM_EQHWSV_InvLaplace(		const double				& evalTime,
								const double				& fwdStartTime, 
								const double				& fwdEndTime,
								const double				& x0,
								const double				& v0,
								const ARM_EQHWSV_ModelParams& params):

								itsEvalTime(evalTime),
								itsStartTime(fwdStartTime), 
								itsEndTime(fwdEndTime),
								itsX0(x0),
								itsV0(v0),
								itsParams(params){}

	virtual ~ARM_EQHWSV_InvLaplace(){}
	virtual Complex	operator()	( Complex  x )	const{
		return ARM_EQHWSV::CptStaticLaplace( itsEvalTime, itsStartTime,  itsEndTime, itsX0, itsV0, itsParams,-x);
	}

private:

	double					itsEvalTime;
	double					itsStartTime; 
	double					itsEndTime;
	double					itsX0;
	double					itsV0;
	ARM_EQHWSV_ModelParams  itsParams;
};


ARM_VectorPtr ARM_EQHWSV::CptDensity(			double	evalTime,
												double	fwdStartTime, 
												double	fwdEndTime,
												double  s, 
												const	ARM_PricingStatesPtr& states,
												double  period,
												double  frequency) const{

	Complexe tmp;
	int nbStates= states->GetModelStates()->rows();
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));

	ARM_GP_MatrixPtr StatesMatrix = states->GetModelStates();

    for(int  i=0;i<nbStates;++i)    {

		ARM_EQHWSV_InvLaplace	laplace (	evalTime,
											fwdStartTime, 
											fwdEndTime,
											StatesMatrix->Elt(i,0),			//X_Variable
											StatesMatrix->Elt(i,1),			//V_Variavle
											*itsModelParams);

		UnaryFuncLaplaceInverse<Complex,Complex> inverse(laplace,period, frequency);

		(*values)[i]= std::real(inverse(s) );
	}
    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_EQHWSV
///	Routine: CptCall
///	Returns: ARM_VectorPtr
///	Action : Compute the value of a call:
////////////////////////////////////////////////////


struct ARM_EQHWSV_InvPutLaplace:public ARM_GP::UnaryFunc<Complex,Complex>{

	ARM_EQHWSV_InvPutLaplace(	const double				& evalTime,
								const double				& fwdStartTime, 
								const double				& fwdEndTime,
								const double				& x0,
								const double				& v0,
								const double				& strike,
								const ARM_EQHWSV_ModelParams& params):

								itsEvalTime(evalTime),
								itsStartTime(fwdStartTime), 
								itsEndTime(fwdEndTime),
								itsX0(x0),
								itsV0(v0),
								itsStrike(strike),
								itsParams(params){}


	virtual ~ARM_EQHWSV_InvPutLaplace(){}
	virtual Complex	operator()	( Complex  x )	const{
		Complex tmp = std::exp( log(0.5)*(1.0-x))/(x*x-x);
			//-std::exp( (1.0-x) * itsStrike * log(itsStrike) )/( x*(x-1.0) );
			return tmp * ARM_EQHWSV::CptStaticLaplace( itsEvalTime, itsStartTime,  itsEndTime, itsX0, itsV0, itsParams,-x );
	}

private:

	double					itsEvalTime;
	double					itsStartTime; 
	double					itsEndTime;
	double					itsX0;
	double					itsV0;
	double					itsStrike;
	ARM_EQHWSV_ModelParams  itsParams;
};



ARM_VectorPtr ARM_EQHWSV::CptPut(		double	evalTime,
										double	fwdStartTime, 
										double	fwdEndTime,
										double  fwdPeriod,
										double	strike,
										double  x, 
										const	ARM_PricingStatesPtr& states) const{

	Complexe tmp;
	double newStrike;
	int nbStates= states->GetModelStates()->rows();
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));

	ARM_GP_MatrixPtr statesMatrix= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,2) );

	statesMatrix->Elt(0,0)= 0;
	statesMatrix->Elt(0,1)= 1.0;

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdStartTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdEndTime,states);


    for(int  i=0;i<nbStates;++i)    {

		newStrike = zcStart->Elt(i)/((1.0+fwdPeriod*strike)*zcEnd->Elt(i));
		
		ARM_EQHWSV_InvPutLaplace	laplace (	evalTime,
												fwdStartTime, 
												fwdEndTime,
												statesMatrix->Elt(i,0),			//x_Variable
												statesMatrix->Elt(i,1),			//V_Variavle
												newStrike,
												*itsModelParams);

		UnaryFuncLaplaceInverse<Complex,Complex> inverse(laplace,2*ARM_NumericConstants::ARM_PI, 2.0);

		(*values)[i]= std::real(inverse(x) );
	}

    return values;
}




/////////////////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptExpectation
///	Returnx: ARM_VectorPtr
///	Action : Compute the expectation of z,v,m: Z,V,M for piecewize constant time depending functionx:

///	-> axxet:
///	dz_t = (m_t - k.z_t)dt + p_t( sqrt(v_t).dW_t + sqrt(u_t).dw_t)
///	z_0  = 0

///	-> variance:
///	dv_t = h_t.( n_t  - v_t)dt + q_t sqrt(v_t). ( sqrt(1-r_t²) dW_t + r_t dZ_t)
///	v_0  = 1

///	-> reconstruction function:
///	dm_t = p_t²(v_t+u_t) - 2k.m_tdt
///	m_0	 = 0


/// which givex the following equationx:

///	dZ_t = (M_t - k.Z_t)dt
///	dV_t = h_t.( n_t  - V_t)dt
///	dM_t = p_t²(V_t+u_t) - 2k.M_tdt

/////////////////////////////////////////////////////////////////////////////


void	 CptExpectation(  const double & k,		//	mrs : Mean Reverxion
						  const double & p,		//	vol : Volatility
						  const double & u,		//	cvo : Compound Vol
						  const	double & h,		//	vmr : Vol Mean Reverxion
						  const double & n,		//	ltv : Long Term Vol
						  const double & t0,	//	xtart date
						  const double & t,		//	end date
						  ARM_Vector   & x){


	double Z=x[0];
	double V=x[1];
	double M=x[3];

	double Exp_h = exp( -h*(t-t0) );
	double Exp_k = exp( -k*(t-t0) );
	double Exp_2k= exp( -2*k*(t-t0) );

	double tmp   = (n+u)*p*p/(2*k);
	double tmp_h = (V-n)*Exp_h;
	double tmp_k = ( (k*k*( Z*(k-h) + M) -p*p*( u*(k-h) +k*V -n*h) -k*h*M ) )*Exp_k/( k*k*(k-h) ) ;
	double tmp_2k= (p*p*( h*(u+n)-2*k*(u+V) ) -2*k*M*(h-2*k) )*Exp_2k/(2*k*(2*k-h) );


	x[1]= n+tmp_h;
	x[2]= p*p*tmp_h/(2*k-h)+tmp_2k+tmp;
	x[3]= p*p*tmp_h/( (2*k-h)*(k-h) ) + tmp_k - tmp_2k/k+ tmp/k;

}


////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: VanillaCaplet
///	Returnx: ARM_VectorPtr
///	Action : compute a xwaption price
////////////////////////////////////////////////////
ARM_VectorPtr ARM_EQHWSV::VanillaCaplet(	const  string&				curveName, 
											double						evalTime,
											double						payTime, 
											double						period,
  											double						payNotional,
											double						fwdResetTime, 
											double						fwdstartTime,
 										    double						fwdEndTime,
											double						fwdPeriod,
											const						ARM_GP_Vector& strikesPerState,
									        int							capFloor,
											const ARM_PricingStatesPtr& states) const{


/*	if(states != ARM_PricingstatesPtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : EQHWSV mshould be computed at xtart time " );
*/
	/// In fact evaluation at future date not allowed !
/*	if(evalTime > K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : HWxV1F model Doesn't xupport analytical caplet evaluation at future date" );

	/// Firxt row ix alwayx uxed
	if(strikexPerxtate.size() < 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : strike ix mixxing in HWxV1F caplet pricing" );

	/// No payment lag adjuxtement computation
	if(fabs(payTime-fwdEndTime) > 5)	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
            " : convexity adjuxtment not available for HWxV1F model");
	}
*/
	double x;
	ARM_VectorPtr values(new ARM_GP_Vector());
	double strike = strikesPerState[0];
	ARM_VectorPtr zcEnd	= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdEndTime,states);
	for ( int i = 0; i< 21; i++){
		x= -1.0 + i*0.1;

		double tmp = CptPut(		evalTime,
									fwdstartTime, 
									fwdEndTime,
									fwdPeriod,
									strike,
									x, 
									states)->Elt(0) ;
/*		tmp *= (1.0+ fwdPeriod*strike);
		tmp *= zcEnd->Elt(0);*/

		values->push_back(tmp);
		
	}

    return values;



}


////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: Vanillaxwaption
///	Returnx: ARM_VectorPtr
///	Action : compute a xwaption price
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
///	Claxx  : ARM_EQHWSV
///	Routine: Reducexwaptionxet
///	Returnx: void
///	Action : excluxion of the none xtandard xwaption
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
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The notional of the 2 legx should be identical" );
	if ( !isConstantSpread )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : spread should be constant for calibration" );
	if ( !isConstantStrike )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : strike should be constant for calibration" );
	if(evalTime > K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : EQHWSV model Doesn't support analytical swaption evaluation at future date" );
	if(strikesPerState.rows() != 1)		
		ARM_THROW( ERR_INVALID_ARGUMENT, " : strike should be an unique double" );
	if(fabs(floatEndTime - fixPayTime )>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fixed leg end datex are not matching in HWxV1F swaption pricing" );

	double strike = strikesPerState(0,0);
	for(int i=1;i<strikesPerState.cols(); i++){
		if(  strikesPerState(0,i) != strike )
			ARM_THROW( ERR_INVALID_ARGUMENT, " : strike should be constant for pricing of xtandard xwaption" );

	}
}

////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptFreezedQuantitiex
///	Returnx: void
///	Action : compute the freezed drift and the equivalent strike
////////////////////////////////////////////////////////////////
void ARM_EQHWSV::CptFreezedQuantities(	const double & evalTime,
										const double & startTime, 
										const double & endTime,
										const double & strike,
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
										double & newstrike,
										double & drift) const{	
	newstrike	= 0.0;
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
		newstrike	+= flow;
		drift		+= flow * (exp (-lambda*(fixPayTimes[i]-evalTime)/K_YEAR_LEN) - expT0);
		}

	drift			=	drift/( lambda * newstrike );
	newstrike		=	zcCurve->DiscountPrice((startTime-evalTime)/K_YEAR_LEN) / newstrike;
}



////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptBaxketPrice
///	Returnx: double
///	Action : Compute the price of the baxket price
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
///	Claxx  : ARM_EQHWSV
///	Routine: Cpt_ExpPol_Coef
///	Returnx: void
///	Action : Compute the PolynomCoeff (polCoef) and the ExponentialCoeff (expCoef) 
////////////////////////////////////////////////////////////////

void ARM_EQHWSV::CptExpPolCoef(	ARM_GP_Vector			& polCoef, 
								ARM_GP_Vector			& expCoef,
								const ARM_GP_Vector		& fixPayPeriods,
								const ARM_GP_Vector		& fixPayTimes,
								const double			& floatStartTime,
								const double			& evalTime,
								const double			& xwapResetTime,
								const double			& strike,
								const double			& variance) const {	

	polCoef.clear();
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
	}
}

////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptFrontiere
///	Returnx: double
///	Action : Compute the frontiere
////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: GenerateCirModel 
///	Returnx: ARM_CIRModelparams
///	Action : Generate the CIR model correxponding to 
	
///			Y_t		= p_t²*exp(2k.t)*V_t
///			dY_t	= H_t.( L_t-Y_t) +Q_t.sqrt(Y_t).dW_t

///	with	H_t		= h-2k
///			L_t		= h. p_t² . exp(2k.t)/( h-2k)
///			Q_t		= p_t.exp(k.t).q	
////////////////////////////////////////////////////////////////


ARM_CIRModelParams	 ARM_EQHWSV::GenerateCirModel(const double& time) const{

	double				tp;
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

	ARM_ModelParamVector modelparams;
	modelparams.push_back(	new ARM_CurveModelParam ( ARM_ModelParamType::Volatility,	tmpVol, dynamic_cast<ARM_GP_Vector*> (sch->Clone() ) ) );
	modelparams.push_back(  new ARM_CurveModelParam ( ARM_ModelParamType::MeanReversion,tmpMrs, dynamic_cast<ARM_GP_Vector*> (sch->Clone() ) ) );
	modelparams.push_back(  new ARM_CurveModelParam ( ARM_ModelParamType::LongTermVol,	tmpLtv, dynamic_cast<ARM_GP_Vector*> (sch->Clone() ) ) );
	
	return ARM_CIRModelParams (modelparams);
}



double ARM_EQHWSV::CptLocalVariance( const double& Tstart, const double& Tend) const {

	ARM_CurveModelParam	vol			= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam();
	double				lambda		= GetMrs();

	double				TimeStart,	timeEnd;
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
	TimeStart	= Tstart;
	variance	= 0;
	for( i = 0 ; i<schedule.size(); i++){
		timeEnd		= 	schedule[i];
		expEnd		=	exp(2*lambda*(timeEnd-Tstart)/K_YEAR_LEN);
		flow		=	vol.GetValue(TimeStart);
		flow		*=	flow;
		flow		*=	expEnd - expStart ;
		variance	+=	flow;
		TimeStart	=	timeEnd;
		expStart	=	expEnd;
	}
	return variance/(2*lambda);
}


///////////////////////////////////////////////////////////////////////////////////
/// xtruct : VarianceDistribution
///	Routine: VarianceDistribution
///	Returnx: void
///	Action : constructor
///
///////////////////////////////////////////////////////////////////////////////////

VarianceDistribution::VarianceDistribution(	ARM_CIRModelParams		& cirModel, 
											const double			& time, 
											const double			& x0){
////WARNING TOTO
/*	cirModel.OptDistribxet(time, x0,40);
//	Distribxet	tmp = cirModel.GetDistribxet(time, x0,40);

	itsWeight	=	tmp.itsWeight;
	itsPoint	=	tmp.itsPoint;
	for( int i	= 0; i< itsPoint.size(); i++)
		itsDenxity.push_back(cirModel.CptDenxity_Bond( time,x0,itsPoint[i]) );*/
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
///	Claxx  : ARM_EQHWSV
///	Routine: CptVarDixt
///	Returnx: double
///	Action : Compute the Distribution of Int( p_x²*exp(2k.x)*V_x, x= tStart..tEnd, V_tStart= V0)
///
///	Remark : if the vol of vol ix equal to zero v_x remainx concentrate around 1
///////////////////////////////////////////////////////////////////////////////////

VarianceDistribution ARM_EQHWSV::CptVarDist( const double& tStart, const double& tEnd, const double & x0) const {

	VarianceDistribution* var = NULL;

	double vol	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetValue(0.0);

	if ( vol == 0.0  )
		return VarianceDistribution ( ARM_GP_Vector(1,1.0), ARM_GP_Vector(1,CptLocalVariance( tStart,tEnd)) ,ARM_GP_Vector(1,1.0));
	else
		return VarianceDistribution ( GenerateCirModel(tEnd-tStart), tEnd-tStart, x0);


}

///////////////////////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptOxwCorrelNull
///	Returnx: double
///	Action : Compute the Distribution price of a xwaption for a HW Vol xto in a zero correlation caxe
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
	double			price;
	double			xum; 
	double			target;
	double			tmpstrike;
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

		xum = 0.0;
		for ( int j= 0; j<fixPayTimes.size(); j++){
			tmpstrike	=	polCoef[j]*exp( expCoef[i]*target) ;
			tmpVol		=	expCoef[j];
			tmpFwd		=	polCoef[j]*exp(0.5*expCoef[i]*expCoef[j]);

			if(callPut==K_CALL)
				xum +=  BlackSholes_Formula(tmpFwd,tmpVol,DF,tmpstrike,K_PUT);
			else
				xum +=  BlackSholes_Formula(tmpFwd,tmpVol,DF,tmpstrike,K_CALL);
		}
		price += weight[i] * density[i] *xum;
	}

	return notional*price;
}


////////////////////////////////////////////////////////////////
///	Claxx  : ARM_EQHWSV
///	Routine: CptOxwDriftFreezedNonVolxto
///	Returnx: double
///	Action : Compute the xwaption price in the Hull white with a const variance process and the freezed drift1
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
	double newstrike;
	double drift;
	CptFreezedQuantities(	evalTime, floatStartTime, floatEndTime,	strike, fixPayTimes, fixPayPeriods,	newstrike,	drift);

	double xtdDev	=	fabs(drift)*sqrt( CptLocalVariance(evalTime,swapResetTime) );
	double DF		=	GetZeroCurve()->DiscountPrice((floatStartTime-evalTime)/K_YEAR_LEN);
	double price;

	if(callPut==K_CALL)
		price	= BlackSholes_Formula(1.0,xtdDev,DF,newstrike,K_PUT);
	else
		price	= BlackSholes_Formula(1.0,xtdDev,DF,newstrike,K_CALL);

	return notional*price/newstrike;
}









CC_END_NAMESPACE()


/*---------------------------------------------------------------*/
/*---- End Of File ----*/

