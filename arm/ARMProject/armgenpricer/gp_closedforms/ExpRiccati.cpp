/*!
 *
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date June 2007
 */
 

#include "gpclosedforms/ExpRiccati.h"
#include "gpclosedforms/whittaker.h"


using namespace std;


CC_BEGIN_NAMESPACE(ARM)

////////////////////////////////////////////////////
///	Class  : ARM_ExpRiccati
///	Routine: ARM_ExpRiccati
///	Returns: double
///	Action : constructor
////////////////////////////////////////////////////
ARM_ExpRiccati::ARM_ExpRiccati (		const Complexe&			alpha,
										const Complexe&			beta0,
										const Complexe&			beta1,
										const Complexe&			gamma0,
										const Complexe&			gamma1,
										const Complexe&			gamma2,
										const Complexe&			lambda,
										const Complexe&			x0,
										const Complexe&			t0):	ARM_ODEFunc ( ),	
																		itsLambda	(	lambda	),
																		itsAlpha	(	alpha	),
																		itsBeta0	(	beta0	),
																		itsBeta1	(	beta1	){





	itsMu = (beta1*lambda-beta0*beta1+2.0*alpha*gamma1)/2.0;
	itsNu = sqrt(beta0*beta0-4.0*alpha*gamma0)/2.0;
	itsEta= sqrt(beta1*beta1-4.0*alpha*gamma2);

	Complexe tmp   = (itsEta+itsBeta1)*exp(itsLambda*t0);
	tmp -= itsLambda-itsBeta0;
	tmp /= 2.0*itsAlpha;
	tmp -= itsMu/(itsAlpha*itsEta);
	tmp += x0;

	Complexe Exp0= itsEta*exp(itsLambda*t0)/itsLambda;
	Complexe WM0 = Hypergeometric_Whittaker_M(itsMu/(itsEta*itsLambda),itsNu/itsLambda,Exp0);
	Complexe LM0 = Hypergeometric_Whittaker_M(itsMu/(itsEta*itsLambda)+1.0,itsNu/itsLambda,Exp0);
	Complexe WW0 = Hypergeometric_Whittaker_W(itsMu/(itsEta*itsLambda),itsNu/itsLambda,Exp0);
	Complexe LW0 = Hypergeometric_Whittaker_W(itsMu/(itsEta*itsLambda)+1.0,itsNu/itsLambda,Exp0);


	itsConstant  =  itsAlpha*tmp*WM0 +(itsMu/itsEta+itsNu+itsLambda/2.0)*LM0;
	itsConstant  /=-itsAlpha*tmp*WW0 + itsLambda*LW0;
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpRiccati
///	Routine: operator ()
///	Returns: double
///	Action : return the primitive of the solution between t and T
////////////////////////////////////////////////////

Complexe ARM_ExpRiccati::operator () (double t, double T) const{

	Complexe Expt= itsEta*exp(itsLambda*t)/itsLambda;
	Complexe ExpT= itsEta*exp(itsLambda*T)/itsLambda;

	Complexe WMt = Hypergeometric_Whittaker_M(itsMu/(itsEta*itsLambda),itsNu/itsLambda,Expt);
	Complexe WWt = Hypergeometric_Whittaker_W(itsMu/(itsEta*itsLambda),itsNu/itsLambda,Expt);
	Complexe WMT = Hypergeometric_Whittaker_M(itsMu/(itsEta*itsLambda),itsNu/itsLambda,ExpT);
	Complexe WWT = Hypergeometric_Whittaker_W(itsMu/(itsEta*itsLambda),itsNu/itsLambda,ExpT);

	
	Complexe tmpt = itsBeta1*exp(itsLambda*t)/itsLambda;
	tmpt += ( itsBeta0 - itsLambda)*t;
	tmpt /= 2.0;
	tmpt += log( itsConstant*WWt+WMt);

	Complexe tmpT = itsBeta1*exp(itsLambda*T)/itsLambda;
	tmpT += ( itsBeta0 - itsLambda)*T;
	tmpT /= 2.0;
	tmpT += log( itsConstant*WWT+WMT);


	return	(tmpt-tmpT)/itsAlpha;
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpRiccati
///	Routine: operator []
///	Returns: double
///	Action : return the solution at time t
////////////////////////////////////////////////////
Complexe ARM_ExpRiccati::operator [] (double t) const{

	Complexe Expt= itsEta*exp(itsLambda*t)/itsLambda;

	Complexe WMt = Hypergeometric_Whittaker_M(itsMu/(itsEta*itsLambda),itsNu/itsLambda,Expt);
	Complexe LMt = Hypergeometric_Whittaker_M(itsMu/(itsEta*itsLambda)+1.0,itsNu/itsLambda,Expt);
	Complexe WWt = Hypergeometric_Whittaker_W(itsMu/(itsEta*itsLambda),itsNu/itsLambda,Expt);
	Complexe LWt = Hypergeometric_Whittaker_W(itsMu/(itsEta*itsLambda)+1.0,itsNu/itsLambda,Expt);
	

	Complexe tmp = -(itsEta+itsBeta1)*exp(itsLambda*t);	
	tmp += (itsLambda-itsBeta0);
	tmp /= 2.0*itsAlpha;
	tmp	+= itsMu/(itsAlpha*itsEta);

	Complexe res = itsLambda*itsConstant*LWt;
	res -= (itsMu/itsEta+itsNu+itsLambda/2.0)*LMt;
	res /= itsConstant*WWt+WMt;
	res /= itsAlpha;

	return tmp + res;

}

void ARM_ExpRiccati::derivs(double x, std::vector<double>* yt, std::vector<double>* dyt) const	{
		ARMTHROW(ERR_INVALID_ARGUMENT, "method ARM_ExpRiccati::derivs is not implemented" );
	}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/