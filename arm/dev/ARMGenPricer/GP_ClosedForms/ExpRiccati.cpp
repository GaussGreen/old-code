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
ARM_ExpRiccati::ARM_ExpRiccati (	const Complexe&			alpha,
									const Complexe&			beta,
									const Complexe&			delta,
									const Complexe&			lambda,
									const Complexe&			x0,
									const Complexe&			x1,
									const Complexe&			x2,
									const Complexe&			y0,
									const Complexe&			t0):ARM_ODEFunc(),	
																itsLambda	(lambda),
																itsDelta	(delta),
																itsAlpha	(alpha/delta),
																itsY0		(delta*y0/lambda),
																itsT0		(t0),
																itsX		(x0*alpha/(2.0*lambda) ){

	Complexe tmp;
	Complexe tmpLambda = lambda/delta;

	tmp			= beta/delta;
	itsBeta		= sqrt(tmp);//

	tmp			= 1.0-4.0*beta/(itsAlpha*itsAlpha*delta);
	itsEta		= itsAlpha*sqrt(tmp)/tmpLambda;//
	itsSin		= sqrt(tmp); //

	itsGamma	= itsX * sqrt(1.0-4.0*(x1/x0)*(x2/x0)*itsBeta*itsBeta/(itsAlpha*itsAlpha) );
	itsKappa	= 2.0*(x1+x2)*itsBeta/(itsAlpha*x0);
	itsKappa	= (1.0+itsX*(2.0-2.0*itsBeta*itsKappa/itsAlpha))/(2.0*itsSin);


	Complexe Exp0= itsEta*exp(-itsLambda*itsT0);
	Complexe WM0 = Hypergeometric_Whittaker_M(itsKappa,itsGamma,Exp0);
	Complexe LM0 = Hypergeometric_Whittaker_M(itsKappa+1.0,itsGamma,Exp0);
	Complexe WW0 = Hypergeometric_Whittaker_W(itsKappa,itsGamma,Exp0);
	Complexe LW0 = Hypergeometric_Whittaker_W(itsKappa+1.0,itsGamma,Exp0);

	Complexe coef= 0.5*(1.0+1.0/itsSin)*Exp0;
	Complexe eps = itsKappa+itsX+0.5;
	Complexe eta = itsGamma+itsKappa+0.5;

	itsOmega = -LW0-WW0*( eps + itsY0 - coef );
	itsTheta = (itsGamma+itsKappa+0.5)*LM0-WM0*( eps + itsY0 - coef );

}


////////////////////////////////////////////////////
///	Class  : ARM_ExpRiccati
///	Routine: operator ()
///	Returns: double
///	Action : return the primitive of the solution between t and T
////////////////////////////////////////////////////

Complexe ARM_ExpRiccati::operator () (double t, double T) const{

	Complexe Expt= itsEta*exp(-itsLambda*t);
	Complexe ExpT= itsEta*exp(-itsLambda*T);

	Complexe WMt = Hypergeometric_Whittaker_M(itsKappa,itsGamma,Expt);
	Complexe WWt = Hypergeometric_Whittaker_W(itsKappa,itsGamma,Expt);
	Complexe WMT = Hypergeometric_Whittaker_M(itsKappa,itsGamma,ExpT);
	Complexe WWT = Hypergeometric_Whittaker_W(itsKappa,itsGamma,ExpT);

	Complexe tmp =(itsX+0.5)*itsLambda*(t-T);
	tmp			+= 1.0/(2.0*itsSin)*(Expt-ExpT);
	tmp			+= log(itsOmega*WMt-itsTheta*WWt);
	tmp			-= log(itsOmega*WMT-itsTheta*WWT);

	return		tmp/itsDelta;
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpRiccati
///	Routine: operator []
///	Returns: double
///	Action : return the solution at time t
////////////////////////////////////////////////////
Complexe ARM_ExpRiccati::operator [] (double t) const{

	Complexe Expt= itsEta*exp(-itsLambda*t);

	Complexe WMt = Hypergeometric_Whittaker_M(itsKappa,itsGamma,Expt);
	Complexe LMt = Hypergeometric_Whittaker_M(itsKappa+1.0,itsGamma,Expt);
	Complexe WWt = Hypergeometric_Whittaker_W(itsKappa,itsGamma,Expt);
	Complexe LWt = Hypergeometric_Whittaker_W(itsKappa+1.0,itsGamma,Expt);


	Complexe tmp =-(itsKappa+itsX+0.5);
	tmp			+= 0.5*(1.0+1.0/itsSin)*Expt;
	tmp			+= ( (itsGamma+itsKappa+0.5)*itsOmega*LMt+itsTheta*LWt)/(itsOmega*WMt-itsTheta*WWt);

	return itsLambda*tmp/itsDelta;

}

void ARM_ExpRiccati::derivs(double x, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const	{
		ARMTHROW(ERR_INVALID_ARGUMENT, "method ARM_ExpRiccati::derivs is not implemented" );
	}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/