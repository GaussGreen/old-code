/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file normal_heston.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include <cmath>
#include <complex>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal_heston.h"
#include "gpclosedforms/oscillatory_integrals.h"
#include "gpclosedforms/vanille_bs_formula.h"

#include "gpbase/numericconstant.h"

using std::sqrt;
using std::log;
using std::real;
using std::imag;
using std::exp;
using std::norm;
using std::complex;

//using std::complex<double>;
//using std::complex<long double>;


///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ (exp(Jt)-1) dq)
///			dV=(omega-theta*V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
///
/////////////////////////////////////////////////////////////////////////////:


CC_BEGIN_NAMESPACE(ARM)


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

complex<long double> GaussianHestonFondamentalTransform(double rho,double lambdaV,double thetaV,
														double kappaV,double V0, complex<long double> omega, double T)
{
	complex<long double> crho(rho,0.0),clambdaV(lambdaV,0.0),cthetaV(thetaV,0.0),
		ckappaV(kappaV,0.0),cV0(V0,0.0),cI(0.0,1.0),cT(T,0.0),deux(2.0,0),un(1.0,0.0);
	
	complex<long double> zeta,B,A,psip,psim,aux,kappaomega,aux2,ckappa2,ckappa4,expzeta,exp0;
	ckappa2=ckappaV*ckappaV;
	ckappa4=ckappa2*ckappa2;
	kappaomega=ckappaV*omega;
	aux=clambdaV+crho*kappaomega*cI;

	zeta=sqrt(aux*aux+kappaomega*kappaomega);
	if(norm(zeta)<1e-50)
	{
		aux2=clambdaV+cI*crho*omega*ckappaV;
		B=aux2/ckappa2-aux2*ckappa2/(cT*aux2+ckappa4);
		A=cthetaV*clambdaV*(deux*(cT*aux2-cI*ckappa4*arctan(kappaomega*crho*cT/(ckappa4+clambdaV*cT)))+
			ckappa4*log(ckappa4*ckappa4/(ckappa4*ckappa4+deux*ckappa4*clambdaV*cT+clambdaV*clambdaV*cT*cT+crho*crho*cT*cT*kappaomega*kappaomega)));
		return exp(A+B*cV0);
	}
	else
	{
		expzeta=exp(-zeta*cT);
		psip=-aux+zeta;
		psim=aux+zeta;
		A=-cthetaV*clambdaV*(cT*psip+deux*log((psim+expzeta*psip)/(deux*zeta)))/ckappa2;
		B=-omega*omega*(un-expzeta)/(psim+psip*expzeta);
		exp0=A+B*cV0;
		if(real(exp0)<-100.0)
		{
			return 0.0;
		}
		else
		{
			return exp(exp0);
		}
	}
}

// lambdaB ici shifte l'evaluation complexe et donne un resultat equivalent dans la mesure ou il n'ya pas de pole entre 0 et lambdaB
complex<long double> FourierPayoff(complex<long double> omega, complex<long double> k)
{
	complex<long double> I(0.0,1.0);
	return -exp(I*k*omega)/(omega*omega);
}


double NormalHestonBasic(double rho,double lambdaV,double thetaV,
						double kappaV,double V0,  double S0,double k,double T,double lambdaB,
						int callput,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec)
{
	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double rho0;
		double lambdaV0;
		double thetaV0;
		double kappaV0;
		double Vinitial0;
		double Sinitial0;
		double k0;
		double T0;
		double lambdaB0;
		int callput0;
	

	public:
		integralToBeDone(double rho0_,double lambdaV0_,double thetaV0_, double kappaV0_,double Vinitial0_,double Sinitial0_,double k0_,
			double T0_,double lambdaB0_, int callput_, int nbfirst_,int nb_,int NbStage_,int NbOscill_,double speed_ ,double prec_) :
		rho0(rho0_),lambdaV0(lambdaV0_),thetaV0(thetaV0_),kappaV0(kappaV0_),Vinitial0(Vinitial0_),Sinitial0(Sinitial0_),k0(k0_),T0(T0_),
			lambdaB0(lambdaB0_),callput0(callput_),OscillatoryIntegral(nbfirst_,nb_,NbStage_,NbOscill_,speed_,prec_)	/// the controled version
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double omega)
		{
			complex<long double> res1;
			complex<long double> res2;
			if(callput0==K_CALL)
			{
				complex<long double> omegashifted(omega,lambdaB0);
				res1=GaussianHestonFondamentalTransform(rho0,lambdaV0,thetaV0,kappaV0,Vinitial0,omegashifted,T0);
				res2=FourierPayoff(omegashifted,k0-Sinitial0);
			}
			else
			{
				complex<long double> omegashifted(omega,-lambdaB0);
				res1=GaussianHestonFondamentalTransform(rho0,lambdaV0,thetaV0,kappaV0,Vinitial0,omegashifted,T0);
				res2=FourierPayoff(omegashifted,k0-Sinitial0);
				
			}
			return real(res1*res2);
		}
		double mapper(double k)
		{
			return tan(k);
		}
		double inversemapper(double S)
		{
			return atan(S);
		}
		double mapperjacobian(double k)   // should be equal to D[mapper[x],x]
		{
			return 1./(k*k+1);
		}
	};
			
	int legendre_nb_pts=nb;
/// here speed gives the number of period per unit for the integral
	double speed1=fabs(k);
	double speed2=sqrt(V0);
	double speed=(speed1>speed2)? speed1 : speed2 ;
	integralToBeDone IS(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,callput,nbfirst,nb,NbStage,NbOscill,speed,prec);

	double IS_value=IS.value();
	return IS_value/ARM_NumericConstants::ARM_PI;
}

double NormalHestonBasicTry(double rho,double lambdaV,double thetaV,
						double kappaV,double V0,  double S0,double k,double T,double lambdaB,
						int callput, int nb)
{
	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double rho0;
		double lambdaV0;
		double thetaV0;
		double kappaV0;
		double Vinitial0;
		double Sinitial0;
		double k0;
		double T0;
		double lambdaB0;
		int callput0;
	

	public:
		integralToBeDone(double rho0_,double lambdaV0_,double thetaV0_, double kappaV0_,double Vinitial0_,double Sinitial0_,double k0_,
			double T0_,double lambdaB0_, int callput_, int nb_) :
		rho0(rho0_),lambdaV0(lambdaV0_),thetaV0(thetaV0_),kappaV0(kappaV0_),Vinitial0(Vinitial0_),Sinitial0(Sinitial0_),k0(k0_),T0(T0_),
			lambdaB0(lambdaB0_),callput0(callput_),OscillatoryIntegral(0.,nb_,5,1.)	/// the controled version
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double omega)
		{
			complex<long double> res1;
			complex<long double> res2;
			if(callput0==K_CALL)
			{
				complex<long double> omegashifted(omega,lambdaB0);
				res1=GaussianHestonFondamentalTransform(rho0,lambdaV0,thetaV0,kappaV0,Vinitial0,omegashifted,T0);
				res2=FourierPayoff(omegashifted,k0-Sinitial0);
			}
			else
			{
				complex<long double> omegashifted(omega,-lambdaB0);
				res1=GaussianHestonFondamentalTransform(rho0,lambdaV0,thetaV0,kappaV0,Vinitial0,omegashifted,T0);
				res2=FourierPayoff(omegashifted,k0-Sinitial0);
				
			}
			return real(res1*res2);
		}
		double mapper(double k)
		{
			return tan(k);
		}
		double inversemapper(double S)
		{
			return atan(S);
		}
		double mapperjacobian(double k)   // should be equal to D[mapper[x],x]
		{
			return 1./(k*k+1);
		}
	};
			
	int legendre_nb_pts=nb;
/// here speed gives the number of period per unit for the integral
	double speed1=fabs(k);
	double speed2=sqrt(V0);
	double speed=(speed1>speed2)? speed1 : speed2 ;
	integralToBeDone IS(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,callput,legendre_nb_pts);

	double IS_value=IS.value();
	return IS_value/ARM_NumericConstants::ARM_PI;
}


double NormalHeston(double rho,double lambdaV,double thetaV1,
						double kappaV1,double V01,  double S01,double k1,double T,
						int callput,double lambdaB,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec)
{
	if(lambdaB < 0.)
	{
		return NormalHestonTry(rho, lambdaV, thetaV1, kappaV1, V01, S01, k1, T, callput, nb);
	}

	double L,kappaV,thetaV,V0,S0,k;
	L=1.0/sqrt(V01);
	thetaV=L*L*thetaV1;
	kappaV=L*kappaV1;
	V0=1.0;
	S0=L*S01;
	k=L*k1;
	if(callput==K_CALL)
	{
		if(k>=S0)
		{
			return  NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_CALL,nbfirst,nb,NbStage,NbOscill,prec)/L;
		}
		else
		{
			return  (S0-k+NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_PUT,nbfirst,nb,NbStage,NbOscill,prec))/L;
			
		}
	}
	else
	{
		if(k<=S0)
		{
			return  NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_PUT,nbfirst,nb,NbStage,NbOscill,prec)/L;
		}
		else
		{
			return  (k-S0+NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_CALL,nbfirst,nb,NbStage,NbOscill,prec))/L;
			
		}
	}

}

double NormalHestonTry(double rho,double lambdaV,double thetaV1,
						double kappaV1,double V01,  double S01,double k1,double T,
						int callput,int nb, double lambdaB)
{
	double L,kappaV,thetaV,V0,S0,k;
	L=1.0/sqrt(V01);
	thetaV=L*L*thetaV1;
	kappaV=L*kappaV1;
	V0=1.0;
	S0=L*S01;
	k=L*k1;
	if(callput==K_CALL)
	{
		if(k>=S0)
		{
			return  NormalHestonBasicTry(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_CALL,nb)/L;
		}
		else
		{
			return  (S0-k+NormalHestonBasicTry(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_PUT,nb))/L;
			
		}
	}
	else
	{
		if(k<=S0)
		{
			return  NormalHestonBasicTry(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_PUT,nb)/L;
		}
		else
		{
			return  (k-S0+NormalHestonBasicTry(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,
				K_CALL,nb))/L;
			
		}
	}

}

/// Version complexe, necessité par le calcul de bisar digital qui paye float

complex<long double> complexified_GaussianHestonFondamentalTransform(
				complex<long double> crho,complex<long double> clambdaV,complex<long double> cthetaV,
				complex<long double> ckappaV,complex<long double> cV0, complex<long double> omega, complex<long double> cT)
{
	complex<long double> cI(0.0,1.0),deux(2.0,0),un(1.0,0.0);

	complex<long double> zeta,B,A,psip,psim,aux,kappaomega,aux2,ckappa2,ckappa4,expzeta,exp0;
	ckappa2=ckappaV*ckappaV;
	ckappa4=ckappa2*ckappa2;
	kappaomega=ckappaV*omega;
	aux=clambdaV+crho*kappaomega*cI;

	zeta=sqrt(aux*aux+kappaomega*kappaomega);
	if(norm(zeta)<1e-50)
	{
		aux2=clambdaV+cI*crho*omega*ckappaV;
		B=aux2/ckappa2-aux2*ckappa2/(cT*aux2+ckappa4);
		A=cthetaV*clambdaV*(deux*(cT*aux2-cI*ckappa4*arctan(kappaomega*crho*cT/(ckappa4+clambdaV*cT)))+
			ckappa4*log(ckappa4*ckappa4/(ckappa4*ckappa4+deux*ckappa4*clambdaV*cT+clambdaV*clambdaV*cT*cT+crho*crho*cT*cT*kappaomega*kappaomega)));
		return exp(A+B*cV0);
	}
	else
	{
		expzeta=exp(-zeta*cT);
		psip=-aux+zeta;
		psim=aux+zeta;
		A=-cthetaV*clambdaV*(cT*psip+deux*log((psim+expzeta*psip)/(deux*zeta)))/ckappa2;
		B=-omega*omega*(un-expzeta)/(psim+psip*expzeta);
		exp0=A+B*cV0;
		if(real(exp0)<-100.0)
		{
			return 0.0;
		}
		else
		{
			return exp(exp0);
		}
	}
}




double Complexified_NormalHestonBasic(complex<double> rho,complex<double> lambdaV,complex<double> thetaV,
						complex<double> kappaV,complex<double> V0,  complex<double> S0,complex<double> k,complex<double> T,
						int callput,double lambdaB,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec)
{

	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		complex<double> rho0;
		complex<double> lambdaV0;
		complex<double> thetaV0;
		complex<double> kappaV0;
		complex<double> Vinitial0;
		complex<double> Sinitial0;
		complex<double> k0;
		complex<double> T0;
		double lambdaB0;
		int callput0;
	

	public:
		integralToBeDone(complex<double> rho0_,complex<double> lambdaV0_,complex<double> thetaV0_, complex<double> kappaV0_,
			complex<double> Vinitial0_,complex<double> Sinitial0_,complex<double> k0_,complex<double> T0_,
			double lambdaB0_, int callput_, int nbfirst_,int nb_,int NbStage_,int NbOscill_,double speed_ ,double prec_) :
		rho0(rho0_),lambdaV0(lambdaV0_),thetaV0(thetaV0_),kappaV0(kappaV0_),Vinitial0(Vinitial0_),Sinitial0(Sinitial0_),k0(k0_),T0(T0_),
			lambdaB0(lambdaB0_),callput0(callput_),OscillatoryIntegral(nbfirst_,nb_,NbStage_,NbOscill_,speed_,prec_)	/// the controled version
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double omega)
		{
			complex<long double> res1;
			complex<long double> res2;
			if(callput0==K_CALL)
			{
				complex<long double> omegashifted(omega,lambdaB0);
				res1=complexified_GaussianHestonFondamentalTransform(rho0,lambdaV0,thetaV0,kappaV0,Vinitial0,omegashifted,T0);
				res2=FourierPayoff(omegashifted,k0-Sinitial0);
			}
			else
			{
				complex<long double> omegashifted(omega,-lambdaB0);
				res1=complexified_GaussianHestonFondamentalTransform(rho0,lambdaV0,thetaV0,kappaV0,Vinitial0,omegashifted,T0);
				res2=FourierPayoff(omegashifted,k0-Sinitial0);
				
			}
			return real(res1*res2);
		}
		double mapper(double k)
		{
			return tan(k);
		}
		double inversemapper(double S)
		{
			return atan(S);
		}
		double mapperjacobian(double k)   // should be equal to D[mapper[x],x]
		{
			return 1./(k*k+1);
		}
	};
			
	int legendre_nb_pts=nb;
/// here speed gives the number of period per unit for the integral
	double speed1=sqrt(norm(k));
	double speed2=sqrt(norm(V0));
	double speed=(speed1>speed2)? speed1 : speed2 ;
	integralToBeDone IS(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,lambdaB,callput,nbfirst,nb,NbStage,NbOscill,speed,prec);

	double IS_value=IS.value();
	return IS_value/ARM_NumericConstants::ARM_PI;
}




double Complexified_NormalHeston(complex<double> rho,complex<double> lambdaV,complex<double> thetaV1,
						complex<double> kappaV1,complex<double> V01,  complex<double> S01,complex<double> k1,
						complex<double> T,
						int callput,double lambdaB,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec)
{
	complex<double> Lc,kappaV,thetaV,V0,S0,k,un(1.0,0.0);
	double L;
	L=1.0/sqrt(sqrt(norm(V01)));
	Lc=L;
	thetaV=Lc*Lc*thetaV1;
	kappaV=Lc*kappaV1;
	V0=un;
	S0=Lc*S01;
	k=Lc*k1;

	if(callput==K_CALL)
	{
		if(real(k-S0)>=0)
		{
			return  Complexified_NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,
				K_CALL,lambdaB,nbfirst,nb,NbStage,NbOscill,prec)/L;
		}
		else
		{
			return  (real(S0-k)+Complexified_NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,
				K_PUT,lambdaB,nbfirst,nb,NbStage,NbOscill,prec))/L;
			
		}
	}
	else
	{
		if(real(k-S0)>=0)
		{
			return  Complexified_NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,
				K_PUT,lambdaB,nbfirst,nb,NbStage,NbOscill,prec)/L;
		}
		else
		{
			return  (real(k-S0)+Complexified_NormalHestonBasic(rho,lambdaV,thetaV,kappaV,V0,S0,k,T,
				K_CALL,lambdaB,nbfirst,nb,NbStage,NbOscill,prec))/L;
			
		}
	}

}

complex<long double> SuperNormalHestonFondamentalTransform(double rho1,double lambda1,double theta1,double kappa1,double V01,
														   double rho2,double lambda2,double theta2,double kappa2,double V02,
														   complex<long double> omega, double T)
{
	complex<long double> crho1(rho1,0.), clambda1(lambda1,0.), ctheta1(theta1,0.), ckappa1(kappa1,0.), cV01(V01,0.);
	complex<long double> crho2(rho2,0.), clambda2(lambda2,0.), ctheta2(theta2,0.), ckappa2(kappa2,0.), cV02(V02,0.);
	complex<long double> cT(T,0.), I(0.,1.), un(1.,0.), deux(2.,0.);

	complex<long double> aux1 = clambda1 + crho1*ckappa1*I*omega;
	complex<long double> zeta1 = sqrt(aux1*aux1 + ckappa1*ckappa1*omega*omega);
	complex<long double> p1 = zeta1-aux1;
	complex<long double> m1 = zeta1+aux1;
	complex<long double> r1 = m1+p1*exp(-zeta1*cT);
	complex<long double> B1 = -omega*omega*(un-exp(-zeta1*cT))/r1;

	complex<long double> aux2 = clambda2 + crho2*ckappa2*I*omega;
	complex<long double> zeta2 = sqrt(aux2*aux2 + ckappa2*ckappa2*omega*omega);
	complex<long double> p2 = zeta2-aux2;
	complex<long double> m2 = zeta2+aux2;
	complex<long double> r2 = m2+p2*exp(-zeta2*cT);
	complex<long double> B2 = -omega*omega*(un-exp(-zeta2*cT))/r2;

	complex<long double> A1 = -ctheta1*clambda1*(cT*p1+deux*log(r1/deux/zeta1))/ckappa1/ckappa1;
	complex<long double> A2 = -ctheta2*clambda2*(cT*p2+deux*log(r2/deux/zeta2))/ckappa2/ckappa2;

	complex<long double> res = A1+A2+B1*cV01+B2*cV02;

	if(real(res)<-100.)
	{
		return 0.;
	}
	else
	{
		return exp(res);
	}
}

double SuperNormalHeston(double rho1, double Kappa1, double Theta1, double Nu1, double V01,
						 double rho2, double Kappa2, double Theta2, double Nu2, double V02,
						 double S0, double k, double T, int callput, int nb)
{
	if(callput == K_CALL && k < S0)
	{
		return S0 - k + SuperNormalHeston(rho1, Kappa1, Theta1, Nu1, V01, 
										  rho2, Kappa2, Theta2, Nu2, V02,
										  S0, k, T, K_PUT, nb);
	}
	else if(callput == K_PUT && k > S0)
	{
		return k - S0 + SuperNormalHeston(rho1, Kappa1, Theta1, Nu1, V01, 
										  rho2, Kappa2, Theta2, Nu2, V02,
										  S0, k, T, K_CALL, nb);
	}

	double L = 1. / sqrt(V01);
	// on passe aux notations OC
	double kappa1 = Nu1 * L;
	double kappa2 = Nu2 * L;
	double theta1 = Theta1 * L * L;
	double theta2 = Theta2 * L * L;
	double lambda1 = Kappa1;
	double lambda2 = Kappa2;
	double NormS0 = S0 * L;
	double NormK = k * L;
	double V1 = 1.;
	double V2 = V02 * L * L;
	double lambdaB = 0.1;

	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double rho1;
		double lambda1;
		double theta1;
		double kappa1;
		double V01;
		double rho2;
		double lambda2;
		double theta2;
		double kappa2;
		double V02;
		double S0;
		double k;
		double T;
		int callput;
		double lambdaB;

	public:
		integralToBeDone(double rho1_,double lambda1_,double theta1_, double kappa1_,double V01_,
			double rho2_,double lambda2_,double theta2_,double kappa2_,double V02_,
			double S0_,double k_,
			double T_,int callput_, int nb_, double lambdaB_) :
			rho1(rho1_), lambda1(lambda1_), theta1(theta1_), kappa1(kappa1_), V01(V01_),
			rho2(rho2_), lambda2(lambda2_), theta2(theta2_), kappa2(kappa2_), V02(V02_),
			S0(S0_), k(k_), T(T_),callput(callput_),
			lambdaB(lambdaB_), OscillatoryIntegral(0.,nb_,5,1.)	/// the controled version
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double omega)
		{
			complex<long double> res1;
			complex<long double> res2;
			if(callput==K_CALL)
			{
				complex<long double> omegashifted(omega,lambdaB);
				res1=SuperNormalHestonFondamentalTransform(rho1,lambda1,theta1,kappa1,V01,
														   rho2,lambda2,theta2,kappa2,V02,
														   omegashifted,T);
				res2=FourierPayoff(omegashifted,k-S0);
			}
			else
			{
				complex<long double> omegashifted(omega,-lambdaB);
				res1=SuperNormalHestonFondamentalTransform(rho1,lambda1,theta1,kappa1,V01,
														   rho2,lambda2,theta2,kappa2,V02,
														   omegashifted,T);
				res2=FourierPayoff(omegashifted,k-S0);
				
			}
			return real(res1*res2);
		}
		double mapper(double k)
		{
			return tan(k);
		}
		double inversemapper(double S)
		{
			return atan(S);
		}
		double mapperjacobian(double k)   // should be equal to D[mapper[x],x]
		{
			return 1./(k*k+1);
		}
	};
			
	integralToBeDone IS(rho1,lambda1,theta1,kappa1,V1,rho2,lambda2,theta2,kappa2,V2,NormS0,NormK,T,callput,nb,lambdaB);
	double IS_value=IS.value();
	return IS_value/ARM_NumericConstants::ARM_PI/L;
}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
