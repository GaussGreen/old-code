/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.cpp
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
#include "gpclosedforms/heston.h"
#include "gpclosedforms/oscillatory_integrals.h"
#include "gpclosedforms/vanille_bs_formula.h"

#include "gpbase/numericconstant.h"

using std::sqrt;
using std::log;
using std::real;
using std::imag;
using std::exp;
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

#define ARM_CF_EPS 1.0e-14

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////



double GeneralizedHeston_ImplicitVol(double F,double K,double sig, double t,double omega,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,int nb)

{
	double opt=GeneralizedHeston(F,K,sig,t,omega,theta,ksi,rho,muJ,sigmaJ,lambda,nb);
	return ARM_CF_BS_Formula::callimplicit_totalvolatility(F, 1.,K, 1,opt,ARM_CF_EPS)/sqrt(t);
}


///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=theta*(lgtvol^2-V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///			V(0)=sig^2
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ

double GeneralizedHeston(double F,double K,double sig, double t,double lgtvol,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,int nb)
{
	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double F;
		double K;
		double V0;
		double T;
		double omega;
		double theta;
		double ksi;
		double rho;
		double muJ;
		double sigmaJ;
		double lambda;
		int nb;
	public:
		integralToBeDone(double F_,double K_,double V0_, double T_,double omega_,double theta_,double ksi_,
			double rho_, double muJ_, double sigmaJ_, double lambda_,vector<double> blist, int ptnb ) :
		F(F_),K(K_),V0(V0_),T(T_),omega(omega_),theta(theta_),ksi(ksi_),rho(rho_),muJ(muJ_),sigmaJ(sigmaJ_),lambda(lambda_),
		OscillatoryIntegral(blist,ptnb)
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double k)
		{
			return GeneralizedHeston_Integrand(F,K,k,V0,T,omega,theta,ksi,rho,muJ,sigmaJ,lambda);
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
	vector<double> blist(8);	// list of explicit boundaries
	blist[0]=0.;
	blist[1]=1.5;
	blist[2]=1.55;
	blist[3]=1.557;
	blist[4]=1.562;
	blist[5]=1.564;
	blist[6]=1.56545;
	blist[7]=1.568;

	int legendre_nb_pts=nb;
	integralToBeDone IS(F,K,sig*sig,t,lgtvol*lgtvol*theta,theta,ksi,rho,muJ,sigmaJ,lambda,blist,legendre_nb_pts); 
	double IS_value=IS.value();
	return F-K*IS_value/ARM_NumericConstants::ARM_PI;
}



long double GeneralizedHeston_Integrand(double F,double K,double k, double V0, double t,double omega,double theta,double ksi,
										double rho, double muJ, double sigmaJ, double lambda)
{
	complex<long double> psiplus,psimoins,zeta,X,z,A,B,L,w,expmoinszetatau;
	complex<long double> I(0,1.), Fc(F,0.), Kc(K,0.), kc(k,0.),rhoc(rho,0.),ksic(ksi,0.),thetac(theta,0.),omegac(omega,0.);
	complex<long double> tau(t,0.),deux(2.,0),un(1.,0),muJc(muJ,0.),sigmaJc(sigmaJ,0.),V0c(V0,0.),lambdatau(lambda*t,0.);
	z=-kc-I/deux;
	X=log(Fc/Kc);
	zeta=sqrt((thetac-I*z*rhoc*ksic)*(thetac-I*z*rhoc*ksic)+ksic*ksic*(I*z+z*z));
	if (real(zeta*tau)>700) expmoinszetatau=0; else expmoinszetatau=exp(-zeta*tau);
	psiplus=zeta-thetac+I*z*rhoc*ksic;
	psimoins=zeta+thetac-I*z*rhoc*ksic;
	A=-omegac/(ksic*ksic) * (psiplus*tau+deux*log((psimoins+psiplus*expmoinszetatau)/(deux*zeta)));
	B=-(I*z+z*z)*(un-expmoinszetatau)/(psimoins+psiplus*expmoinszetatau);
	complex<long double> expmujplusigmaj2sur2(exp(muJ+sigmaJ*sigmaJ/2.),0.);
	if(real(muJc*z*I-sigmaJc*sigmaJc*z*z/deux)<-22000.)
	{
		L=-un-I*z*(expmujplusigmaj2sur2-un);
	}
	else
	{
		L=exp(muJc*z*I-sigmaJc*sigmaJc*z*z/deux)-un-I*z*(expmujplusigmaj2sur2-un);
	}
	w=X*I*z+A+B*V0c+L*lambdatau;
//	cout<<"\n k="<<k<<" L="<<L<<" z="<<z<<" muJc*z*I-sigmaJc*sigmaJc*z*z/deux="<<muJc*z*I-sigmaJc*sigmaJc*z*z/deux<<" w="<<w;
	return exp(real(w))*cos(imag(w))/(k*k+0.25);
}


double GeneralizedHeston2(double F,double K,double sig, double t,double lgtvol,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,
						 int nb1, int nb,int NbStage, int NbOscill,double prec)
{
	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double F;
		double K;
		double V0;
		double T;
		double omega;
		double theta;
		double ksi;
		double rho;
		double muJ;
		double sigmaJ;
		double lambda;
		int nb;
	public:
		integralToBeDone(double F_,double K_,double V0_, double T_,double omega_,double theta_,double ksi_,
			double rho_, double muJ_, double sigmaJ_, double lambda_, 
			int nb1_,int nb_,int NbStage_,int NbOscill_,double speed_ ,double prec_ ) :
		F(F_),K(K_),V0(V0_),T(T_),omega(omega_),theta(theta_),ksi(ksi_),rho(rho_),muJ(muJ_),sigmaJ(sigmaJ_),lambda(lambda_),
		OscillatoryIntegral(nb1_,nb_,NbStage_,NbOscill_,speed_,prec_)	/// the controled version
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double k)
		{
			return GeneralizedHeston_Integrand(F,K,k,V0,T,omega,theta,ksi,rho,muJ,sigmaJ,lambda);
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
	integralToBeDone IS(F,K,sig*sig,t,lgtvol*lgtvol*theta,theta,ksi,rho,muJ,sigmaJ,lambda,
		nb1, nb, NbStage, NbOscill,fabs(log(F/K)),prec); 
	double IS_value=IS.value();
	return F-K*IS_value/ARM_NumericConstants::ARM_PI;
}


///**********************************************************************************************
/// essai sur une approach autoadaptative, ne marche pas encore

long double GeneralizedHeston_Integrand_OscillatoryComponent(double F,double K,double k, double V0, double t,double omega,double theta,double ksi,
										double rho, double muJ, double sigmaJ, double lambda)
{
	complex<long double> psiplus,psimoins,zeta,X,z,A,B,L,w,expmoinszetatau;
	complex<long double> I(0,1.), Fc(F,0.), Kc(K,0.), kc(k,0.),rhoc(rho,0.),ksic(ksi,0.),thetac(theta,0.),omegac(omega,0.);
	complex<long double> tau(t,0.),deux(2.,0),un(1.,0),muJc(muJ,0.),sigmaJc(sigmaJ,0.),V0c(V0,0.),lambdatau(lambda*t,0.);
	z=-kc-I/deux;
	X=log(Fc/Kc);
	zeta=sqrt((thetac-I*z*rhoc*ksic)*(thetac-I*z*rhoc*ksic)+ksic*ksic*(I*z+z*z));
	expmoinszetatau=exp(-zeta*tau);
	psiplus=zeta-thetac+I*z*rhoc*ksic;
	psimoins=zeta+thetac-I*z*rhoc*ksic;
	A=-omegac/(ksic*ksic) * (psiplus*tau+deux*log((psimoins+psiplus*expmoinszetatau)/(deux*zeta)));
	B=-(I*z+z*z)*(un-expmoinszetatau)/(psimoins+psiplus*expmoinszetatau);
	complex<long double> expmujplusigmaj2sur2(exp(muJ+sigmaJ*sigmaJ/2.),0.);
	if(real(muJc*z*I-sigmaJc*sigmaJc*z*z/deux)<22000.)
	{
		L=-un-I*z*(expmujplusigmaj2sur2-un);
	}
	else
	{
		L=exp(muJc*z*I-sigmaJc*sigmaJc*z*z/deux)-un-I*z*(expmujplusigmaj2sur2-un);
	}
	w=X*I*z+A+B*V0c+L*lambdatau;
	return imag(w);
}

// here is the autoadaptatif version 
double GeneralizedHeston2(double F,double K,double sig, double t,double omega,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,int nb)
{
	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double F;
		double K;
		double V0;
		double T;
		double omega;
		double theta;
		double ksi;
		double rho;
		double muJ;
		double sigmaJ;
		double lambda;
		int nb;
	public:
		integralToBeDone(double F_,double K_,double V0_, double T_,double omega_,double theta_,double ksi_,
			double rho_, double muJ_, double sigmaJ_, double lambda_,int cpntnb, int ptnb ) :
		F(F_),K(K_),V0(V0_),T(T_),omega(omega_),theta(theta_),ksi(ksi_),rho(rho_),muJ(muJ_),sigmaJ(sigmaJ_),lambda(lambda_),
		OscillatoryIntegral(cpntnb,0.,ptnb)
		{}
		double oscillatorycomponant(double k)
		{
			return GeneralizedHeston_Integrand_OscillatoryComponent(F,K,k,V0,T,omega,theta,ksi,rho,muJ,sigmaJ,lambda);
		}
		double integrand(double k)
		{
			return GeneralizedHeston_Integrand(F,K,k,V0,T,omega,theta,ksi,rho,muJ,sigmaJ,lambda);
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
	integralToBeDone IS(F,K,sig*sig,t,omega,theta,ksi,rho,muJ,sigmaJ,lambda,8,legendre_nb_pts); 
	double IS_value=IS.value();
	return F-K*IS_value/ARM_NumericConstants::ARM_PI;
//	return 0;
}


#undef ARM_CF_EPS

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
