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
#include "firsttoinc.h"
#include "gpbase/port.h"


#include <cmath>
#include <complex>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/simple_heston.h"
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



CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-14

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Process:  dSt=Sqrt[Vt] dW1t
///				dVt=kappa(theta-Vt) dt + epsolon Sqrt[Vt] dW2t
///				dW1t dW2t = rho dt
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

long double Simple_Heston_Integrand(double F,double K,double k, double V0, double t,double kappa,double theta,double ksi,
										double rho)
{
	complex<long double> psiplus,psimoins,zeta,X,z,A,B,L,w,expmoinszetatau;
	complex<long double> I(0,1.), Fc(F,0.), Kc(K,0.), kc(k,0.),rhoc(rho,0.),ksic(ksi,0.),kappac(kappa,0.),thetac(theta,0.);
	complex<long double> tau(t,0.),deux(2.,0),un(1.,0),V0c(V0,0.);
	z=kc+I/deux;
	X=log(Fc/Kc);
	zeta=sqrt((kappac+I*z*rhoc*ksic)*(kappac+I*z*rhoc*ksic)+ksic*ksic*(z*z-I*z));
	if (real(zeta*tau)>700) expmoinszetatau=0; else expmoinszetatau=exp(-zeta*tau);
	psiplus=zeta-kappac-I*z*rhoc*ksic;
	psimoins=zeta+kappac+I*z*rhoc*ksic;
	A=-thetac*kappac/(ksic*ksic) * (psiplus*tau+deux*log((psimoins+psiplus*expmoinszetatau)/(deux*zeta)));
	B=-(z*z-I*z)*(un-expmoinszetatau)/(psimoins+psiplus*expmoinszetatau);
	w=-X*I*z+A+B*V0c;
	return exp(real(w))*cos(imag(w))/(k*k+0.25);
}




double Shifted_Heston(double S,double K,double V, double t,double kappa,double theta,double ksi,
						 double rho, double m,int nb1, int nb,int NbStage, int NbOscill,double prec)
{
	if (fabs(m/S)<1e-7)
	{
		return Simple_Heston(S,K,t,V,kappa,theta,ksi,rho,nb1,nb,NbStage,NbOscill,prec);
	}
	else
	{
		if (m>0)
		{
			double Ka=(1.-m)*S+m*K;
			double Va=m*m*V;
			double ksia=m*ksi;
			double thetaa=theta*m*m;
			return Simple_Heston(S,Ka,t,Va,kappa,thetaa,ksia,rho,nb1,nb,NbStage,NbOscill,prec)/m;
			
		}
		else
		{
			double Ka=(1.-m)*S+m*K;
			double Va=m*m*V;
			double ksia=-m*ksi;
			double thetaa=theta*m*m;
			return ((Ka-S)+Simple_Heston(S,Ka,t,Va,kappa,thetaa,ksia,rho,nb1,nb,NbStage,NbOscill,prec))/(-m);
		}
	}
	
}



double Simple_Heston(double F,double K,double t, double V0,double kappa,double theta,double ksi,
						 double rho, int nb1, int nb,int NbStage, int NbOscill,double prec)
{
	class integralToBeDone : public OscillatoryIntegral
	{
	private:
		double F;
		double K;
		double V0;
		double T;
		double kappa;
		double theta;
		double ksi;
		double rho;

	public:
		integralToBeDone(double F_,double K_,double V0_, double T_,double kappa_,double theta_,double ksi_,
			double rho_,  int nb1_,int nb_,int NbStage_,int NbOscill_,double speed_ ,double prec_) :
		F(F_),K(K_),V0(V0_),T(T_),kappa(kappa_),theta(theta_),ksi(ksi_),rho(rho_),
			OscillatoryIntegral(nb1_,nb_,NbStage_,NbOscill_,speed_,prec_)	/// the controled version
		{}
		double oscillatorycomponant(double k)
		{
			return 0;
		}
		double integrand(double k)
		{
			return Simple_Heston_Integrand(F,K,k,V0,T,kappa,theta,ksi,rho);
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
/// here log(F/K) gives the number of period per unit for the integral
	integralToBeDone IS(F,K,V0,t,kappa,theta,ksi,rho,nb1, nb, NbStage, NbOscill,fabs(log(F/K)),prec); 
//	integralToBeDone IS(F,K,V0,t,kappa,theta,ksi,rho,nb1, nb, NbStage, NbOscill,0.01,prec); 
	double IS_value=IS.value();
	return F-K*IS_value/ARM_NumericConstants::ARM_PI;
}




long double Simple_Heston_Integrand_OscillatoryComponent(double F,double K,double k, double V0, double t,double omega,double theta,double ksi,
										double rho)
{
	return 0;
}




#undef ARM_CF_EPS

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
