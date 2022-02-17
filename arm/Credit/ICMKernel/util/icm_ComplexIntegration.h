//		© Fixed Income R&D 
//		Florent Serre - 7/10/2003
//

//		History
//
//		October 2003			création

#ifndef _F1_ComplexIntegration_H_
#define _F1_ComplexIntegration_H_

#include "icm_functors.h"
#include <math.h>

// Compute integral
//     /- b
//	  / 
//   /    f(x) exp(-x²)dx
//  /
//-/ a

template <class F> class ComplexIntegration_t
{
public:
	ComplexIntegration_t(F f):m_f(f){}

	virtual ~ComplexIntegration_t() {}

	// Simpson Method
	void SimpsonMethod(double x_min,
					   double x_max,
					   int nbSteps,
					   double& reRes,
					   double& imRes); // real part and imaginary part

	// Adaptative Simpson Method
	void AdaptativeSimpsonMethod(double x_min,
								 double x_max,
								 int nbSteps,
								 double& reRes,
								 double& imRes,
								 double eps=1.e-4);

	void	GaussLegendreMethod(double x_min,
								 double x_max,
								 int nbSteps,
								 double& reRes,
								 double& imRes,
								 double eps=1.e-11);

	void	GaussLaguerreMethod(int degree,
								double alpha,
								double& reRes,
								double& imRes,
								int NbIterMax=10,
								double eps=1.e-11);

	void	QuasiMonteCarlo(int nbIter,
							double& reRes,
							double& imRes,
							double eps=1.e-11);


	void Trapeze(double a,
				double b,
				int n,
				double& reRes,
				double& imRes);
protected:
	void CoefficientLegendre(double x1,
							double x2,
							std::vector<double>& x,
							std::vector<double>& w,
							int n,
							double eps=1.e-11);

	void CoefficientLaguerre(double alpha,
						     std::vector<double>& x,
							 std::vector<double>& w,
							 int NbIterMax=10,
							 double eps=1.e-11);

public:
	static bool test(std::string& errStr);
private :
	F	m_f;

	// for Legendre && Laguerre
	static std::vector<double> m_x;
	static std::vector<double> m_w;
	static double	m_x_min; // = alpha for laguerre
	static double	m_x_max;
	static	int	m_size;
};
//----------------------------------------------------------------------------
//	Helper. 
template <class F> inline ComplexIntegration_t<F> ComplexIntegration(F f)
{ return ComplexIntegration_t<F>(f) ; }
//----------------------------------------------------------------------------
template <class F>
double
ComplexIntegration_t<F>::m_x_min = 0.;
//----------------------------------------------------------------------------
template <class F> 
double
ComplexIntegration_t<F>::m_x_max = 0.;
//----------------------------------------------------------------------------
template <class F> 
int
ComplexIntegration_t<F>::m_size = 0;
//----------------------------------------------------------------------------
template <class F> 
std::vector<double>
ComplexIntegration_t<F>::m_x = std::vector<double>(0);
//----------------------------------------------------------------------------
template <class F> 
std::vector<double>
ComplexIntegration_t<F>::m_w = std::vector<double>(0);
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::SimpsonMethod(double x_min,
									   double x_max,
									   int nbSteps,
									   double& reRes,
									   double& imRes)
{
	/* checks */
	EP_ASSERT (nbSteps>0);
	EP_ASSERT (x_min<x_max);

	/* compute */
	double  lag_x    = (x_max - x_min) / double(nbSteps);

	reRes = 0.;
	imRes = 0.;

	for (int i(0); i<=nbSteps; i++)
	{
		double x = x_min + i * lag_x;
		double factor = (i==0||i==nbSteps) ? 0.5 * lag_x : lag_x ;
		double* res = m_f(x);
		reRes += factor * res[0];
		imRes += factor * res[1];
	}

}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::Trapeze(double a,
								double b,
								int n,
								double& reRes,
								double& imRes)
{
	if (n == 1) {
		double* res_a = m_f(a);
		double* res_b = m_f(b);
		reRes = 0.5*(b-a)*(res_a[0] + res_b[0]);
		imRes = 0.5*(b-a)*(res_a[1] + res_b[1]);
	}
	else {
		for (int it=1,j=1;j<n-1;j++) it <<= 1;
		double tnm = it;
		double del = (b-a)/tnm;
		double x = a + 0.5*del;
		double reSum(0.);
		double imSum(0.);

		for (j=1;j<=it;j++,x+=del) {
			double* res = m_f(x);
			reSum += res[0];
			imSum += res[1];
		}
		reRes = 0.5*(reRes+(b-a)*reSum/tnm);
		imRes = 0.5*(imRes+(b-a)*imSum/tnm);
	}
}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::AdaptativeSimpsonMethod(double x_min,
												 double x_max,
												 int	nbSteps,
												 double& reRes,
												 double& imRes,
												 double eps)
{
	reRes = 0.;
	imRes = 0.;

	double olds = -1.0e30;

	for (int j=1;j<=nbSteps;j++) {

		Trapeze(x_min, x_max, j, reRes, imRes);
		// Check CVG
		if (fabs((reRes*reRes+imRes*imRes)-olds) < eps*fabs(olds))			return;
		else if (reRes*reRes+imRes*imRes == 0.0 && olds == 0.0 && j > 6)	return;
		olds = reRes*reRes+imRes*imRes;

	}

	EP_THROWEX("Too many steps in routine integrate");
}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::CoefficientLegendre(double x1,
											 double x2,
											 std::vector<double>& x,
											 std::vector<double>& w,
											 int n,
											 double eps)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m = (n+1)/2;
	xm = 0.5*(x2+x1);
	xl = 0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		}
		while (fabs(z-z1) > eps);
		//EP_LOG << "m=" << m << "n=" << n << "   i-1="<< i-1 << "   n+1-i+1" << n+1-i+1 << std::endl;
		x[i-1] = xm-xl*z;
		x[n-i] = xm+xl*z;
		w[i-1] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n-i] = w[i-1];
	}
}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::GaussLegendreMethod(double x_min,
											 double x_max,
											 int	 size,
											 double& reRes,
											 double& imRes,
											 double eps)
{
	if((x_min!=m_x_min)||(x_max!=m_x_max)||(m_size!=size)) {
		EP_ASSERT(size>0)
		m_size = size;
		m_x_max = x_max;
		m_x_min = x_min;
		m_x.resize(m_size);
		m_w.resize(m_size);	
		CoefficientLegendre(m_x_min, m_x_max, m_x, m_w, m_size); // calcul des coeff
	}

	reRes = 0.;
	imRes = 0.;

	for (int i=0;i<m_size;i++) {
		double* res = m_f(m_x[i]);
		reRes += m_w[i]*res[0];
		imRes += m_w[i]*res[1];
	}
}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::CoefficientLaguerre(double alpha,
											 std::vector<double>& x,
											 std::vector<double>& w,
											 int NbIter,
											 double eps)
{
	double ai, p1(0.), p2(0.), p3(0.), pp(0.), z(0.), z1; //High precision is a good idea for this routine.
	
	for (int i=1;i<=m_size;i++) { //Loop over the desired roots.
		
		if (i == 1) { //Initial guess for the smallest root.
			z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*m_size+1.8*alpha);
		} 
		else if (i == 2) { //Initial guess for the second root.
			z += (15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*m_size);
		}
		else { //Initial guess for the other roots.
			ai = i-2.;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alpha/(1.0+3.5*ai))*(z-x[i-2-1])/(1.0+0.3*alpha);
		}

		for (int its=1;its<=NbIter;its++) { //Refinement by Newton’s method.
			p1 = 1.;
			p2 = 0.;
			for (int j=1;j<=m_size;j++) { //Loop up the recurrence relation to get the Laguerre polynomial evaluated at z.
				p3=p2;
				p2=p1;
				p1=((2.*j-1.+alpha-z)*p2-(j-1.+alpha)*p3)/((double)j);
			}
			//p1 is now the desired Laguerre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp = (m_size*p1-(m_size+alpha)*p2)/z;
			z1 = z;
			z = z1-p1/pp; //Newton’s formula.
			if (fabs(z-z1) <= eps) break;
		}
	
		if (its > NbIter) EP_THROWEX("too many iterations in gauss laguerre");
		x[i-1] = z; //Store the root and the weight.
		w[i-1] = -exp(ep::MathSrv::gammln(alpha+m_size)-ep::MathSrv::gammln((double)m_size))/(pp*m_size*p2);
		
	}
	
}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::GaussLaguerreMethod(int	 degree,
											 double  alpha,
											 double& reRes,
											 double& imRes,
											 int NbIterMax,
											 double eps)
{
	if((m_size!=degree)||(alpha!=m_x_min)) {
		EP_ASSERT(degree>0)
		m_size = degree;
		m_x_min = alpha;
		m_x.resize(m_size);
		m_w.resize(m_size);	
		CoefficientLaguerre(alpha, m_x, m_w, NbIterMax, eps); // calcul des coeff
	}
	
	reRes = 0.;
	imRes = 0.;

	for (int i=0;i<m_size;i++) {
		double* res = m_f(m_x[i]);
		reRes += m_w[i]*res[0];
		imRes += m_w[i]*res[1];
	}
}
//----------------------------------------------------------------------------
template <class F>
void
ComplexIntegration_t<F>::QuasiMonteCarlo(int nbIter,
										double& reRes,
										double& imRes,
										double eps)
{
	EP_ASSERT(nbIter>0)

	reRes = 0.;
	imRes = 0.;

	EpGrayFaure::sobolseqinit (1);
	ep::SimpleVector<double> u(1);

	for(int k=1;k<=nbIter;k++) {

		EpGrayFaure::sobolseq(1, u);
		double x = ep::MathSrv::invCumNorm (u[0]);
		double* res = m_f(x);
		reRes += (res[0]-reRes)/((double)k);
		imRes += (res[1]-imRes)/((double)k);
	}
}
//----------------------------------------------------------------------------
// For test purpose
class BEx {
	public: 
	BEx(){c = new double[2]; c[0] = 1.; c[1] = 1.;};
	virtual ~BEx(){delete [] c; c=NULL;};
	double* f (double t) {c[0] = exp(-t*t); c[1] = exp(-t*t); return c;}
	double* f2 (double t) {c[0] = 1.; c[1] = 1.; return c;}
	private:
	double* c;
};
//----------------------------------------------------------------------------
template <class F>
bool
ComplexIntegration_t<F>::test(std::string& errStr)
{
	bool ret(true);

	errStr = "Testing ComplexIntegration";

	BEx b;
	double reRes, imRes;

	ComplexIntegration(f1::mem_call(&BEx::f, b)).SimpsonMethod(-10,10, 1000, reRes, imRes);
	SDMTEST((reRes>1.77245385)&&(reRes<1.77245386))
	SDMTEST((imRes>1.77245385)&&(imRes<1.77245386))

	ComplexIntegration(f1::mem_call(&BEx::f2, b)).GaussLaguerreMethod(10, 1., reRes, imRes);
	SDMTEST((reRes>0.999999)&&(reRes<1.00000001))
	SDMTEST((imRes>0.999999)&&(imRes<1.00000001))

	return ret;
}
#endif	// _F1_ComplexIntegration_H_


