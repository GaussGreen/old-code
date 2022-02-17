

#ifndef _F1_COMPLEXHERMITEINTEGRATION_H_
#define _F1_COMPLEXHERMITEINTEGRATION_H_

#include "icm_functors.h"
#include <math.h>
//#include "SDMTrace.h"

// Compute integral
//     /- +inf
//	  / 
//   /    f(x) exp(-x²)dx
//  /
//-/-inf
// with an hermite approximation at degree 20

/*********************************************************************************/
/*! \class  ICM_Cache_Distrib icm_ComplexHermiteIntegration.h "icm_ComplexHermiteIntegration.h"
 *  \author F Dezormes
 *	\version 1.0
 *	\date   March 2004
 *	\file   icm_ComplexHermiteIntegration.h
 *	\brief  for hermite integral calculations
/***********************************************************************************/

template <class F> class ComplexHermiteIntegration_t
{
public:
	ComplexHermiteIntegration_t(F f):m_f(f){}
	virtual ~ComplexHermiteIntegration_t() {}
	void Integrate(double& reRes, double& imRes); // real part and imaginary part

	static bool test(std::string& errStr);
private :
	F	m_f;
	static double HermiteCoeff[20][2];
};
//----------------------------------------------------------------------------
//	Helper. 
template <class F> inline ComplexHermiteIntegration_t<F> ComplexHermiteIntegration(F f)
{ return ComplexHermiteIntegration_t<F>(f) ; }
//----------------------------------------------------------------------------
template <class F> double ComplexHermiteIntegration_t<F>::HermiteCoeff[20][2] = { 
						{0.2453407083009 , 4.622436696006*pow(10.,-1)},
						{0.7374737285454 , 2.866755053628*pow(10.,-1)},
						{1.2340762153953 , 1.090172060200*pow(10.,-1)},
						{1.7385377121166 , 2.481052088746*pow(10.,-2)},
						{2.2549740020893 , 3.243773342238*pow(10.,-3)},
						{2.7888060584281 , 2.283386360163*pow(10.,-4)},
						{3.3478545673832 , 7.802556478532*pow(10.,-6)},
						{3.9447640401156 , 1.086069370769*pow(10.,-7)},
						{4.6036824495507 , 4.399340992273*pow(10.,-10)},
						{5.3874808900112 , 2.229393645534*pow(10.,-13)},
						{-0.2453407083009 , 4.622436696006*pow(10.,-1)},
						{-0.7374737285454 , 2.866755053628*pow(10.,-1)},
						{-1.2340762153953 , 1.090172060200*pow(10.,-1)},
						{-1.7385377121166 , 2.481052088746*pow(10.,-2)},
						{-2.2549740020893 , 3.243773342238*pow(10.,-3)},
						{-2.7888060584281 , 2.283386360163*pow(10.,-4)},
						{-3.3478545673832 , 7.802556478532*pow(10.,-6)},
						{-3.9447640401156 , 1.086069370769*pow(10.,-7)},
						{-4.6036824495507 , 4.399340992273*pow(10.,-10)},
						{-5.3874808900112 , 2.229393645534*pow(10.,-13)} };
//----------------------------------------------------------------------------
template <class F>
void
ComplexHermiteIntegration_t<F>::Integrate(double& reRes, double& imRes)
{
	reRes = 0.;
	imRes = 0.;
	for(int i=0;i<20;i++) {
		double* res = m_f(HermiteCoeff[i][0]);
		reRes += res[0]*HermiteCoeff[i][1];
		imRes += res[1]*HermiteCoeff[i][1];
	}
}
//----------------------------------------------------------------------------
// For test purpose
class B {
	public: 
	B(){c = new double[2]; c[0] = 1.; c[1] = 1.;};
	~B(){delete [] c; c=NULL;};
	double* f (double) {return c;}
	private:
	double* c;
};
//----------------------------------------------------------------------------
template <class F>
bool
ComplexHermiteIntegration_t<F>::test(std::string& errStr)
{
	bool ret(true);

	errStr = "Testing ComplexHermiteIntegration";

	B b;
	double reRes, imRes;

	ComplexHermiteIntegration(f1::mem_call(&B::f, b)).Integrate(reRes, imRes);

//	SDMTEST((reRes>1.77245385)&&(reRes<1.77245386))
//	SDMTEST((imRes>1.77245385)&&(imRes<1.77245386))

	return ret;
}
#endif	// _F1_COMPLEXHERMITEINTEGRATION_H_


