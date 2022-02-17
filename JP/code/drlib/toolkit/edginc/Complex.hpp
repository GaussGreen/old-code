//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Complex.hpp
//
//   Description : 
//
//   Author      : Regis Guichard
//
//   Date        : 23 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_COMPLEX_HPP
#define EDR_COMPLEX_HPP

#if defined(_MSC_VER)
// disable warning
#pragma warning(disable : 4275)
#endif

#include <complex>
#include <cmath>
#include "edginc/imsl.h"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL CComplex: public complex<double>{
public:
    friend TOOLKIT_DLL CComplex exp(const CComplex& z);
    friend TOOLKIT_DLL CComplex sqrt(const CComplex& z);
    friend TOOLKIT_DLL CComplex log(const CComplex& z);

    explicit CComplex(double re, double im): complex<double>(re, im){}
    explicit CComplex(){}
    CComplex(double re): complex<double>(re){}
    CComplex(const complex<double>& rhs): complex<double>(rhs){}
	CComplex& operator=(double re){
        return (*this = CComplex(re, 0.0));
    }

    static double absSquare(const CComplex& z){
        return (z.real() * z.real() + z.imag() * z.imag());
    }
	
	//GAD 19/02/2006
    
    //computes the argument of a complex
    static double argument(const CComplex& z){
        double res = atan( z.imag() / z.real() ); 
        return res;
    }

    //end of GAD 19/02/2006

    bool isReal() const;
    bool isImag() const;

private:
    explicit CComplex(const d_complex& rhs): complex<double>(rhs.re, rhs.im){}
};
typedef CComplex Complex;

TOOLKIT_DLL CComplex exp(const CComplex& z);
TOOLKIT_DLL CComplex sqrt(const CComplex& z);
TOOLKIT_DLL CComplex log(const CComplex& z);

using ::log;
using ::exp;
using ::sqrt;


typedef vector<CComplex> ComplexArray;

DRLIB_END_NAMESPACE

#endif
