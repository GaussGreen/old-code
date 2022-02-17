//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Complex.cpp
//
//   Description : 
//
//   Author      : Regis Guichard
//
//   Date        : 03 Dec 02
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

bool Complex::isReal() const{
    return Maths::isZero(this->imag());
}

bool Complex::isImag() const{
    return Maths::isZero(this->real());
}

CComplex exp(const CComplex& z){
    d_complex imsl_z = {z.real(), z.imag()};
    return CComplex(imsl_z_exp(imsl_z));
}

CComplex sqrt(const CComplex& z){
    d_complex imsl_z = {z.real(), z.imag()};
    return CComplex(imsl_z_sqrt(imsl_z));
}

CComplex log(const CComplex& z){
    d_complex imsl_z = {z.real(), z.imag()};
    return CComplex(imsl_z_log(imsl_z));
}


DRLIB_END_NAMESPACE
