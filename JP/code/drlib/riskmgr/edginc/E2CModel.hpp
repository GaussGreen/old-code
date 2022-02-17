//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IE2CModel.hpp
//
//   Description : Interface that must be implemented if an Instrument 
//                 
//   Author      : Andre Segger
//
//   Date        : 07 November 2002	
//
//
//----------------------------------------------------------------------------
#ifndef I_E2C_MODEL_HPP
#define I_E2C_MODEL_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL IE2CModel {
public:
    static CClassConstSP const TYPE;

    virtual ~IE2CModel() {}

    virtual void setParSpreadPricing(const bool setPar) = 0;

protected:
    IE2CModel();

};

DRLIB_END_NAMESPACE
#endif
