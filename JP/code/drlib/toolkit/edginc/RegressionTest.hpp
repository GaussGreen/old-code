//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskMgrInterface.hpp
//
//   Description : Defines interface to RiskMgr
//
//   Author      : Andrew J Swain
//
//   Date        : 19 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_REGRESSIONTEST_HPP
#define EDR_REGRESSIONTEST_HPP

#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL IRegressionTest{
public:
    virtual ~IRegressionTest(){}

    virtual IObjectSP runTest() const = 0;
protected:
};


DRLIB_END_NAMESPACE

#endif
