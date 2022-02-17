//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestDVF.hpp
//
//   Description : Deterministic Vol Func request
//
//   Author      : JNJ
//
//   Date        : 06 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VOL_REQUEST_DVF_HPP
#define EDR_VOL_REQUEST_DVF_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** Captures how to interpolate a Deterministic Volatility Function */
class MARKET_DLL CVolRequestDVF: public CVolRequest {
public:
    static CClassConstSP const TYPE;

    virtual DateTime getStartDate() const = 0;
    virtual bool isFwdStarting() const = 0;	
    virtual double getStrikeTweakUnscaled() const = 0;
    virtual double getTimeTweakUnscaled() const = 0;
    virtual double getProbDensRatioMin() const = 0;

    virtual bool getUseUnscaledTimeTweak() const = 0;
	virtual bool getUseNextStepDerivs() const = 0;
	virtual bool getUseMidPoint() const = 0;

    virtual ~CVolRequestDVF();
protected:
    CVolRequestDVF(const CClassConstSP& clazz);
private:
    CVolRequestDVF(const CVolRequestDVF &rhs);
    CVolRequestDVF& operator=(const CVolRequestDVF& rhs);
};

// smart pointers for CVolRequestDVF
typedef smartConstPtr<CVolRequestDVF> CVolRequestDVFConstSP;
typedef smartPtr<CVolRequestDVF> CVolRequestDVFSP;

// array of vol request
typedef array<CVolRequestDVFSP, CVolRequestDVF> CVolRequestDVFArray;

// smart pointers for CVolRequestDVFArray
typedef smartConstPtr<CVolRequestDVFArray> CVolRequestDVFArrayConstSP;
typedef smartPtr<CVolRequestDVFArray> CVolRequestDVFArraySP;

DRLIB_END_NAMESPACE

#endif
