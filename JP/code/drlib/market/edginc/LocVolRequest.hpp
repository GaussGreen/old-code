//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LocVolRequest.hpp
//
//   Description : Local vol request 
//
//   Author      : Regis Guichard
//
//   Date        : 24 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef LOC_VOLREQUEST_HPP
#define LOC_VOLREQUEST_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequestDVF.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL LocVolRequest: public CVolRequestDVF {
public:
    static CClassConstSP const TYPE;
    friend class LocVolRequestHelper;
    
    LocVolRequest(DateTime        startDate,
                  bool            fwdStarting,
                  bool            useUnscaledTimeTweak,
                  bool            useNextStepDerivs,
                  bool            useMidPoint,
                  double          strikeTweakUnscaled,
                  double          timeTweakUnscaled,
                  double          probDensRatioMin);

    LocVolRequest(DateTime        startDate,
                  bool            fwdStarting,
                  bool            useUnscaledTimeTweak,
                  bool            useNextStepDerivs,
                  bool            useMidPoint,
                  double          strikeTweakUnscaled,
                  double          timeTweakUnscaled,
                  double          probDensRatioMin,
                  string          speedup);

    /** Returns the start date for forward starting volatility 
        interpolation. */
    virtual DateTime getStartDate() const;
    virtual bool isFwdStarting() const;	
    virtual double getStrikeTweakUnscaled() const;
    virtual double getTimeTweakUnscaled() const;
    virtual double getProbDensRatioMin() const;

    virtual bool getUseUnscaledTimeTweak() const;
    virtual bool getUseNextStepDerivs() const;
    virtual bool getUseMidPoint() const;

    virtual string getSpeedup() const; 

protected:
    DateTime                startDate;
    bool                    fwdStarting;
    double                  strikeTweakUnscaled;
    double                  timeTweakUnscaled;
    double                  probDensRatioMin;

    bool                    useUnscaledTimeTweak;
    bool                    useNextStepDerivs;
    bool                    useMidPoint;

    string                  speedup;
private:
    LocVolRequest();
    LocVolRequest(const LocVolRequest &rhs);
    LocVolRequest& operator=(const LocVolRequest& rhs);
};

typedef smartConstPtr<LocVolRequest> LocVolRequestConstSP;
typedef smartPtr<LocVolRequest> LocVolRequestSP;

DRLIB_END_NAMESPACE

#endif
