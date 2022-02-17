//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDSVolRequestSimpleEuropean.hpp
//
//   Description : Vol request for a european option of the same type as the
//                 benchmark instruments in the vol matrix/cube
//
//   Author      : Charles Morcom
//
//   Date        : 20 December 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CDSVOLREQUESTSIMPLEEUROPEAN_HPP
#define EDR_CDSVOLREQUESTSIMPLEEUROPEAN_HPP

#include "edginc/VolRequest.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CDSVolRequestSimpleEuropean)

/** Captures how to interpolate a volatility, i.e., it captures both
    instrument specific data needed to calculate volatilities/variance
    together with the methodology to be used for the interpolation.
    The methodology is captured by the type of the VolRequest used */
class MARKET_DLL CDSVolRequestSimpleEuropean: public CVolRequest {
public:
    static CClassConstSP const TYPE;
    //virtual double getOptionPrice();
    //virtual double getImpliedVolatility();
    virtual ~CDSVolRequestSimpleEuropean();
    CDSVolRequestSimpleEuropean(
        double strike, 
        double forward, 
        bool strikeIsSpread, 
        bool isCall, 
        DateTime exerciseDate, 
        DateTime ulMaturity,
        bool isMultiQ);
    double getStrike() const;
    double getForward() const;
    bool getIsCall() const;
    bool isMultiQ() const;
    const DateTime& getExerciseDate() const;
    const DateTime& getULMaturity() const;
protected:
    CDSVolRequestSimpleEuropean();
private:
    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    /**Strike of option*/
    double strike;
    /**Forward that is relevant to strike*/
    double forward;
    /**True if call, else put*/
    bool isCall;
    /**Option exercise date*/
    DateTime exerciseDate;
    /**Underlying maturity*/
    DateTime ulMaturity;
    /**True if want multi-q computations. Else BS impled vol*/
    bool multiQ;
    /* TODO: add something which specifies the underlying more completely? No
       shouldn't need that, as included in volsurface/cube, assuming that the
       option is the same type.*/

    static void load(CClassSP& clazz);
    friend class CDSVolRequestSimpleEuropeanHelper;
    friend class ICDSVol;
};

DRLIB_END_NAMESPACE

#endif
