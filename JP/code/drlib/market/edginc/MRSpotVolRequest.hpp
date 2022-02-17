//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : MRSpotVolRequest.hpp
//
//   Description : A way of asking for mean reverting 'spot' vols. 
//                 Here 'spot' means in the traditional IR sense
//
//   Author      : Mark A Robson
//
//   Date        : 14 Dec 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SPOT_VOLREQUEST_HPP
#define EDR_SPOT_VOLREQUEST_HPP

#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** A way of asking for 'spot' vols. 
    Here 'spot' means in the traditional IR sense */
class MARKET_DLL MRSpotVolRequest: public CVolRequest {
public:
    static CClassConstSP const TYPE;
    //// initially this is going to be used for a flat cds vol curve so we
    //// can defer for now what interpolation data we need
    MRSpotVolRequest();

private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
    MRSpotVolRequest(const MRSpotVolRequest &rhs);
    MRSpotVolRequest& operator=(const MRSpotVolRequest& rhs);
};

typedef smartConstPtr<MRSpotVolRequest> MRSpotVolRequestConstSP;
typedef smartPtr<MRSpotVolRequest> MRSpotVolRequestSP;

DRLIB_END_NAMESPACE

#endif
