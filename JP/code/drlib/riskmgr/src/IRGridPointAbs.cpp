//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridPointAbs.hpp
//
//   Description : Identifies a point on an IR vol surface
//
//   Author      : Mark A Robson
//
//   Date        : 17 June 2005
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/IRGridPointAbs.hpp"

DRLIB_BEGIN_NAMESPACE
IRGridPointAbs::~IRGridPointAbs(){}

IRGridPointAbs::IRGridPointAbs(ExpiryConstSP   expiry,
                               const DateTime& swaptionEnd,
                               ExpiryConstSP   tenor,
                               const DateTime& swapMaturity):
    expiry(expiry), swaptionEnd(swaptionEnd),
    tenor(tenor), swapMaturity(swapMaturity) {}

/** Returns an IRGridPoint from this object (IRGridPoint is the type
    returned to the client) */
IRGridPointSP IRGridPointAbs::toGridPoint() const{
    return IRGridPointSP(new IRGridPoint(tenor, expiry));
}

/** returns the expiry associated with this result */
ExpiryConstSP IRGridPointAbs::getExpiry() const{
    return expiry;
}

/** returns the tenor associated with this result */
ExpiryConstSP IRGridPointAbs::getTenor() const{
    return tenor;
}

/** Modify gridPts so that all points whose expiries are after endDate,
    except for the very first, are removed. This is useful when determining
    when to stop tweaking */
void IRGridPointAbs::trim(IRGridPointAbsArray& gridPts, /* (M) */
                          const DateTime&      endDate){
    DateTime firstExpiryAfterEndDate;
    for (unsigned int i = 0; i < gridPts.size(); i++){
        const DateTime& expiry = gridPts[i]->swaptionEnd;
        if (expiry > endDate && 
            (firstExpiryAfterEndDate.empty() || 
             expiry < firstExpiryAfterEndDate)){
            firstExpiryAfterEndDate = expiry;
        }
    }
    if (!firstExpiryAfterEndDate.empty()){
        // junk gridPts whose expiry is after firstExpiryAfterEndDate
        for (IRGridPointAbsArray::iterator iter = gridPts.begin();
             iter != gridPts.end(); /* ++ in loop */){
            const DateTime& expiry = (*iter)->swaptionEnd;
            if (expiry > firstExpiryAfterEndDate){
                iter = gridPts.erase(iter);
            } else {
                ++iter;
            }
        }
    }
}

/** Finds the greatest swapMaturity date in the supplied list of points. The
    supplied array must not be empty */
DateTime IRGridPointAbs::maxSwapMaturity(const IRGridPointAbsArray& gridPts){
    DateTime maxMat(gridPts.front()->swapMaturity);
    for (unsigned int i = 1; i < gridPts.size(); i++){
        const DateTime& swapEnd = gridPts[i]->swapMaturity;
        if (swapEnd > maxMat){
            maxMat = swapEnd;
        }
    }
    return maxMat;
}


DRLIB_END_NAMESPACE
