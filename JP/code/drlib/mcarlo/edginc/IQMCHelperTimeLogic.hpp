//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : IQMCHelperTimeLogic.hpp
//
//   Description :Class to manage DateTime for assets
//
//   Author      : Vladimir Grebinskiy
//
//   Date        :
//
//----------------------------------------------------------------------------

#ifndef EDR_IQMCTIMELOGIC_HPP
#define EDR_IQMCTIMELOGIC_HPP

#include "edginc/DateTime.hpp"
//#include <set>
//#include <algorithm>

#include "edginc/TemplateIdx.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


// Encapsulates IQMC asset time logic
/* This interface captures the current use of state variable generators: we should be able to collect dates from them and then answer the questions like "what is the index of such and such date?"
In addition, some scenarios require us to be careful about memory usage. So the current implementation trades speed efficiency for memory optimization. (The problem is that different DateTimeArrays can be stored in tens of thousands of assets and if each of them contains thousands of elements we have a real problem, as sizeof(DateTime) is 24 bytes on 32 bit and 40 bytes on 64 bit Linux.
*/

class IQMCHelperTimeLogic : public virtual VirtualDestructorBase
{
    public:

    static const int npos = SpotIdx::npos;
    // Add block of Dates
    virtual void   addAggregatedDates(
                        const DateTimeArray& dfDates,
                        const DateTimeArray& edfReqDates,
                        const DateTimeArray& edfFwdDates) = 0;
    virtual void   addAggregatedDates(
                        DateTimeArrayConstSP dfDates,
                        DateTimeArrayConstSP edfReqDates,
                        DateTimeArrayConstSP edfFwdDates) = 0;

    virtual DateTimeArray   getDFDates() const = 0;     ///< return union of all DF requests
    virtual DateTimeArray   getReqEDFDates() const = 0; ///< return union of all requested Dates
    virtual DateTimeArray   getFwdEDFDates() const = 0; ///< return union of all Forward dates

    virtual DateTimeArray   getAssetDates() const {
        return DateTime::merge(getDFDates(), getReqEDFDates());
    }

    virtual SpotIdx             getDFIdx    (const DateTime& date) = 0; ///< returns Idx of a spot date
    virtual FwdIdx              getFwdEDFIdx(const DateTime& date) = 0; ///< returns Idx of a forward date
    virtual SpotIdx             getReqEDFIdx(const DateTime& date) = 0; ///< returns Idx of a requested date
    virtual SpotIdx             getReqEDFIdx(FwdIdx fwdIdx) = 0;        ///< translate between forward and requested indices
    virtual DateTime            getReqDate(SpotIdx idx) = 0; ///< return date such that getReqEDFIdx(date) == idx
    
    virtual void                cleanup() {} ///< does not change external state; mainly needed to save memory in some memory critical cases (ex: CR asset)
    virtual void                debug() {} ///< Auxilary hook
    virtual void                trim(const DateTime& maxDiffMat, const DateTime& maxCurveMat) {} // more memory tricks
    virtual ~IQMCHelperTimeLogic() {}
};

DECLARE(IQMCHelperTimeLogic);

DRLIB_END_NAMESPACE
#endif
