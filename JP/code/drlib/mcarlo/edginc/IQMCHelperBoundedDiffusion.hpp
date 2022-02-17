//----------------------------------------------------------------------------
//
//   Group       : xAsset SRM
//
//   Filename    : IQMCHelperBoundedDiffusion.hpp
//
//   Description : Interface that IQMCDiffusibleAssets need to implement to support MaxDiffusion and MaxCurveMaturity
//

//
//----------------------------------------------------------------------------

#ifndef EDR_IQMCBOUNDEDDIFFUSION_HPP
#define EDR_IQMCBOUNDEDDIFFUSION_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Each interested asset can support calculation of MaxDiff(usion)Date.
    It is required for SimpathI assets.
    This class allows one to update and retrieve dates that are needed for diffusion and for Util functions (i.e. last {t_i} and last of {T_j}
 */

class IQMCHelperBoundedDiffusion {
    public:

        /** Returns currently known MaxDiffusionDate */
        virtual DateTime getMaxDiffDate() const = 0;
        /** Returns currently known MaxCurveMaturity */
        virtual DateTime getMaxCurveMat() const = 0;

        virtual void updateMaxDiffDate(const DateTime& date) = 0;
        virtual void updateCurveMat   (const DateTime& date) = 0;
        virtual bool empty() const = 0; ///< true iff updateMaxDiffDate was never called
        virtual ~IQMCHelperBoundedDiffusion() {}
};


DRLIB_END_NAMESPACE
#endif
