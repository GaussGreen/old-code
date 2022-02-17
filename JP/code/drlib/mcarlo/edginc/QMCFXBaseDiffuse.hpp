//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCFXBaseDiffuse.hpp
//
//   Description :
//
//   Date        : Fri Apr 21 2006
//
//   Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>
//
//
//----------------------------------------------------------------------------

#ifndef QMCFXBaseDiffuse_HPP
#define QMCFXBaseDiffuse_HPP

#include <vector>

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

// We have to avoid circular dependency of headers, as IR refers FX and FX refers to IR;
// so we extract interface that IR calls from FX
class QMCFXBaseDiffuse : public IQMCDiffusibleFX
{
    public:
        /** Prepare to simulate a new path. Returns the sigmaFX for the first
        step */
        virtual double begin(IQMCRNGManagerSP rngMgr) = 0;

    /** Diffuses the fx to the next date. Returns the sigmaFX for the next
        date. idx will be 0 to move from the first date to the second. */
        virtual double moveToNextDate(double forLnMONEY,
                                      int    simDateIdx) = 0;

        /** [Must be] called at the end of the path */
        virtual void end() = 0;

        /** retrieves information for quanto adjustment of non-IR assets (say: credit, etc)*/
        virtual const std::vector<double>& getSigmaFX() const = 0;

        /** retrieves full path for 'struck' currency adjustment of non-IR assets (say: equity, etc)*/
        virtual const std::vector<double>& requestFullSpotPath() = 0;
};

DECLARE(QMCFXBaseDiffuse);


DRLIB_END_NAMESPACE

#endif
