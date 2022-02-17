//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCHelperVectorDiffusionDate
//
//   Description : Abstract class partially implementing  IQMCHelperMaxDiffusionDate via keeping a vector of dependent assets
//
//   Date        : Thu Apr 20 2006
//
//   Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>
//
//----------------------------------------------------------------------------

#ifndef QMCHelperVectorDiffusionDate_HPP
#define QMCHelperVectorDiffusionDate_HPP

#include "edginc/IQMCHelperMaxDiffusionDate.hpp"
#include <vector>

DRLIB_BEGIN_NAMESPACE


// Interface to propagate MaxDiffusion date in an acyclic graph of assets
// IQMCDiffusible assets need to implement updateDependent() method (which requires knowledge of dependent assets)

class QMCHelperVectorDiffusionDate : public  IQMCHelperMaxDiffusionDate {
    private:
        std::vector<IQMCHelperMaxDiffusionDate*> dependent; // list of dependent objects, possibly NULLs (which are skipped)
    protected:
        virtual void        updateDependent(const DateTime& maturity) const;

    public:
        QMCHelperVectorDiffusionDate(); // zero dependent : IR
        QMCHelperVectorDiffusionDate(IQMCHelperMaxDiffusionDate *p);// 1 dependent CR, EQ
        QMCHelperVectorDiffusionDate(IQMCHelperMaxDiffusionDate *p1, IQMCHelperMaxDiffusionDate *p2);
        /** append a new dependent to the list of dates to be notified */
        void add(IQMCHelperMaxDiffusionDate *dep);
};

DRLIB_END_NAMESPACE

#endif
