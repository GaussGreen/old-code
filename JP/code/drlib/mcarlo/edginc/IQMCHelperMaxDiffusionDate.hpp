//----------------------------------------------------------------------------
//
//   Group       : xAsset SRM
//
//   Filename    : IQMCDiffusionMaturity.hpp
//
//   Description : all assets share the same timeLine, but we only need dates that are listed in the Generators
//

//
//----------------------------------------------------------------------------

#ifndef IQMCHelperMaxDiffusionDate_HPP
#define IQMCHelperMaxDiffusionDate_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

// Interface to propagate MaxDiffusion date in an acyclic graph of assets
// IQMCDiffusible assets need to implement updateDependent() method (which requires knowledge of dependent assets)

class IQMCHelperMaxDiffusionDate {
    DateTime    date;
    bool        initialized;
    protected:
        virtual void        updateDependent(const DateTime& maturity)  const = 0; // call updateMaturity on dependent
    public:
        DateTime            currentMaturity(void) const {return date;} // returns currently known max maturity
        void                updateMaturity (const DateTime& maturity) {
            if (! initialized || currentMaturity() < maturity) {
                initialized = true;
                updateDependent(date=maturity);
            }
        }

        IQMCHelperMaxDiffusionDate() : date(0, 0), initialized(false) {}
        bool    empty() const {return !initialized;} // true iff updateMaturity was never called

        virtual ~IQMCHelperMaxDiffusionDate() {}

};
DRLIB_END_NAMESPACE
#endif //EDR_IQMCDIFFUSIBLEASSETBASE_HPP

