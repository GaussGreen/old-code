//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EndDateCollector.cpp
//
//   Description : End date collector class
//
//   Author      : Andrew McCleery
//
//   Date        : 28 May 2004
//
//
//   $Log:
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EndDateCollector.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ObjectIteration.hpp"

DRLIB_BEGIN_NAMESPACE

EndDateCollector::EndDateCollector(const DateTime&    endDate,
                                   const Sensitivity* sensitivity):
    endDate(endDate),  sensitivity(sensitivity) {}

DateTime EndDateCollector::getMaxEndDate(IObjectConstSP obj) const{
    // class for handling call back
    class Action: public ObjectIteration::IActionConst{
    public:
        DateTime           endDate;
        const Sensitivity* sensitivity;
        Action(const DateTime& endDate, const Sensitivity* sensitivity):
            endDate(endDate), sensitivity(sensitivity){}

        // called by ObjectIteration
        bool invoke(const ObjectIteration::State& state, IObjectConstSP obj) {
            // Cast matched component object to LastSensDate and
            // compare endDate to max so far
            const LastSensDate* lsd = dynamic_cast<const LastSensDate *>(obj.get());
            DateTime newEndDate = lsd->endDate(sensitivity);
            if (newEndDate > endDate) {
                endDate = newEndDate;
            }
            // Keep recursing through sub-components
            return true;
        }
    };
    Action action(endDate, sensitivity);
    // Build an instance of the class that drives the recursion
    ObjectIteration iteration(LastSensDate::TYPE);
    // Recurse over components of obj
    iteration.recurse(action, obj);

    return action.endDate;
}

DRLIB_END_NAMESPACE
