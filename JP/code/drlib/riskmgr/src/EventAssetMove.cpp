//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EventAssetMove.cpp
//
//   Description : EventAssetMove class
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EventAssetMove.hpp"

DRLIB_BEGIN_NAMESPACE

EventAssetMove::EventAssetMove(){};

// constructor
void EventAssetMove::SetEvents(HolidayConstSP hol, 
                DateTimeArrayConstSP eventDates,
                DoubleArrayConstSP moveSizes,
                vector<TMove> moveTypes)
{
    Hol = hol;
    EventDates = eventDates;
    MoveSizes = moveSizes;
    MoveTypes = moveTypes;
}
// destructor
EventAssetMove::~EventAssetMove()
{
}

/** precise time point placement rule */
DateTimeArraySP EventAssetMove::getCritDate(int noExerciseWindow, bool isCallTypeAdjustment) const
{
    // rule for call: if jumpUp add pt at EOD
    //                else add pt at previous business EOD

    // rule for put: if jumpDown add pt at EOD
    //                else add pt at previous business EOD

    DateTimeArraySP critDates(new DateTimeArray(EventDates->size()));

    for (int i=0; i<EventDates->size(); i++)
    {
        if ( ((MoveTypes[i]==DOLLAR_DIV || MoveTypes[i]==YIELD)
            && isCallTypeAdjustment) ||
            (!(MoveTypes[i]==DOLLAR_DIV || MoveTypes[i]==YIELD) 
            && !isCallTypeAdjustment))
        {
            (*critDates)[i] = DateTime(Hol->addBusinessDays((*EventDates)[i],-1-noExerciseWindow).getDate(), DateTime::END_OF_DAY_TIME);
        }
        else
        {
            (*critDates)[i] = DateTime((*EventDates)[i].getDate(), DateTime::END_OF_DAY_TIME);
        }
    }
    /*

    bool placeBefore;

    placeBefore =
      (  (*MoveSizes)[i] >= 0.0 &&  isCallTypeAdjustment) || 
      (  (*MoveSizes)[i]  < 0.0 && !isCallTypeAdjustment) )

    if ( placeBefore ) {
          (*critDates)[i] = DateTime(Hol->addBusinessDays((*EventDates)[i],-1-noExerciseWindow).getDate(), DateTime::END_OF_DAY_TIME);
    } else {
          (*critDates)[i] = DateTime((*EventDates)[i].getDate(), DateTime::END_OF_DAY_TIME);
    }
    */

    DateTime::removeDuplicates(*critDates, false);

    return critDates;
}


/** get events dates, sizes and types */
void EventAssetMove::getEvents(DateTimeArrayConstSP& dates, DoubleArrayConstSP& sizes, vector<TMove>& types) const
{
    dates = EventDates;
    sizes = MoveSizes;
    types = MoveTypes;
}

DRLIB_END_NAMESPACE
