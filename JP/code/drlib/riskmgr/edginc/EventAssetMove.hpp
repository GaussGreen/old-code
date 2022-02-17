//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EventAssetMove.hpp
//
//   Description : EventAssetMove class
//
//
//----------------------------------------------------------------------------

#ifndef EDG_EVENT_ASSET_MOVE_HPP
#define EDG_EVENT_ASSET_MOVE_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL EventAssetMove
{
public:

    typedef enum{DOLLAR_DIV, YIELD, OTHER} TMove; // what else ? 

    EventAssetMove();

    void SetEvents(HolidayConstSP hol, 
                DateTimeArrayConstSP eventDates,
                DoubleArrayConstSP moveSizes,
                vector<TMove> moveType);

    virtual ~EventAssetMove();

    /** precise time point placement rule */
    virtual DateTimeArraySP getCritDate(int noExerciseWindow, bool isCallTypeAdjustment) const;

    /** access function: get event dates, sizes and types */
    virtual void getEvents(DateTimeArrayConstSP& dates, DoubleArrayConstSP& sizes, vector<TMove>& types) const;

private:
    HolidayConstSP              Hol;
    DateTimeArrayConstSP        EventDates;
    vector<TMove>          MoveTypes; 
    DoubleArrayConstSP     MoveSizes;
};

DRLIB_END_NAMESPACE
#endif
