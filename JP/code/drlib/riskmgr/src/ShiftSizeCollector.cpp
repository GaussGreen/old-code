//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ShiftSizeCollector.cpp
//
//   Description : Delta Shift size collector class
//
//   Author      : André Segger
//
//   Date        : 20 Jun 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
/** Invoked when TestCollect is 'loaded' */
void ShiftSizeCollector::load(CClassSP& clazz){
    REGISTER(ShiftSizeCollector, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICollector);
}

CClassConstSP const ShiftSizeCollector::TYPE = CClass::registerClassLoadMethod(
    "ShiftSizeCollector", typeid(ShiftSizeCollector), load);

ShiftSizeCollector::ShiftSizeCollector(const Delta*          sensControl,
                                       const double          strike,
                                       const TAdjustmentType adjType):
        CObject(TYPE), shiftSize(0.0), strike(strike), adjType(adjType)
{
    if ( adjType != FWD_START_ADJUSTMENT      &&
         adjType != SPOT_START_ADJUSTMENT )
    {
        throw ModelException(
            "ShiftSizeCollector::ShiftSizeCollector",
            "Delta shift size ajustment type has to be either forward or"
            "spot starting");
    }

    sensCtrl  = DeltaConstSP(DeltaConstSP::attachToRef(sensControl));
    shiftSize = sensControl->getShiftSize();
}

/** Vols should call this with their asset's current spot price 
    through the IVolDeltaShiftSize interface. 
    If the information is not available
    then this routine should not be called. This method then
    alters the shift size in view of the information. */
void ShiftSizeCollector::adjustShiftSize(const DoubleArraySP& strikeList,
                                         const string&        assetName,
                                         const double         spot)
{
    OutputNameSP name(new OutputName(assetName));


//    if (adjType == FWD_START_ADJUSTMENT && 
//        sensCtrl->marketDataNameMatches(name)){

    if (adjType == FWD_START_ADJUSTMENT) {
        double       loBound, hiBound;
        double       newShift;
        int          numLevels     = strikeList->size();
        int          i;
        
        // check that we have a valid strike list
        if ( numLevels < 1 ) {
            throw ModelException(
                "ShiftSizeCollector::shiftSizeValidate",
                "Empty strike list passed to shift size validation");
        }
        
        // find bounding points in list of levels
        if (numLevels == 1) {
            loBound = hiBound = spot;
        } else {
            // handle case where spot is above end of levels
            if ( spot > (*strikeList)[numLevels - 1] ) {
                loBound = (*strikeList)[numLevels - 1];
                hiBound = spot;
            }
            // ... and where spot is below the end of levels
            else if ( spot < (*strikeList)[0] )
            {
                loBound = spot;
                hiBound = (*strikeList)[0];
            } else  {
                // loop round to find the upper bounding column
                for (i = 1; spot > (*strikeList)[i]; i++) {
                    // empty loop
                }
                // Check that spot is not a point in the levels, if 
                //it is then it would be strikeList[i] 
                if ( Maths::equals((*strikeList)[i], spot)) {
                    loBound = hiBound = spot;
                } else {
                    loBound = (*strikeList)[i-1];
                    hiBound = (*strikeList)[i];
                }
            }
        }
        
        if ( Maths::equals(loBound, spot)) {
            if ( !Maths::equals(hiBound, spot)) {
                
                newShift = Maths::max( (hiBound - spot)/(2.0*spot),
                                       Delta::MINIMUM_SHIFT);
                if (shiftSize < newShift) {
                    newShift = 0.0; // Don't want to use the new shift
                }
            } else {
                // If we're on a curve, or the spot is a strike don't change 
                // shift unless specified shift is less than minimum shift
                if (shiftSize < Delta::MINIMUM_SHIFT) {
                    newShift = Delta::MINIMUM_SHIFT;
                } else {
                    // Don't want to use the new shift
                    newShift = 0.0;
                }
            }
        } else {
            // Either we're in the surface, or off the top of it
            if ( !Maths::equals(hiBound,spot) ) {
                newShift = Maths::max(Delta::MINIMUM_SHIFT,
                                      Maths::min((spot - loBound)/(2.0*spot),
                                                 (hiBound - spot)/(2.0*spot)));
                if (shiftSize < newShift) {
                    // Don't want to use the new shift
                    newShift = 0.0;
                }
            } else {
                newShift = Maths::max( (spot - loBound)/(2.0*spot),
                                       Delta::MINIMUM_SHIFT);
                if (shiftSize < newShift) {
                    // Don't want to use the new shift
                    newShift = 0.0;
                }
            }
        }           
        
        if (newShift >= Delta::MINIMUM_SHIFT &&  newShift < shiftSize) {
            shiftSize = newShift;
        }
    }
}

double ShiftSizeCollector::getShiftSize() const
{
    return shiftSize;
}

void ShiftSizeCollector::setShiftSize(double newShiftSize)
{
    shiftSize = newShiftSize;
}


ShiftSizeCollector::TAdjustmentType
ShiftSizeCollector::getAdjustmentType() const
{
    return adjType;
}

/** returns the strike*/
double ShiftSizeCollector::getStrike() const
{
    return strike;
}

/** returns the original sens control */
DeltaConstSP ShiftSizeCollector::getSensControl() const
{
    return sensCtrl;
}

DRLIB_END_NAMESPACE 
