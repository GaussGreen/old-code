//---------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : ResetSchedule.hpp
//
//   Description : A reset schedule class
//
//   Author      : André Segger
//
//   Date        : 07 January 2003
//
//
//----------------------------------------------------------------------------

#ifndef RESETSCHEDULE_HPP
#define RESETSCHEDULE_HPP

#include <string>
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** A schedule of dates & levels you can interpolate on */

class MARKET_DLL ResetSchedule : public CObject {
public:
    static CClassConstSP const TYPE;

    ResetSchedule(const DateTime&      valueDate,
                  const DateTimeArray& dates,
                  const DoubleArray&   minResetStrikes,
                  const DoubleArray&   maxResetStrikes);

    virtual ~ResetSchedule();

    /** how long is the schedule ? */
    int length() const;

    /** validation */
    virtual void validatePop2Object();

    /** return date list (deep copy) */
    DateTimeArray getDates() const;

    /** returns reference to date array */
    const DateTimeArray &getDateArray() const;

    /** return min reset strike list (deep copy) */
    DoubleArray getMinResetStrikes() const;

    /** return max reset strike list (deep copy) */
    DoubleArray getMaxResetStrikes() const;

    /** return the parity at which the bond gets reset (deep copy) */
    DoubleArray getParity() const;

    /** get reset information for a particular date */
    bool hasReset(const DateTime& date) const;

    // Get conversion ratio following reset based on spot and reset parameters
    // faceValue, maxStrike, minStrike, spot quantities passed must be in consistent ccy (payoff or underlying)
    // Ratio = FaceValue / (spot * resetPremium) subject to minStrike<=spot*resetPremium<=maxStrike
    static double getResetRatio(double faceValue,    
                                double maxStrike,    // Max conversion price following reset 
                                double minStrike,    // Min conversion price following reset 
                                double resetPremium,   // Multiple of spot to give new conversion price after reset (otherwise known as parity)
                                double spot);

    /** get current conversion ratio based on the initial conv ratio and a history of fixings */
    double getCurrentConversionRatio(const double&   initialConvRatio,
                                     const DateTime& valueDate,
                                     const double&   faceValue);

    double getMaxResetStrike(const DateTime& date) const;

    double getMinResetStrike(const DateTime& date) const;

    double getParity(const DateTime& date) const;

    // string getResetType() const { return resetType;}

    string getResetType(const DateTime& date);

    bool getReceiveIntrinsic(const DateTime& date);

    bool isFlat(const DateTime& date);

    CStringArraySP getResetTypes() const { return resetTypes;}

    void rollDate(const DateTime& oldValueDate, const DateTime& newValueDate, 
                  const double newSpot,         const bool isCcyStruck, const double newFX);

    void setValueDate(const DateTime& valueDate);

    void preProcessSchedule(const bool isCcyStruck, const double& fxSpot, DoubleArraySP fwdFXs, 
                            const double& initialConvPrice);

    void scaleLevels(const double scaleFactor);

    bool isFwdStart();

    bool isStarted(DateTime valDate);

    DateTime getStartDate(DateTime valDate);

    // check it's already started or not.  If it's started, return "altLvl".
    double getStartLevel(DateTime valDate, double altLvl);


private:
    friend class ResetScheduleHelper;

    ResetSchedule();
    ResetSchedule(const ResetSchedule &rhs);
    ResetSchedule& operator=(const ResetSchedule& rhs);

    DateTimeArraySP resetDates;
    DoubleArraySP   minResetStrike;
    DoubleArraySP   maxResetStrike;
    DoubleArraySP   parity;
    DoubleArraySP   resetLevel;
    DoubleArraySP   resetFX;
    CStringArraySP  resetTypes;
    string          resetType;
    DateTime        valueDate;

    bool            fwdStarting;
    CashFlowArray   initLevel;

    BoolArraySP     receiveIntrinsic;

};

typedef smartConstPtr<ResetSchedule> ResetScheduleConstSP;
typedef smartPtr<ResetSchedule> ResetScheduleSP;

DRLIB_END_NAMESPACE
#endif


