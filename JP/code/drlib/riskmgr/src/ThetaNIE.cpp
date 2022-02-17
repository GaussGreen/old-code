//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ThetaNIE.cpp
//
//   Description : Theta Net Interested Earned
//
//   Author      : Stephen Hope
//
//   Date        : 8 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/ThetaNIE.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Results.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for theta */
const string ThetaNIE::NAME = "THETA_NIE";


/** for reflection */
ThetaNIE::ThetaNIE(): Theta(TYPE, NAME), 
    timeOffsetNIE(DateTime::END_OF_DAY_TIME){}

ThetaNIE::ThetaNIE(int offset, HolidaySP hols):
    Theta(offset, hols, TYPE, NAME),
    timeOffsetNIE(DateTime::END_OF_DAY_TIME){}

/** If timeOffsetNIE = END_OF_DAY_TIME return same date at EOD,
    If timeOffsetNIE = BEFORE_EX_DIV_TIME return date + offset @ BEX */
DateTime ThetaNIE::rollDate(const DateTime &date) const
{
    DateTime rollDate;

    if (timeOffsetNIE == DateTime::END_OF_DAY_TIME)
    {
        // just roll to the end of today
        rollDate = DateTime(date.getDate(), DateTime::END_OF_DAY_TIME);
    }
    else if (timeOffsetNIE == DateTime::BEFORE_EX_DIV_TIME)
    {
        rollDate = hols->addBusinessDays(date, offset);
        // adjust the time
        rollDate = DateTime(rollDate.getDate(), DateTime::BEFORE_EX_DIV_TIME);
    }

    return rollDate;
}

/** calculates given sensitivity - invoked by calculateSens */
void ThetaNIE::calculate(TweakGroup*  tweakGroup,
                         CResults*    results)
{
    OutputNameConstSP outputName(new OutputName(getSensOutputName()));
    try
    {
        if (!results->exists(Results::INSTRUMENT_PACKET, outputName)) 
        {
            double divisor = this->divisor();
            timeOffsetNIE = DateTime::END_OF_DAY_TIME;
            double origPrice = getSensPrice(results, tweakGroup->getInstrument(), tweakGroup->getModel(), getControl());

            double shiftedPriceEOD = shiftAndPrice(tweakGroup, origPrice);

            timeOffsetNIE = DateTime::BEFORE_EX_DIV_TIME;

            double shiftedPricePSOD = shiftAndPrice(tweakGroup, origPrice);

            double thetaNIE = (shiftedPricePSOD - shiftedPriceEOD)/divisor;
            results->storeScalarGreek(thetaNIE,
                                      Results::INSTRUMENT_PACKET, 
                                      outputName);
        }
    }
    catch (exception& e)
    {
        results->storeGreek(IObjectSP(new Untweakable(e)),
                            Results::INSTRUMENT_PACKET, 
                            outputName);
    }
}

/** Factory class dictates what methods of building this
    sensitivity are supported. Here we support a default */
class ThetaNIEFactory: public SensitivityFactory::IDefault {
public:
    virtual Sensitivity* createDefault(){
        HolidaySP hols(Holiday::weekendsOnly());
        return new ThetaNIE(Theta::DEFAULT_SHIFT, hols);
    }
};

/** Invoked when ThetaNIE is 'loaded' */
void ThetaNIE::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ThetaNIE, clazz);
    SUPERCLASS(Theta);
    EMPTY_SHELL_METHOD(defaultThetaNIE);
    FIELD(timeOffsetNIE, "ThetaNIE offset");
    FIELD_MAKE_TRANSIENT(timeOffsetNIE);
    // register how to build our sensitivity
    SensitivityFactory::addSens(ThetaNIE::NAME, 
                                new ThetaNIEFactory(), 
                                new ThetaNIE(),
                                ThetaNIE::Shift::TYPE);
}

IObject* ThetaNIE::defaultThetaNIE(){
    return new ThetaNIE();
}

CClassConstSP const ThetaNIE::TYPE = CClass::registerClassLoadMethod(
    "ThetaNIE", typeid(ThetaNIE), load);


DRLIB_END_NAMESPACE
