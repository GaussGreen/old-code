//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergySwapScheduleBase.cpp
//
//   Description : Energy Swap instrument
//
//   Author      : Sean Chen
//
//   Date        : Aug. 29, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/EnergySwapScheduleBase.hpp"
#include "edginc/EnergyStreamSchedule.hpp"

#include <algorithm>
#include <string>
using namespace std;

DRLIB_BEGIN_NAMESPACE


void EnergySwapScheduleBase::validatePop2Object()
{
    static const string method("EnergySwapScheduleBase::validatePop2Object");

      
    // ABSOLUTE
    if (nearbyLabel.size())
    {
        nearbyAbsLabel = EnergyContractLabel(nearbyLabel);
        nearbyRel = 1;
    }
    else
        nearbyAbsLabel = EnergyContractLabel(0,0);

}


double EnergySwapScheduleBase::calculatePV(const EnergyFuturesCurveSP& energyFuturesCurve, 
                                   const YieldCurveSP& yieldCurve, 
                                   double rate,
                                   double notionalAmount, 
                                   const string& notionalType,
                                   const string& dealDerection,
                                   double pastAverage, 
                                   bool pastAverageInclToday)
{
    const string method = "EnergySwapScheduleBase::calculatePV()";

    double value;

    // Fixed leg...
    value = fixedEnergyStreamScheduleSP->getPV(energyFuturesCurve,
                                               yieldCurve, 
                                               rate,
                                               notionalAmount,
                                               notionalType,
                                               dealDerection);

    // Floaing leg...
    value += floatingEnergyStreamScheduleSP->getPV(energyFuturesCurve,
                                                   yieldCurve, 
                                                   notionalAmount, 
                                                   notionalType,
                                                   dealDerection,
                                                   pastAverage, 
                                                   pastAverageInclToday);

    
    return value;

}

DateTime EnergySwapScheduleBase::getValueDate() const
{
    return valueDate;
}

DateTime EnergySwapScheduleBase::getStartDate() const
{
    return startDate;
}

DateTime EnergySwapScheduleBase::getEndDate() const
{
    return endDate;
}

EnergyUnderlyerConstSP EnergySwapScheduleBase::getEnergyUnderlyer() const
{
    return energyUnderlyer.getSP();
}

EnergySwapScheduleBase::EnergySwapScheduleBase(CClassConstSP clazz) : 
      CObject(clazz),
      avgDaysB4End(0),settleDays(5),nearbyRel(1)
{
}

EnergySwapScheduleBase::EnergySwapScheduleBase():
      CObject(TYPE),
      avgDaysB4End(0),settleDays(5),nearbyRel(1)
{    
}

EnergySwapScheduleBase::~EnergySwapScheduleBase()
{    
}

EnergyStreamScheduleFloatingConstSP EnergySwapScheduleBase::getFloatingEnergyStreamSchedule() const
{
	return floatingEnergyStreamScheduleSP;
}

EnergyStreamScheduleFixedConstSP EnergySwapScheduleBase::getFixedEnergyStreamSchedule() const
{
	return fixedEnergyStreamScheduleSP;
}

ObjectArraySP EnergySwapScheduleBase::getFixedLegDetails() const 
{

    static const string method = "EnergySwapScheduleBase::getFixedLegDetails";
  
    StringArraySP paymentDates( new StringArray );
    StringArraySP couponStartDates( new StringArray );
    StringArraySP couponEndDates( new StringArray );
    StringArraySP notionalStartDates( new StringArray );
    StringArraySP notionalEndDates( new StringArray );
    IntArraySP numFixings( new IntArray );
    StringArraySP fixingDates( new StringArray );
    StringArraySP fixingLabels( new StringArray );
    DoubleArraySP fixingRates( new DoubleArray );

    // 6 columns
    ObjectArraySP output = ObjectArraySP(new ObjectArray(6));
       
    fixedEnergyStreamScheduleSP->getScheduleDetails(*paymentDates,
                                         *couponStartDates,
                                         *couponEndDates,
                                         *notionalStartDates,
                                         *notionalEndDates,
                                         *numFixings,
                                         *fixingDates,
                                         *fixingLabels,
                                         *fixingRates);
            
    (*output)[0]=IObjectSP(paymentDates);
    (*output)[1]=IObjectSP(couponStartDates);
    (*output)[2]=IObjectSP(couponEndDates);
    (*output)[3]=IObjectSP(notionalStartDates);
     (*output)[4]=IObjectSP(notionalEndDates);
     (*output)[5]=IObjectSP(fixingRates);
 
     return output;
 
}

ObjectArraySP EnergySwapScheduleBase::getFloatingLegDetails() const 
{

    static const string method = "EnergySwapScheduleBase::getFloatingdLegDetails";
  
    StringArraySP paymentDates( new StringArray );
    StringArraySP couponStartDates( new StringArray );
    StringArraySP couponEndDates( new StringArray );
    StringArraySP notionalStartDates( new StringArray );
    StringArraySP notionalEndDates( new StringArray );
    IntArraySP numFixings( new IntArray );
    StringArraySP fixingDates( new StringArray );
    StringArraySP fixingLabels( new StringArray );
    DoubleArraySP fixingRates( new DoubleArray );

    // 9 columns
    ObjectArraySP output = ObjectArraySP(new ObjectArray(9));
       
    floatingEnergyStreamScheduleSP->getScheduleDetails(*paymentDates,
                                         *couponStartDates,
                                         *couponEndDates,
                                         *notionalStartDates,
                                         *notionalEndDates,
                                         *numFixings,
                                         *fixingDates,
                                         *fixingLabels,
                                         *fixingRates);
            
    (*output)[0]=IObjectSP(paymentDates);
    (*output)[1]=IObjectSP(couponStartDates);
    (*output)[2]=IObjectSP(couponEndDates);
    (*output)[3]=IObjectSP(notionalStartDates);
     (*output)[4]=IObjectSP(notionalEndDates);
     (*output)[5]=IObjectSP(numFixings);
     (*output)[6]=IObjectSP(fixingDates);
     (*output)[7]=IObjectSP(fixingLabels);
    (*output)[8]=IObjectSP(fixingRates);
           
     return output;
 
}


class EnergySwapScheduleBaseHelper
{

public:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); 
        REGISTER(EnergySwapScheduleBase, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSwapScheduleBase);

        FIELD(valueDate,          "Value Date");
        FIELD(energyUnderlyer,    "Underlyer Wrapper");
        FIELD(dealPeriod,         "Deal Period"); // handle to EnergyDealPeriod
    
        FIELD(avgPeriod,       "Averaging Period"); // nD,M,W,Y
        //FIELD_MAKE_OPTIONAL(avgPeriod); // default, whole coupon period
        FIELD(avgFrequency,       "Averaging Frequency"); // nD,M,W,Y
        //FIELD_MAKE_OPTIONAL(avgFrequency); // default, whole coupon period
        FIELD(avgDaysB4End,   "Average Period Convention"); //END, DAY1B4END, DAY2B4END
        FIELD_MAKE_OPTIONAL(avgDaysB4End); // default is 0
        
        FIELD(settleDays,         "Days to Settlement");  // default is 5 days
        FIELD_MAKE_OPTIONAL(settleDays);
        //RELATIVE
        FIELD(nearbyRel,             "Relative Nearby Contract"); // Relative to current contract
        FIELD_MAKE_OPTIONAL(nearbyRel);  // default is 1. front month.
        // ABSOLUTE
        FIELD(nearbyLabel,         "Absolute Nearby Contract"); // overwrite 'nearby', eg, Sep08.
        FIELD_MAKE_OPTIONAL(nearbyLabel);  // default, using nearby

        // START END
        FIELD(startDate,           "Deal Start Date");
        FIELD_MAKE_TRANSIENT(startDate);
        FIELD(endDate,             "Deal End Date");
        FIELD_MAKE_TRANSIENT(endDate);
        FIELD(nearbyAbsLabel,         "Internal Absolute Nearby Contract"); // overwrite 'nearby', eg, Sep08.
        FIELD_MAKE_TRANSIENT(nearbyAbsLabel);  // default, using nearby
        
        FIELD(fixedEnergyStreamScheduleSP,   "fixed stream");
        FIELD_MAKE_TRANSIENT(fixedEnergyStreamScheduleSP);
        FIELD(floatingEnergyStreamScheduleSP,   "floating stream");
        FIELD_MAKE_TRANSIENT(floatingEnergyStreamScheduleSP);


    }

    static IObject* defaultSwapScheduleBase()
    {
        return new EnergySwapScheduleBase();
    }
};

CClassConstSP const EnergySwapScheduleBase::TYPE = CClass::registerClassLoadMethod(
    "EnergySwapScheduleBase", typeid(EnergySwapScheduleBase), EnergySwapScheduleBaseHelper::load);

class GetEnergySwapScheduleDetailsAddin : public CObject 
{

public:

    static CClassConstSP const TYPE;

    // addin parameters
    EnergySwapScheduleBaseSP swapScheduleSP; // a handle to EnergySwapScheduleBase object
    string streamType;


    static IObjectSP getStreamDetails(GetEnergySwapScheduleDetailsAddin *params)
    {
        static const string method = "GetEnergySwapScheduleDetailsAddin::getStreamDetails";
        
        /**
        vector<DateTime> paymentDates;
        vector<DateTime> couponStartDates;
        vector<DateTime> couponEndDates;
        vector<DateTime> notionalStartDates;
        vector<DateTime> notionalEndDates;
        vector<int> numFixings;
        vector<string> fixingDatesString;
        ***/

        StringArraySP paymentDates( new StringArray );
        StringArraySP couponStartDates( new StringArray );
        StringArraySP couponEndDates( new StringArray );
        StringArraySP notionalStartDates( new StringArray );
        StringArraySP notionalEndDates( new StringArray );
        IntArraySP numFixings( new IntArray );
        StringArraySP fixingDates( new StringArray );
        StringArraySP fixingLabels( new StringArray );
        // No rates available now, just to fool the interface
        DoubleArraySP fixingRates( new DoubleArray ); 


        // 9 columns
        ObjectArraySP output = ObjectArraySP(new ObjectArray(9));
       
        if (EnergySwapScheduleBase::TYPE->isInstance(params->swapScheduleSP))
        {
            EnergySwapScheduleBase& aEnergySwapSchedule = dynamic_cast<EnergySwapScheduleBase&>(*params->swapScheduleSP);

            string tmpString = params->streamType;
            transform(tmpString.begin(), tmpString.end(), tmpString.begin(), (int(*)(int))tolower);

            EnergyStreamScheduleSP theStream;
            if (tmpString == "fixed" )
                theStream = aEnergySwapSchedule.fixedEnergyStreamScheduleSP;
            else
                theStream = aEnergySwapSchedule.floatingEnergyStreamScheduleSP;
            
            theStream->getScheduleDetails(*paymentDates,
                                         *couponStartDates,
                                         *couponEndDates,
                                         *notionalStartDates,
                                         *notionalEndDates,
                                         *numFixings,
                                         *fixingDates,
                                         *fixingLabels,
                                         *fixingRates);
            
            (*output)[0]=IObjectSP(paymentDates);
            (*output)[1]=IObjectSP(couponStartDates);
            (*output)[2]=IObjectSP(couponEndDates);
            (*output)[3]=IObjectSP(notionalStartDates);
            (*output)[4]=IObjectSP(notionalEndDates);
            (*output)[5]=IObjectSP(numFixings);
            (*output)[6]=IObjectSP(fixingDates);
            (*output)[7]=IObjectSP(fixingLabels);

            
            return output;
        }
        else
        {
            throw ModelException("Objecy is not an Energy Underlyer");
            
        }  
    }
    
    GetEnergySwapScheduleDetailsAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        REGISTER(GetEnergySwapScheduleDetailsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetEnergySwapScheduleDetailsAddin);
        // order of registration effects order of parameters in addin function
        FIELD(swapScheduleSP, "Handle to EnergySwapSchedule object");
        FIELD(streamType, "FIXED or FLOAING");
        
        Addin::registerClassObjectMethod("GET_ENERGY_SWAP_SCHEDULE_DETAILS",
                                         Addin::MARKET,
                                         "Get energy swap stream details",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getStreamDetails);
    }
    
    static IObject* defaultGetEnergySwapScheduleDetailsAddin()    
    {
        return new GetEnergySwapScheduleDetailsAddin();
    }
};

CClassConstSP const GetEnergySwapScheduleDetailsAddin::TYPE = CClass::registerClassLoadMethod(
    "GetEnergySwapScheduleDetailsAddin", typeid(GetEnergySwapScheduleDetailsAddin), GetEnergySwapScheduleDetailsAddin::load);


bool  EnergySwapScheduleBaseLoad() { return (EnergySwapScheduleBase::TYPE != 0);   }

DRLIB_END_NAMESPACE

