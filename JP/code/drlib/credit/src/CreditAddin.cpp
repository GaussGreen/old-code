//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditAddin.cpp
//
//   Description : Credit addin utilities
//
//   Author      : Ning Shen
//
//   Date        : 26 Aug 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditEngine.hpp"

DRLIB_BEGIN_NAMESPACE

////////////////////////////////////////////////
/** net two path value grids */
class NetPathValues : public CObject
{
public:
    CreditPathValuesOutArray   pathValues;
	CMarketDataSP		 market;            /** market cache $unregistered */ 
    MarginAcctSP         marginAcct;        /** margin collateral account */
	
    double               percentile;        /** percentile used for peak computation */
    string               creditCcyName;      /** name of currency we wish to compute exposure in */

	CreditSpreadCurveSP	 creditSpreads;		/** optional field: required if cvr is to be computed. */
	bool				 outputGrids;		/** if true, returns all values on all paths */
	
	
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz){
        REGISTER(NetPathValues, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultNetPathValues);
        clazz->setPublic(); 
        FIELD(pathValues, "path value arrays");
		FIELD(marginAcct,"margin collateral account");
        FIELD(percentile, "percentile used for peak computation");
        FIELD(creditCcyName, "name of currency we wish to compute exposure in");
		FIELD(creditSpreads, "optional field: required if cvr is to be computed.");
		FIELD(outputGrids, "if true, returns all values on all paths");

        Addin::registerClassObjectMethod("NET_PATH_VALUES",
                                         Addin::RISK, // CREDIT ??
                                         "net array of path value grids, computing aggregate peak and average",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)netPathValues);
    }
	
    static IObject* defaultNetPathValues(){
        return new NetPathValues();
    }

    NetPathValues() : CObject(TYPE){};

    static IObjectSP netPathValues(NetPathValues* data)
    {
        if (data->pathValues.size() == 0)
            return IObjectSP(   );

        CreditPathValuesOutSP result(CreditPathValuesOutSP::dynamicCast((IObjectSP)(data->pathValues[0]->clone())));
        if (data->pathValues.size() == 1)
            return result;

        for (int i=0; i<data->pathValues.size(); i++)
            result->add(data->pathValues[i].get());

		if (!!result)
		{
			result->calcPeakAverageStdAndDRE(data->market, 
									   data->marginAcct, 
									   data->percentile, 
									   data->creditCcyName,
									   data->creditSpreads);
		}

        return result;
    }
};

CClassConstSP const NetPathValues::TYPE = CClass::registerClassLoadMethod("NetPathValues", typeid(NetPathValues), load);

#if 0 /** don't need anymore */
////////////////////////////////////////////////

////////////////////////////////////////////////
/** calc peak and average exposures */
class PeakAndAverage : public CObject
{
public:
    CreditPathValuesOutSP   pathValues;
    double                  peakPercentile; // 0.975 default

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz){
        REGISTER(PeakAndAverage, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPeakAndAverage);
        clazz->setPublic(); 
        FIELD(pathValues, "path values");
        FIELD(peakPercentile, "percentile for peak exposure");
        FIELD_MAKE_OPTIONAL(peakPercentile);

        Addin::registerClassObjectMethod("PEAK_AND_AVERAGE",
                                         Addin::RISK, // CREDIT ??
                                         "calc peak and average exposures",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)peakAverage);
    }
	
    static IObject* defaultPeakAndAverage(){
        return new PeakAndAverage();
    }

    PeakAndAverage() : CObject(TYPE), peakPercentile(0.975) {};

    static IObjectSP peakAverage(PeakAndAverage* data)
    {
        CreditPathValuesOutSP result(CreditPathValuesOutSP::dynamicCast((IObjectSP)(data->pathValues->clone())));
        result->calcPeakAndAverage(data->peakPercentile);
	    return result;
    }
};
CClassConstSP const PeakAndAverage::TYPE = CClass::registerClassLoadMethod("PeakAndAverage", typeid(PeakAndAverage), load);
#endif

////////////////////////////////////////////////

////////////////////////////////////////////////
/** load addin types */
bool CreditAddinLoad()
{
    bool success =
        NetPathValues::TYPE && true;
//        && PeakAndAverage::TYPE;


    return success;

}


DRLIB_END_NAMESPACE

