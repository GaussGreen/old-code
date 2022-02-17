//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : MQQuasiIRVol.cpp
//
//   Description : Base class of IRVol for MultiQ 
//
//   Author      : Keith Jia
//
//   Date        : June 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/MQQuasiIRVol.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/MultiQDistribution.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/Maths.hpp"
#include "edginc/TimeMetric.hpp"


DRLIB_BEGIN_NAMESPACE

/* Make sure the class links */
bool  MQQuasiIRVolLoad() {
    return (MQQuasiIRVol::TYPE != 0);
}

//-----------------------------------------------------
// MarketObject methods
//-----------------------------------------------------
string MQQuasiIRVol::getName() const
{
    return name;
}

void MQQuasiIRVol::getMarket(const IModel* model, const MarketData *market)
{
    static const string method = "MQQuasiIRVol::getMarket";
    
    try
    {
        market->GetReferenceDate(valueDate);
        hols.getData(model, market);
        //hols.setObject(MarketObjectSP(Holiday::weekendsOnly()));

        MaturityPeriod period("1D");
        Holiday* holiday = hols.get();

        if(!metric.get())
        {
            metric = TimeMetricSP(new TimeMetric(1., holiday));
        }


        //spotDate = bdc->adjust(valueDate, holiday);

        
        spotDate = DateFwdThenAdjust(valueDate,
                                     period,
                                     spotOffset,
                                     *bdc,
                                     *holiday);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

//-----------------------------------------------------
// MQQuasiIRVol methods
//-----------------------------------------------------
MQQuasiIRVol::MQQuasiIRVol(CClassConstSP clazz): IRVolCommon(clazz)
{}


MQQuasiIRVol::~MQQuasiIRVol()
{}

TimeMetricConstSP MQQuasiIRVol::GetTimeMetric() const
{
    return TimeMetricConstSP(metric);
}

HolidayConstSP MQQuasiIRVol::getHoliday() const
{
    return hols.getSP();
}

/** interp vols to get the data required to create the mq, and create mq.
 */ 
MultiQDistribution* MQQuasiIRVol::CreateMQ(ExpiryConstSP inputTenor,
                                           double parFwdRate, 
                                           const DateTime& inputResetDate) const
{
    static const string method = "MQQuasiIRVol::CreateMQ";

    try{
        tenor = inputTenor;
        //interp vols, prepare for MQ creation
        //get atmBSVol, inputVols, inputStrikes, optionTimeToExpiry
        int idx = 0;
        int jdx = 0;

        //tenorDate = reset + tenor
        //atmTenorDates = reset + tenors
        resetDate = inputResetDate;
        DateTime tenorDate = bdc->adjust(tenor->toDate(resetDate), hols.get());
        DateTimeArraySP atmTenorDates = DateTimeArraySP(new DateTimeArray(atmTenors->size()));

        for(int i=0; i<atmTenors->size(); i++){
            (*atmTenorDates)[i] = bdc->adjust((*atmTenors)[i]->toDate(resetDate), hols.get());
        }

        //find tenor idx location
        if(tenorDate > atmTenorDates->front())
        {
            idx = tenorDate.findUpper(*atmTenorDates) - 1;
        }

        //find reset jdx location
        if(resetDate > atmMatDates->front())
        {
            jdx = tenorDate.findUpper(*atmMatDates) - 1;
        }

        //---------- MQ parameters ------------
        //double optionTimetoExpiry;  //yearfrac from valueDate to resetDate 
        //double atmBSVol;            //interpolated from market data atmVols
        //doubleArraySP inputStrikes; //interpolated from strikeBoundaries
        //doubleArraySP inputVols;    //interpolated from impliedVols

        double optionTimeToExpiry = metric->yearFrac(valueDate, resetDate);

        //linear interp -- atmBSVol
        double x = metric->yearFrac( (*atmTenorDates)[idx], tenorDate)/
            metric->yearFrac( (*atmTenorDates)[idx], (*atmTenorDates)[idx+1] );
        double y = metric->yearFrac( (*atmMatDates)[jdx], resetDate )/
            metric->yearFrac( (*atmMatDates)[jdx], (*atmMatDates)[jdx+1] );

        atmBSVol = (*atmVols)[idx][jdx] * (1-x)*(1-y) + (*atmVols)[idx+1][jdx] * x*(1-y)
            + (*atmVols)[idx][jdx+1] * (1-x)*y + (*atmVols)[idx+1][jdx+1] * x*y;


        //otm strikes and inputVols

        //find the tenor.  Exception will be thrown if not found.
        int k = tenor->search(otmTenors.get());
        
        int numCols = strikeBoundaries->size();
        DoubleArraySP inputStrikes = DoubleArraySP(new DoubleArray(*strikeBoundaries));
        if(isInputQs)
        {
            //MultiQDistribution makes it this way
            inputStrikes->erase(inputStrikes->end()-1);
        }

        DoubleArraySP inputVols = DoubleArraySP(new DoubleArray(numCols));
        CDoubleMatrixSP volSurface = (*impliedVols)[k];

        jdx = 0;
        if(resetDate > otmMatDates->front())
        {
            jdx = tenorDate.findUpper(*otmMatDates) - 1;
        }
        
        y = metric->yearFrac((*otmMatDates)[jdx], resetDate)/
            metric->yearFrac((*otmMatDates)[jdx], (*otmMatDates)[jdx+1]);

        for(int i = 0; i < numCols; i++)
        {
            (*inputVols)[i] = (*volSurface)[i][jdx] * (1-y) 
                + (*volSurface)[i][jdx+1] * y;
        }


        //crate MQ
        string s = "bsImpliedVol";
        if(isInputQs) {
            s = "explicitQ";
        }

        MultiQDistribution *mymq = 
            new MultiQDistribution(optionTimeToExpiry,
                                   parFwdRate,
                                   atmBSVol,
                                   *inputVols,
                                   *inputStrikes,
                                   s);

        return (mymq);
    } catch(exception& e) {
            throw ModelException(e, method);
    }
}

void MQQuasiIRVol::validatePop2Object()
{
    static const string method = "MQQuasiIRVol::validatePop2Object";

    try
    {
        dcc.reset(DayCountConventionFactory::make(swapDCC));

        bdc.reset(BadDayConventionFactory::make(swapBDC));

        if(atmTenors->empty()) {
            throw ModelException(method, name + " atmTenors are empty");
        }

        if(atmMatDates->empty()) {
            throw ModelException(method, name + " atmMaturities are empty");
        }

        if(atmVols->empty()) {
            throw ModelException(method, name + " atmVols are empty");
        }

        if(otmTenors->empty()) {
            throw ModelException(method, name + " otmTenors are empty");
        }

        if(otmMatDates->empty()) {
            throw ModelException(method, name + " otmMaturities are empty");
        }

        if(impliedVols->empty()) {
            throw ModelException(method, name + " impliedVols are empty");
        }

        if(strikeBoundaries->empty()) 
        {
            if(isInputQs)
            {
                /** Since currently for 2Q, Aladdin implicitly assumes the deltas 
                    are simply 0, 0.5, 1.  No deltas are explicitly supplied.  In this case, 
                    strikeBoundaries may be empty.  Now check if it is 2Q.
                **/
                CDoubleMatrixSP volSurface = (*impliedVols)[0];
                int surfaceCols = volSurface->numCols();
                if(surfaceCols == 2)
                {
                    strikeBoundaries = DoubleArraySP(new DoubleArray(2));
                    (*strikeBoundaries)[0] = 0.5;
                    (*strikeBoundaries)[1] = 1;
                } else {
                    throw ModelException(method, name + " strikeBoundaries are empty");
                }
            } else 
            {
                throw ModelException(method, name + " strikeBoundaries are empty");
            }
        }
        

        //atmVols
        if(atmVols->numCols() != atmTenors->size()) {
            throw ModelException(method, name + " atmVols has different numCols " 
                                 + Format::toString(atmVols->numCols()) + " than atmTenor size "
                                 + Format::toString(atmTenors->size()));
        }

        if(atmVols->numRows() != atmMatDates->size()) {
            throw ModelException(method, name + " atmVols has different numRows " 
                                 + Format::toString(atmVols->numRows()) + " than atmMaturities size "
                                 + Format::toString(atmMatDates->size()));
        }

        //otmVols
        if(otmTenors->size() != impliedVols->size()) {
            throw ModelException(method, name + " impliedVols has different z-axis size " 
                                 + Format::toString(impliedVols->size()) + " than tenors size "
                                 + Format::toString(otmTenors->size()));
        }

        int otmMatSize = otmMatDates->size();
        int strikeSize = strikeBoundaries->size();

        for(int i = 0; i < impliedVols->size(); i++)
        {
            CDoubleMatrixSP volSurface = (*impliedVols)[i];
            int surfaceRows = volSurface->numRows();
            int surfaceCols = volSurface->numCols();
            if(otmMatSize != surfaceRows) {
                throw ModelException
                    (method, name + Format::toString(i) + "th impliedVols slice has different numRows " 
                     + Format::toString(surfaceRows) + " than otmMatDates size "
                     + Format::toString(otmMatSize));
            }
            
            if(strikeSize != surfaceCols) {
                throw ModelException
                    (method, name + Format::toString(i) + "th impliedVols slice has different numCols " 
                     + Format::toString(surfaceCols) + " than strikeBoundaries "
                     + Format::toString(strikeSize));
            }
        }

    } catch(exception& e) {
            throw ModelException(e, method);
    }
}

//-----------------------------------------------------
// Invoked when class is loaded
//-----------------------------------------------------

void MQQuasiIRVol::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(MQQuasiIRVol, clazz);
    SUPERCLASS(IRVolCommon);
    FIELD(name, "name");
    FIELD(valueDate, "today");
    FIELD(spotOffset, "offset days from base to spot");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(metric,"");
    FIELD_MAKE_OPTIONAL(metric);
    FIELD(atmTenors, "At the money vols tenors");
    FIELD(atmMatDates, "At the money vols maturities");
    FIELD(atmVols, "At the money vols");
    FIELD(otmTenors, "At the money vols tenors. z-axis for implieVols cube");
    FIELD(otmMatDates, "Off the money strike maturities. y-axis for impliedVols cube");
    FIELD(strikeBoundaries, "strike ratios or deltas. x-axis for impliedVols cube");
    FIELD(impliedVols, "vols or q's cube");
    FIELD(isInputQs, "true if deltas and Qs are inputs");  
    FIELD(swapDCC, "");
    FIELD(swapBDC, "");
    FIELD(hols, "Holidays");

    //FIELD_NO_DESC(mq);
    //FIELD_MAKE_TRANSIENT(mq);
    FIELD_NO_DESC(dcc);
    FIELD_MAKE_TRANSIENT(dcc);
    FIELD_NO_DESC(bdc);
    FIELD_MAKE_TRANSIENT(bdc);

}

//-----------------------------------------------------
// Static variables
//-----------------------------------------------------

CClassConstSP const MQQuasiIRVol::TYPE = CClass::registerClassLoadMethod(
    "MQQuasiIRVol", typeid(MQQuasiIRVol), load);

DEFINE_TEMPLATE_TYPE(MQQuasiIRVolWrapper);

DRLIB_END_NAMESPACE

