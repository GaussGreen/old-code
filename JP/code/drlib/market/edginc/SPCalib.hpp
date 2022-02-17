//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPCalib.hpp
//
//   Description : smile params and flat vol for basis spread
//
//   Date        : 
//
//----------------------------------------------------------------------------

#ifndef SPCALIB_HPP
#define SPCALIB_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/VolRequest.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/MRSpotVolProcessed.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class BootstrappedBasisIndexCurve;

/** SP Vol Parameters consisting of only qSpread, FlatSpotVol, meanRev */
class MARKET_DLL SPCalib: 
    public MarketObject {
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const { return name; }

    class Processed;
    friend class Processed;

    class MARKET_DLL Processed : public MRSpotVolProcessed {
    public:
        static CClassConstSP const TYPE;

        virtual string getName() const { return vol->name; }
        double getQSpread () const { return vol->qSpread; }
        double getFlatVol () const { return vol->flatSpotVol; }
        virtual double meanReversion() const { return vol->meanRev; }

        /** Returns spot vols between pairs of dates ie between
        [initialStartDate, subsequentDates[0]], 
        [subsequentDates[0], subsequentDates[1], ... */
        virtual void spotVol(const DateTime&     initialStartDate,
                            const DateTimeArray& subsequentDates,
                            DoubleArray&         spotvol) const;

        /** calculates the trading time between two dates */
        virtual double calcTradingTime(const DateTime &date1, 
                                        const DateTime &date2) const;

        /** retrieve time measure for the vol - no time metric here */
        virtual TimeMetricConstSP GetTimeMetric() const;

        /** Returns the dates (possible none!) that the spot vol is
            defined on */
        virtual DateTimeArraySP getSpotVolDates() const {
            // no dates!
            return DateTimeArraySP(new DateTimeArray());
        }

        Processed(const SPCalib* vol) :  
            MRSpotVolProcessed(TYPE), vol(vol) {}

    private:
        static IObject* defaultConstructor(){ return new Processed(); }
        static void load(CClassSP& clazz);
        Processed() : MRSpotVolProcessed(TYPE) {}

        smartConstPtr<SPCalib> vol;        
    };


    static bool recognizedVolRequest(const CVolRequest* volRequest);

    /** Returns this object containing smile parameters */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;

    static void load(CClassSP& clazz);

    /* How to ask for the spot diffusion parameters */
    class MARKET_DLL SPCalibRequest: public CVolRequest {
    public:
        static CClassConstSP const TYPE;
        SPCalibRequest() : CVolRequest(TYPE) {}
    private:
        static IObject* defaultConstructor(){ return new SPCalibRequest(); }
        static void load(CClassSP& clazz);
    };


private:
    static IObject* defaultConstructor(){
        return new SPCalib();
    }

    SPCalib();
    SPCalib(const SPCalib& rhs);
    SPCalib& operator=(const SPCalib& rhs);
    // fields 
    string name;    // name of this market data object
    double qSpread;
    double flatSpotVol;
    double meanRev;
};

DECLARE(SPCalib);

// support for wrapper class
typedef MarketWrapper<SPCalib> SPCalibWrapper;

DECLARE(SPCalibWrapper);


DRLIB_END_NAMESPACE

#endif
