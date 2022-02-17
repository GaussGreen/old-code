//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Extendible.cpp
//
//   Description : Extendible Instrument (generic 1 factor)
//
//   Author      : Stephen Hope
//
//   Date        : 03 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/PayStream.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Average.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/Maths.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Future.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/CashSettleDate.hpp"

DRLIB_BEGIN_NAMESPACE

class AverageSpot;


class Extendible: public Generic1Factor,
                  virtual public CClosedFormLN::IIntoProduct,
                  virtual public LastSensDate{
public:
    static CClassConstSP const TYPE;  
    static const double MINIMUM_SPREAD;
    

    virtual void Validate(){
        static const string method = "Extendible::Validate";
        
        // Call Generic1Factor validate()
        validate();

        /* ensure that we are pricing for one contract as fixed 
           notional does not make any sense for this model */
        if (!oneContract)
        {
            throw ModelException(method,
                                 "Extendible cannot be priced for fixed notional");
        }

        DateTime lodate;
        DateTime equityMaturity;

        avgOut->getBoundingDates(lodate, equityMaturity);


        // fwdStarting options cannot start after the first avgOut date
        if (fwdStarting && startDate.isGreater(lodate)) {
                throw ModelException(method,
                                     "first average out date (" + 
                                     lodate.toString() + 
                                     ") must be >= start date (" + 
                                     startDate.toString() + 
                                     ") for a forward starting option");
        }

        // check that the equity maturity is before or on the first libor refix date
        if (extensionSchedule->length() > 0)
        {
            if (equityMaturity.isGreater(extensionSchedule->firstDate()))
            {
                throw ModelException(method, 
                                     "last average out date ("+
                                     equityMaturity.toString()+
                                     ") is not before the first libor refix date ("+
                                     extensionSchedule->firstDate().toString()+")");
            }
        }
        else
        {
            throw ModelException(method,
                                 "Extension schedule is empty !");
        }


        // ensure lowKstrike is not positive
        if (Maths::isPositive(lowKspread))
        {
            throw ModelException(method,
                                 "lowKspread ("+ 
                                 Format::toString(lowKspread)+
                                 ") must be 0 or -ve");
        }
        
        // ensure highKstrike is not negative
        if (Maths::isNegative(highKspread))
        {
            throw ModelException(method,
                                 "highKspread ("+ 
                                 Format::toString(highKspread)+
                                 ") must be 0 or +ve");
        }
        
        if ((highKspread - lowKspread) < MINIMUM_SPREAD)
        {
            ErrorHandler::writeMsg("The difference between highKspread and "
                                   "lowKspread must be >= to the DR recommended "
                                   "minimum spread! ("+
                                   Format::toString(MINIMUM_SPREAD)+
                                   "). highKspread will default to ("+
                                   Format::toString(MINIMUM_SPREAD/2)+
                                   ") and lowKspread will default to ("+
                                   Format::toString(-MINIMUM_SPREAD/2)+")");
            lowKspread = -MINIMUM_SPREAD/2;
            highKspread = MINIMUM_SPREAD/2;
        }
        
        // Check that the extension values are descending
        DoubleArray extStrikes = extensionSchedule->getValues();
        for (int i = 0; i < extensionSchedule->length()-1; i++)
        {
            if (extStrikes[i+1] >= extStrikes[i])
            {
                throw ModelException(method,
                                     "extension schedule values must be in descending order!");
            } 
        }

        // check that the last extensionSchedule value is zero
        if (!(Maths::isZero(extensionSchedule->lastValue())))
        {
            throw ModelException(method, 
                                 "The last value in the extensionSchedule "
                                 "must = zero.");
        }

        AssetUtil::assetCrossValidate(asset.get(),
                                      fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);
        
        // check that underlying is not a future
        if ( Future::TYPE->isInstance(asset.get()))
        {
            throw ModelException(method,
                                 "Average options on Futures are not allowed");
        }

        // check that underlying is not an fx asset
        if ( FXAsset::TYPE->isInstance(asset.get()))
        {
            throw ModelException(method,
                                 "Options on FX assets are not allowed yet");
        }
        
        // check that there are at least 2 extension dates
        if (extensionSchedule->length() < 2)
        {
            throw ModelException(method,
                                 "There must be at least two dates in the extension schedule.");
        }

        // check that the floating rate reset interval is not greater than the extension level intervals
        DateTimeArray extDates = extensionSchedule->getDates();
        DateTime firstExtDate = extDates[0];
        MaturityPeriod resetInterval(liborReset);
        if ((resetInterval.toDate(firstExtDate)).isGreater(extDates[1]))
        {
            throw ModelException(method,
                                 "The floating rate reset interval must not be greater than the "
                                 "interval between extension dates.");
        }
        
         /* While we're here create an instrument settlement that settles
            on the maturity date (=last avgOut date) since this is the only
            settlement that makes sense for this product */
        settle = InstrumentSettlementSP(new CashSettleDate(equityMaturity));
        
        // generate a single libor pay stream
        generatePayStream();
    }

private:
    friend class ExtendibleHelper;
    friend class ExtendibleClosedForm;

    Extendible():Generic1Factor(TYPE),lowKspread(0.0), highKspread(0.0),
       spread(0.0){}; 
    Extendible(const Extendible& rhs);
    Extendible& operator=(const Extendible& rhs);

    /** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** Generate the payStream object */
    void generatePayStream(){
        static const string method = "Extendible::generatePayStream";
        
        try
        {

            /* Validation has already checked that the interval between dates
               in the extension schedule is greater than the floatRate reset interval */
            int periodsPerExtIntvl = 1;
            
            DateTimeArray extDates = extensionSchedule->getDates();
            DateTime refixDate = extDates[0];
            DateTime payDate = extDates[1];
            MaturityPeriod resetIntvl(liborReset);
        
            while ((refixDate = resetIntvl.toDate(refixDate)).isLess(payDate))
            {
                periodsPerExtIntvl++;
            }
            
            // can now generate the refixDates
            DateTimeArraySP refixDates(new DateTimeArray(0));
            int numExtDates = extensionSchedule->length();
            int i;
            for (i = 0; i < numExtDates-1; i++)
            {
                refixDates->push_back(extDates[i]);
                for (int j = 0; j < periodsPerExtIntvl-1; j++)
                {
                    refixDates->push_back(resetIntvl.toDate((*refixDates)[refixDates->size()-1]));
                }
            }
            
            // generate the accrualDates
            DateTimeArraySP accrualDates(new DateTimeArray(refixDates->size()+1));
            for (i = 0; i < refixDates->size(); i++)
            {
                (*accrualDates)[i] = (*refixDates)[i];
            }
            (*accrualDates)[accrualDates->size()-1] = extDates[numExtDates-1];
            
            // generate the payment dates
            DateTimeArraySP payDates(new DateTimeArray(accrualDates->size()-1));
            for (i = 1; i < accrualDates->size(); i++)
            {
                (*payDates)[i-1] = (*accrualDates)[i];
            }

            /* generate the compound bools. Assume every pay date is
               a 'real' pay date. */
            CBoolArraySP compound(new BoolArray(refixDates->size()));
            for (i = 0; i < refixDates->size(); i++)
            {
                (*compound)[i] = false;
            }

            /* generate the notionals (all equal to the G1F notional)
               This is used as the Libor notional (the underlying is oneContract remeber) */
            DoubleArraySP notionals(new DoubleArray(refixDates->size()));
            for (i = 0; i < refixDates->size(); i++)
            {
                (*notionals)[i] = notional;
            }
            
            // generate the refixLevels (all equal zero)
            DoubleArraySP refixLevels(new DoubleArray(refixDates->size()));
            
            // generate the spreads (all equal)
            DoubleArraySP spreads(new DoubleArray(refixDates->size()));
            for (i = 0; i < refixDates->size(); i++)
            {
                (*spreads)[i] = spread;
            }
            
            // generate the weights (all one)
            DoubleArraySP weights(new DoubleArray(refixDates->size()));
            for (i = 0; i < refixDates->size(); i++)
            {
                (*weights)[i] = 1.0;
            }
            
            // Now generate the FloatRate object
            string liborBDC = "None";
            DayCountConventionConstSP rateDayCountConv(DayCountConventionFactory::make(liborDCC));
            BadDayConventionConstSP rateBadDayConv(BadDayConventionFactory::make(liborBDC));
            MaturityPeriodSP resetInterval(new MaturityPeriod(liborReset));
            // Extendible only allows resetInterval = matInterval
            MaturityPeriodSP matInterval(new MaturityPeriod(liborReset));
            
            HolidaySP liborHols(Holiday::noHolidays());
            
            FloatRateSP floatRate(new FloatRate(discount.get(),
                                                liborHols.get(),
                                                rateDayCountConv.get(),
                                                rateBadDayConv.get(),
                                                resetInterval.get(),
                                                matInterval.get(),
                                                false,  // not multiPeriod
                                                0,
                                                0));
            

            DayCountConventionConstSP accrualDCC(DayCountConventionFactory::make(liborDCC));
            
            // construnct the PayStream object
            payStream = PayStreamSP(new PayStream(accrualDCC.get(),
                                                  false, // compoundFlat
                                                  floatRate.get(),
                                                  (*notionals),
                                                  (*refixDates),
                                                  (*accrualDates),
                                                  (*payDates),
                                                  (*refixLevels),
                                                  (*spreads),
                                                  (*weights),
                                                  (*compound)));
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    double cashFlowPV(const DateTime& date,
                      const CashFlowArray& cfl)const{
        static const string method = "Extendible::cashFlowPV";
        double value = 0.0;
        try
        {
            for (int idx = 0; idx < cfl.size(); idx++)
            {
                DateTime cfDate = cfl[idx].date;
                if (cfDate.isGreater(date))
                {
                    // pv the cash flow
                    value += cfl[idx].amount * discount->pv(date,
                                                            cfDate);
                }
            }
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
        
        return value;
    }

    double priceSpread(double& lowStrike,
                       double& highStrike,
                       const DateTime& maturityDate)const{
        static const string method = "Extendible::priceSpread";
        double spreadValue = 0.0;

        try
        {
            spreadValue = Average::priceSpotSpread(valueDate,
                                                   startDate,
                                                   maturityDate,
                                                   fwdStarting,
                                                   false,  // put
                                                   true,   // oneContract
                                                   notional, 
                                                   initialSpot,
                                                   lowStrike,
                                                   highStrike,
                                                   settle.get(),
                                                   asset.get(),
                                                   discount.get(),
                                                   avgOut.get());
            
            return spreadValue;
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
        
    }
    
    /** Indicates whether VEGA_MATRIX is sensible for this instrument.*/
    bool avoidVegaMatrix(const IModel*model){
        return false;
    }

    /** Returns all strikes the Extendible is sensitive to  */
    DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                      const IModel*      model){
        static const string method = "Extendible::getSensitiveStrikes";
        
        try
        {
            DoubleArraySP sensStrikes(new DoubleArray(0));
            SensitiveStrikeDescriptor sensStrikeDesc;
            sensStrikeDesc.forwardOnly = false;
            DateTime instStartDate = fwdStarting? startDate : valueDate;
            DateTime maturityDate;
            DateTime loDate; // not used
            avgOut->getBoundingDates(loDate, maturityDate);

            // Now get the Extendible specific sensitive strikes
            DoubleArray extStrikes = extensionSchedule->getValues();
            DoubleArraySP extSensStrikes(new DoubleArray(extStrikes.size()*2 - 2));

            // Adjust the strike list with the low and high spreads
            int i;
            for (i = 0; i < extStrikes.size()-1; i++)
            {
                (*extSensStrikes)[2*i] = extStrikes[i]+lowKspread;
                (*extSensStrikes)[2*i+1] = extStrikes[i]+highKspread;

                if (!fwdStarting)
                {
                    (*extSensStrikes)[2*i] *= initialSpot;
                    (*extSensStrikes)[2*i+1] *= initialSpot;
                }
            }

            // for each strike create a vol request to fire at the asset
            for (i = 0; i < extSensStrikes->size(); i++)
            {
                LinearStrikeVolRequest volRequest((*extSensStrikes)[i],
                                                  instStartDate,
                                                  maturityDate,
                                                  fwdStarting);
                
                // Generic does the rest
                getSensStrikes(outputName,
                              &volRequest,
                               sensStrikeDesc,
                               sensStrikes);
            }

            return sensStrikes;
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

    /** when to stop tweaking */
    DateTime endDate(const Sensitivity* sensControl) const{
        return extensionSchedule->lastDate();
    }

    bool sensShift(Theta* shift) {
        static const string method = "Extendible::sensShift";
        try  {
            // roll the parent (updates value date etc)
            Generic1Factor::sensShift(shift);
            // fill in any samples
            avgOut->roll(shift->getUtil(valueDate), 0 /* iAsset */,
                         asset.get());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
        return true; // our components have theta type sensitivity
    }
            
/**  a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    bool priceDeadInstrument(CControl* control, CResults* results) const{
        
        bool deadInstrument = false;
        
        DateTime equityMaturity;
        DateTime loDate; // not used
        avgOut->getBoundingDates(loDate, equityMaturity);

        if (valueDate.isGreaterOrEqual(equityMaturity))
        {
            deadInstrument = true;
            results->storePrice(0.0, discount->getCcy());
            addRequests(control, results, 0.0, equityMaturity);
        }

        return deadInstrument;
    }

    
    
    /** Price the Extendible instrument */
    void price(Control* control, CResults* results)const{
        static const string method = "Extendible::price";

        try
        {
            
            // valueDate >= equityMaturityDate is taken care of here
            if(priceDeadInstrument(control, results)){
                return; // dead instrument priced
            }
            
            DateTimeArray extDates = extensionSchedule->getDates();
            int numLibors = extDates.size()-1;
            DoubleArray liborValues(numLibors);
            CashFlowArrayConstSP cfl;
            DateTime lastPayDate;
            DateTime equityMaturity;
            DateTime loDate; // not used
            avgOut->getBoundingDates(loDate, equityMaturity);
            int i;
            
            /* Price the libor stream for each extension date.
               The libor streams extend from each date in the 
               extension schedule back to the first date in the schedule */
            for (i = 1; i < extDates.size(); i++)
            {
                lastPayDate = extDates[i];
                
                // generate cash flow list for this liborStream
                cfl = payStream->cashFlowList(equityMaturity,
                                              lastPayDate);
                
                // PV the cash flow list to the equity maturity date
                liborValues[i-1] = cashFlowPV(equityMaturity,
                                          *cfl);
            }
            
            // now price the put spreads using the AverageSpot pricer
            DoubleArray spreadValues(numLibors);
            DoubleArray digitals(numLibors);
            DoubleArray extLevels = extensionSchedule->getValues();
            double highK, lowK;
            
            for (i = 0; i < numLibors; i++)
            {
                highK = extLevels[i] + highKspread;
                lowK = extLevels[i] + lowKspread;
                
                if (!fwdStarting)
                {
                    highK *= initialSpot;
                    lowK *= initialSpot;
                }
                
                /* lowK and highK are passed by reference therefore :
                   If fwdStarting = TRUE, lowK  and highK will be converted
                   to absolute values in Average::priceSpotSpread so no need
                   to calculate the asset fwd price at the start date again here */
                double spreadValue = priceSpread(lowK,highK, equityMaturity); 
                
            spreadValue /= (highK - lowK);
            
            spreadValues[i] = spreadValue;
            }
            
            /* calculate the digitals from the spreads
               e.g the first libor must be multiplied by 
               the first spread - second spread */
            for (i = 0; i < numLibors-1; i++)
            {
                digitals[i] = spreadValues[i] - spreadValues[i+1];
            }
            /* the last digital is just the last spread value
               since the value of a put with zero strike is zero */
            digitals[numLibors-1] = spreadValues[numLibors-1];
            
            // calculate the total price
            double value = 0.0;
            for (i = 0; i < numLibors; i++)
            {
                value += digitals[i] * liborValues[i];
            }
            
            // record the total price in the output 
            results->storePrice(value, discount->getCcy());
            if (control && control->isPricing())
            {
                addRequests(control, results, value, equityMaturity);
            }
            
        }
        catch (exception& e)
        {
            throw ModelException(e, method);
        }
    }

private:
    
    double                    lowKspread;          // low spread to apply to put strike
    double                    highKspread;         // high spread to apply to put strike
    double                    spread;              // libor spread
    string                    liborDCC;            // day count conv for float rate
    string                    liborReset;          // reset interval for libor stream
    ScheduleSP                extensionSchedule;   // schedule of libor extensions
    SampleListSP              avgOut;              // averaging out samples

    // transient
    PayStreamSP               payStream;           // generated internally
    InstrumentSettlementSP    settle;              // ignore input and set to last avgOut date
};

const double Extendible::MINIMUM_SPREAD = 0.001;

class ExtendibleHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Extendible, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate); 
        EMPTY_SHELL_METHOD(defaultExtendible);
        FIELD(lowKspread, "low strike spread");
        FIELD(highKspread, "high strike spread");
        FIELD(spread, "libor spread");
        FIELD(liborDCC, "libor day count convention");
        FIELD(liborReset, "libor reset interval");
        FIELD(extensionSchedule, "libor extension schedule");
        FIELD(avgOut, "averaging out dates");
        FIELD(payStream, "payment stream");
        FIELD(settle, "override of instrument settlement");

        // hide from dd interface but tweak pay stream
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(payStream); 
        FIELD_MAKE_TRANSIENT(settle);
    }

    static IObject* defaultExtendible(){
        return new Extendible();
    }
};

CClassConstSP const Extendible::TYPE = CClass::registerClassLoadMethod(
    "Extendible", typeid(Extendible), ExtendibleHelper::load);

/** private class */
class ExtendibleClosedForm: public CClosedFormLN::IProduct{
private:
    const Extendible* extb; // a reference

public:
    ExtendibleClosedForm(const Extendible* extb): extb(extb){}

    void price(CClosedFormLN* model,
               Control*    control, 
               CResults*   results)const{
        extb->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* Extendible::createProduct(CClosedFormLN* model) const
{
    return new ExtendibleClosedForm(this);
}

/* for class loading */
bool ExtendibleLoad() {
    return (Extendible::TYPE != 0);
}


DRLIB_END_NAMESPACE

