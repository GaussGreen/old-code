//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : EnergyUnderlyer.hpp
//
//   Description : An Energy 'asset'
//                 E.g. WTI_NYMEX (Oil), NG_NYMEX (for natural gas).
//                 We chose not to use the word asset as in energy there is 
//                 generally no spot price and the market is futures based
//                 
//                 This class defines the date rules for the underlyer - its
//                 futures expiry ruel and if it as options on these futures
//                 the option expiry.
//
//                 The equivalent class in FXLIB was DRCommodityIndex but as
//                 with Asset CommodityIndex already has a meaning in QLIB
//                 for things like GSCI DJ-AIG
//
//   Author      : Sean Chen
//
//   Date        : April 1, 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_ENERGY_UNDERLYER_HPP
#define EDR_ENERGY_UNDERLYER_HPP

#include "edginc/MarketFactor.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/EnergyContractLabel.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL EnergyUnderlyer : public MarketObject,
                        virtual public IMarketFactor
{
    public:
    
        static CClassConstSP const TYPE;

        friend class EnergyUnderlyerHelper;
        
        virtual ~EnergyUnderlyer();

        static const int INDEX;
        static const int OPTION;

      
        string getName() const;
        void getMarket(const IModel* model, const MarketData* market);
        string getPricingCurrency() const;
        double getNotional() const;
        string getUnderlyerType() const;
        string getParentName() const;
        HolidayConstSP getHoliday() const;

        // functionality
        DateTime expiryDate(const string& benchmark) const;
        DateTime expiryDate(int month, int year) const;
        DateTime expiryDate(const EnergyContractLabel& label) const;
        DateTime expiryDate(const EnergyContractLabel& label, int contractType) const;  // index and option

        void expiryDateList(
            vector<DateTime>* dates, 
            const EnergyContractLabel& startBenchmark, 
            const EnergyContractLabel& endBenchmark,
            int contractType,
            int contractIncrement = 1) const;

        void expiryDateList(
            vector<DateTime>* dates,
            const DateTime& startDate,
            int numContracts,
            int contractType,
            int contractIncrement = 1) const;

        void expiryDateList(
            vector<DateTime>* dates,
            const EnergyContractLabel& startLabel,
            int numContracts,
            int contractType,
            int contractIncrement = 1) const;

        void expiryDateList(
            vector<DateTime>* dates,
            const DateTime& startDate,
            const DateTime& endDate,
            int contractType,
            int contractIncrement = 1) const;

        void energyContractLabelList(
            vector<EnergyContractLabel>* labelList,
            const EnergyContractLabel& startLabel,
            const EnergyContractLabel& endLabel,
            int contractIncrement = 1) const;
            
        DateTime optionExpiryDate(const string& theBenchmark) const;
        DateTime optionExpiryDate(int theMonth, int theYear) const;
        DateTime optionExpiryDate(const EnergyContractLabel& theLabel) const;

        const YieldCurve* getYieldCurve() const;

        // determine the contract for the given fixing date
        void calculateContract(const DateTime& fixingDate, EnergyContractLabel* label, int contractType) const;
        EnergyContractLabel calculateContractLabel(const DateTime& fixingDate, int contractType) const;

        int calculateContractNumber(const DateTime& today, const EnergyContractLabel& label, int contractType) const;
        int calculateContractNumber(const DateTime& today, const DateTime& fixingDate, int contractType) const;
        EnergyContractLabel calculateContractLabel(const DateTime& today, int contractNumber, int contractType) const;
        void calculateContractMonthBoundaryDates(
            const DateTime& fixingDate, 
            int contractType,
            DateTime* startDate, 
            DateTime* endDate) const;

        EnergyContractLabel addContract( const EnergyContractLabel& originalLabel, int numContracts = 1) const;  // alter input parameter and return

        bool isHoliday(const DateTime& date) const;
        bool isHoliday(int day, int month, int year) const;
        bool isBusinessDay(const DateTime& date) const;
        bool isBusinessDay(int day, int month, int year) const;
        bool hasOptionOnUnderlyer() const;
        bool isExpiryDate(const DateTime& fixingDate, int contractType) const;

        int generateContractSchedule(
            vector<DateTime>* initialDates, 
            vector<DateTime>* finalDates, 
            vector<EnergyContractLabel>* labels,
            const DateTime& startDate,
            const DateTime& endDate,
            int contract,
            int rollOffset) const;

		DateTime busDayAdjustFollowing(const DateTime& date) const;
        DateTime busDayAdjustPrior(const DateTime& date) const;
		DateTime busDayAdjust(const DateTime& date) const;
		DateTime addMonthsToDate(const DateTime& date, int numMonths) const;
		DateTime addIBORMonthsToDate(const DateTime& date, int numMonths) const;
		DateTime addIBORMonthsToDateNoAdj(const DateTime& date, int numMonths) const;
		DateTime addBusDaysToDate(const DateTime& date, int numDays) const;
		DateTime addDaysToDateAndAdjust(const DateTime& date, int numDays) const;

		DateTime addContractPeriod(const DateTime& date, int numContracts, int contractType) const;

    private:

        EnergyUnderlyer();
    
        string            name;
        string            underlyerType;
        HolidayWrapper            holidayWrapper;
        YieldCurveWrapper        pricingCurrency;
        string                units;
        double                notional;
        string            expiryRule;
        string        expiryMonthRule;
        bool                optionOnUnderlyer;  // true/false flag
        int                optionExpiryRule;
        int                settlementDaysRule;
        string            parentName;
};

typedef MarketWrapper<EnergyUnderlyer> EnergyUnderlyerWrapper;
typedef smartPtr<EnergyUnderlyer> EnergyUnderlyerSP;
typedef smartConstPtr<EnergyUnderlyer> EnergyUnderlyerConstSP;

DRLIB_END_NAMESPACE

#endif

