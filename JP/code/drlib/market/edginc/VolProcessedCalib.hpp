//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : VolProcessedCalib.hpp
//
//   Description : Vol Processed Object for Calibration
//
//   Author      : Anwar E Sidat
//
//   Date        : 06-Nov-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_VolProcessedCalib_HPP
#define QLIB_VolProcessedCalib_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/VolProcessed.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(Expiry);

/** Vol processed object for calibration. */
class MARKET_DLL VolProcessedCalib : public CObject,
                                     public virtual CVolProcessed
{
public:
    static CClassConstSP const TYPE;
    
    VolProcessedCalib(
        DateTimeArraySP swapStartDates,
        DateTimeArraySP swapMaturityDates,
        DoubleArraySP   vols,
        ExpiryArraySP   selectedTenors,
        IntArraySP      selectedTenorIndices,
        ExpiryArraySP   selectedExpiries,
        IntArraySP      selectedExpiryIndices,
        DateTimeArraySP selectedSwapMaturityDates);

    /** identifies the market data name of the volatility */
    virtual string getName() const;

    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const;

    /** retrieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const;

    VolProcessedCalib(const CClassConstSP& clazz);
    virtual ~VolProcessedCalib();
    
    /** Processed vol data */
    ExpiryArraySP    getSelectedTenors() const { return m_selectedTenors; }
    IntArraySP       getSelectedTenorIndices() const { return m_selectedTenorIndices; }
    ExpiryArraySP    getSelectedExpiries() const { return m_selectedExpiries; }
    IntArraySP       getSelectedExpiryIndices() const { return m_selectedExpiryIndices; }
    DateTimeArraySP  getSelectedSwapMaturityDates() const { return m_selectedSwapMaturityDates; }
    DoubleArraySP    getVols() const { return m_vols; }
    DateTimeArraySP  getSwapStartDates() const { return m_swapStartDates; }
    DateTimeArraySP  getSwapMaturityDates() const { return m_swapMaturityDates; }

protected:
    string       name;
    TimeMetricSP metric;

private:
    static void load(CClassSP& clazz);
    VolProcessedCalib() : CObject(TYPE) {}

    DateTimeArraySP m_swapStartDates; 
    DateTimeArraySP m_swapMaturityDates; 
    DoubleArraySP   m_vols; 
    ExpiryArraySP   m_selectedTenors;
    IntArraySP      m_selectedTenorIndices;
    ExpiryArraySP   m_selectedExpiries;
    IntArraySP      m_selectedExpiryIndices;
    DateTimeArraySP m_selectedSwapMaturityDates; 
};

typedef smartConstPtr<VolProcessedCalib> VolProcessedCalibConstSP;
typedef smartPtr<VolProcessedCalib> VolProcessedCalibSP;
#ifndef QLIB_VolProcessedCalib_CPP
    EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<VolProcessedCalib>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<VolProcessedCalib>);
#else
    INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<VolProcessedCalib>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<VolProcessedCalib>);
#endif

DRLIB_END_NAMESPACE

#endif
