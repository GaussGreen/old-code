#include "drrate.h"
#include "drexception.h"


map <upper_string, DRRate*, less <upper_string>, MYALLOC(DRRate*) > DRRate::m_allRates;

DRRate LIBOR1M ("LIBOR1M", 12, 1, GTO_ACT_360, GTO_CURVE_DISCOUNT, 0);
DRRate LIBOR3M ("LIBOR3M", 4, 3, GTO_ACT_360, GTO_CURVE_DISCOUNT, 0);
DRRate LIBOR1Y ("LIBOR1Y", 1, 12, GTO_ACT_360, GTO_CURVE_DISCOUNT, 0);
DRRate LIBOR2Y ("LIBOR2Y", 2, 24, GTO_B30_360, GTO_CURVE_DISCOUNT, 0);
DRRate LIBOR10Y ("LIBOR10Y", 2, 120, GTO_B30_360, GTO_CURVE_DISCOUNT, 0);

DRRate CMT1M ("CMT1M", 12, 1, GTO_ACT_365F, GTO_CURVE_INDEX1, 0);
DRRate CMT3M ("CMT3M", 4, 3, GTO_ACT_365F, GTO_CURVE_INDEX1, 0);
DRRate CMT1Y ("CMT1Y", 2, 12, GTO_B30_360, GTO_CURVE_INDEX1, 0);
DRRate CMT2Y ("CMT2Y", 2, 24, GTO_B30_360, GTO_CURVE_INDEX1, 0);
DRRate CMT10Y ("CMT10Y", 2, 120, GTO_B30_360, GTO_CURVE_INDEX1, 0);

DRRate REFIRATE ("REFIRATE", 0,0,0,0,0);	// mbs specific refirate 

DRRate::DRRate (string name, long cpnsPerYear, long matInMonths, long dayCountConv, 
		long curveIndex, long numSettleDays)
		:m_cpnsPerYear(cpnsPerYear), m_matInMonths(matInMonths), m_name(name),
		m_dayCountConv(dayCountConv), m_curveIndex(curveIndex), m_numSettleDays(numSettleDays)
{
	if (matInMonths != 0) {
		SetFloatArray(cpnsPerYear, matInMonths, dayCountConv, 
			curveIndex, numSettleDays);
	}
	else {
		m_floatRateArray = NULL;
	}
	m_allRates[m_name] = this;
}

/***************************************************************************
 *  SetSingleFloatArray
 *  Constructor for TFloatRateArray to define a single
 *  floating rate (i.e., weight=1.0)
	void SetFloatArray (long cpnsPerYear, long matInMonths, long dayCountConv, long curveIndex, 
		long numSettleDays);

 ***************************************************************************/
void DRRate::SetFloatArray (long cpnsPerYear, long matInMonths, long dayCountConv, 
		long curveIndex, long numSettleDays)
{
    double     defaultWgt = 1.0;
    double     defaultSprd = 0.;
    long       cpnIntInMonths;
	TDateAdjIntvl adjInterval;
	
    // call simple GTO constructor 
    if ((m_floatRateArray = GtoFloatRateArrayMakeSimple
		(&cpnsPerYear,		// payments/yr 
		&matInMonths,		// mat in months 
		&dayCountConv,
		&defaultWgt,		// weight 
		1)) IS NULL)		// # index rates 
		throw DRException("Failed to alloc TFloatRateArray");
	
	/* for simple rates, we also force the cpn interval 
     * to be no longer than the maturity interval (even if 
     * GTO allows this), since we need this to be true
     * for some tree pricing code 
     */
    /* for zero coupon, always force pay pd to mat pd */
    if( cpnsPerYear IS 0 )
    {
        m_floatRateArray->defs[0].payInterval = m_floatRateArray->defs[0].matInterval;
    }
    /* for cpn-bearing, check if pay pd > mat pd */
    else
    {
        cpnIntInMonths = 12 / cpnsPerYear;
        if( cpnIntInMonths > matInMonths )
        {
            m_floatRateArray->defs[0].payInterval = m_floatRateArray->defs[0].matInterval;
        }
    }
    /* also set hidden member data */
	GTO_SET_ADJ_INTERVAL_DAYS (adjInterval, numSettleDays)
    m_floatRateArray->defs[0].spotOffset = adjInterval;
    m_floatRateArray->defs[0].spread = defaultSprd;
    m_floatRateArray->curveIndices[0] = curveIndex;
    m_floatRateArray->rateInfo = NULL;
}

const DRRate* TranslateRate (DRString r)
{
	upper_string rate = r;

	if (rate == "LIBOR1M") return &LIBOR1M;
	else if (rate == "LIBOR3M") return &LIBOR3M;
	else if (rate == "LIBOR1Y") return &LIBOR1Y;
	else if (rate == "LIBOR2Y") return &LIBOR2Y;
	else if (rate == "LIBOR10Y") return &LIBOR10Y;
	else if (rate == "CMT1M") return &CMT1M;
	else if (rate == "CMT3M") return &CMT3M;
	else if (rate == "CMT1Y") return &CMT1Y;
	else if (rate == "CMT2Y") return &CMT2Y;
	else if (rate == "CMT10Y") return &CMT10Y;
	else if (rate == "REFIRATE") return &REFIRATE;
	else
		throw DRException ("Invalid rate name: ") <<  rate;

	return &LIBOR1M;
}

DRString TranslateRate (const DRRate& rate)
{
	if (rate == LIBOR1M) return DRString("LIBOR1M");
	else if (rate == LIBOR3M) return DRString("LIBOR3M");
	else if (rate == LIBOR1Y) return DRString("LIBOR1Y");
	else if (rate == LIBOR2Y) return DRString("LIBOR2Y");
	else if (rate == LIBOR10Y) return DRString("LIBOR10Y");
	else if (rate == CMT1M) return DRString("CMT1M");
	else if (rate == CMT3M) return DRString("CMT3M");
	else if (rate == CMT1Y) return DRString("CMT1Y");
	else if (rate == CMT2Y) return DRString("CMT2Y");
	else if (rate == CMT10Y) return DRString("CMT10Y");
	else if (rate == REFIRATE) return DRString("REFIRATE");
	else
		throw DRException ("Invalid rate");

	return DRString("LIBOR1M");
}

ostream& operator<< (ostream& s, const DRRate& rate)
{
	s << rate.m_name;
	return s;
}

/*
istream& operator>> (istream& s, DRRate* rate)
{
	upper_string name;
	s >> name;

	map <upper_string, DRRate*>::iterator iter = DRRate::m_allRates.find(name);
	if (iter == DRRate::m_allRates.end())
		throw DRException ("Attempt to find a nonexistant rate ") << name;

	rate = (*iter).second;
	return s;
}
*/
