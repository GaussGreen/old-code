///****************************************************************************/
///*      CHFInterface class: interface to CHF prepayment model        .      */
///****************************************************************************/
///*      CHFIntfc.cpp                                                        */
///****************************************************************************/

#include "chfintfc.h"
#include <sstream>

ostringstream CHFErrorstream; 

#define DUMMY_TRANCHE_ID "X"
#define DUMMY_TRANCHE_ID_LEN 1

extern CForecast ForecastObj;

CHFInterface::CHFInterface(
        long today,
		PoolType poolType, 
		double wac, 
		double originalTerm, 
		double originalBalance, 
		string dataDir) {

        m_today = today;
        m_poolType = poolType;
		m_wac = wac;
		m_originalTerm = originalTerm;
		m_originalBalance = originalBalance;
		m_CHFDataDirectory = dataDir;

        m_cForecast=NULL;

        if (CHFINTFCDDG) {

            fprintf (stderr, "\nCHFInterface::CHFInterface\n");
            fprintf (stderr, "\ttoday %ld\n", m_today);
            fprintf (stderr, "\tpoolType %d\n", m_poolType);
            fprintf (stderr, "\twac %lf\n", m_wac);
            fprintf (stderr, "\toriginalTerm %lf\n", m_originalTerm);
            fprintf (stderr, "\toriginalBalance %lf\n", m_originalBalance);

            fprintf (stderr, "\n\n");
            fprintf (stderr, "path\tt\twala\tbal\tmtg_0\tmtg_1\tmtg_2\tswap2_2\tswap10_2\tcpr\n");
        }
}

CHFInterface::~CHFInterface() {
    if (m_cForecast) {
        m_cForecast->Unload();
        //delete m_cForecast;
        m_cForecast=NULL;
    }
}

double CHFInterface::calculateCPR(
        int p,
        int t,             // months
		double wala,       // months
		double remBalance,
        int bufSz,
        double* mtg,       
	    double* swap2Y, 
		double* swap10Y   
        )
{
  
    // set up CHF_RATE_STRUCTURE
   
    CHF_RATE_STRUCTURE *rates =  m_cForecast->m_stFwdRates;

    int i,j=0;
    
    /* logic
    
      rates[t-1] = SRM3[0] // current
      rates[t-2] = SRM3[1] // 1M lag, etc
      rates[t-3] = SRM3[2]
      rates[t-4] = SRM3[3]
      
        also note, at time point t (months), CHF looks at rates[t-1] 
        for the current rate. 
        
    */
    for (j= 0, i = t-1; j< bufSz; j++, i--) {
        if (i>=0) {
            rates[i].dSwap_2y = swap2Y[j];
            rates[i].dSwap_10y = swap10Y[j];
            rates[i].dMortgage_rate = mtg[j];
        }
    }

    // tells the CForecast object that CHF_RATE_STRUCTURE has changed

    m_cForecast->m_bReloadRates = 0.; 

	// set up CHFPREPAYMENTMODELSTRUCT
    
	m_cForecast->m_lpPrepayParams->simulation_year = t > 12 ? t/12 -1 : 0;

    m_cForecast->m_lpPrepayParams->simulation_month 
        = t - m_cForecast->m_lpPrepayParams->simulation_year *12 ;

    switch (m_poolType) {
        
    case FNMA30:
        m_cForecast->m_lpPrepayParams->loan_type = 0;  
        break;
        
    case GNMA30:
        m_cForecast->m_lpPrepayParams->loan_type = 1;  
        break;
        
    case FNMA15:
        m_cForecast->m_lpPrepayParams->loan_type = 2 ;
        break;
        
    default:
        CHFErrorstream << "Pool type not supported by CHF prepay model"
            << std::endl;
        throw CHFErrorstream.str();
        
    }
	
	m_cForecast->m_lpPrepayParams->gross_wac = m_wac;
	m_cForecast->m_lpPrepayParams->original_term = m_originalTerm;
	m_cForecast->m_lpPrepayParams->remain_term = m_originalTerm - wala;
	m_cForecast->m_lpPrepayParams->wala = wala;
	m_cForecast->m_lpPrepayParams->remain_balance = remBalance;

	// call prepay model
    double cpr =0;
    cpr = m_cForecast->CalcPrepay();
	
    if (CHFINTFCDDG) {
        fprintf (stderr, "%d\t%d\t%3.2lf\t%3.1lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\t%1.4lf\n",
            p,t,wala,remBalance,mtg[0],mtg[1],mtg[2],swap2Y[2],swap10Y[2],cpr);
     }   

    return cpr;

}

void CHFInterface::initialize() {

    // Set pointer to CForecast object
    // We want to used this global copy, since it is referenced in some
    // of the member functions.
   
    // **NOTE: talk to CHF about  the ForecastObj global ...
    m_cForecast = &ForecastObj;  

    m_cForecast->Load();

    // initialize the CHF rate structure

    for (int i =0; i<CHF_MAX_FORWARD_RATES+1; i++) {
        memset((void*)&(m_cForecast->m_stFwdRates[i]),0,sizeof(CHF_RATE_STRUCTURE));
    }

	// deal with original balance. CHF code expects original balance to be an
	// underlying CTranche object that is referenced by a string ID.

	CTranche *tranche = new CTranche(m_originalBalance);  // dealloc in CForecast.unload()
	m_cForecast->m_pTrancheDataObj->m_TrancheMap[(LPCSTR)DUMMY_TRANCHE_ID] = tranche; 
	m_cForecast->pTranche = tranche;
	strncpy(m_cForecast->m_lpPrepayParams->szTranche_id,DUMMY_TRANCHE_ID,DUMMY_TRANCHE_ID_LEN);   

    // Set m_cForecast->m_dtValuationMonth; Observe that this a CTime object, 
    // not just an integer month.
    long year = m_today/10000;
    long month = (m_today - year*10000)/100;
    long day = m_today - year*10000 - month*100;
    m_cForecast->m_dtValuationMonth = CTime(year,month,day,0,0,0);

	if (!m_cForecast->m_bLoaded) {
        CHFErrorstream << "CForecast object failed to load" << std::endl;
        throw CHFErrorstream.str();
	}

}

// Placeholder, not currently used
double CHFInterface::getHistoricalRate(int i) {
	return 0.0;

}