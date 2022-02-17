/* DRIntfc_c.cpp */
#pragma warning (disable : 4786 )

#define MAXRATES 400 

#include "DRIntfc_c.h"

#include <frmpp.h>   /* DR prepay model */
#include "ppbase.h"
#include <assert.h>
#include <sstream>
#include <string>
#include <primus_rates.h>

#if defined (_MSC_VER)  // Microsoft compiler

    #define P_DIR "//facetwin_p/kapital/work/drwrapper/eurombs/froasdir/refirates/"

#elif defined (__SUNPRO_CC) // Sun compiler

    #define P_DIR  "/work/drwrapper/eurombs/froasdir/refirates/"

#else // if defined (_LINUX)

    #define P_DIR  "/work/drwrapper/eurombs/froasdir/refirates/"

#endif

#define P_CC_RATES "/ppm_regress_historical_cc_rates.dat"
#define P_RATES "/ppm_regress_historical_rates.dat"

/**/
namespace {

// accumulate alib error msg in a static string
static std::string s_alibErrorMsgAccumulator;


extern "C" {
    
    static TBoolean AlibErrorAccumulate(char *string, void *callBackData)
    {
        // typedef TBoolean (TErrCallBackFunc)(char *string, void *callBackData);
        s_alibErrorMsgAccumulator += string;
        return FALSE;  // don't log the error to file
    }
}


// Empty accumulator and return accumlated
static inline std::string EmptyAlibErrorAccumulator()
{
    std::string tmp = s_alibErrorMsgAccumulator;
    s_alibErrorMsgAccumulator = std::string();
    return tmp;
}


// Accumulate alib error msg and append them to exceptions
// Approach is to leverage Alib error handling to capture Alib
// and Q3 errors.
static inline int RegisterAlibErrorToException()

{
    /* Register Alib error callback */
    GtoErrMsgAddCallback(AlibErrorAccumulate, false, NULL);
    GtoErrMsgOn();
  
    return 0;
}

// Make the registration happen

class AlibErrorAccumulate
{

public:
        AlibErrorAccumulate() {
                RegisterAlibErrorToException();
        }

} _some_var;

}
/**/

using namespace std;

static ostringstream DRErrorstream; 

class DRHandle {

public:

    TMbsRateEnv *rates;
    TMbsDeal *cohort;
    TMbsPrepayDefaults *model;
    long today;
   
    std::vector<long> dates;
    std::vector<double> crates;
    std::vector<double> points;

    void loadRates(char*path, long origDate);
    void loadRatesx(char*path);
    void initRates();
    
    DRHandle(TMbsRateEnv *r, TMbsDeal* c, TMbsPrepayDefaults* p, long t) {
        rates = r;
        cohort = c;
        model = p;
        today= t;
    }

    ~DRHandle() {
        free(rates->adjAmortIndexRates);
        free(rates->amortIndexRates);
        free(rates->fhCommitPts);
        free(rates->histAmortIndexRates);
        free(rates->histFhCommitPts);
        free(rates->adjHistAmortIndexRates);
        free(rates->histAmortIndexDates);
        delete rates;
        delete cohort->mbsPrepayAssump;
        delete cohort;
        delete model;
    }

    TDate lastHistDate;
    double lastHistRate;

} ;


void DRHandle::loadRates(char* path, long origDate) {

    DRDate _today = toDate(today);

    DRDate _origDate = toDate(origDate);

    DRDate startDate = _origDate;
    startDate.SetDay(1);
    DRDateIn x(cohort->mbsPrepayAssump->wala + 9, 'M', FALSE);
    startDate -= x;

    string ccFile(path);
    if (ccFile == "DEFAULT")
        ccFile = P_DIR;
    ccFile+=P_CC_RATES;
   
    string rateFile(path);
    if (rateFile == "DEFAULT")
        rateFile = P_DIR;
    rateFile+=P_RATES;


    try {
        
        PrimusHistRates primus(_today+ONE_DAY,ccFile.c_str(),rateFile.c_str());
        
        int term = 30;
        if (cohort->mbsTerm == MBS_PP_MBSTERM_15)
            term = 15;
        
        lastHistDate = primus.GetLastHistRateDate();
        lastHistRate = primus.GetMtgRate(lastHistDate, term);
        
        
        DRDate tempDate = startDate;
        while (tempDate <= _today && tempDate < lastHistDate) {
            
            //TDate dt = tempDate;
            dates.push_back(tempDate.GetYear()*100+tempDate.GetMonth());
            crates.push_back(primus.GetMtgRate(tempDate, term));
            points.push_back(primus.GetPoint(tempDate, term)); 
            //        std::cerr << "load rates dt " << tempDate.GetYear()*100+tempDate.GetMonth()
            //            << " r "
            //            << primus.GetMtgRate(tempDate, term) << " pt "
            //            << primus.GetPoint(tempDate, term) << std::endl;
            tempDate += ONE_MONTH;
            
        }
        
    }
    catch (...) {
        DRErrorstream << "Failed to load primus rates " << std::endl;
        DRErrorstream << "CC File " << ccFile  << std::endl;
        DRErrorstream << "Rate File " << rateFile << std::endl;
        throw DRErrorstream.str();
    }
}

void DRHandle::loadRatesx(char* path) {

    std::string txt;
    std::ifstream from(path);
    if (!from) {
        DRErrorstream << "can't open hist data file";
        throw DRErrorstream.str();
    }

    // read past comments
    while (from)
        if (from.peek() == '#' )
            getline(from, txt);
        else 
            break;

    //read data    
    long date;
    double rate,point;
    while ( from >> date >> rate >> point ) { 
        dates.push_back(date);
        crates.push_back(rate);
        points.push_back(point); 
//        std::cerr << " hist " << date << " r " << rate << " pts " << point
//            << std::endl;
    }

    if (dates.size() == 0) {
        DRErrorstream << "no historical rates loaded";
        throw DRErrorstream.str();
    }
    

   from.close();
}

/* Notes: DRFRM model is quite robust with respect to historical rate
 * date. If data is not current, DRFRM code fills in missing data 
 * based on last know set of rates.  Secondly, DRFRM model ignores points
 * that are > today 
 */

void DRHandle::initRates() {
    int i=0;
    
    // first simulated rate will be next month; DRFRM model
    // will automatically use historical rates up through 
    // month previous to this date.
    // bmc:  we need a #define for this day-of-month thing....
    // 
    TDate bom;
   
    TDateOf(today,&bom);
    SetDOM(1,TRUE,&bom);
    //note:if we change prepay start date to be sim dt+1, then we need
    // to make this bom + 2
    NxtMth(bom,1,&(rates->amortIndexStartDate));

    // declared as TDate, but really expected to be YYYYMMDD:
    rates->histAmortIndexStartDate = dates[0]*100+1;
    rates->numHistAmortIndexRates = dates.size();
    // declared as TDate, but expected to be YYYYMM:
    rates->histAmortIndexDates 
        = (TDate*)calloc(dates.size(),sizeof(TDate));  
    rates->histAmortIndexRates 
        = (double*)calloc(crates.size(),sizeof(double));
    rates->histFhCommitPts 
        = (double*)calloc(points.size(),sizeof(double));
    rates->adjHistAmortIndexRates 
        = (double*)calloc(crates.size(),sizeof(double));

     // alloc non-hist rate structures

    rates->numAmortIndexRates = 0;
    rates->amortIndexRates 
        = (double*)calloc(MAXRATES,sizeof(double));
    rates->adjAmortIndexRates 
        = (double*)calloc(MAXRATES,sizeof(double));
    rates->fhCommitPts 
        = (double*)calloc(MAXRATES,sizeof(double));

    if (! rates->histAmortIndexDates  ||
        ! rates->histAmortIndexRates  ||
        ! rates->histFhCommitPts ||
        ! rates->adjHistAmortIndexRates ||
        ! rates->amortIndexRates ||
        ! rates->adjAmortIndexRates ||
        ! rates->fhCommitPts) {

        DRErrorstream << "Allocation error";
        throw DRErrorstream.str();
    }

    for (i=0; i<dates.size(); i++)
        rates->histAmortIndexDates[i] = dates[i];

    for (i=0; i<crates.size(); i++)
        rates->histAmortIndexRates[i]=crates[i];
    
    for (i=0; i<points.size(); i++) 
        rates->histFhCommitPts[i]=points[i];
    
    // DRFRM function that applies points to the mortgage rate
    // to produce the points-adjusted mortgage rate
    if (PointAdjustCommitRate(rates->numHistAmortIndexRates,
			    rates->amortIndexType,
			    rates->histAmortIndexRates,
			    rates->histFhCommitPts,
			    model->desiredRefiPoints,
			    model->ioMultiplier,
			    rates->adjHistAmortIndexRates) != SUCCESS)
    {
      DRErrorstream << "hist point adjustment failed";
      throw DRErrorstream.str();
    }

    // Determine current points

    long idx;
    TDate histSt;
    TDateOf(rates->histAmortIndexStartDate,&histSt);
    MonsDiff(histSt,bom,&idx);

    if (idx < 0 || idx > points.size()) {
        DRErrorstream << "Error determining current points";
        throw DRErrorstream.str();
    }
    
    double currPts =rates->histFhCommitPts[idx-1];
    for (int k=0; k<MAXRATES; k++)
        rates->fhCommitPts[k]=currPts;

}

void* DR_create_handle(
    PP_POOL* pool, 
    PP_DRFRM* params,
    char* ratePath, 
    long today)
{

	TMbsRateEnv *rates=NULL;
	TMbsPrepayDefaults* ppm=NULL;
    TMbsPrepayAssump* assump=NULL;
    TMbsDeal* cohort=NULL;
    DRHandle* handle=NULL;

    int i;

    try {
        
        /* cohort */
        cohort = new TMbsDeal;
        
        cohort->mbsPrepayAssump = NULL;
        
        switch (pool->agency) {
            
        case PP_FNMA:
            cohort->mbsAgency = MBS_AGENCY_FNMA;
            break;
            
        case PP_FHLMC:
            cohort->mbsAgency = MBS_AGENCY_FHLMC;
            break;
            
        case PP_GOLD:
            cohort->mbsAgency = MBS_AGENCY_GOLD;
            break;
            
        case PP_GNMAI:
            cohort->mbsAgency = MBS_AGENCY_GNMAI;
            break;
            
        case PP_GNMAII:
            cohort->mbsAgency = MBS_AGENCY_GNMAII;
            break;
            
        case PP_WHOLE:
            cohort->mbsAgency = MBS_AGENCY_WHOLE;
            break;
            
        default:
            DRErrorstream << "DR_create_handle: bad agency type" 
                << std::endl;
            throw DRErrorstream.str();
            
        }
        
        if (pool->WALA+pool->WARM <= 180) {
            cohort->mbsTerm = MBS_PP_MBSTERM_15;
        }
        else {
            cohort->mbsTerm = MBS_PP_MBSTERM_30;
        }
        
        assump = new TMbsPrepayAssump;
        
        cohort->mbsPrepayAssump = assump;
        
        TDate today_t;
        TDateOf(today,&today_t);
        
        assump->startDate = today_t;
        assump->numAmortMons = 1;  // always do 1 month of cpr calcs at a time
        assump->grossCpn = pool->WAC;
        assump->warm = pool->WARM;
        assump->wala = pool->WALA;
        assump->inclSchedAmort = 0;
        assump->amortForm = MBS_PP_SPD_CPR;
        assump->turnoverTwk = params->turnoverTwk;
        assump->incentPtTwk = params->incentPtTwk;
        assump->totSteepTwk = params->totSteepTwk;
        assump->baseSteepTwk = params->baseSteepTwk;
        assump->ageRampTwk = params->ageRampTwk;

        assump->prepayModel = MBS_PP_MODEL_FIX_MGRP;
        assump->inclBullet = 0;
        assump->useNewWacInfo = 0;
        assump->bulletDate = 0;
        assump->grossMargin = 0;
        assump->prepayVector = NULL;
        assump->stdSmmSpdTwk = 1.0;
        assump->stdSmmFlrTwk = 1.0;
        assump->stdSteepTwk = 1.0;
        assump->stdIncentPtTwk  = 1.0;
        for (i=0; i<12; i++)
            assump->seasonality[i] = 0;
        for (i=0; i<8; i++)
            assump->mgrpScurveLogist[i] = 0;
        for (i=0; i<4; i++)
            assump->mgrpSeasLogist[i] = 0;
        
        /* defaults */
        ppm = new TMbsPrepayDefaults;
        
        ppm->desiredRefiPoints = params->RefiPt;
        ppm->ioMultiplier = params->IOMult;
        memcpy(ppm->seasonality,params->SeasMlt,sizeof(ppm->seasonality));
        ppm->numGroups = params->NbGrp;
        ppm->gpAddlSmm = params->AddSMM;
        ppm->critRat = params->CritRt;
        memcpy(ppm->scurveLogist,params->SCrvPrm,sizeof(ppm->scurveLogist));
        memcpy(ppm->seasLogist,params->SLogPrm,sizeof(ppm->seasLogist));
        memcpy(ppm->groupPs,params->GrpPrm,sizeof(ppm->groupPs));
        memcpy(ppm->ltAges,params->LTAPrm,sizeof(ppm->ltAges));
        memcpy(ppm->mgrpLagWgts,params->RfLgWt,sizeof(ppm->mgrpLagWgts));

        ppm->maxLag = 0;
        ppm->smmLogistRight = 0;
        ppm->smmLogistRight = 0;
        ppm->smmLogistWidth = 0;
        ppm->smmLogistInflec = 0;
        ppm->seasLogistRight = 0;
        ppm->seasLogistLeft = 0;
        ppm->seasLogistWidth = 0;
        ppm->seasLogistInflec =0;
        ppm->histFh30Idx =0;
        ppm->histFh15Idx =0;
        ppm->histCmt10Idx  = 0;
        ppm->refiIncShift = 0;
        ppm->newWacRefiIncShift = 0;
        ppm->absoluteRateEffect = 0;

        /* rate env */
        rates = new TMbsRateEnv;
        
        if (cohort->mbsTerm == MBS_PP_MBSTERM_30) {
            rates->amortIndexType = MBS_PP_REFI_TYPE_FH30;
        }
        else {
            rates->amortIndexType = MBS_PP_REFI_TYPE_FH15;
        }
        
        rates->histAmortIndexDates = NULL;
        rates->histAmortIndexRates = NULL;
        rates->histFhCommitPts = NULL;
        rates->adjHistAmortIndexRates = NULL;
        
        rates->amortIndexRates = NULL;
        rates->adjAmortIndexRates = NULL;
        rates->fhCommitPts = NULL;
        
        handle = new DRHandle(rates,cohort,ppm,today);
        
        // load historical data
        handle->loadRates(ratePath, pool->origDate);
        // put it into the DRFRM structure
        handle->initRates();
    }

    catch (...) {
        DRErrorstream << "DR_create_handle failed" << std::endl;
        handle=NULL;
    }
   
    return (void*)handle;
}

void DR_destroy_handle (
    void* handle) 
{
      
    DRHandle* h = (DRHandle*)handle;
    delete h;
}

int DR_Cpr(
           void* handle, 
           long simDate,
           PP_POOL* pool,
           int sz,
           double* r,
           double* cpr)
{
    int status = FAILURE;
    int i;
    *cpr = 0.;
    double q,blnd;
    long dom;
    DRHandle* h = (DRHandle*)handle;
    
    TDate simDate_t;
    TDateOf(simDate,&simDate_t);
    
    // update sim date
    h->cohort->mbsPrepayAssump->startDate = simDate_t;
     
    
    //NxtMth(simDate_t,1,&(h->cohort->mbsPrepayAssump->startDate));
    
    // update pool data
    h->cohort->mbsPrepayAssump->grossCpn = pool->WAC;
    h->cohort->mbsPrepayAssump->warm = pool->WARM;
    h->cohort->mbsPrepayAssump->wala = pool->WALA;
    
    // update rates
    if (sz > MAXRATES) {
        DRErrorstream << "DR_Cpr: too many rates" << std::endl;
        goto RETURN;
    }
    
    h->rates->numAmortIndexRates = sz;
    for (i = 0; i<sz; i++) 
        h->rates->amortIndexRates[i] = r[i];
    
    // apply point adjustment
    if( PointAdjustCommitRate(h->rates->numAmortIndexRates,
        h->rates->amortIndexType,
        h->rates->amortIndexRates,
        h->rates->fhCommitPts,
        h->model->desiredRefiPoints,
        h->model->ioMultiplier,
        h->rates->adjAmortIndexRates) != SUCCESS ) {
        
        DRErrorstream << "DR_Cpr: point adjust failure" << std::endl;
        goto RETURN;
    }
    
    /* blending logic that old pricer does*/
    if ( h->rates->numAmortIndexRates > 0 
        && h->rates->numHistAmortIndexRates > 0) {
        
        GetDOM(h->lastHistDate,&dom);
        blnd = dom/31.0;
        q = h->rates->adjHistAmortIndexRates[h->rates->numHistAmortIndexRates-1];
        h->rates->adjHistAmortIndexRates[h->rates->numHistAmortIndexRates-1] =
            q*blnd + 
            h->rates->adjAmortIndexRates[0]*(1-blnd) ;
    }

    double c[1];
    
    if (frm_mgrp_prepays(h->cohort,h->rates,h->model, c) != SUCCESS) {
        DRErrorstream << EmptyAlibErrorAccumulator() << std::endl;
        
        DRErrorstream << "DR_Cpr: frm_mgrp_prepays(...) failed" 
            << std::endl;
        goto RETURN;
    }
    
    *cpr = c[0];
    status = SUCCESS;
    
RETURN:    
    

    return status;;
}

void DR_get_error(char* s, int sz) {
    int i;
    string e = DRErrorstream.str();
    for (i=0; i<sz; s[i]=(char)NULL, i++);
    for (i=0; i<sz && i<e.size(); i++)
        s[i] = e[i];

    DRErrorstream.str("");
}
