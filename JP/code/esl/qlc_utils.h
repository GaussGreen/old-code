#ifndef QLC_UTILS_H
#define QLC_UTILS_H
// ESL QLib Conversion Utils
// These functions are used to simplify the printing of QLib objects in XML format.
// This is used to port regression tests to QLib.

#include <stdio.h>
#include <eslhead.h>

#ifdef  __cplusplus
extern "C" {
#endif

/* buffer size needed for "str" in qlc_toDate and qlc_toInterval */
#define qlc_strBufferLen 20

char*       qlc_toDate(long date, char* str);             /* ex: 18-Jan-2006 */
void        qlc_pBenchmarkDate(FILE *f, const char *tag, long d);

char*       qlc_toInterval(long date1, long date2, char* str); /* ex: 1M, 5Y */
const char* qlc_freq2mat(char freq);   /* ex: ESL_FREQ_SEMI_ANNUAL 'S' -> 6M */
const char *qlc_freq2matStr(const char *freq);

const char* qlc_dcc(char dcc); /* ex: ESL_DCC_30_360 '3' -> B30360 */

void qlc_pString(FILE *f, const char *tag, const char *val);
void qlc_pBool(FILE *f, const char *tag, int val);
void qlc_pInt(FILE *f, const char *tag, int val);
void qlc_pDouble(FILE *f, const char *tag, double val);
void qlc_pDoubleArray(FILE *f, const char *tag, int len, double *val);
void qlc_pDoubleArrayFac(FILE *f, const char *tag, int len, double *val, double factor);
void qlc_pStringArray(FILE *f, const char *tag, int len, const char **array);
void qlc_pDateArray(FILE *f, const char *tag, int len, long *array);
void qlc_pBenchmarkDateArray(FILE *f, const char *tag, int len, long *array);
void qlc_pDoubleMatrix(FILE *f, const char *tag, int nbRows, int nbCols, double *vals);
void qlc_pDcc(FILE *f, const char *tag, char dcc);
void qlc_pDate(FILE *f, const char *tag, long d, int EOD=0);
void qlc_pFlexDates(FILE *f, const char *tag, int nbDates, long *dates);

void qlc_pSchedule(FILE *f, const char *tag, int nb, long *dates, double *val);

void qlc_printMarketBegin(FILE *f, long date, int incValDatePricingEvents);
void qlc_printMarketEnd(FILE *f);

void qlc_printRiskMgrBegin(FILE *f);
void qlc_printRiskMgrEnd(FILE *f);

void qlc_printInstrumentBegin(FILE *f);
void qlc_printInstrumentEnd(FILE *f);

void qlc_printControl(FILE *f);
void qlc_printControl2(FILE *f, int nbOutputRequest, const char **outputRequest);

#ifdef  __cplusplus
}
/* begin of C++ section */
#include <string>

namespace qlc {
	
inline std::string toString(int i) {
    char buf[20];
    sprintf(buf,"%d",i);
    return buf;
}

inline std::string toInterval(long date1, long date2) { /* ex: 1M, 5Y */
    char buf[qlc_strBufferLen];
    qlc_toInterval(date1, date2, buf);
    return buf;
}

inline void pString(FILE *f, const std::string& tag, const std::string& val) {
    qlc_pString(f, tag.c_str(), val.c_str());
}
inline void pBenchmarkDate(FILE *f, const std::string& tag, long val) {
    qlc_pBenchmarkDate(f, tag.c_str(), val);
}
inline void pDouble(FILE *f, const std::string& tag, double val) {
    qlc_pDouble(f, tag.c_str(), val);
}
inline void pDate(FILE *f, const std::string& tag, long val) {
    qlc_pDate(f, tag.c_str(), val);
}


class BaseCounter {
private:
	static int counter;
	int ref;
    std::string localTag;
protected:
    bool printRefOk(FILE *f, const std::string& tag, const std::string& type);
    void printClose(FILE *f);
public:
    virtual void print(FILE *f, const std::string& tag) = 0;
	BaseCounter() : ref(0) {}
    virtual ~BaseCounter() {}
};

struct SimpleHols : public BaseCounter {
	std::string name;
	SimpleHols(const std::string& name) : name(name) {}
    void print(FILE *f, const std::string& tag);
};

struct SimpleMetric : public BaseCounter {
	std::string name;
	std::string hols;
	SimpleMetric(const std::string& name, const std::string& hols) : name(name), hols(hols) {}
	void print(FILE *f, const std::string& tag);
};

struct SwapVol : public BaseCounter {
	std::string name;
	SimpleMetric *metric;
	SimpleHols   *hols;
	SWAPVOL_DATA *sv;
	SwapVol(const std::string& name, SimpleMetric *metric, SimpleHols *hols, SWAPVOL_DATA *sv) 
	: name(name), metric(metric), hols(hols), sv(sv) {}
	SwapVol() {}
	void print(FILE *f, const std::string& tag);
};

struct BaseVol : public BaseCounter {
	std::string name;
	SimpleMetric *metric;
	SimpleHols   *hols;
	BASEVOL_DATA *bv;
    long baseDate;
	BaseVol(const std::string& name, SimpleMetric *metric, SimpleHols *hols, 
        BASEVOL_DATA *bv, long baseDate) 
	: name(name), metric(metric), hols(hols), bv(bv), baseDate(baseDate) {}
	BaseVol() {}
	void print(FILE *f, const std::string& tag);
};

struct VolPair : public BaseCounter {
	std::string name;
	BaseVol    *bv;
	SwapVol    *sv;
	VolPair(const std::string& name, BaseVol *bv, SwapVol *sv)
	: name(name), bv(bv), sv(sv) {}
	VolPair() {}
	void print(FILE *f, const std::string& tag);
};

struct ZeroCurve : public BaseCounter {
	std::string name;
	T_CURVE     *tc;
	ZeroCurve(const std::string& name, T_CURVE *tc) : name(name), tc(tc) {}
	void print(FILE *f, const std::string& tag);
};

struct UntweakableYC : public BaseCounter {
	std::string name;
	std::string ccy;
	T_CURVE    *tc;
	VolPair    *vp;
	UntweakableYC(const std::string& name, const std::string&ccy, T_CURVE *tc, VolPair *vp) : name(name), ccy(ccy), tc(tc), vp(vp) {}
	UntweakableYC() {}
	void print(FILE *f, const std::string& tag);
};  

struct IrCalib : public BaseCounter {
	std::string name;
	std::string smile;
	std::string model;
	IrCalib(const std::string& name, const std::string& smile, const std::string& model) 
	: name(name), smile(smile), model(model) {}
	void print(FILE *f, const std::string& tag);
};

struct IrCalibSmile2Q : public BaseCounter {
	std::string name;
	std::string style;
	double qLeft;
	double qRight;
	double fwdShift;
    int cetNbIter;
	IrCalibSmile2Q(const std::string& name, const std::string& style, double qLeft, double qRight, 
                   double fwdShift, int cetNbIter) 
	: name(name), style(style), qLeft(qLeft), qRight(qRight), fwdShift(fwdShift), cetNbIter(cetNbIter) {}
	void print(FILE *f, const std::string& tag);		
};

struct IrCalibModelXFL : public BaseCounter {
	std::string name;
	std::string style;
	double meanReversion;
	double weight;
    int ppy;
    double nbStdDevs;
    int nbFactors;
	IrCalibModelXFL(const std::string& name, const std::string& style, double meanReversion, 
	    double weight, int ppy, double nbStdDevs, int nbFactors) 
	: name(name), style(style), meanReversion(meanReversion), weight(weight), ppy(ppy), 
	    nbStdDevs(nbStdDevs), nbFactors(nbFactors) {}
	void print(FILE *f, const std::string& tag);		
};

struct Correlation : public BaseCounter {
	std::string name;
	std::string asset1;
	std::string asset2; 
	double     correl;
	Correlation(const std::string& name, const std::string& asset1, const std::string& asset2, double correl) 
	: name(name), asset1(asset1), asset2(asset2), correl(correl) {}
	void print(FILE *f, const std::string& tag);		
};

struct CcyIRParams : public BaseCounter {
	std::string volCalibIndex;
	std::string irCalibName; 
	std::string smileSet; 
	std::string modelSet; 
	std::string curveToDiffuse; 
	std::string curveToDiscount; 
	CcyIRParams(const std::string& volCalibIndex, const std::string& irCalibName, const std::string& smileSet,
	            const std::string& modelSet, const std::string& curveToDiffuse, const std::string& curveToDiscount) 
	: volCalibIndex(volCalibIndex), irCalibName(irCalibName), smileSet(smileSet), 
	  modelSet(modelSet), curveToDiffuse(curveToDiffuse), curveToDiscount(curveToDiscount) {}

	void print(FILE *f, const std::string& tag);		
};
	
struct IndexSpecIR : public BaseCounter {
	std::string name;
    std::string factor;
	std::string tenor;// ESL frequency or '3M' style
    std::string freq; // ESL frequency or '3M' style
    char dcc; // ESL dcc
    std::string history;
    std::string zeroBankMethod;
    int nbExplicitZeroDates;
    long *explicitZeroDates;
    int nbZerosAtOneTime;

    void setTenorFreqDcc(const std::string &tenor, T_CURVE const *t_curve);
	IndexSpecIR() : nbExplicitZeroDates(0), explicitZeroDates(0), nbZerosAtOneTime(-1) {}
	void print(FILE *f, const std::string& tag);		
};

struct IndexSpecEQ : public BaseCounter {
    IndexSpecEQ() {}
	void print(FILE *f, const std::string& tag);		
};

struct AssetHistory : public BaseCounter {
    std::string name;
    std::string assetName;
    int nb;
    long *dates;
    double *values;
    AssetHistory() : nb(0), dates(0), values(0) {}
    void print(FILE *f, const std::string& tag);		
};

struct OptionSchedDates : public BaseCounter {
    long *notifDate;
    long *exerciseDate;
    int nbDates;
	OptionSchedDates() : notifDate(0), exerciseDate(0) {}
	void print(FILE *f, const std::string& tag);		
};

struct KOption : public BaseCounter {
    const char *longOrShort;
    const char* optionType;
	const char* exerType;
	const char* smoothing;
    OptionSchedDates *sched;
    double *strikes;
    std::string discount;

	KOption() : longOrShort("OPTION_LONG"), optionType(0), exerType(0), smoothing(0), sched(0) {}
	void print(FILE *f, const std::string& tag);		
};

struct CouponSchedDates : public BaseCounter {
    long *accStart;
    long *accEnd;
    long *reset;
    long *resetEff;
    long *pay;
    int nbDates;
	CouponSchedDates() : accStart(0), accEnd(0), reset(0), resetEff(0), pay(0), nbDates(0) {}
	void print(FILE *f, const std::string& tag);		
};

struct KFloatLeg : public BaseCounter {
    CouponSchedDates *sched;
    bool payInitialPrincipal;
    bool payPrincipal;
    BaseCounter *index;

    // for each coupon - begin
    double *notionals;
    double *weights;
    double *spreads;    
    double *dcfs;
    // for each coupon - end

    const char* rateType;    // simple or continuous compounding rate
    char dcc; // day count convention

    void print(FILE *f, const std::string& tag);

    KFloatLeg() : sched(0), payInitialPrincipal(false), payPrincipal(false), index(0), 
        notionals(0), weights(0), spreads(0), dcfs(0), rateType(0), dcc(0) {}
};

struct KOEvent : public BaseCounter {
    BaseCounter *und;
    long *barDates;    /* lo, hi */
    double *barVal[2]; /* lo, hi */
    double idxWeight;
    int nb;

    void print(FILE *f, const std::string& tag);
    KOEvent() : und(0), idxWeight(1.), nb(0) {
        for (int i=0; i<2; ++i) {
            barDates=0;
            barVal[i]=0;
        }
    }
};

struct KSum : public BaseCounter {
    int nb;
    BaseCounter **listK;
    double *weights;
    double constant;
    bool recordOutputName;
    void print(FILE *f, const std::string& tag);
    KSum(double constant=0.) : nb(0), listK(0), weights(0), constant(constant), recordOutputName(true) {}
};


struct KKnockOut : public BaseCounter {
    KOEvent *koEv[2];
    std::string smoothing;
    std::string barrierType;
    BaseCounter *upUnd;
    BaseCounter *midUnd;
    BaseCounter *downUnd;

    void print(FILE *f, const std::string& tag);
    KKnockOut() : upUnd(0), midUnd(0), downUnd(0) 
    {
        koEv[0]=0; 
        koEv[1]=0;
    }
};

struct KRib : public BaseCounter {
    BaseCounter *obsUnd;
    long *obsDates;
    int nbObs;
    long *obsEffDates;

    void print(FILE *f, const std::string& tag);
    KRib() : obsUnd(0), obsDates(0), nbObs(0), obsEffDates(0) {}
};

struct KRib2 : public BaseCounter {
    BaseCounter *obsUnd;
    long *obsDates;
    int nbObs;

    void print(FILE *f, const std::string& tag);
    KRib2() : obsUnd(0), obsDates(0), nbObs(0) {}
};

struct KRibFloatLeg : public BaseCounter {
    CouponSchedDates *sched;
    double *notionals;
    double *spreads;
    double *weights;
    double *dcfs;
    double *capRates;
    double *floorRates;
    std::string rateType;
    bool payInitialPrincipal;
    bool payPrincipal;
    char dcc;
    BaseCounter *index;
    KRib *rib;

    void print(FILE *f, const std::string& tag);
    KRibFloatLeg() : notionals(0), spreads(0), weights(0), dcfs(0), 
        payInitialPrincipal(false), payPrincipal(false), dcc(0), index(0), rib(0) {}
};

struct KRibTurbo_IndexLeg : public BaseCounter {
    int nbDates;
    double *notionals;
    double *spreads;
    double *weights;
    bool payInitialPrincipal;
    bool payPrincipal;
    std::string rateType;
    char dcc;
    double *dcfs;
    BaseCounter *index;
    std::string discount;

    void print(FILE *f, const std::string& tag);
    KRibTurbo_IndexLeg() : nbDates(0), notionals(0), spreads(0), weights(0), 
        payInitialPrincipal(false), payPrincipal(false), dcc(0), 
        dcfs(0), index(0) {}
};

struct KRibTurbo : public BaseCounter {
    KRibTurbo_IndexLeg *legs;
    int nbLegs;
    CouponSchedDates *sched;
    KRib2 *rib;
    double *floorRates;
    double *capRates;
    double *capFloorNotionals;
    char capFloorDcc;
    double *capFloorDcfs;

    std::string stubType;
    std::string discount;

    void print(FILE *f, const std::string& tag);
    KRibTurbo() : legs(0), nbLegs(0), sched(0), rib(0), floorRates(0), 
        capRates(0), capFloorNotionals(0), capFloorDcc(0), capFloorDcfs(0) {}
};

void printOneCcyMarket(
	FILE *f,
	T_CURVE t_curve[3],/* (I) Structure of zero curve data */
	char   *baseVolFilename,
	char   *swapVolFilename,
    MKTVOL_DATA *mktvol_data,
    const char *ccyStr,
    int ppy,
    int nbSigmaMax,
    int nbFactors,
    UntweakableYC *yc /* I/O */
    );


} /* end of namespace qlc */

/* end of C++ section */
#endif

#endif
