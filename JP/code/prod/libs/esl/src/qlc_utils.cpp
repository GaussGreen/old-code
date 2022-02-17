#include "qlc_utils.h"
#include "esl_market.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>

#ifdef ESL_NEW_CURVE
#include "irx/zerocurve.h"
#endif

using namespace std;
namespace qlc {

/*******************************************************/
int BaseCounter::counter=0;

bool BaseCounter::printRefOk(FILE *f, const string& tag, const std::string& type) {
    localTag = tag;
    if (ref) {
	    fprintf(f,"<%s REF='%d' />\n",tag.c_str(),ref);
	    return true;
    }
	ref=++counter; 
    fprintf(f,"<%s TYPE='%s' REF='%d'>\n", tag.c_str(), type.c_str(), ref);
	return false; 
}

void BaseCounter::printClose(FILE *f) {
    fprintf(f,"</%s>\n", localTag.c_str());
}

/*******************************************************/

void SimpleHols::print(FILE *f, const string& tag) {
    if (printRefOk(f,tag, "Holiday")) return;
    pString(f, "name", name);
    qlc_pBool(f, "useWeekends", true);
    qlc_pDateArray(f, "holidays", 0, 0);
    printClose(f);
}

/*******************************************************/

void SimpleMetric::print(FILE *f, const string& tag) {
    if (printRefOk(f,tag, "TimeMetric")) return;
    pDouble(f, "nonTradTimeFrac", 1);
    pString(f, "marketHols", hols);
	printClose(f);
 }

/*******************************************************/

void SwapVol::print(FILE *f, const string& tag) {
	int i,j,rows,cols;
	if (printRefOk(f,tag, "IRVol")) return;
    pString(f, "name", name);
	metric->print(f, "metric");

	rows = sv->NbSwaptionExpiries;
	cols = sv->NbSwapTenors;
	fprintf(f,"<matrix TYPE='DoubleMatrix' cols='%d' rows='%d'>\n",cols, rows);
	for (i=0; i<cols; ++i) {
		fprintf(f,"<Column%d>\n",i);
		for (j=0; j<rows; ++j) fprintf(f,"%.20f\n",sv->VolMatrix[j][i]);
		fprintf(f,"</Column%d>\n",i);
	}
	fprintf(f,"</matrix>\n");
	
	fprintf(f,"<expiry TYPE='ExpiryArray' length='%d'>\n",rows);
	for (i=0; i<rows; ++i) {
		fprintf(f,"<Item%d TYPE='MaturityPeriod'>%ldM</Item%d>\n",i,sv->SwaptionExpiries[i],i);
	}
	fprintf(f,"</expiry>\n");
	fprintf(f,"<spotOffset TYPE='Int'>0</spotOffset>\n");
	
	fprintf(f,"<tenor TYPE='ExpiryArray' length='%d'>\n",cols);
	for (i=0; i<cols; ++i) {
		fprintf(f,"<Item%d TYPE='MaturityPeriod'>%ldM</Item%d>\n",i,sv->SwapTenors[i],i);
	}
	fprintf(f,"</tenor>\n");
    qlc_pDate(f, "baseDate", 0);
    pString(f, "params", "");
    pString(f, "hols", hols->name);
	fprintf(f,"<bdc TYPE='BadDayNone'/>\n");
	fprintf(f,"<swapDayCount TYPE='Null'/>\n");
	fprintf(f,"<swapFrequency TYPE='Null'/>\n");
	printClose(f);
}

/*******************************************************/

void BaseVol::print(FILE *f, const string& tag) {
	int i;

	if (printRefOk(f,tag, "IRVol")) return;
    pString(f, "name", name);
	metric->print(f, "metric");

    int volShift=0;
    while (volShift<bv->NbVols-1 && bv->VolDates[volShift]<baseDate) ++volShift;

    fprintf(f,"<matrix TYPE='DoubleMatrix' cols='1' rows='%d'>\n",bv->NbVols-volShift);
	fprintf(f,"<Column0>\n");
	for (i=0; i<bv->NbVols-volShift; ++i) fprintf(f,"%.20f\n",bv->Vols[i+volShift]);
	fprintf(f,"</Column0>\n");
	fprintf(f,"</matrix>\n");
	
	fprintf(f,"<expiry TYPE='ExpiryArray' length='%d'>\n",bv->NbVols-volShift);
	for (i=0; i<bv->NbVols-volShift; ++i) {
        pBenchmarkDate(f, "Item"+toString(i), bv->VolDates[i+volShift]);
	}
	fprintf(f,"</expiry>\n");
	fprintf(f,"<spotOffset TYPE='Int'>0</spotOffset>\n"
	          "<tenor TYPE='ExpiryArray' length='1'>\n"
	          "<Item0 TYPE='MaturityPeriod'>%s</Item0>\n",qlc_freq2mat(bv->Frequency));
	fprintf(f,"</tenor>\n");
    qlc_pDate(f, "baseDate", 0);
    pString(f, "params", "");
    pString(f, "hols", hols->name);
	fprintf(f,"<bdc TYPE='BadDayNone'/>\n"
              "<swapDayCount TYPE='Null'/>\n"
              "<swapFrequency TYPE='Null'/>\n");
    printClose(f);
}

/*******************************************************/

void ZeroCurve::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag, "FourPlusIZeroCurve")) return;
#ifdef ESL_NEW_CURVE
    int i;  // VC98 does not handle scoping of variables defined in for loops correctly.

    qlc_pDate(f, "valueDate", YMDDateFromIRDate(tc->baseDate));

    fprintf(f, "<dates TYPE='DateTimeArray' length='%d'>\n", tc->numItems-1);
    for (i = 0; i < tc->numItems; i++) {
        if (tc->startDates[i] == tc->baseDate)
            continue;
        pDate(f, "Item"+toString(i), YMDDateFromIRDate(tc->startDates[i]));
    }
    fprintf(f, "</dates>\n");

    fprintf(f, "<rates TYPE='DoubleArray' length='%d'>\n", tc->numItems-1);
    for (i = 0; i < tc->numItems; i++) {
        if (tc->startDates[i] == tc->baseDate)
            continue;

        double rate;
        if (irxZeroRate(tc,
                        tc->startDates[i],
                        //strcmp(tc->SwapDCC, "360") == 0 ? IRX_ACT_360 : strcmp(tc->SwapDCC, "ACT") == 0 ? IRX_ACT_ACT : IRX_ACT_365F,
                        IRX_ACT_365F,
                        IRX_ANNUAL_RATE,
                        &rate) != SUCCESS)
            break;
        
        pDouble(f, "Item"+toString(i), rate);
    }
    fprintf(f, "</rates>\n");

    qlc_pInt(f, "basis", 1);  // Annually compounded rates ?
	fprintf(f,"<dayCountConv TYPE='Actual365'/>\n");
#else
    qlc_pDate(f, "valueDate", tc->Today);
    qlc_pDateArray(f, "dates", tc->NbZero, tc->ZeroDate);
    qlc_pDoubleArray(f, "rates", tc->NbZero, tc->Zero);
    qlc_pInt(f, "basis", (tc->SwapFreq == 'A' ? 1 : 2)); // not used by QLib IrConverter anyway
	fprintf(f,"<dayCountConv TYPE='%s'/>\n",
		(!strcmp(tc->SwapDCC,"ACT") ? "B30360":
         !strcmp(tc->SwapDCC,"360") ? "Actual360"
         : "Actual365"));
#endif
	printClose(f);
}

/*******************************************************/

void VolPair::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag, "IRVolPair")) return;
    pString(f, "name", name);
	bv->print(f,"baseVol");
	sv->print(f,"swapVol");
	printClose(f);
}
		   
/*******************************************************/

void UntweakableYC::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag, "UntweakableYC")) return;
    pString(f, "name", name);
    pString(f, "ccy", ccy);
    qlc_pDate(f, "today", tc->Today);
    qlc_pInt(f, "spotOffset", 0);
	fprintf(f,"<moneyMarketDayCount TYPE='%s'/>\n", (tc->MMB==360 ? "B30360" : "Actual365"));
	fprintf(f,"<swapDayCount TYPE='%s'/>\n",
      (!strcmp(tc->SwapDCC,"ACT") ? "B30360":
       !strcmp(tc->SwapDCC,"360") ? "Actual360"
       : "Actual365"));
    qlc_pInt(f, "swapFrequency", (tc->SwapFreq=='A'?1:2));
	ZeroCurve zc(ccy, tc);
	zc.print(f,"zeroCurve");
    qlc_pDate(f, "valueDate", tc->Today);
	if (vp) vp->print(f, "irVol");
	printClose(f);
}

/*******************************************************/

void IrCalib::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag, "IRCalib")) return;
    pString(f, "name", name);
    pString(f, "smileCalib", smile);
    pString(f, "modelCalib", model);
	printClose(f);
}

/*******************************************************/

void IrCalibSmile2Q::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag, "IRCalib::Smile2Q")) return;
    pString(f, "name", name);
    const char *p = style.c_str();
    qlc_pStringArray(f, "style", 1, &p);

	fprintf(f,"<params TYPE='DoubleArrayArray' length='1'>\n");
    double styleParams[]={1.-qLeft, 1.-qRight, fwdShift, cetNbIter};
    qlc_pDoubleArray(f, "Item0", 4, styleParams);
	fprintf(f,"</params>\n");
	printClose(f);
}

/*******************************************************/

void IrCalibModelXFL::print(FILE *f, const string& tag) {
    const char *(paramLabel[]) = {"meanReversion", "factorWeight", "ppy", "nbStdDevs", 
        "backBone", "nbStateVariables" ,"stateVarStdDevs" };
    double param[] = {meanReversion, weight, ppy, nbStdDevs,
        0, 0, 0 };

    const char *typeStr[]={"IRCalib::Model1FL", "IRCalib::Model2FL", "IRCalib::Model3FL"};
    if (nbFactors<0 || nbFactors>3) {
        fprintf(f, "************ error nbFactors=%d\n",nbFactors);
        return;
    }
    if (printRefOk(f,tag,typeStr[nbFactors-1])) return;
    pString(f, "name", name);
    const char *p = style.c_str();
    qlc_pStringArray(f, "style", 1, &p);
    qlc_pStringArray(f, "paramLabel", 7, paramLabel);
	fprintf(f,"<params TYPE='DoubleArrayArray' length='1'>\n");
    qlc_pDoubleArray(f, "Item0", 7, param);
    fprintf(f,"</params>\n");
	printClose(f);
}

/*******************************************************/

void Correlation::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"Correlation")) return;
    pString(f, "name", name);
    pString(f, "asset1", asset1);
    pString(f, "asset2", asset2);
	fprintf(f,"<correlation TYPE='Double'>%.20f</correlation>\n",correl);
	fprintf(f,"<corrExpiry TYPE='MaturityTimePeriod'>3Y SOD</corrExpiry>\n");
	printClose(f);
}

/*******************************************************/

void CcyIRParams::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"RateTree::CcyIRParams")) return;
    pString(f, "volCalibIndex", volCalibIndex);
    pString(f, "irCalib", irCalibName);
    pString(f, "smileSet", smileSet);
    pString(f, "modelSet", modelSet);
    pString(f, "curveToDiffuse", curveToDiffuse);
    pString(f, "curveToDiscount", curveToDiscount);
	printClose(f);
}

/*******************************************************/

void IndexSpecIR::setTenorFreqDcc(const string &indexName, T_CURVE const *t_curve) {
    char LiborDCC = (t_curve->MMB == 360 ? '0' : '5');
    char SwapDCC;
            
    if (!strcmp (t_curve->SwapDCC, "ACT"))
        SwapDCC = '3';
    else if (!strcmp (t_curve->SwapDCC, "360"))
        SwapDCC = '0';
    else 
        SwapDCC = '5';

    assert(indexName.size()!=0);
    int n = atoi(indexName.c_str());
    char l = indexName[indexName.size()-1];
    if (l=='m') {
        tenor = toString(n)+"M";
        freq = tenor;
        dcc = LiborDCC;
    }
    if (l=='y') {
        tenor = toString(n)+"Y";
        freq = qlc_freq2mat(t_curve->SwapFreq);
        dcc = SwapDCC;
    }
}

void IndexSpecIR::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"IndexSpecIR")) return;
    pString(f, "name", name);
    pString(f, "factor", factor);
	fprintf(f, "<tenor TYPE='MaturityPeriod'>%s</tenor>\n", qlc_freq2matStr(tenor.c_str()));
	fprintf(f, "<frequency TYPE='MaturityPeriod'>%s</frequency>\n", qlc_freq2matStr(freq.c_str()));
    qlc_pDcc(f, "dcc", dcc);
	if (!history.empty()) pString(f, "history", history);
    if (zeroBankMethod.size()) pString(f, "zeroBankMethod", zeroBankMethod);
    if (nbExplicitZeroDates) qlc_pDateArray(f, "explicitZeroDates", nbExplicitZeroDates, explicitZeroDates);
    qlc_pInt(f, "nbZerosAtOneTime", nbZerosAtOneTime);
    printClose(f);
}

void AssetHistory::print(FILE *f, const string& tag) {
    if (printRefOk(f,tag,"AssetHistory")) return;
    pString(f, "name", name);
    pString(f, "assetName", assetName);
    qlc_pDateArray(f, "dates", nb, dates);
    qlc_pDoubleArray(f, "values", nb, values);
    printClose(f);
}

void IndexSpecEQ::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"IndexSpecEQ")) return;
    pString(f, "name", "JPM");
    pString(f, "asset", "JPM");
	printClose(f);
}

void OptionSchedDates::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"OptionSchedDates")) return;
    fprintf(f, "<initializer TYPE='Null' />\n"); 
    qlc_pFlexDates(f, "notifDate", nbDates, notifDate);
    qlc_pFlexDates(f, "exerciseDate", nbDates, exerciseDate);
    printClose(f);
}    

void CouponSchedDates::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"CouponSchedDates")) return;
    qlc_pFlexDates(f, "accStart", nbDates, accStart);
    qlc_pFlexDates(f, "accEnd", nbDates, accEnd);
    qlc_pFlexDates(f, "reset", nbDates, reset);
    qlc_pFlexDates(f, "resetEff", nbDates, resetEff);
    qlc_pFlexDates(f, "pay", nbDates, pay);
    printClose(f);
}    

void KOption::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KOption")) return;
    fprintf(f, "<und TYPE='Null' />\n");
    pString(f, "longOrShort", longOrShort);
    pString(f, "optionType", optionType);
    pString(f, "exerType", exerType);
    pString(f, "smoothing", smoothing);
    sched->print(f ,"sched");
    qlc_pDoubleArray(f, "strikes", sched->nbDates, strikes);
    // American and interpolation styles not supported at this stage - todo
    if (discount.size()) pString(f, "discount", discount.c_str());
    printClose(f);
}    

void KFloatLeg::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KFloatLeg")) return;
    sched->print(f, "sched");
    qlc_pBool(f, "payInitialPrincipal", payInitialPrincipal);
    qlc_pBool(f, "payPrincipal", payPrincipal);
    if (index) index->print(f, "index");
    qlc_pDoubleArray(f, "notionals", sched->nbDates, notionals);
    if (weights) qlc_pDoubleArray(f, "weights", sched->nbDates, weights);
    if (weights) qlc_pDoubleArray(f, "spreads", sched->nbDates, spreads);
    if (dcfs) qlc_pDoubleArray(f, "dcfs", sched->nbDates, dcfs);
    pString(f, "rateType", rateType);
    qlc_pDcc(f, "dcc", dcc);
    printClose(f);
}

void KRibFloatLeg::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KRibFloatLeg")) return;
    sched->print(f, "sched");
    qlc_pBool(f, "payInitialPrincipal", payInitialPrincipal);
    qlc_pBool(f, "payPrincipal", payPrincipal);
    index->print(f, "index");
    qlc_pDoubleArray(f, "notionals", sched->nbDates, notionals);
    qlc_pDoubleArray(f, "weights", sched->nbDates, weights);
    qlc_pDoubleArray(f, "spreads", sched->nbDates, spreads);
    if (capRates) qlc_pDoubleArray(f, "capRates", sched->nbDates, capRates);
    if (floorRates) qlc_pDoubleArray(f, "floorRates", sched->nbDates, floorRates);


    if (dcfs) qlc_pDoubleArray(f, "dcfs", sched->nbDates, dcfs);
    pString(f, "rateType", rateType);
    if (dcc) qlc_pDcc(f, "dcc", dcc);
    rib->print(f, "rib");
    printClose(f);
}

void KRibTurbo_IndexLeg::print(FILE *f, const string& tag) {
    if (printRefOk(f,tag,"KRibTurbo::IndexLeg")) return;
    qlc_pDoubleArray(f, "notionals", nbDates, notionals);
    if (spreads) qlc_pDoubleArray(f, "spreads", nbDates, spreads);
    if (weights) qlc_pDoubleArray(f, "weights", nbDates, weights);
    qlc_pBool(f, "payInitialPrincipal", payInitialPrincipal);
    qlc_pBool(f, "payPrincipal", payPrincipal);
    pString(f, "rateType", rateType);
    if (dcfs) qlc_pDoubleArray(f, "dcfs", nbDates, dcfs);
    else qlc_pDcc(f, "dcc", dcc);

    fprintf(f,"<indexes TYPE='IProdCreatorArray' length='%d'>\n",1);
    for (int i=0; i<1; ++i) {
        index->print(f, "Item"+toString(i));
    }
    fprintf(f,"</indexes>\n");
    if (discount.size()) pString(f, "discount", discount);
    printClose(f);
}

void KRibTurbo::print(FILE *f, const string& tag) {
    if (printRefOk(f,tag,"KRibTurbo")) return;

    fprintf(f,"<legs TYPE='KRibTurbo::IndexLegArray' length='%d'>\n",nbLegs);
    for (int i=0; i<nbLegs; ++i) {
        legs[i].print(f, "Item"+toString(i));
    }
    fprintf(f,"</legs>\n");
    sched->print(f, "sched");
    if (rib) rib->print(f, "rib");
    if (floorRates) qlc_pDoubleArray(f, "floorRates", sched->nbDates, floorRates);
    if (capRates) qlc_pDoubleArray(f, "capRates", sched->nbDates, capRates);
    if (capFloorNotionals) qlc_pDoubleArray(f, "capFloorNotionals", sched->nbDates, capFloorNotionals);
    if (capFloorDcfs) qlc_pDoubleArray(f, "capFloorDcfs", sched->nbDates, capFloorDcfs);
    else if (capFloorNotionals) qlc_pDcc(f, "capFloorDcc", capFloorDcc);
    pString(f, "stubType", stubType);
    if (discount.size()) pString(f, "discount", discount);
    printClose(f);
}

void KRib::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KRib")) return;
    obsUnd->print(f, "obsUnd");;
    qlc_pFlexDates(f, "obsDates", nbObs, obsDates);
    qlc_pFlexDates(f, "obsEffDates", nbObs, obsEffDates);
    printClose(f);
}

void KRib2::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KRib2")) return;
    fprintf(f,"<obsUnds TYPE='IProdCreatorArray' length='%d'>\n",1);
    for (int i=0; i<1; ++i) {
        obsUnd->print(f, "Item"+toString(i));;
    }
    fprintf(f,"</obsUnds>\n");
    qlc_pFlexDates(f, "obsDates", nbObs, obsDates);
    printClose(f);
}

void KSum::print(FILE *f, const string& tag) {
    if (printRefOk(f,tag,"KSum")) return;
    fprintf(f, "<listK TYPE='IProdCreatorArray' length='%d'>\n", nb);
    for (int i=0; i<nb; ++i) {
        listK[i]->print(f, "Item"+toString(i));
    }
    fprintf(f, "</listK>\n");
    qlc_pDoubleArray(f, "weights", nb, weights);
    pDouble(f, "constant", constant);
    qlc_pBool(f, "recordOutputName", recordOutputName);
    printClose(f);
}

void KKnockOut::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KKnockOut")) return;
    fprintf(f, "<knockouts TYPE='KKnockOut::KOEventArray' length='%d'>\n", (koEv[1] ? 2 : 1));
    koEv[0]->print(f, "Item0");
    if (koEv[1]) koEv[1]->print(f, "Item1");
    fprintf(f, "</knockouts>\n");
    pString(f, "smoothing", smoothing);
    if (barrierType.size()) pString(f, "barrierType", barrierType);
    upUnd->print(f, "upUnd");
    midUnd->print(f, "midUnd");
    downUnd->print(f, "downUnd");
    printClose(f);
}

void KOEvent::print(FILE *f, const string& tag) {
	if (printRefOk(f,tag,"KKnockOut::KOEvent")) return;
    und->print(f, "und");
    qlc_pSchedule(f, "loBarrier", nb, barDates, barVal[0]);
    qlc_pSchedule(f, "hiBarrier", nb, barDates, barVal[1]);
    pDouble(f, "idxWeight", idxWeight);
    printClose(f);
}

void printOneCcyMarket(
	FILE *f,
	T_CURVE *t_curve,/* (I) Structure of zero curve data */
	char   *baseVolFilename,
	char   *swapVolFilename,
    MKTVOL_DATA *mktvol_data,
    const char *ccyStr,
    int ppy,
    int nbSigmaMax,
    int nbFactors,
    UntweakableYC yc[3] /* I/O */
    )
{
    SWAPVOL_DATA svData;
	BASEVOL_DATA bvData;
	if (EslReadVolsW(&bvData, &svData, baseVolFilename, swapVolFilename)!=SUCCESS) {
	        fprintf(f,"************************** ERROR, qlc_printOneCcyMarket\n");
	        return;
	}

    string ccy(ccyStr);
	SimpleHols hols("SIMPLE_"+ccy);
	SimpleMetric metric("SIMPLE_"+ccy, "SIMPLE_"+ccy);
	BaseVol bv(ccy+"-BASEVOL", &metric, &hols, &bvData, t_curve->ValueDate);
	SwapVol sv(ccy+"-SWAPVOL", &metric, &hols, &svData);
	VolPair vp(ccy, &bv, &sv);

	hols.print(f, "entry");
	vp.print(f, "entry");
    for (int j=0; j<3; ++j) { /* curve 0,1,2 */
		yc[j] = UntweakableYC(ccy+toString(j), ccy, &t_curve[j], &vp);
		yc[j].print(f, "entry");
	}

	IrCalib(ccy, ccy, ccy).print(f, "entry");
	IrCalibSmile2Q(ccy,"local",mktvol_data->QLeft,mktvol_data->QRight,mktvol_data->FwdShift,
        mktvol_data->CetNbIter).print(f, "entry");
	IrCalibModelXFL(ccy,"local",mktvol_data->Beta[0],
        mktvol_data->Alpha[0], ppy, nbSigmaMax, nbFactors).print(f, "entry");
}

} /* end of namespace qlc */

using namespace qlc;

/*******************************************************/

char *qlc_toDate(long date, char* str) {
    int day, mon, year;
    char month[][10] = { "???", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

    day = date % 100;
    mon = date % 10000 / 100;
    year = date / 10000;

    sprintf(str, "%02d-%s-%d", day, month[mon], year);
    return str;
}

/*******************************************************/

char *qlc_toInterval(long date1, long date2, char* str) {
	int d,m,y;

	d = (date2 % 100) - (date1 % 100);
	m = (date2 % 10000 / 100) - (date1 % 10000 / 100);
	y = (date2 / 10000) - (date1 / 10000);
	if (d<0) { d+=30; --m; }
	if (d>=30) { d-=30; ++m; }
	if (m<0) { m+=12; --y; }
	if (m>=12) { m-=12; ++y; }
	if (y && m==0) { sprintf(str,"%dY",y); return str; } else m+=12*y;
	if (m && d==0) { sprintf(str,"%dM",m); return str; } else d+=30*m;
	if (d%7==0) { sprintf(str,"%dW",d/7); return str; }
	sprintf(str,"%dD",d);
	return str;
}

/*******************************************************/

const char *qlc_freq2matStr(const char *freq) {
    if (!freq || !*freq) return 0;
    if (strlen(freq)>1) return freq;
    return qlc_freq2mat(freq[0]);
}

const char *qlc_freq2mat(char freq) {
	switch (freq) {
		case ESL_FREQ_ANNUAL: return "12M";
		case ESL_FREQ_SEMI_ANNUAL: return "6M";
		case ESL_FREQ_QUARTERLY: return "3M";
		case ESL_FREQ_MONTHLY: return "1M";
		case ESL_FREQ_WEEKLY: return "1W";
		case ESL_FREQ_DAILY: return "1D";
		default: return "?unknown";
	}
}
		
const char *qlc_dcc(char dcc) {
    switch (dcc) {
        case ESL_DCC_30_360 /* '3' */: return "B30360";
        case ESL_DCC_ACT_360 /* '0' */: return "Actual360";
        case ESL_DCC_ACT_365 /* '5' */: return "Actual365";
        case ESL_DCC_ACT_ACT /* 'A' */: return "ActualActual";
        default: return "?unknown";
    }
}


void qlc_pString(FILE *f, const char *tag, const char *val) {
    fprintf(f, "<%s TYPE='String'>%s</%s>\n",tag, val, tag);
}

void qlc_pBool(FILE *f, const char *tag, int val) {
    fprintf(f, "<%s TYPE='String'>%s</%s>\n",tag, (val ? "Y" : "N"), tag);
}

void qlc_pInt(FILE *f, const char *tag, int val) {
    fprintf(f, "<%s TYPE='Int'>%d</%s>\n",tag, val, tag);
}

void qlc_pDouble(FILE *f, const char *tag, double val) {
    fprintf(f, "<%s TYPE='Double'>%30.18f</%s>\n",tag, val, tag);
}

void qlc_pDoubleArrayFac(FILE *f, const char *tag, int len, double *val, double factor) {
    fprintf(f, "<%s TYPE='DoubleArray' length='%d'>\n", tag, len);
    for (int i = 0; i < len; i++) {
        pDouble(f, "Item"+toString(i), val[i]*factor);
    }
    fprintf(f, "</%s>\n", tag);
}

void qlc_pDoubleArray(FILE *f, const char *tag, int len, double *val) {
    qlc_pDoubleArrayFac(f, tag, len, val, 1.0);
}

void qlc_pStringArray(FILE *f, const char *tag, int len, const char **array) {
    fprintf(f, "<%s TYPE='StringArray' length='%d'>\n", tag, len);
    for (int i = 0; i < len; i++) {
        pString(f, "Item"+toString(i), array[i]);
    }
    fprintf(f, "</%s>\n", tag);
}

void qlc_pDateArray(FILE *f, const char *tag, int len, long *array) {
    fprintf(f, "<%s TYPE='DateTimeArray' length='%d'>\n", tag, len);
    for (int i = 0; i < len; i++) {
        pDate(f, "Item"+toString(i), array[i]);
    }
    fprintf(f, "</%s>\n", tag);
}

void qlc_pBenchmarkDateArray(FILE *f, const char *tag, int len, long *array) {
    fprintf(f, "<%s TYPE='ExpiryArray' length='%d'>\n", tag, len);
    for (int i = 0; i < len; i++) {
        pBenchmarkDate(f, "Item"+toString(i), array[i]);
    }
    fprintf(f, "</%s>\n", tag);
}

void qlc_pDoubleMatrix(FILE *f, const char *tag, int nbRows, int nbCols, double *vals) {
    fprintf(f, "<%s TYPE='DoubleMatrix' cols='%d' rows='%d'>\n", tag, nbCols, nbRows);
    for (int j=0; j<nbCols; ++j) {
        fprintf(f, "<Column%d>", j);
        for (int i=0; i<nbRows; ++i) {
            fprintf(f, "%30.18f ", vals[i*nbCols+j]);
        }
        fprintf(f, "</Column%d>", j);
    }
    fprintf(f, "</%s>", tag);
}

void qlc_pDcc(FILE *f, const char *tag, char dcc) {
    fprintf(f,"<%s TYPE='%s' />\n", tag, qlc_dcc(dcc));
}

void qlc_pDate(FILE *f, const char *tag, long d, int EOD) {
    char dateStr[qlc_strBufferLen];
    if (d==0) {
        fprintf(f, "<%s TYPE='DateTime'/>\n",tag);
        return;
    }
    fprintf(f, "<%s TYPE='DateTime'>%s %s</%s>\n",tag, qlc_toDate(d, dateStr), (EOD?"EOD":"SOD"), tag);
}

void qlc_pBenchmarkDate(FILE *f, const char *tag, long d) {
    char dateStr[qlc_strBufferLen];
    if (d==0) {
        fprintf(f, "<%s TYPE='BenchmarkDate'/>\n",tag);
        return;
    }
    fprintf(f, "<%s TYPE='BenchmarkDate'>%s SOD</%s>\n",tag, qlc_toDate(d, dateStr), tag);
}

void qlc_pFlexDates(FILE *f, const char *tag, int nbDates, long *dates) {
    fprintf(f, "<%s TYPE='FlexDates'>\n",tag);
    qlc_pDateArray(f, "datesExpl", nbDates, dates);
    fprintf(f, "<modifier TYPE='Null' /> \n");
    fprintf(f, "</%s>\n",tag);
}

void qlc_pSchedule(FILE *f, const char *tag, int nb, long *dates, double *val) {
    fprintf(f, "<%s TYPE='Schedule'>\n", tag);
    qlc_pDateArray(f, "dates", nb, dates);
    qlc_pDoubleArray(f, "values", nb, val);
    pString(f, "interp", "L");
    fprintf(f, "</%s>\n", tag);
}

/*******************************************************/

void qlc_printMarketBegin(FILE *f, long date, int incValDatePricingEvents) {
	fprintf(f,"<market TYPE='MarketData'>\n");
    qlc_pDate(f, "RefDate", date, !incValDatePricingEvents);
	fprintf(f,"<DataCache>\n");
}

void qlc_printMarketEnd(FILE *f) {
	fprintf(f,"</DataCache>\n");
    fprintf(f,"</market>\n");
}

void qlc_printRiskMgrBegin(FILE *f) {
	fprintf(f,
"<!-- WARNING: The format of this XML does NOT constitute a public interface.\n"
"Quantitative Research & Development reserve the right to change the format at will and do\n"
"not support any applications that attempt to parse, interpret or in any way\n"
"manipulate this XML other than via publicly supported QR&D interfaces.\n"
"Anyone assuming the style, format or content of this XML does so entirely at their own risk. -->\n"
"<!-- QLIB Version -Rates QLC library- -->\n"
"<RISK-MGR TYPE='MultiRiskMgrInterface'>\n"
	);
}

void qlc_printRiskMgrEnd(FILE *f) {
	fprintf(f,"</RISK-MGR>\n");
}

void qlc_printInstrumentBegin(FILE *f)
{
	fprintf(f,
"<insts TYPE='ArrayInstrumentCollection'>\n"
"<instruments TYPE='InstrumentArray' length='1'>\n"
	);
}

void qlc_printInstrumentEnd(FILE *f)
{
	fprintf(f,
"</instruments>\n"
"</insts>\n"
	);
}

void qlc_printControl2(FILE *f, int nbOutputRequest, const char **outputRequest)
{
    int i;
	fprintf(f,"<ctrl TYPE='Control'>\n");
    fprintf(f,"<sens TYPE='SensitivityArray' length='0' />\n");
    fprintf(f,"<outputRequests TYPE='OutputRequestArray' length='%d'>\n", nbOutputRequest);
    for (i=0; i<nbOutputRequest; ++i) {
        fprintf(f,"<Item%d TYPE='OutputRequest'>\n", i);
        pString(f, "requestName", outputRequest[i]);
	    fprintf(f,"</Item%d>\n", i);
    }
    fprintf(f,"</outputRequests>\n");
    qlc_pBool(f, "writeToFile", false);
    pString(f, "fileName", "input.xml");
    qlc_pBool(f, "scaleResults", false);
    pString(f, "assetPriceSource", "ASSET_THEO");
    fprintf(f,"</ctrl>\n");
}

void qlc_printControl(FILE *f) {
    const char *array[]={"UNDERLYING", "FWD_UNDERLYING", "FLOAT_LEG_VALUE", "FIXED_LEG_VALUE", "OPTION_PRICE", 0};
    qlc_printControl2(f, 5, array);
}

/*******************************************************/
