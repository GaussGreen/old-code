#include "qlc_hyb3.h"
#include "hyb3_market.h"
#include "esl_market.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>

using namespace std;

namespace qlc {
	
void FxVol::print(FILE *f, const std::string &tag) {
	int i;
	char dateStr[qlc_strBufferLen];
	
	if (printRefOk(f,tag, "SRMFX::Vol")) return;
	
    qlc_pString(f, "name", name);
	/* compVol */
	fprintf(f,"<compVolExpiry TYPE='ExpiryArray' length='%d'>\n",fxVol->NbBaseVols);
	for (i=0; i<fxVol->NbBaseVols; ++i) {
		fprintf(f,"<Item%d TYPE='MaturityPeriod'>%s</Item%d>\n",i,qlc_toInterval(fxVol->ValueDate,fxVol->BaseVolDates[i],dateStr),i);
	}
	fprintf(f,"</compVolExpiry>\n");
	fprintf(f,"<compVolMatExpiry TYPE='Null'/>\n");
    qlc_pDoubleArrayFac(f, "compVol", fxVol->NbBaseVols, fxVol->BaseVols,1./100.);

	/* spotVol */
	fprintf(f,"<spotVolExpiry TYPE='Null'/>\n");
    qlc_pDoubleArray(f, "spotVol", 0, 0);
	
	/* smile */
	fprintf(f,"<smileExpiry TYPE='ExpiryArray' length='%d'>\n",fx->NbSmilePt);
	for (i=0; i<fx->NbSmilePt; ++i) {
		fprintf(f,"<Item%d TYPE='BenchmarkDate'>%s SOD</Item%d>\n",i,qlc_toDate(fx->SmileDate[i],dateStr),i);
	}
	fprintf(f,"</smileExpiry>\n");

    qlc_pDoubleArray(f, "smileA1", fx->NbSmilePt, fx->a1);
    qlc_pDoubleArray(f, "smileA2", fx->NbSmilePt, fx->a2);
    qlc_pDoubleArray(f, "smileA3", fx->NbSmilePt, fx->a3);

	printClose(f);
}

/*******************************************************/

void FxAsset::print(FILE *f, const std::string &tag) {
    if (printRefOk(f,tag, "FXAsset")) return;
    qlc_pString(f, "name", name);
    qlc_pDate(f, "today", 0);
    qlc_pString(f, "holidays", "");
	yc1->print(f, "baseCcy");
	yc2->print(f, "riskCcy");
    pDouble(f, "spotFX", fvh->fxVol->FXSpotRate);
	fvh->print(f,"fxVol");
    qlc_pBool(f, "coherentGreeks", false);
    /*if (nbPastFixings)
    {
        int i;
        char dateStr[MAXNBDATE];
        fprintf(f, "<history TYPE='AssetHistory'>\n");
        fprintf(f, "<name TYPE='String'>%s</name>\n", name);
        fprintf(f, "<assetName TYPE='String'>%s</assetName>\n", name);
        fprintf(f, "<dates TYPE='DateTimeArray' length='%d' REF='dates.DateTimeArray.5'>\n",
                nbPastFixings);
        for (i = 0; i < nbPastFixings; i++)
        {
            fprintf(f, "<Item%d TYPE='DateTime'>%s SOD</Item%d>\n",
                    i,
                    qlc_toDate(pastDate[i],dateStr),
                    i);
        }
        fprintf(f, "</dates>\n");
        fprintf(f, "<values TYPE='DoubleArray' length='%d' REF='values.DoubleArray.5'>\n",
                nbPastFixings);
        for (i = 0; i < nbPastFixings; i++)
        {
            fprintf(f, "<Item%d TYPE='Double'>%lf</Item%d>\n", 
                i,
                pastRate[i],
                i);
        }
        fprintf(f, "</values>\n");
        fprintf(f, "</history>\n");
    }*/
	printClose(f);
}

} /* end of namespace qlc */

using namespace qlc;

int qlc_convDivType(int t) {
    switch (t) {
        case 'D': return 0; // AMOUNT
        case 'Y': return 1; // PERCENT
        case 'C': return 2; // CONTINUOUS
        default: assert(0); return -1;
    }
}

void qlc_printSimpleEquity(FILE *f, const char *tag, EQ_DATA *eq_data, long valueDate) {
    fprintf(f, "<%s TYPE='SimpleEquity'>\n", tag);

    fprintf(f, "<equity TYPE='Equity'>\n");
    pString(f, "name", "JPM");
    pDouble(f, "stockPrice", eq_data->Spot);

    switch (eq_data->SettleType) {
        case 'R':
            fprintf(f,"<settlement TYPE='RollingSettlement'>\n");
            qlc_pInt(f, "period", eq_data->NbSettle);
            pString(f, "hols", "SIMPLE_DOM");
            fprintf(f,"</settlement>\n");
            break;

        case 'F':
            fprintf(f,"<settlement TYPE='FixedSettlement'>\n");
            qlc_pDateArray(f, "firstTradingDates", eq_data->NbSettle, eq_data->LastTrading);
            qlc_pDateArray(f, "settlementDates", eq_data->NbSettle, eq_data->SettleDate);
            pString(f, "marketHols", "SIMPLE_DOM");
            fprintf(f,"</settlement>\n");
            break;

        default: assert(0);
    }

    pString(f, "yc", "DOM1");
    qlc_pDate(f, "valueDate", valueDate /*!!!*/);
    qlc_pDate(f, "stockDate", valueDate /*!!!*/);
    fprintf(f, "<divList TYPE='DividendList'>\n");
    fprintf(f, "<divArray TYPE='DividendArray' length='%d'>\n", eq_data->NbFwd);
    for (int i=0; i<eq_data->NbFwd; ++i) {
        fprintf(f, "<Item%d TYPE='Dividend'>\n", i);
        qlc_pDate(f, "exDivDate", eq_data->FwdDate[i]);
        pDouble(f, "divAmount", eq_data->Fwd[i]);
        qlc_pInt(f, "divType", qlc_convDivType(eq_data->FwdType[i]));
        qlc_pDate(f, "payDivDate", eq_data->FwdDate[i]); // just because payDivDate is mandatory
        fprintf(f, "</Item%d>\n", i);
    }
    fprintf(f, "</divArray>\n");
    fprintf(f, "</divList>\n");

    fprintf(f, "<borrowCurve TYPE='BorrowCurve::Interface'>\n");
    qlc_pDate(f, "baseDate", valueDate /*!!!*/);
    int skipBorrow=0;
    while (skipBorrow<eq_data->NbBorrow-1 && eq_data->BorrowDate[skipBorrow]<valueDate) ++skipBorrow;
    qlc_pBenchmarkDateArray(f, "expiries", eq_data->NbBorrow-skipBorrow, eq_data->BorrowDate+skipBorrow);
    qlc_pDoubleArray(f, "rates", eq_data->NbBorrow-skipBorrow, eq_data->Borrow+skipBorrow);
    pString(f, "name", "JPM");
    qlc_pInt(f, "basis", 0);
    pString(f, "dayCount", "30/360");
    fprintf(f, "</borrowCurve>\n");
    fprintf(f, "</equity>\n");
    pString(f, "vol", "JPM");
    pString(f, "name", "JPM");
    fprintf(f, "</%s>\n", tag);
}

void qlc_pSRMEQVol(FILE *f, const char *tag, EQ_DATA *eq_data) {
    char dateStr[qlc_strBufferLen];
    fprintf(f, "<%s TYPE='SRMEQ::Vol'>\n", tag);
    pString(f, "name", "JPM");

    fprintf(f, "<compVolBMs TYPE='StringArray' length='%d'>\n", eq_data->NbVol);
    for (int i = 0; i < eq_data->NbVol; i++) {
        pString(f, "Item"+toString(i), string(qlc_toDate(eq_data->VolDate[i+1], dateStr))+" SOD");
        // toInterval(eq_data->ValueDate, eq_data->VolDate[i]));
    }
    fprintf(f, "</compVolBMs>\n");
    qlc_pDoubleArray(f, "compVol", eq_data->NbVol, eq_data->Vol+1);
    {
        const int dateLen=30;
        char *dates;
        char **datesPtr;
        dates = new char[dateLen*eq_data->NbSmilePt];
        datesPtr = new char*[eq_data->NbSmilePt];
        for (int i=0; i<eq_data->NbSmilePt; ++i) {
            datesPtr[i] = dates+dateLen*i;
            qlc_toDate(eq_data->SmileDate[i], datesPtr[i]);
            strcat(datesPtr[i], " SOD");
        }
        qlc_pStringArray(f, "smileBMs", eq_data->NbSmilePt, (const char **)datesPtr);
        delete[] datesPtr;
        delete[] dates;
    }
    qlc_pDoubleArray(f, "smileA1", eq_data->NbSmilePt, eq_data->a1);
    qlc_pDoubleArray(f, "smileA2", eq_data->NbSmilePt, eq_data->a2);
    qlc_pDoubleArray(f, "smileA3", eq_data->NbSmilePt, eq_data->a3);
    fprintf(f, "</%s>\n", tag);
}

int qlc_printHyb3Market(
	FILE *f,
	T_CURVE t_curve[2][3],/* (I) Structure of zero curve data */
    char   *baseVolFilenameFor,
    char   *swapVolFilenameFor,
    char   *baseVolFilenameDom,
    char   *swapVolFilenameDom,
    char   *FXVolFilename,
    FX_DATA *fx_data,
    MKTVOL_DATA *mktvol_data,
    int    nbPastFX,
    long   pastFXDates[],
    double pastFXRates[],
    HYB3_TREE_DATA *tree_data)
{
    UntweakableYC yc[2][3];

    printOneCcyMarket(f, t_curve[FOR], baseVolFilenameFor, swapVolFilenameFor, 
        &mktvol_data[FOR], "FOR", tree_data->Ppy, tree_data->NbSigmaMax, 1, yc[FOR]);

    printOneCcyMarket(f, t_curve[1], baseVolFilenameDom, swapVolFilenameDom, 
        &mktvol_data[DOM], "DOM", tree_data->Ppy, tree_data->NbSigmaMax, 1, yc[DOM]);

	FXVOLATILITY_DATA fxVolData;
	FXSMILE_DATA smileData;
    char blank[] = "";

    if (HYB3ReadFXEnvW(&fxVolData, &smileData, FXVolFilename, blank)!=SUCCESS) return FAILURE;

	FxVol fxVol("FOR_DOM", &fxVolData, fx_data);
	FxAsset fxAsset("DOMFOR", &yc[DOM][1], &yc[FOR][1], &fxVol);
    fxAsset.nbPastFixings = nbPastFX;
    fxAsset.pastDate = pastFXDates;
    fxAsset.pastRate = pastFXRates;

	fxVol.print(f,"entry");
	fxAsset.print(f, "entry");
	
	Correlation("FOR_DOMFOR","FOR1","DOMFOR", fx_data->Rho[1][0]).print(f, "entry");
	Correlation("DOM_DOMFOR","DOM1","DOMFOR", fx_data->Rho[2][0]).print(f, "entry");
	Correlation("FOR_DOM","FOR1","DOM1", fx_data->Rho[0][0]).print(f, "entry");
	
	return SUCCESS;
}


int qlc_printHyb3EQMarket(
	FILE *f,
	T_CURVE *t_curve,/* (I) Structure of zero curve data */
    char   *baseVolFilename,
    char   *swapVolFilename,
    MKTVOL_DATA *mktvol_data,
    HYB3_TREE_DATA *tree_data,
    EQ_DATA *eq_data)
{
    UntweakableYC yc[3];

    printOneCcyMarket(f, t_curve, baseVolFilename, swapVolFilename, 
        mktvol_data, "DOM", tree_data->Ppy, tree_data->NbSigmaMax, 1, yc);
    qlc_printSimpleEquity(f, "entry", eq_data, t_curve->ValueDate);
    qlc_pSRMEQVol(f, "entry", eq_data);

    Correlation eqIrCorr("DOM_JPM","DOM0", "JPM", eq_data->Rho[0][0]);
    eqIrCorr.print(f, "entry");

	return SUCCESS;
}

int qlc_printHyb3FXModel(
    FILE *f,
    HYB3_TREE_DATA *tree_data,
    const char *termPRNFileName)
{
    fprintf(f, "<model TYPE='Hyb3FX'>\n");
    qlc_pInt(f, "ppy", tree_data->Ppy);
    qlc_pInt(f, "maxStdDeviations", tree_data->NbSigmaMax);
    qlc_pString(f, "FXSmileParams", "");
    qlc_pString(f, "FXIRCorrelations", "");
    qlc_pString(f, "zeroInterpStyle", "LINEAR");
    qlc_pString(f, "productInfoFileName", "");
    CcyIRParams(tree_data->Index[FOR], "FOR", "local", "local", 
        "FOR"+toString(tree_data->CvDiff[FOR]), 
        "FOR"+toString(tree_data->CvDisc[FOR])).print(f,"forIRParams");
    CcyIRParams(tree_data->Index[DOM], "DOM", "local", "local", 
        "DOM"+toString(tree_data->CvDiff[DOM]), 
        "DOM"+toString(tree_data->CvDisc[DOM])).print(f,"domIRParams");
    if (termPRNFileName) qlc_pString(f, "treeDataDebugFile", termPRNFileName);
    fprintf(f, "</model>\n");

    return SUCCESS;
}

int qlc_printHyb3EQModel(
    FILE *f,
    HYB3_TREE_DATA *tree_data,
    const char *termPRNFileName, 
    long matDate,
    const char *eqSpotVolOverride
    )
{
    fprintf(f, "<model TYPE='Hyb3EQ'>\n");
    qlc_pInt(f, "ppy", tree_data->Ppy);
    qlc_pInt(f, "maxStdDeviations", tree_data->NbSigmaMax);
    qlc_pString(f, "zeroInterpStyle", "LINEAR");
    qlc_pString(f, "productInfoFileName", "");
    CcyIRParams(tree_data->Index[0], "DOM", "local", "local", 
        "DOM"+toString(tree_data->CvDiff[0]), "DOM"+toString(tree_data->CvDisc[0])).print(f,"IRParams");

    if (termPRNFileName) qlc_pString(f, "treeDataDebugFile", termPRNFileName);
    if (matDate) qlc_pDate(f,"lastDividendDate", matDate);
    qlc_pBool(f, "incValDatePricingEvents", true);
    if (strcmp(eqSpotVolOverride,"last")==0) {
        strcpy((char*)eqSpotVolOverride, "nil"); // does not matter for most tests
        // printf("eqSpotVolOverride==last not supported\n");
        // exit(0);
    }
    if (strcmp(eqSpotVolOverride,"nil")) 
        qlc_pDouble(f, "eqSpotVolOverride", atof(eqSpotVolOverride)/100.);

    fprintf(f, "</model>\n");

    return SUCCESS;
}
