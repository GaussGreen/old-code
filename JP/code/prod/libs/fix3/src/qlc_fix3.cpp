#include "qlc_fix3.h"
#include "esl_market.h"
#include <string.h>
#include <stdio.h>

using namespace qlc;
using namespace std;

void qlc_printFix3Market(
	FILE *f,
	T_CURVE t_curve[3],/* (I) Structure of zero curve data */
	char   *baseVolFilename,
	char   *swapVolFilename,
	FIX3_TREE_DATA *tree_data,
    MKTVOL_DATA *mktvol_data
) {
    UntweakableYC yc[3];
    printOneCcyMarket(f, t_curve, baseVolFilename, swapVolFilename, 
        mktvol_data, "DOM", tree_data->Ppy, tree_data->NbSigmaMax, tree_data->NbFactor, yc);
}

void qlc_printFix3Param(
    FILE *f,
    char *volCalibIndex,
    FIX3_TREE_DATA *tree_data,
    char *termPRNFileName,
    int cetSmoothing,
    int useFix3CB,
    int incValDatePricingEvents)
{
    const char *modelName = (useFix3CB ? "Fix3CB" : "Fix3");
    fprintf(f, "<model TYPE='%s'>\n", modelName);
    qlc_pInt(f, "nbFactors", tree_data->NbFactor);
    CcyIRParams ccyIRParams(volCalibIndex, "DOM", "local", "local", 
        "DOM"+toString(tree_data->CvDiff), "DOM"+toString(tree_data->CvDisc));
    ccyIRParams.print(f, "IRParams");
    if (termPRNFileName) qlc_pString(f, "treeDataDebugFile", termPRNFileName);
    qlc_pString(f, "zeroInterpStyle", "FLAT_FWD");
    qlc_pBool(f, "cetSmoothing", cetSmoothing);
    qlc_pBool(f, "incValDatePricingEvents", incValDatePricingEvents);
    fprintf(f, "</model>\n");
}
