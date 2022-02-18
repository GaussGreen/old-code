
#include "TargetNoteCaller.h"

#include "LGMSVClosedFormApprox.h"
#include "LGMSVUtil.h"
#include "TargetNote.h"
#include "TargetNoteCalib.h"
#include "TargetNoteProdStruct.h"
#include "TargetNoteSV.h"
#include "TargetNoteUtil.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srtaccess.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"

/* ----------------------------------------------------------------------------------------------------------------
 */
/* Calling routine for Target Note */

char* TargetNoteCaller(TARN_Struct* in_tarn)
{
    /* Variable decalaration */
    char*       err = 0;
    TARN_AUX    AUX, *aux      = &AUX;
    TARN_Struct tarnObj, *tarn = &tarnObj;

    /* copy the input struct */
    smessage("copy_TARN_Struct");
    copy_TARN_Struct(in_tarn, tarn);

    /* Reverse the eod flags because we got the convention reversed */
    tarn->market.eodFixFlag = !(tarn->market.eodFixFlag);
    tarn->market.eodPayFlag = !(tarn->market.eodPayFlag);

    /* initialize the auxilliary structure */
    smessage("init_TARN_AUX");
    if (err = init_TARN_AUX(tarn, aux))
        goto FREE_RETURN;

    /* calculate any historical coupons */
    smessage("TARN_calcHistory");
    if (err = TARN_calcHistory(tarn, aux))
        goto FREE_RETURN;

    /* if the deal is knocked out (but does not have a final payment)or has finished, record the
     * remaining PV and return */
    if ((tarn->output.iIsKnockedOut &&
         (tarn->deal.couponType != TN_FINAL || aux->i1stCpn == (tarn->deal.nCouponDates - 1))) ||
        tarn->output.iIsExpired)
    {
        tarn->output.dmLGM2F_PV[0][0]   = aux->dHistFundPV;
        tarn->output.dmLGM2F_PV[1][0]   = aux->dHistCpnPV;
        tarn->output.dmLGM1F_PV[0][0]   = aux->dHistFundPV;
        tarn->output.dmLGM1F_PV[1][0]   = aux->dHistCpnPV;
        tarn->output.dmLGM1FSV_PV[0][0] = aux->dHistFundPV;
        tarn->output.dmLGM1FSV_PV[1][0] = aux->dHistCpnPV;
        tarn->output.dmLGMSV_PV[0][0]   = aux->dHistFundPV;
        tarn->output.dmLGMSV_PV[1][0]   = aux->dHistCpnPV;
        goto FREE_RETURN;
    }

    /* Calibrate underlyings (if necessary) */
    smessage("TARN_calib");
    if (err = TARN_calib(tarn, aux))
        goto FREE_RETURN;

    /* calculate LGM2FCV price */
    if (tarn->pricing.iPrice2FCV)
    {
        smessage("2F");
        if (err =
                TargetNoteMC_2F(tarn->output.szTARN_LGM2F_UND, tarn, aux, tarn->output.dmLGM2F_PV))
            goto FREE_RETURN;

        smessage("1F");
        if (err =
                TargetNoteMC_2F(tarn->output.szTARN_LGM1F_UND, tarn, aux, tarn->output.dmLGM1F_PV))
            goto FREE_RETURN;

        smessage("1FSV");
        if (err = TargetNoteMC_SV(
                tarn->output.szTARN_LGM1FSV_UND, tarn, aux, tarn->output.dmLGM1FSV_PV))
            goto FREE_RETURN;
    }

    /* calculate LGMSV price */
    if (tarn->pricing.iPriceSV)
    {
        smessage("SV");
        if (err =
                TargetNoteMC_SV(tarn->output.szTARN_LGMSV_UND, tarn, aux, tarn->output.dmLGMSV_PV))
            goto FREE_RETURN;
    }

FREE_RETURN:
    /* copy the output struct */
    smessage("copy_Output");
    copy_TARN_Output_Struct(&tarn->output, &in_tarn->output);

    /* free the memory */
    smessage("free_TRAN_Struct");
    free_TARN_Struct(tarn);
    smessage("free_TRAN_Aux");
    free_TARN_AUX(aux);
    return err;
}

/* Calling routine for Target Note */
/* ----------------------------------------------------------------------------------------------------------------
 */
