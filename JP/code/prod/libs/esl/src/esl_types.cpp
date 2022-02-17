#include "esl_types.h"


#define FTAG "%-#22.15g"

void Print_MKTVOL_DATA(FILE *f, MKTVOL_DATA *m) {
    int i,j;
    int NbVolUsed=0;

    if (!f) {
        fprintf(stderr,"Print_MKTVOL_DATA: file not ready\n");
        return;
    }
    fprintf(f,"Contents of MKTVOL_DATA\n\n");
    if (!m) {
        fprintf(f,"NULL\n");
        return;
    }
    fprintf(f,"BaseDate: %ld\n", m->BaseDate);
    fprintf(f,"NbCetVol: %d\n", m->NbCetVol);
    fprintf(f,"CetNbIter: %d\n", m->CetNbIter);
    fprintf(f,"CetVegaError: "FTAG"\n", m->CetVegaError);
    fprintf(f,"VolUnit: %d\n", m->VolUnit);
    if (m->NbVol) {
        fprintf(f, "\n  i  VolDate   Vol                   VolUsed  VolType  SwapSt   SwapMat\n");
        for (i=0; i<m->NbVol; ++i) {
            fprintf(f, "%3d  %8ld  " FTAG "  %5d  %7d  %ld %ld\n" , 
            i, m->VolDate[i], m->Vol[i], m->VolUsed[i], (int)m->VolType[i], m->SwapSt[i], m->SwapMat[i]);
            if (m->VolUsed[i])
                ++NbVolUsed;
        }
    }
    if (m->NbBmkMr) {
        fprintf(f, "\n   i  ");
        for (j=0; j<6; ++j) fprintf(f, "Aweight[%d]              ",j);
        fprintf(f, "\n");
        for (i=0; i<m->NbBmkMr; ++i) {
            fprintf(f, "%4d  " , i);
            for (j=0; j<6; ++j) {
                fprintf(f, FTAG "  ", (m->Aweight[j] ? m->Aweight[j][i] : -1e99));
            }
            fprintf(f,"\n");
        }
    }
    fprintf(f,"\n");
    fprintf(f,"Freq:         %c\n", m->Freq);
    fprintf(f,"DCC:          %c\n", m->DCC);
    fprintf(f,"SkipFlag:     %d\n", m->SkipFlag);
    fprintf(f,"CalibFlag:    %d\n", m->CalibFlag);
    fprintf(f,"FilterSpotVolFlag: %d\n", m->FilterSpotVolFlag);
    fprintf(f,"SmoothingFlag:%c\n", m->SmoothingFlag);
    fprintf(f,"TraceFlag:    %c\n", m->TraceFlag);
    fprintf(f,"ModelChoice:  %d\n", m->ModelChoice);
    fprintf(f,"NbFactor:     %d\n", m->NbFactor);
    fprintf(f,"QLeft:        "FTAG"\n", m->QLeft);
    fprintf(f,"QRight:       "FTAG"\n", m->QRight);
    fprintf(f,"FwdShift:     "FTAG"\n", m->FwdShift);
    fprintf(f,"Alpha[0,1,2]: "FTAG" "FTAG" "FTAG" \n", m->Alpha[0], m->Alpha[1], m->Alpha[2]);
    fprintf(f,"Beta[0,1,2]:  "FTAG" "FTAG" "FTAG" \n", m->Beta[0], m->Beta[1], m->Beta[2]);
    fprintf(f,"Rho[0,1,2]:   "FTAG" "FTAG" "FTAG" \n", m->Rho[0], m->Rho[1], m->Rho[2]);
    fprintf(f,"Bbq:          "FTAG"\n", m->Bbq);
    fprintf(f,"VolNorm:      "FTAG"\n", m->VolNorm);
    fprintf(f,"VolLogn:      "FTAG"\n", m->VolLogn);
    fprintf(f,"NbTDInp:      %d\n", m->NbTDInp);
    fprintf(f,"NbBmkMr:      %d\n", m->NbBmkMr);

    if (m->NbTDInp) {
        fprintf(f,"\n  i  TDInpDate BmkDate\n");
        for (i=0; i<m->NbTDInp; ++i) {
            fprintf(f, "%3d   %8ld %8ld\n", i, m->TDInpDate[i], m->BmkDate[i]);
        }
        fprintf(f,"\n  i  TDInpDate BetaTD[0]               BetaTD[1]               BetaTD[2]\n");
        for (i=0; i<m->NbTDInp; ++i) {
            fprintf(f, "%3d  %8ld  "FTAG"  "FTAG"  "FTAG"\n", i,m->TDInpDate[i], m->BetaTD[0][i], m->BetaTD[1][i], m->BetaTD[2][i]);
        }
        fprintf(f,"\n  i  TDInpDate BetaBmk[0]              BetaBmk[1]              BetaBmk[2]\n");
        for (i=0; i<m->NbTDInp; ++i) {
            fprintf(f, "%3d  %8ld  "FTAG"  "FTAG"  "FTAG"\n", i,m->TDInpDate[i], m->BetaBmk[0][i], m->BetaBmk[1][i], m->BetaBmk[2][i]);
        }
        fprintf(f,"\n  i  TDInpDate AlphaTD[0]              AlphaTD[1]              AlphaTD[2]\n");
        for (i=0; i<m->NbTDInp; ++i) {
            fprintf(f, "%3d  %8ld  "FTAG"  "FTAG"  "FTAG"\n", i,m->TDInpDate[i], m->AlphaTD[0][i], m->AlphaTD[1][i], m->AlphaTD[2][i]);
        }
        fprintf(f,"\n  i  TDInpDate RhoTD[0]                RhoTD[1]                RhoTD[2]\n");
        for (i=0; i<m->NbTDInp; ++i) {
            fprintf(f, "%3d  %8ld  "FTAG"  "FTAG"  "FTAG"\n", i,m->TDInpDate[i], m->RhoTD[0][i], m->RhoTD[1][i], m->RhoTD[2][i]);
        }
    }
    fprintf(f,"\n");
    fprintf(f,"  i  SmileDate QLeftTD                 QRightTD                FwdShiftTD\n");
    for (i=0; i<m->NbSmileDates; ++i) {
        fprintf(f,"%3d  %8ld  "FTAG"  "FTAG"  "FTAG"\n", i, m->SmileDate[i], m->QLeftTD[i], m->QRightTD[i], m->FwdShiftTD[i]);
    }
    fprintf(f,"\n");
    fprintf(f,"Amap:         "FTAG"\n", m->Amap);
    fprintf(f,"Bmap:         "FTAG"\n", m->Bmap);
    fprintf(f,"Afac:         "FTAG"\n", m->Afac);
    fprintf(f,"Bfac:         "FTAG"\n", m->Bfac);
    fprintf(f,"Cfac:         "FTAG"\n", m->Cfac);
    fprintf(f,"Dfac:         "FTAG"\n", m->Dfac);
    fprintf(f,"\n");
}

void Print_T_CURVE(FILE *f, T_CURVE *tc) {
    fprintf(f, "Print_T_CURVE\n\n");
#ifdef ESL_NEW_CURVE
    fprintf(f, "Type: _IrxTZeroCurve\n");
#else
    fprintf(f,"Type: _T_CURVE\n");
    fprintf(f,"Today: %ld\n", tc->Today);
    fprintf(f,"SpotDays: %d\n", tc->SpotDays);
    fprintf(f,"ValueDate: %ld\n", tc->ValueDate);
    fprintf(f,"SwapFreq: %c\n", tc->SwapFreq);
    fprintf(f,"SwapDCC: %s\n", tc->SwapDCC);
    fprintf(f,"MMB: %d\n", tc->MMB);
    fprintf(f,"NbZero: %d\n", tc->NbZero);
    fprintf(f,"  i ZeroDate Zero\n");
    {
        int i;
        for (i=0; i<tc->NbZero; ++i) {
            fprintf(f,"%3d %8ld " FTAG "\n", i, tc->ZeroDate[i], tc->Zero[i]);
        }
    }
    fprintf(f,"InterpType: %ld\n", tc->InterpType);
#endif
}

void Esl_CbkInit(CLAIM_BANK *CBK)      /* (I/O) the target claim bank   */
{
    if (CBK != NULL)
    {
        CBK->Slices = NULL;
        CBK->EvDates = NULL;
        CBK->ErDates = NULL;
        CBK->TotNbSlices = 0;
        CBK->NbActiveSlices = 0;
        CBK->MaxErDate = -1;
        CBK->Locked = FALSE;
        CBK->auxSlice = NULL;
    }
}
