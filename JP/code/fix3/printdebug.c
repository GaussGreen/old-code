#include "fix123head.h"

void Fix3_Print_FIX3_TREE_DATA(FILE *f, FIX3_TREE_DATA *t) {
#define FTAG "%-#22.15g"
    int i,j;
    if (!f) {
        fprintf(stderr,"Cannot open file to print FIX3_TREE_DATA\n");
        return;
    }
    fprintf(f,"Contents of FIX3_TREE_DATA\n\n");
    if (!t) {
        fprintf(f,"NULL\n");
        return;
    }
/*  // check relevance of the following lines
    fprintf(f, "\n i  Type  CritDate  CritType  NbZeros | SuppDate[0,1,2] | Value[0,1,2,3,4]\n");
    for (i=0; i<NBCRITDATE; ++i) {
        if (t->CritDate[i] == NULL) {
            fprintf(f, "\tNULL\n");
            continue;
        }
        fprintf(f, "%2d  %4d  %8ld         %c  %7d ", i, t->CritDate[i]->Type, t->CritDate[i]->CritDate,
            t->CritType[i], t->NbZeros[i]);
        fprintf(f, "| ");
        for (j=0; j<3; ++j) {
            fprintf(f, "%ld ", t->CritDate[i]->SuppDate[j]);
        }
        fprintf(f, "          | ");
        for (j=0; j<5; ++j) {
            fprintf(f, FTAG, t->CritDate[i]->Value[j]);
        }
        fprintf(f,"\n");
        fflush(f);
    }
    fprintf(f, "\n   i    TPDate | critical dates #\n");
    for (i=0; i<t->NbTP; ++i) {
        for (j=0; j<NBCRITDATE; ++j) {
            if (t->TPtype[j][i])
                break;
        }
        if (j==NBCRITDATE)
            continue;
        fprintf(f, "%4d  %8ld | ", i, t->TPDate[i]);
        for (j=0; j<NBCRITDATE; ++j) {
            if (t->TPtype[j][i]) {
                fprintf(f, "%d ",j);
            }
        }
        fprintf(f, "\n");
    }
*/
    if (t->NbEDevDates) {
        fprintf(f, "   i  EDevDate\n");
        for (i=0; i<t->NbEDevDates; ++i) {
            fprintf(f, "%4d %9ld\n", i, t->EDevDate[i]);
        }
    }
    else fprintf(f, "\n");

    fprintf(f, "NbTP: %d\n\n",t->NbTP);

    for (j=0; j<3; ++j) {
        fprintf(f, "   i   TPDate Length                 ");
        fprintf(f, "ZeroCoupon%d            ZeroRate%d              FwdRate%d               ", j, j, j);
        fprintf(f,"\n");
        fflush(f);
        for (i=0; i<t->NbTP; ++i) {
            fprintf(f, "%4d %8ld " FTAG " ", i, t->TPDate[i], t->Length[i]);
            fprintf(f, FTAG " " FTAG " " FTAG " ", t->ZeroCoupon[j][i], t->ZeroRate[j][i], t->FwdRate[j][i]);
            fprintf(f, "\n");
        }
        fprintf(f,"\n");
    }
    fprintf(f, "\nPpy: %d\n", t->Ppy);
    for (i=0; i<3; ++i) {
        fprintf(f, "PpyCet[%d]: %d\n",i, t->PpyCet[i]);
    }
    fprintf(f, "JumpPpy: %d\n", t->JumpPpy);
    fprintf(f, "NbDailyPts: %d\n", t->NbDailyPts);
    fprintf(f, "\n");
    fprintf(f, "CvDiff: %d\n", t->CvDiff);
    fprintf(f, "CvIdx1: %d\n", t->CvIdx1);
    fprintf(f, "CvIdx2: %d\n", t->CvIdx2);
    fprintf(f, "CvDisc: %d\n", t->CvDisc);
    fprintf(f, "\n");
    fprintf(f, "NbFactor: %d\n", t->NbFactor);

    fprintf(f, "\n   i   TPDate Length                 ");
    for (j=0; j<6; ++j) fprintf(f, "Aweight[%d]              ",j);
    fprintf(f, "\n");
    for (i=0; i<t->NbTP; ++i) {
        fprintf(f, "%4d %8ld " FTAG " ", i, t->TPDate[i], t->Length[i]);
        for (j=0; j<6; ++j) {
            fprintf(f, FTAG "  ", (t->Aweight[j] ? t->Aweight[j][i] : -1e99));
        }
        fprintf(f,"\n");
    }
    fprintf(f,"\n");
    fprintf(f,"NbSigmaMax: %d\n", t->NbSigmaMax);
    fprintf(f,"\n  i Width HalfWidth\n");
    for (i=0; i<3; ++i)
        fprintf(f,"%3d %5d %9d\n", i, t->Width[i], t->HalfWidth[i]);

    fprintf(f, "\n   i   TPDate ZCenter                LengthJ                  Top1  Bottom1  OutTop1  OutBottom1\n");
    for (i=0; i<t->NbTP; ++i) {
        fprintf(f, "%4d %8ld " FTAG " " FTAG " %6d   %6d   %6d      %6d\n", i, t->TPDate[i], t->ZCenter[i], t->LengthJ[i], t->Top1[i], t->Bottom1[i], t->OutTop1[i], t->OutBottom1[i]);
    }
    fprintf(f, "\n");
}
