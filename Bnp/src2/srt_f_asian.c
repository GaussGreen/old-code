/******************************************************************************/
/*                                                                            */
/*       MODULE      :	SRT_F_TMEOPTION.C
 */
/*		 AUTHOR		 :	EF (modif. from asian equity of CY)                   */
/*       CREATED     :	April 1998
 */
/*                                                                            */
/******************************************************************************/

#include "opfnctns.h"
#include "srt_h_all.h"

/* -------------------------------------------------------------------------- */
/* solve the system in (F,d,Sig) knowing (M1,M2,M3) s.t.
                                F+d=M1
                                F^2*(exp(sig^2*T)-1)=M2-M1^2
                                F^3*(exp(3*sig^2*T)-1)+3*(M2-M1^2)*(M1-F)=M3-M1^3 */

//void find_d(double* ans, double spot, double A1, double A2, double A3)
//{
//    double aux, auxder, v2, v3, x, y, forward_, disp_, tau_;
//    int    i;
//
//    v2 = A2 - A1 * A1;
//    v3 = A3 - 3.0 * v2 * A1 - A1 * A1 * A1;
//    y  = A2 / A1;
//
//    for (i = 0; i < 10; i++)
//    {
//        aux    = y * y * y + 3.0 * v2 * y - v3;
//        auxder = 3.0 * (y * y + v2);
//        y      = y - aux / auxder;
//    }
//
//    x        = 1.0 + y * y / v2;
//    tau_     = sqrt((log(1.0 + y * y / v2)));
//    forward_ = y / (x - 1.0);
//    disp_    = A1 - forward_;
//
//    ans[0] = disp_ * spot;
//    ans[1] = forward_ * spot;
//    ans[2] = tau_;
//}

/* -------------------------------------------------------------------------- */
/* Asian price */
//double eval_asian(
//    int            ndate,
//    double*        date,
//    double*        fwds,
//    double*        vols,
//    double         spot,
//    double         maturity,
//    double         strike,
//    int            npast_fix,
//    double         avg_cur,
//    double         disc,
//    SrtCallPutType call_put)
//{
//    long    i;
//    double *forwards, *fforwards, *ffforwards;
//    double  Mom1, Mom2, Mom3, extra, aux1, aux2;
//    double  ffMom1, fMom1, fMom2, sffMom1, sfMom1, sfMom2, wf;
//    double  ans, *sans;
//
//    /*----------------------------*/
//    /* computing the first 3 moments */
//    forwards   = (double*)calloc(ndate, sizeof(double));
//    fforwards  = (double*)calloc(ndate, sizeof(double));
//    ffforwards = (double*)calloc(ndate, sizeof(double));
//    sans       = (double*)calloc(4, sizeof(double));
//
//    Mom1 = fMom1 = ffMom1 = 0.0;
//    Mom2 = fMom2 = 0.0;
//    Mom3         = 0.0;
//    for (i = 0; i < ndate; i++)
//    {
//        forwards[i]   = fwds[i] / (npast_fix + ndate) / spot;
//        fforwards[i]  = forwards[i] * exp(vols[i] * vols[i] * date[i]);
//        ffforwards[i] = forwards[i] * exp(2 * vols[i] * vols[i] * date[i]);
//        sfMom1        = fMom1;
//        sffMom1       = ffMom1;
//        sfMom2        = fMom2;
//        Mom1 += forwards[i];
//        fMom1 += fforwards[i];
//        ffMom1 += ffforwards[i];
//        Mom2 += 2.0 * sfMom1 * forwards[i];
//        Mom2 += forwards[i] * fforwards[i];
//        fMom2 += 2.0 * sffMom1 * fforwards[i];
//        fMom2 += fforwards[i] * ffforwards[i];
//        Mom3 += 3.0 * sfMom2 * forwards[i];
//        Mom3 += 3.0 * sffMom1 * forwards[i] * fforwards[i];
//        wf = fforwards[i];
//        Mom3 += (wf * wf * wf);
//    }
//    aux1  = Mom1;
//    aux2  = Mom2;
//    extra = (avg_cur * npast_fix) / (npast_fix + ndate) / spot;
//    Mom1 += extra;
//    Mom2 += extra * (2.0 * aux1 + extra);
//    Mom3 += extra * (3.0 * aux2 + extra * (extra + 3.0 * aux1));
//
//    /* solve the system in (F,d,Sig) knowing (M1,M2,M3) */
//    find_d(sans, spot, Mom1, Mom2, Mom3);
//    /* compute the price of the asian */
//    ans = srt_f_optblksch(sans[1], strike - sans[0], sans[2], maturity, disc, call_put, PREMIUM);
//
//    free(sans);
//    free(forwards);
//    free(fforwards);
//    free(ffforwards);
//
//    return (ans);
//}

/* -------------------------------------------------------------------------- */

//Err srt_f_asian(
//    int            nforwards,
//    double*        forwards_date,
//    double*        forwards,
//    int            nfixing,
//    double*        fixing_date,
//    int            nvol,
//    double*        vol_date,
//    double*        volat,
//    double         spot,
//    double         maturity,
//    double         strike,
//    int            npast_fix,
//    double         avg_cur,
//    double         disc,
//    SrtCallPutType call_put,
//    SrtGreekType   greek,
//    double*        answer)
//{
//    Err     err = NULL;
//    double  eps = 0.01, x;
//    double *fwds, *vols;
//    int     i, nf, nv;
//
//    /* interpolating data on fixing dates */
//    fwds = (double*)calloc(nfixing, sizeof(double));
//    vols = (double*)calloc(nfixing, sizeof(double));
//
//    nf = nv = 0;
//    for (i = 0; i < nfixing; i++)
//    {
//        while ((nf < nforwards) && (forwards_date[nf] < fixing_date[i]))
//            nf++;
//        if (nf == 0)
//            fwds[i] = forwards[0];
//        else if (nf == nforwards)
//            fwds[i] = forwards[nforwards - 1];
//        else
//        {
//            x = (fixing_date[i] - forwards_date[nf - 1]) /
//                (forwards_date[nf] - forwards_date[nf - 1]);
//            fwds[i] = x * forwards[nf] + (1.0 - x) * forwards[nf - 1];
//        }
//
//        while ((nv < nvol) && (vol_date[nv] < fixing_date[i]))
//            nv++;
//        if (nv == 0)
//            vols[i] = volat[0];
//        else if (nf == nforwards)
//            vols[i] = volat[nforwards - 1];
//        else
//        {
//            x       = (fixing_date[i] - vol_date[nv - 1]) / (vol_date[nv] - vol_date[nv - 1]);
//            vols[i] = x * volat[nv] + (1.0 - x) * volat[nv - 1];
//        }
//    }
//
//    /*------------------------------------------*/
//    switch (greek)
//    {
//    case PREMIUM:
//
//        *answer = eval_asian(
//            nfixing,
//            fixing_date,
//            fwds,
//            vols,
//            spot,
//            maturity,
//            strike,
//            npast_fix,
//            avg_cur,
//            disc,
//            call_put);
//        break;
//
//    default:
//        *answer = 0.0;
//        break;
//    }
//
//    free(fwds);
//    free(vols);
//
//    return (err);
//}

/*--------------------------------- End of File -------------------------------------*/