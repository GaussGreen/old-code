/* =======================================================================
   FILE NAME: swp_h_cmt_capfloor.h
   ======================================================================= */

#ifndef SWP_H_CMT_CAPFLOOR_H
#define SWP_H_CMT_CAPFLOOR_H

Err swp_f_CMT_capfloor(
    long        start,
    long        end_or_nfp,
    String      cap_freq_str,
    String      cap_basis_str,
    double      strike,
    double      flatvol,
    String      vol_type_str,
    String      rec_pay_str,
    String      greek_str,
    SrtCurvePtr crvptr,
    double*     result);

Err swp_f_CMT_capfloor_impvol(
    long        start,
    long        end_or_nfp,
    String      cap_freq_str,
    String      cap_basis_str,
    double      strike,
    double      premium,
    String      rec_pay_str,
    SrtCurvePtr crvptr,
    double*     imp_vol);

Err swp_f_CMCapFloor(
    long        TheoEnd,
    long        numperiodcap,
    double      Underlying,
    char*       CMTFrequency,
    char*       CMTBasis,
    char*       CMTRefRateCode,
    double      Strike,
    char*       CapFrequency,
    char*       CapBasis,
    char*       szCMTVolName,
    char*       szSwaptionVolName,
    char*       szYieldCurveName,
    char*       cRefRateCode,
    char*       cCapFloor,
    char*       method,
    SRT_Boolean AdjForSpread,
    int         NumStrikesInVol,
    double*     Strikes,
    double**    result);

Err swp_f_CMCapFloorNY(
    long   TheoEnd,
    long   numperiodcap,
    double Underlying,
    char*  cCMTFrequency,
    char*  cCMTBasis,
    char*  cCMTRefRateCode,
    double Strike,
    char*  cCapFrequency,
    char*  cCapBasis,
    char*  szCMTVolName,
    char*  szSwaptionVolName,
    char*  szYieldCurveName,
    char*  cRefRateCode,
    char*  cCapFloor,
    char*
        method, /*BS: BS on CMS Full Smile + Fwd Spread - FS: Full Smile on Fwd Swap + Fwd Spread -
                   FVBS: BS on CMS Flat Vol + Fwd Spread - FV: Flat Vol on Fwd Swap + Fwd Spread*/
    SRT_Boolean AdjForSpread,
    int         NumStrikesInVol,
    double*     Strikes,
    double**    result);

Err swp_f_GetCMVols(
    char*    szYieldCurveName,
    char*    szSwaptionVolName,
    char*    cRefRateCode,
    double   Underlying,
    char*    cCMTFrequency,
    char*    cCMTBasis,
    char*    cCMTRefRateCode,
    int      NumStrikes,
    double*  Strikes,
    int      NumMats,
    double*  Mats,
    double*  alpha,
    double*  beta,
    double*  rho,
    char*    Method,
    double*  chi,
    double   dLambda,
    double   dSpreadLocVol,
    double   dCorrel,
    double** CMTVols,
    double** SabrComp);

#endif
