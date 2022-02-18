#include "AmortMidatADIRecst.h"

//////////////////////////////////////////////////////////////////////
///  Precompute constants for discount factors,  given E(R1) and E(R2)
/// DF(t,T) = exp(*pdConstant - X1 * *pdCoeff1 - X2 * *pdCoeff2)
static void _df_constant_(
    double dLAM1_0_t,
    double dLAM2_0_t,
    double dLAM1_0_T,
    double dLAM2_0_T,
    double dEILam1_t,
    double dEILam2_t,
    double dPHI1_t,
    double dPHI2_t,
    double dPHICross_t,
    double dExp_R1,
    double dExp_R2,
    double dRho,
    // results
    double* pdConstant,
    double* pdCoeff1,
    double* pdCoeff2)
{
    /// LAMBDA
    const double dLAM1_t_T = (dLAM1_0_T - dLAM1_0_t) * dEILam1_t;
    const double dLAM2_t_T = (dLAM2_0_T - dLAM2_0_t) * dEILam2_t;

    // XI
    const double dXI1Sqrd     = dLAM1_t_T * dLAM1_t_T * dPHI1_t;
    const double dXI2Sqrd     = dLAM2_t_T * dLAM2_t_T * dPHI2_t;
    const double dXICrossSqrd = dLAM1_t_T * dLAM2_t_T * dPHICross_t;

    /// Bond Variance
    const double dHalfBondVar = 0.5 * (dXI1Sqrd + dXI2Sqrd) + dRho * dXICrossSqrd;

    // assign
    *pdCoeff1   = dLAM1_t_T;
    *pdCoeff2   = dLAM2_t_T;
    *pdConstant = -dHalfBondVar - *pdCoeff1 * dExp_R1 - *pdCoeff2 * dExp_R2;

    _ASSERTE(dLAM1_t_T > 0. && dLAM2_t_T > 0.);
    _ASSERTE(
        (dPHI1_t > 0. && dPHI2_t > 0. && dPHICross_t > 0.) ||     // regular 2f case
        (dPHI1_t == 0. && dPHI2_t == 0. && dPHICross_t == 0.) ||  // initial point at today
        (dPHI1_t > 0. && dPHI2_t == 0. &&
         dPHICross_t == 0.)  // 1f case with 2nd dimension suppressed
    );
}

void _df_constant(
    double dDF_0_t,
    double dDF_0_T,
    double dLAM1_0_t,
    double dLAM2_0_t,
    double dLAM1_0_T,
    double dLAM2_0_T,
    double dEILam1_t,
    double dEILam2_t,
    double dPHI1_t,
    double dPHI2_t,
    double dPHICross_t,
    double dExp_R1,
    double dExp_R2,
    double dRho,
    // results
    double* pdConstant,
    double* pdCoeff1,
    double* pdCoeff2)
{
    _ASSERTE(dDF_0_t >= dDF_0_T);
    if (dDF_0_t == dDF_0_T)
    {
        *pdConstant = 0.;
        *pdCoeff1   = 0.;
        *pdCoeff2   = 0.;
    }
    else
    {
        _df_constant_(
            dLAM1_0_t,
            dLAM2_0_t,
            dLAM1_0_T,
            dLAM2_0_T,
            dEILam1_t,
            dEILam2_t,
            dPHI1_t,
            dPHI2_t,
            dPHICross_t,
            dExp_R1,
            dExp_R2,
            dRho,
            pdConstant,
            pdCoeff1,
            pdCoeff2);
    }
}

static double __DF_nojump_(
    double          dX1,
    double          dX3,
    double          dDF_0_t,
    double          dDF_0_T,
    double          dLAM1_0_t,
    double          dLAM2_0_t,
    double          dLAM1_0_T,
    double          dLAM2_0_T,
    double          dEILam1_t,
    double          dEILam2_t,
    double          dPHI1_t,
    double          dPHI2_t,
    double          dPHICross_t,
    const _PCQ_Seq* pNumeraire,
    double          dRho,
    double          dAlphaRho,
    double          dAlpha_Sqrt1mRhoSqrd)
{
    /// LAMBDA
    const double dLAM1_t_Nume = (*pNumeraire->m_pdLAMZero1_Begin - dLAM1_0_t) * dEILam1_t;
    const double dLAM2_t_Nume = (*pNumeraire->m_pdLAMZero2_Begin - dLAM2_0_t) * dEILam2_t;
    const double dLAM1_t_T    = (dLAM1_0_T - dLAM1_0_t) * dEILam1_t;
    const double dLAM2_t_T    = (dLAM2_0_T - dLAM2_0_t) * dEILam2_t;

    // XI
    const double dXI1Sqrd     = dLAM1_t_T * dLAM1_t_T * dPHI1_t;
    const double dXI2Sqrd     = dLAM2_t_T * dLAM2_t_T * dPHI2_t;
    const double dXICrossSqrd = dLAM1_t_T * dLAM2_t_T * dPHICross_t;

    /// Bond Variance
    const double dBondVar1 = dXI1Sqrd + dRho * dXICrossSqrd;
    const double dBondVar2 = dXI2Sqrd + dRho * dXICrossSqrd;

    // Exponent
    const double dX2 = dX3 * dAlpha_Sqrt1mRhoSqrd + dAlphaRho * dX1;  // back out X2 from X1 and X3
    double       dExponent = -0.5 * (dBondVar1 + dBondVar2);
    dExponent += dBondVar1 * dLAM1_t_Nume / dLAM1_t_T;
    dExponent += dBondVar2 * dLAM2_t_Nume / dLAM2_t_T;
    dExponent -= dX1 * dLAM1_t_T;
    dExponent -= dX2 * dLAM2_t_T;

    _ASSERTE(dLAM1_t_T > 0 && dLAM2_t_T > 0.);
    _ASSERTE(
        (dPHI1_t > 0. && dPHI2_t > 0. && dPHICross_t > 0.) ||     // regular 2f case
        (dPHI1_t == 0. && dPHI2_t == 0. && dPHICross_t == 0.) ||  // initial point at today
        (dPHI1_t > 0. && dPHI2_t == 0. &&
         dPHICross_t == 0.)  // 1f case with 2nd dimension suppressed
    );

    return exp(dExponent) * dDF_0_T / dDF_0_t;
}

static double __DF_(
    double dX1,
    double dX3,
    double dDF_0_t,
    double dDF_0_T,
    double dLAM1_0_t,
    double dLAM2_0_t,
    double dLAM1_0_T,
    double dLAM2_0_T,
    double dEILam1_t,
    double dEILam2_t,
    double dPHI1_t,
    double dPHI2_t,
    double dPHICross_t,
    double dExp_R1,  // contains measure-related information
    double dExp_R2,  // contains measure-related information
    double dRho,
    double dAlphaRho,
    double dAlpha_Sqrt1mRhoSqrd)
{
    /// LAMBDA
    const double dLAM1_t_T = (dLAM1_0_T - dLAM1_0_t) * dEILam1_t;
    const double dLAM2_t_T = (dLAM2_0_T - dLAM2_0_t) * dEILam2_t;

    // XI
    const double dXI1Sqrd     = dLAM1_t_T * dLAM1_t_T * dPHI1_t;
    const double dXI2Sqrd     = dLAM2_t_T * dLAM2_t_T * dPHI2_t;
    const double dXICrossSqrd = dLAM1_t_T * dLAM2_t_T * dPHICross_t;

    /// Bond Variance
    const double dHalfBondVar = 0.5 * (dXI1Sqrd + dXI2Sqrd) + dRho * dXICrossSqrd;

    // Exponent
    const double dR1 = dExp_R1 + dX1;
    const double dX2 = dX3 * dAlpha_Sqrt1mRhoSqrd + dAlphaRho * dX1;  // back out X2 from X1 and X3
    const double dR2 = dExp_R2 + dX2;
    double       dExponent = -dHalfBondVar;

    dExponent -= dR1 * dLAM1_t_T;
    dExponent -= dR2 * dLAM2_t_T;

    _ASSERTE(dLAM1_t_T > 0 && dLAM2_t_T > 0.);
    _ASSERTE(
        (dPHI1_t > 0. && dPHI2_t > 0. && dPHICross_t > 0.) ||
        (dPHI1_t == 0. && dPHI2_t == 0. && dPHICross_t == 0.));

    return exp(dExponent) * dDF_0_T / dDF_0_t;
}

double _DF_nojump_(
    double          dX1,
    double          dX3,
    double          dDF_0_t,
    double          dDF_0_T,
    double          dLAM1_0_t,
    double          dLAM2_0_t,
    double          dLAM1_0_T,
    double          dLAM2_0_T,
    double          dEILam1_t,
    double          dEILam2_t,
    double          dPHI1_t,
    double          dPHI2_t,
    double          dPHICross_t,
    const _PCQ_Seq* pNumeraire,
    double          dRho,
    double          dAlphaRho,
    double          dAlpha_Sqrt1mRhoSqrd)
{
    // const double dScale = 1.e-5;
    _ASSERTE(dDF_0_t >= dDF_0_T);
    if (dDF_0_t == dDF_0_T)
        return 1.;
    return __DF_nojump_(
        dX1,
        dX3,
        dDF_0_t,
        dDF_0_T,
        dLAM1_0_t,
        dLAM2_0_t,
        dLAM1_0_T,
        dLAM2_0_T,
        dEILam1_t,
        dEILam2_t,
        dPHI1_t,
        dPHI2_t,
        dPHICross_t,
        pNumeraire,
        dRho,
        dAlphaRho,
        dAlpha_Sqrt1mRhoSqrd);
}

double _DF_(
    double dX1,
    double dX3,
    double dDF_0_t,
    double dDF_0_T,
    double dLAM1_0_t,
    double dLAM2_0_t,
    double dLAM1_0_T,
    double dLAM2_0_T,
    double dEILam1_t,
    double dEILam2_t,
    double dPHI1_t,
    double dPHI2_t,
    double dPHICross_t,
    double dExp_R1,
    double dExp_R2,
    double dRho,
    double dAlphaRho,
    double dAlpha_Sqrt1mRhoSqrd)
{
    // const double dScale = 1.e-5;
    _ASSERTE(dDF_0_t >= dDF_0_T);
    if (dDF_0_t == dDF_0_T)
        return 1.;
    return __DF_(
        dX1,
        dX3,
        dDF_0_t,
        dDF_0_T,
        dLAM1_0_t,
        dLAM2_0_t,
        dLAM1_0_T,
        dLAM2_0_T,
        dEILam1_t,
        dEILam2_t,
        dPHI1_t,
        dPHI2_t,
        dPHICross_t,
        dExp_R1,
        dExp_R2,
        dRho,
        dAlphaRho,
        dAlpha_Sqrt1mRhoSqrd);
}

//-------------------------------------------------------------------------------------------------------
//	Description: Local helper for _fill_ExpR
//
//  Returns:  the expected value of R(t) under the measure induced by ZeroBond(t,T)
//-------------------------------------------------------------------------------------------------------
static double _ExpR(
    double dLAM_0_t,
    double dLAM_0_T,
    double dAuxLAM_0_t,
    double dAuxLAM_0_T,
    double dEILam1_t,
    double dEIAuxLam1_t,
    double dPHI_t,
    double dPHICross_t,
    double dRho)
{
    const double dLAM1_t_T    = (dLAM_0_T - dLAM_0_t) * dEILam1_t;
    const double dAuxLAM1_t_T = (dAuxLAM_0_T - dAuxLAM_0_t) * dEIAuxLam1_t;
    return -dLAM1_t_T * dPHI_t - dRho * dAuxLAM1_t_T * dPHICross_t;
}

//-------------------------------------------------------------------------------------------------------
//	Description:
//
//  Returns:  the expected value of R1,2(t) under the measure induced by ZeroBond(t,T)
//-------------------------------------------------------------------------------------------------------
void _fill_ExpR(
    const double*   pdLAM1_0_T_Begin,
    const double*   pdLAM1_0_T_End,
    const double*   pdLAM2_0_T_Begin,
    const double*   pdEILam1_T_Begin,
    const double*   pdEILam2_T_Begin,
    const double*   pdPHI1_T_Begin,
    const double*   pdPHI2_T_Begin,
    const double*   pdPHICross_T_Begin,
    double          dRho,
    const _PCQ_Seq* pNumeraire,
    // results
    double* pdExp_R1_Begin,
    double* pdExp_R2_Begin)
{
    // R1(0) = 0,R2(0) = 0.
    double dExp_R1 = 0.;
    double dExp_R2 = 0.;

    // variables at time today (t=0.)
    double dLAM1_0_Tm1   = 0.;
    double dLAM2_0_Tm1   = 0.;
    double dEILam1_Tm1   = 1.;
    double dEILam2_Tm1   = 1.;
    double dPHI1_Tm1     = 0.;
    double dPHI2_Tm1     = 0.;
    double dPHICross_Tm1 = 0.;

    for (; pdLAM1_0_T_Begin < pdLAM1_0_T_End; ++pdLAM1_0_T_Begin,
                                              ++pdLAM2_0_T_Begin,
                                              ++pdEILam1_T_Begin,
                                              ++pdEILam2_T_Begin,
                                              ++pdPHI1_T_Begin,
                                              ++pdPHI2_T_Begin,
                                              ++pdPHICross_T_Begin,
                                              ++pdExp_R1_Begin,
                                              ++pdExp_R2_Begin)
    {
        const double dLAM1_0_NumMAt =
            pNumeraire ? *pNumeraire->m_pdLAMZero1_Begin : *pdLAM1_0_T_Begin;
        const double dLAM2_0_NumMAt =
            pNumeraire ? *pNumeraire->m_pdLAMZero2_Begin : *pdLAM2_0_T_Begin;

        /// expected R1(T-1) under the measure induced by ZeroBond(t,T)
        /// or used-passed in numeraire
        const double dExp_R1Tm1 = _ExpR(
            dLAM1_0_Tm1,
            dLAM1_0_NumMAt,
            dLAM2_0_Tm1,
            dLAM2_0_NumMAt,
            dEILam1_Tm1,
            dEILam2_Tm1,
            dPHI1_Tm1,
            dPHICross_Tm1,
            dRho);

        /// expected R2(T-1) under the measure induced by ZeroBond(t,T)
        /// or used-passed in numeraire
        const double dExp_R2Tm1 = _ExpR(
            dLAM2_0_Tm1,
            dLAM2_0_NumMAt,
            dLAM1_0_Tm1,
            dLAM1_0_NumMAt,
            dEILam2_Tm1,
            dEILam1_Tm1,
            dPHI2_Tm1,
            dPHICross_Tm1,
            dRho);

        /// expected R1(T) under the measure induced by ZeroBond(t,T)
        /// or used-passed in numeraire
        const double dExp_R1T = _ExpR(
            *pdLAM1_0_T_Begin,
            dLAM1_0_NumMAt,
            *pdLAM2_0_T_Begin,
            dLAM2_0_NumMAt,
            *pdEILam1_T_Begin,
            *pdEILam2_T_Begin,
            *pdPHI1_T_Begin,
            *pdPHICross_T_Begin,
            dRho);

        /// expected R2(T) under the measure induced by ZeroBond(t,T)
        /// or used-passed in numeraire
        const double dExp_R2T = _ExpR(
            *pdLAM2_0_T_Begin,
            dLAM2_0_NumMAt,
            *pdLAM1_0_T_Begin,
            dLAM1_0_NumMAt,
            *pdEILam2_T_Begin,
            *pdEILam1_T_Begin,
            *pdPHI2_T_Begin,
            *pdPHICross_T_Begin,
            dRho);

        // Use the fact that E(J) [ R(T)| F(T-1) ]
        // = R(T-1) + E(ZC(T)) [ R(T)| F(0) ] - E(ZC(T)) [ R(T-1)| F(0) ]
        *pdExp_R1_Begin = dExp_R1 + dExp_R1T - dExp_R1Tm1;
        *pdExp_R2_Begin = dExp_R2 + dExp_R2T - dExp_R2Tm1;

        dLAM1_0_Tm1   = *pdLAM1_0_T_Begin;
        dLAM2_0_Tm1   = *pdLAM2_0_T_Begin;
        dEILam1_Tm1   = *pdEILam1_T_Begin;
        dEILam2_Tm1   = *pdEILam2_T_Begin;
        dPHI1_Tm1     = *pdPHI1_T_Begin;
        dPHI2_Tm1     = *pdPHI2_T_Begin;
        dPHICross_Tm1 = *pdPHICross_T_Begin;
        dExp_R1       = *pdExp_R1_Begin;
        dExp_R2       = *pdExp_R2_Begin;
    }
}

//////////////////////////////////////////////////////////////////////
///  Compute discount factors,  given E(R1(t)) and E(R2(t)) at slice t
void _df_singleslice(
    // grid info
    const double* pdX1_Begin,
    const double* pdX1_End,
    const double* pdX3_Begin,
    const double* pdX3_End,
    double        dDF_0_t,
    double        dLAM1_0_t,
    double        dLAM2_0_t,
    double        dExpIntglam1_t,
    double        dExpIntglam2_t,
    double        dPHI1_t,
    double        dPHI2_t,
    double        dPHICross_t,
    double        dDF_0_T,
    double        dLAM1_0_T,
    double        dLAM2_0_T,
    double        dExp_R1,
    double        dExp_R2,
    double        dRho,
    double        dAlphaRho,
    double        dAlpha_Sqrt1mRhoSqrd,
    /// result
    double** ppdDF,
    int      nR_Begin,
    int      nC_Begin)
{
    const int nR_End = nR_Begin + pdX3_End - pdX3_Begin;
    const int nC_End = nC_Begin + pdX1_End - pdX1_Begin;

    const double dInv_DF_0_t = 1. / dDF_0_t;
    const double dBondRatio  = dDF_0_T * dInv_DF_0_t;

    /// Constants needed to utilize the regularity of the grid to compute DF fast
    const double dDeltaX1 = (pdX1_End - pdX1_Begin == 1) ? 0. : pdX1_Begin[1] - pdX1_Begin[0];
    const double dDeltaX3 = (pdX3_End - pdX3_Begin == 1) ? 0. : pdX3_Begin[1] - pdX3_Begin[0];
    const double dDeltaX1_AlphaRho            = dDeltaX1 * dAlphaRho;
    const double dDeltaX3_Alpha_Sqrt1mRhoSqrd = dDeltaX3 * dAlpha_Sqrt1mRhoSqrd;

    const double dX1_Base = *pdX1_Begin - dDeltaX1;
    const double dX3_Base = *pdX3_Begin - dDeltaX3;
    const double dX2_Base = dX3_Base * dAlpha_Sqrt1mRhoSqrd + *pdX1_Begin * dAlphaRho;
    double       dConstant, dCoeff1, dCoeff2;

    int nX1, nX3;

    _ASSERTE(pdX1_End - pdX1_Begin >= 1);
    _ASSERTE(pdX3_End - pdX3_Begin >= 1);

    _df_constant(
        dDF_0_t,
        dDF_0_T,
        dLAM1_0_t,
        dLAM2_0_t,
        dLAM1_0_T,
        dLAM2_0_T,
        dExpIntglam1_t,
        dExpIntglam2_t,
        dPHI1_t,
        dPHI2_t,
        dPHICross_t,
        dExp_R1,
        dExp_R2,
        dRho,
        &dConstant,
        &dCoeff1,
        &dCoeff2);

    {
        const double dExpX1 = exp(-dDeltaX1 * dCoeff1 - dDeltaX1_AlphaRho * dCoeff2);
        const double dExpX3 = exp(-dDeltaX3_Alpha_Sqrt1mRhoSqrd * dCoeff2);
        double       dDF_Base =
            exp(dConstant - (*pdX1_Begin) * dCoeff1 - dX2_Base * dCoeff2) * dBondRatio;

        for (nX1 = nC_Begin; nX1 < nC_End; ++nX1, dDF_Base *= dExpX1)  /// X1 - cols !
        {
            double dDF = dDF_Base * dExpX3;
            for (nX3 = nR_Begin; nX3 < nR_End; ++nX3, dDF *= dExpX3)  /// X3 - rows !
                ppdDF[nX3][nX1] = dDF;
        }
    }
}
