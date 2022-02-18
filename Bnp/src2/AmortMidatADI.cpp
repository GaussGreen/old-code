//#include "AmortMidatADI.h"
//
//#include "AmortMidatADIGrid.h"
//#include "math.h"
//
////-------------------------------------------------------------------------------------------------------
////	Description: common function parameters needed by the binary function PFN_BI
////-------------------------------------------------------------------------------------------------------
// typedef struct _COMM_PFN_BI_
//{
//    const PFN_COEFF_3DPDE
//                 m_FUNC_;  // tertiary function PFN_COEFF_3DPDE(double x,double aux,double
//                 t,void*)
//    const void*  m_pT;     // time slice
//    const double m_dMult;  // for discretization
//    void*        m_pComm;  // common parameters accepted by m_FUNC_
//} _Comm_PFN_BI;
//
////-------------------------------------------------------------------------------------------------------
////	Description: binary functor to update 0.5*drift/delta_x
////
////	Returns :
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _Half_Mu_InvDeltaX_PFN_BI(
//    const void* pX,
//    const void* pAuxX,  // pIn2
//    /// Common parameters used in drift functor
//    const void* pComm,
//    void*       pMu  // pComm, // common function parameters
//)
//{
//    /// cast
//    const _Comm_PFN_BI* p = (const _Comm_PFN_BI*)pComm;
//    // result
//    *(double*)pMu = (p->m_dMult) * p->m_FUNC_(pX, pAuxX, p->m_pT, p->m_pComm);
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description: updates the drift matrix ppdHalf_Mu_InvDeltaX
////						drift is a function of X, AuxX and Time
////
////  Return : ppdHalf_Mu_InvDeltaX = 0.5 * Mu(X,Aux,Time) /DeltaX
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _Update_Mu(
//    /// State X
//    const double* pdX_Begin,
//    const double* pdX_End,
//    // Auxiliary State X
//    const double* pdAuxX_Begin,
//    const double* pdAuxX_End,
//    // Time
//    const double* pdT,
//    /// Common parameters used in drift functor
//    void* pComm,
//    // functor to update the drift : _MU_FUNC_(State X, Aux X, T, Comm)
//    // NB: first argument must be State X, not Auxiliary X
//    PFN_COEFF_3DPDE _MU_FUNC_,
//    /// Results
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppdHalf_Mu_InvDeltaX,  // based
//    int      nX_Begin,
//    int      nAuxX_Begin)
//{
//    // if there is no diffusion in the X direction, return
//    if (pdX_End - pdX_Begin == 1)
//        return;
//
//    else
//    {
//        const double dDelta_X        = pdX_Begin[1] - pdX_Begin[0];
//        const double dHalf_InvDeltaX = 0.5 / dDelta_X;
//
//        _Comm_PFN_BI Comm_BI = {_MU_FUNC_, pdT, dHalf_InvDeltaX, pComm};
//
//        _transform_matrix(
//            pdX_Begin,
//            pdX_End - pdX_Begin,
//            sizeof(*pdX_Begin),
//            pdAuxX_Begin,
//            pdAuxX_End - pdAuxX_Begin,
//            sizeof(*pdAuxX_Begin),
//            sizeof(**ppdHalf_Mu_InvDeltaX),
//            _Half_Mu_InvDeltaX_PFN_BI,
//            &Comm_BI,
//            ppdHalf_Mu_InvDeltaX,
//            nAuxX_Begin,
//            nX_Begin,
//            0);
//    }
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description: updates the local vol and drift of the PDE at time *pdT.
////
////  Return : ppdHalf_MuX_InvDeltaX = 0.5 * Mu(X,Aux,Time)/DeltaX
////				 pdHalf_SigXSqrd_InvDeltaXSqrd = 0.5 *  Sig(Time)/DeltaX^2
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _Update_Slice(
//    // Time
//    const double* pdT,
//    // State X
//    const double* pdX_Begin,
//    const double* pdX_End,
//    // Auxiliary State X
//    const double* pdAuxX_Begin,
//    const double* pdAuxX_End,
//    /// Common parameters used in vol and drift functor
//    void* pComm,
//    /// functors to update local vol : _SIGMA_Sqrd_FUNC_(State X, Aux State X, Time)
//    // NB: first argument must be State X, not Auxiliary X
//    PFN_COEFF_3DPDE _SIGMA_Sqrd_FUNC_,
//    /// functors to update drift : _MU_FUNC_(State X, Aux State X, Time)
//    // NB: first argument must be State X, not Auxiliary X
//    PFN_COEFF_3DPDE _MU_FUNC_,
//    /// Results
//    double* pdHalf_SigSqrd_InvDeltaXSqrd,
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [ nX_Begin, ...)
//    double** ppdHalf_Mu_InvDeltaX,
//    int      nX_Begin,    // index on the matrix
//    int      nAuxX_Begin  // index on the matrix
//)
//{
//    // if there is no diffusion in the X direction, return
//    if (pdX_End - pdX_Begin == 1)
//        return;
//
//    else
//    {
//        const double dDelta_X     = pdX_Begin[1] - pdX_Begin[0];
//        const double dDeltaSqrd_X = dDelta_X * dDelta_X;
//
//        _ASSERTE(pdX_End - pdX_Begin >= 1);
//
//        // if there is no diffusion in the X direction, then no updating is needed!
//        if (dDelta_X == 0.)
//            return;
//
//        // update volatility-related info
//        *pdHalf_SigSqrd_InvDeltaXSqrd = 0.5 * _SIGMA_Sqrd_FUNC_(0, 0, pdT, pComm) / dDeltaSqrd_X;
//
//        // update drift-related info
//        _Update_Mu(
//            pdX_Begin,
//            pdX_End,
//            pdAuxX_Begin,
//            pdAuxX_End,
//            pdT,
//            pComm,
//            _MU_FUNC_,
//            ppdHalf_Mu_InvDeltaX,
//            nX_Begin,
//            nAuxX_Begin);
//    }
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description: decomposes a tridiagonal system
////
////  To do : Misleading argument names used! (will come back ...)
////
////	Returns : Lower, Inverse of the Diagonal, Upper
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _lu_decompose_tridiag_(
//    double        dDiag_Imp_front,     // for boundary
//    double        dOffDiag_Imp_front,  // for boundary
//    double        dDiag_Imp_back,      // for boundary
//    double        dOffDiag_Imp_back,   // for boundary
//    double        dInvHalfDeltaT,
//    double        dHalf_InvDeltaXSqrd_SigSqrd,
//    const double* pdHalf_InvDeltaX_Mu_Begin,
//    const double* pdHalf_InvDeltaX_Mu_End,
//    // results
//    double* pdDInv_Begin,
//    double* pdL_Begin,
//    double* pdU_Begin)
//{
//    const double dInvDeltaXSqrd_SigSqrd = 2. * dHalf_InvDeltaXSqrd_SigSqrd;
//    const double dC                     = dInvHalfDeltaT + dInvDeltaXSqrd_SigSqrd;
//    _ASSERTE(dC != 0.);
//    *pdDInv_Begin = 1. / dC;
//
//    if (pdHalf_InvDeltaX_Mu_End - pdHalf_InvDeltaX_Mu_Begin >= 3)
//    {
//        const double* pdMuQ = pdHalf_InvDeltaX_Mu_Begin + 1;
//        double *      pdU = pdU_Begin + 1, *pdDInv = pdDInv_Begin + 1, *pdL = pdL_Begin + 1;
//
//        // update super diagonal
//        *pdU_Begin = dOffDiag_Imp_front;  // NB : boundary
//        for (; pdMuQ < pdHalf_InvDeltaX_Mu_End - 1; ++pdMuQ, ++pdU)
//            *pdU = -(*pdMuQ) - dHalf_InvDeltaXSqrd_SigSqrd;
//
//        // boundary
//        {
//            // update L and D (main diagonal)
//            _ASSERTE(dDiag_Imp_front != 0.);
//            *pdDInv_Begin = 1. / dDiag_Imp_front;
//            *pdL_Begin =
//                (*(pdHalf_InvDeltaX_Mu_Begin + 1) - dHalf_InvDeltaXSqrd_SigSqrd) /
//                dDiag_Imp_front;
//            _ASSERTE((dC - *pdL_Begin * dOffDiag_Imp_front) != 0.);
//            *(pdDInv_Begin + 1) = 1. / (dC - *pdL_Begin * dOffDiag_Imp_front);
//        }
//
//        for (pdMuQ = pdHalf_InvDeltaX_Mu_Begin + 1; pdMuQ < pdHalf_InvDeltaX_Mu_End - 2;
//             ++pdMuQ, ++pdL, ++pdDInv)
//        {
//            *pdL          = (*(pdMuQ + 1) - dHalf_InvDeltaXSqrd_SigSqrd) * (*pdDInv);
//            *(pdDInv + 1) = dC + (*pdL) * (dHalf_InvDeltaXSqrd_SigSqrd + (*pdMuQ));
//            _ASSERTE((*(pdDInv + 1)) != 0.);
//            *(pdDInv + 1) = 1. / (*(pdDInv + 1));  // invert
//        }
//
//        // boundary
//        {
//            const int nSize = (pdHalf_InvDeltaX_Mu_End - pdHalf_InvDeltaX_Mu_Begin);
//
//            pdL    = pdL_Begin + nSize - 2;
//            pdDInv = pdDInv_Begin + nSize - 2;
//            pdMuQ  = pdHalf_InvDeltaX_Mu_Begin + nSize - 2;
//
//            *pdL          = dOffDiag_Imp_back * (*pdDInv);
//            *(pdDInv + 1) = dDiag_Imp_back + (*pdL) * (dHalf_InvDeltaXSqrd_SigSqrd + *pdMuQ);
//            _ASSERTE((*(pdDInv + 1)) != 0.);
//            *(pdDInv + 1) = 1. / (*(pdDInv + 1));
//        }
//    }
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description: explicit scheme in the State X,traversing from ppd_From to ppd_To
////						at the boundary, the second derivative w.r.t to State X
///is /0
////
////	Returns : ppd_To
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _ADI_SingleSlice_Explicit(
//    // 1/ (Delta_T * 0.5)
//    double dInvHalfDeltaT,
//    // "From matrix" - rows represent State X and cols represent Auxiliary State X
//    const double** ppd_From,
//    int            nX_Begin,     // starting row index of pdd_from
//    int            nX_End,       // end row index of pdd_from, not inclusive
//    int            nAuxX_Begin,  // starting col index of pdd_from
//    int            nAuxX_End,    // end row index of pdd_from, not inclusive
//
//    // drift at time slice from	=  0.5 *  Mu (State X, Aux State X, T)/DeltaX
//    // NB: same size of the TRANSPOSE of ppd_From
//    // Row Indices = [nAuxX_Begin,...nAuxX_End) , Col Indices = [nX_Begin, ...nX_End)
//    const double** ppdHalf_Mu_InvDeltaX_from,
//
//    // local vol at time slice from, states - independent = 0.5 *  Sig(T)/Delta_X^2
//    double dHalf_SigSqrd_InvDeltaXSqrd_from,
//
//    /// Result
//    /// same size of the TRANSPOSE ppd_From
//    /// Rows = [nAuxX_Begin, ... nAuxX_End) ; Cols =  [nX_Begin, ... nX_End)
//    double** ppd_To)
//{
//    const double  dC            = dInvHalfDeltaT - 2. * dHalf_SigSqrd_InvDeltaXSqrd_from;
//    const double *pdCenter_from = ppd_From[nX_Begin], *pdUp_from = 0, *pdDown_from = 0;
//    int           nX, nAuxX;
//
//    _ASSERTE(nX_End - nX_Begin >= 1);
//
//    // if # of nodes in the dominant direction is less than 3, return
//    if (nX_End - nX_Begin == 1)
//    {
//        for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//            ppd_To[nAuxX][nX_Begin] = dC * pdCenter_from[nAuxX];
//    }
//
//    else
//    {
//        for (nX = nX_Begin + 1; nX < nX_End - 1; ++nX)
//        {
//            pdCenter_from = ppd_From[nX];
//            pdUp_from     = ppd_From[nX + 1];
//            pdDown_from   = ppd_From[nX - 1];
//
//            for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//            {
//                double*      p    = &(ppd_To[nAuxX][nX]);
//                const double dMuQ = ppdHalf_Mu_InvDeltaX_from[nAuxX][nX];
//                *p                = dC * pdCenter_from[nAuxX];
//                *p += (dMuQ + dHalf_SigSqrd_InvDeltaXSqrd_from) * pdUp_from[nAuxX];
//                *p += (-dMuQ + dHalf_SigSqrd_InvDeltaXSqrd_from) * pdDown_from[nAuxX];
//            }
//        }
//
//        /// handle the boundary - nRow_Begin
//        {
//            pdCenter_from = ppd_From[nX_Begin];
//            pdUp_from     = ppd_From[nX_Begin + 1];
//            pdDown_from   = ppd_From[nX_Begin - 1];
//
//            for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//            {
//                const double dMuQ = 2. * ppdHalf_Mu_InvDeltaX_from[nAuxX][nX_Begin];
//                ppd_To[nAuxX][nX_Begin] =
//                    (dInvHalfDeltaT - dMuQ) * pdCenter_from[nAuxX] + dMuQ * pdUp_from[nAuxX];
//            }
//        }
//
//        /// handle the boundary - nRow_End-1
//        {
//            pdCenter_from = ppd_From[nX_End - 1];
//            pdDown_from   = ppd_From[nX_End - 2];
//
//            for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//            {
//                const double dMuQ = 2. * ppdHalf_Mu_InvDeltaX_from[nAuxX][nX_End - 1];
//                ppd_To[nAuxX][nX_End - 1] =
//                    (dInvHalfDeltaT + dMuQ) * pdCenter_from[nAuxX] - dMuQ * pdDown_from[nAuxX];
//            }
//        }
//    }
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description: implict scheme in the State X,traversing from ppd_From to ppd_To
////						at the boundary, the second derivative w.r.t to State X
///is /0
////
////  To do : to be optimzed with pointer arithematic
////
////	Returns : ppd_From
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _ADI_SingleSlice_Implicit(
//    // 1/ (Delta_T * 0.5)
//    double dInvHalfDeltaT,
//    // "From matrix" - rows represent State X and cols represent Auxiliary State X
//    // On input : From Matrix
//    // On output : To Matrix
//    double** ppd_From,
//    int      nX_Begin,     // starting row index of pdd_from
//    int      nX_End,       // end row index of pdd_from, not inclusive
//    int      nAuxX_Begin,  // starting col index of pdd_from
//    int      nAuxX_End,    // end row index of pdd_from, not inclusive
//    // drift at time slice from	=  0.5 *  Mu (State X, Aux State X, T) /DeltaX
//    // NB: same size of the TRANSPOSE of ppd_From
//    // Row Indices = [nAuxX_Begin,...nAuxX_End) , Col Indices = [nX_Begin, ...nX_End)
//    const double** ppdHalf_MuX_InvDeltaX_to,
//    // local vol at time slice from, states - independent = 0.5 *  Sig^2(T) /Delta_X^2
//    double dHalf_SigXSqrd_InvDeltaXSqrd_to,
//    /// storage
//    // NB: all are of the same size of the TRANSPOSE of pdd_from
//    /// Rows = [nAuxX_Begin, ... nAuxX_End) ; Cols =  [nX_Begin, ... nX_End)
//    double** ppdStorage,
//    double** ppd_DInv,
//    double** ppd_LSub,
//    double** ppd_USup)
//{
//    const double *pdCenter_from = 0, *pdDown_to = 0;
//    double*       pdCenter_to = 0;
//    const double *pdLSub = 0, *pd_DInv = 0, *pUSup = 0;
//
//    int      nX, nAuxX;
//    double** ppdPlaceHolder = 0;
//
//    // if # of nodes in the dominant direction is less than 3, return
//    //_ASSERTE(nX_End-nX_Begin>=1);
//
//    //// decompose LU
//    for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//    {
//        const double dInvDeltaXR_MuR_front = 2. * ppdHalf_MuX_InvDeltaX_to[nAuxX][nX_Begin];
//        const double dInvDeltaXR_MuR_back  = 2. * ppdHalf_MuX_InvDeltaX_to[nAuxX][nX_End - 1];
//
//        const double dDiag_front    = dInvHalfDeltaT + dInvDeltaXR_MuR_front;
//        const double dOffDiag_front = -dInvDeltaXR_MuR_front;
//        const double dDiag_back     = dInvHalfDeltaT - dInvDeltaXR_MuR_back;
//        const double dOffDiag_back  = dInvDeltaXR_MuR_back;
//
//        _lu_decompose_tridiag_(
//            dDiag_front,
//            dOffDiag_front,
//            dDiag_back,
//            dOffDiag_back,
//            dInvHalfDeltaT,
//            dHalf_SigXSqrd_InvDeltaXSqrd_to,
//            ppdHalf_MuX_InvDeltaX_to[nAuxX] + nX_Begin,
//            ppdHalf_MuX_InvDeltaX_to[nAuxX] + nX_End,
//            ppd_DInv[nAuxX] + nX_Begin,
//            ppd_LSub[nAuxX] + nX_Begin,
//            ppd_USup[nAuxX] + nX_Begin);
//    }
//
//    // Lower Ly = b
//    {
//        /// boundary
//        pdCenter_from = ppd_From[nX_Begin];
//        for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//            ppdStorage[nAuxX][nX_Begin] = pdCenter_from[nAuxX];
//
//        for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//        {
//            pdLSub      = ppd_LSub[nAuxX];
//            pdCenter_to = ppdStorage[nAuxX];
//
//            for (nX = nX_Begin + 1; nX < nX_End; ++nX)
//                pdCenter_to[nX] = ppd_From[nX][nAuxX] - pdLSub[nX - 1] * pdCenter_to[nX - 1];
//        }
//    }
//
//    // SWAP
//    _SWAP_(ppd_From, ppdStorage, ppdPlaceHolder);
//
//    // Lower Ux = y
//    {
//        pdCenter_to = ppdStorage[nX_End - 1];
//        for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//            pdCenter_to[nAuxX] = ppd_From[nAuxX][nX_End - 1] * ppd_DInv[nAuxX][nX_End - 1];
//
//        for (nAuxX = nAuxX_Begin; nAuxX < nAuxX_End; ++nAuxX)
//        {
//            pd_DInv       = ppd_DInv[nAuxX];
//            pUSup         = ppd_USup[nAuxX];
//            pdCenter_from = ppd_From[nAuxX];
//
//            for (nX = nX_End - 2; nX >= nX_Begin; --nX)
//                ppdStorage[nX][nAuxX] =
//                    (pdCenter_from[nX] - pUSup[nX] * ppdStorage[nX + 1][nAuxX]) * pd_DInv[nX];
//        }
//    }
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description: traverses from ppd_From to ppd_To in one single slice
////
////	Returns : ppdU
////
////-------------------------------------------------------------------------------------------------------
// static void __stdcall _ADI_SingleSlice(
//    // states
//    const double* pdX_Begin,  // dominant factor
//    const double* pdX_End,
//    const double* pdAuxX_Begin,  // Auxiliary factor
//    const double* pdAuxX_End,
//    /// On input : "From Matrix"
//    /// On output : "To Matrix"
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    double** ppdU,
//    int      nX_Begin,
//    int      nAuxX_Begin,
//    // 1/ (Delta_T * 0.5)
//    double dInvHalfDeltaT,
//    /// info at slice from
//    double dHalf_SigSqrd_InvDeltaXSqrd_from,
//    double dHalf_SigSqrd_InvDeltaAuxXSqrd_from,
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    const double** ppdHalf_Mu_InvDeltaX_from,
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    const double** ppdHalf_Mu_InvDeltaAuxX_from,
//    // info at nth slice to
//    double dHalf_SigSqrd_InvDeltaXSqrd_to,
//    double dHalf_SigSqrd_InvDeltaAuxXSqrd_to,
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    const double** ppdHalf_Mu_InvDeltaX_to,
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    const double** ppdHalf_Mu_InvDeltaAuxX_to,
//
//    /// Storage
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppdStorage,
//
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppd_AuxXLX,
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppd_AuxXDInvX,
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppd_AuxXUX,
//
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    double** ppd_XLAuxX,
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    double** ppd_XDInvAuxX,
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    double** ppd_XUAuxX)
//{
//    /// size
//    const int nSizeX = pdX_End - pdX_Begin;
//    const int nX_End = nX_Begin + nSizeX;
//
//    const int nSizeAuxX = pdAuxX_End - pdAuxX_Begin;
//    const int nAuxX_End = nAuxX_Begin + nSizeAuxX;
//
//    /// explicit in dominant factor
//    _ADI_SingleSlice_Explicit(
//        dInvHalfDeltaT,
//        ppdU,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        ppdHalf_Mu_InvDeltaX_from,
//        dHalf_SigSqrd_InvDeltaXSqrd_from,
//        ppdStorage);
//
//    /// implicit in secondary factor
//    _ADI_SingleSlice_Implicit(
//        dInvHalfDeltaT,
//        ppdStorage,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        ppdHalf_Mu_InvDeltaAuxX_to,
//        dHalf_SigSqrd_InvDeltaAuxXSqrd_to,
//        ppdU,
//        ppd_XDInvAuxX,
//        ppd_XLAuxX,
//        ppd_XUAuxX);
//
//    // explict in secondary factor
//    _ADI_SingleSlice_Explicit(
//        dInvHalfDeltaT,
//        ppdStorage,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        ppdHalf_Mu_InvDeltaAuxX_from,
//        dHalf_SigSqrd_InvDeltaAuxXSqrd_from,
//        ppdU);
//
//    /// implicit in dominant factor
//    _ADI_SingleSlice_Implicit(
//        dInvHalfDeltaT,
//        ppdU,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        ppdHalf_Mu_InvDeltaX_to,
//        dHalf_SigSqrd_InvDeltaXSqrd_to,
//        ppdStorage,
//        ppd_AuxXDInvX,
//        ppd_AuxXLX,
//        ppd_AuxXUX);
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description:  initializes the storage needed by the call to _ADI_MultipleSlices
////
////	Returns : storage
////-------------------------------------------------------------------------------------------------------
// const char* Init_ADIStorage(_ADI_Storage* pstore, int nSize_X, int nSize_AuxX)
//{
//    pstore->m_pm = pstore->m_pdata = 0;
//    pstore->m_pm                   = calloc(6 * nSize_AuxX + 5 * nSize_X, sizeof(double*));
//    pstore->m_pdata                = calloc(11 * nSize_AuxX * nSize_X, sizeof(double));
//    if (!pstore->m_pm || !pstore->m_pdata)
//        return "Init_ADIStorage(..): mem allocation failed!";
//    return 0;
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description:  free storage
////
////	Returns :
////-------------------------------------------------------------------------------------------------------
// void Free_ADIStorage(const _ADI_Storage* pStorage)
//{
//    if (pStorage->m_pm)
//        free(pStorage->m_pm);
//    if (pStorage->m_pdata)
//        free(pStorage->m_pdata);
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description:   core PDE solving engine
////						 solver traverses from pdT_Begin to pdT_End
////						  this version without mem storage in the interface
////
////	Returns : Grid ppdU
////-------------------------------------------------------------------------------------------------------
// const char* ADI_MultipleSlices(
//    /// slice from pdT_Begin to pdT_End
//    /// IMPORTANT: by default backwards equation is solve
//    /// pdT_Begin is a DESCENDING sequence of time knots
//    const double* pdT_Begin,
//    const double* pdT_End,
//    // orthogonal states
//    const double* pdX_Begin,  // dominant  state
//    const double* pdX_End,
//    const double* pdAuxX_Begin,  // secondary state
//    const double* pdAuxX_End,
//    /// common parameters for updating vol and drift
//    void* pFuncArg,
//    // functor to update vol of the dominant variable
//    // functor must accept the argument seqence (X,AuxX,T)
//    PFN_COEFF_3DPDE _SIGMAX_Sqrd_FUNC_,
//    // functor to update drift of the dominant variable
//    // functor must accept the argument seqence (X,AuxX,T)
//    PFN_COEFF_3DPDE _MUX_FUNC_,
//    // functor to update vol of the secondary variable
//    // functor must accept the argument seqence (AuxX,X,T)
//    PFN_COEFF_3DPDE _SIGMAAuxX_Sqrd_FUNC_,
//    // functor to update drift of the secondary variable
//    // functor must accept the argument seqence (AuxX,X,T)
//    PFN_COEFF_3DPDE _MUAuxX_FUNC_,
//    // On input :  Grid at time pdT_Begin
//    // On output: Grid at time PdT_End-1
//    // Row Indices = [nX_Begin,...] , Col Indices = [nAuxX_Begin, ...]
//    double** ppdU,
//    int      nX_Begin,
//    int      nAuxX_Begin)
//{
//    _ADI_Storage store = {0};
//    const char*  szErr = Init_ADIStorage(&store, pdX_End - pdX_Begin, pdAuxX_End - pdAuxX_Begin);
//
//    if (szErr)
//    {
//        Free_ADIStorage(&store);
//        return szErr;
//    }
//
//    _ADI_MultipleSlices(
//        pdT_Begin,
//        pdT_End,
//        pdX_Begin,
//        pdX_End,
//        pdAuxX_Begin,
//        pdAuxX_End,
//        pFuncArg,
//        _SIGMAX_Sqrd_FUNC_,
//        _MUX_FUNC_,
//        _SIGMAAuxX_Sqrd_FUNC_,
//        _MUAuxX_FUNC_,
//        ppdU,
//        nX_Begin,
//        nAuxX_Begin,
//        &store);
//
//    Free_ADIStorage(&store);
//    return 0;
//}
//
////-------------------------------------------------------------------------------------------------------
////	Description:   core PDE solving engine
////						 solver traverses from pdT_Begin to pdT_End
////						 this version with mem storage in the interface
////
////	Returns : Grid ppdU
////-------------------------------------------------------------------------------------------------------
// void _ADI_MultipleSlices(
//    /// slice from pdT_Begin to pdT_End
//    /// IMPORTANT: by default backwards equation is solve
//    /// pdT_Begin is a DESCENDING sequence of time knots
//    const double* pdT_Begin,
//    const double* pdT_End,
//    // orthogonal states
//    const double* pdX_Begin,  // dominant  state
//    const double* pdX_End,
//    const double* pdAuxX_Begin,  // secondary state
//    const double* pdAuxX_End,
//    /// common parameters for updating vol and drift
//    void* pFuncArg,
//    // functor to update vol of the dominant variable
//    // functor must accept the argument seqence (X,AuxX,T)
//    PFN_COEFF_3DPDE _SIGMAX_Sqrd_FUNC_,
//    // functor to update drift of the dominant variable
//    // functor must accept the argument seqence (X,AuxX,T)
//    PFN_COEFF_3DPDE _MUX_FUNC_,
//    // functor to update vol of the secondary variable
//    // functor must accept the argument seqence (AuxX,X,T)
//    PFN_COEFF_3DPDE _SIGMAAuxX_Sqrd_FUNC_,
//    // functor to update drift of the secondary variable
//    // functor must accept the argument seqence (AuxX,X,T)
//    PFN_COEFF_3DPDE _MUAuxX_FUNC_,
//    // On input :  Grid at time pdT_Begin
//    // On output: Grid at time PdT_End-1
//    // Row Indices = [nX_Begin,...] , Col Indices = [nAuxX_Begin, ...]
//    double**            ppdU,
//    int                 nX_Begin,
//    int                 nAuxX_Begin,
//    const _ADI_Storage* pStorage)
//{
//    // sizes
//    const int nSize_AuxX = pdAuxX_End - pdAuxX_Begin;
//    const int nSize_X    = pdX_End - pdX_Begin;
//    const int nSize_Grid = nSize_AuxX * nSize_X;
//    const int nX_End     = nX_Begin + nSize_X;
//    const int nAuxX_End  = nAuxX_Begin + nSize_AuxX;
//
//    double** pm    = (double**)pStorage->m_pm;
//    double*  pdata = (double*)pStorage->m_pdata;
//
//    // assemble matrices
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppd_Storage = (double**)_alloca_matrix(
//        pm, pdata, nAuxX_Begin, nAuxX_End, nX_Begin, nX_End, sizeof(double));
//
//    /// Row Indices = [nAuxX_Begin,...) , Col Indices = [nX_Begin, ...)
//    double** ppdHalf_Mu1_InvDeltaX1_from = (double**)_alloca_matrix(
//        pm += nSize_AuxX,
//        pdata += nSize_Grid,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        sizeof(double));
//    double** ppdHalf_Mu1_InvDeltaX1_to = (double**)_alloca_matrix(
//        pm += nSize_AuxX,
//        pdata += nSize_Grid,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        sizeof(double));
//    double** ppd_AuxXLX = (double**)_alloca_matrix(
//        pm += nSize_AuxX,
//        pdata += nSize_Grid,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        sizeof(double));
//    double** ppd_AuxXDInvX = (double**)_alloca_matrix(
//        pm += nSize_AuxX,
//        pdata += nSize_Grid,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        sizeof(double));
//    double** ppd_AuxXUX = (double**)_alloca_matrix(
//        pm += nSize_AuxX,
//        pdata += nSize_Grid,
//        nAuxX_Begin,
//        nAuxX_End,
//        nX_Begin,
//        nX_End,
//        sizeof(double));
//
//    /// Row Indices = [nX_Begin,...) , Col Indices = [nAuxX_Begin, ...)
//    double** ppdHalf_Mu2_InvDeltaX2_from = (double**)_alloca_matrix(
//        pm += nSize_AuxX,
//        pdata += nSize_Grid,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        sizeof(double));
//    double** ppdHalf_Mu2_InvDeltaX2_to = (double**)_alloca_matrix(
//        pm += nSize_X,
//        pdata += nSize_Grid,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        sizeof(double));
//    double** ppd_XLAuxX = (double**)_alloca_matrix(
//        pm += nSize_X,
//        pdata += nSize_Grid,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        sizeof(double));
//    double** ppd_XDInvAuxX = (double**)_alloca_matrix(
//        pm += nSize_X,
//        pdata += nSize_Grid,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        sizeof(double));
//    double** ppd_XUAuxX = (double**)_alloca_matrix(
//        pm += nSize_X,
//        pdata += nSize_Grid,
//        nX_Begin,
//        nX_End,
//        nAuxX_Begin,
//        nAuxX_End,
//        sizeof(double));
//
//    double dHalf_Sig1Sqrd_InvDeltaX1Sqrd_from = 0;
//    double dHalf_Sig2Sqrd_InvDeltaX2Sqrd_from = 0;
//    double dHalf_Sig1Sqrd_InvDeltaX1Sqrd_to   = 0;
//    double dHalf_Sig2Sqrd_InvDeltaX2Sqrd_to   = 0;
//
//    // alias
//    double*  pSigXQ_from    = &dHalf_Sig1Sqrd_InvDeltaX1Sqrd_from;
//    double*  pSigAuxXQ_from = &dHalf_Sig2Sqrd_InvDeltaX2Sqrd_from;
//    double** ppMuXQ_from    = ppdHalf_Mu1_InvDeltaX1_from;
//    double** ppMuAuxXQ_from = ppdHalf_Mu2_InvDeltaX2_from;
//    double*  pSigXQ_to      = &dHalf_Sig1Sqrd_InvDeltaX1Sqrd_to;
//    double*  pSigAuxXQ_to   = &dHalf_Sig2Sqrd_InvDeltaX2Sqrd_to;
//    double** ppMuXQ_to      = ppdHalf_Mu1_InvDeltaX1_to;
//    double** ppMuAuxXQ_to   = ppdHalf_Mu2_InvDeltaX2_to;
//
//    const double* pdT = pdT_Begin;
//
//    // Update state info for X1 and X2 at slice from
//    _Update_Slice(
//        pdT,
//        pdX_Begin,
//        pdX_End,
//        pdAuxX_Begin,
//        pdAuxX_End,
//        pFuncArg,
//        _SIGMAX_Sqrd_FUNC_,
//        _MUX_FUNC_,
//        pSigXQ_from,
//        ppMuXQ_from,
//        nX_Begin,
//        nAuxX_Begin);
//
//    _Update_Slice(
//        pdT,
//        pdAuxX_Begin,
//        pdAuxX_End,
//        pdX_Begin,
//        pdX_End,
//        pFuncArg,
//        _SIGMAAuxX_Sqrd_FUNC_,
//        _MUAuxX_FUNC_,
//        pSigAuxXQ_from,
//        ppMuAuxXQ_from,
//        nAuxX_Begin,
//        nX_Begin);
//
//    /// slicing between dT_LB and dT_HB
//    for (pdT = pdT_Begin + 1; pdT < pdT_End; ++pdT)
//    {
//        const double dDeltaT = (*pdT - *(pdT - 1));
//
//        // diffuse only when dDeltaT >0.
//        if (dDeltaT != 0.)
//        {
//            const double dInvHalfDeltaT = fabs(-2. / dDeltaT);
//            _ASSERTE(dInvHalfDeltaT > 0.);
//
//            // Update state info for X1 and X2 at slice from
//            _Update_Slice(
//                pdT,
//                pdX_Begin,
//                pdX_End,
//                pdAuxX_Begin,
//                pdAuxX_End,
//                pFuncArg,
//                _SIGMAX_Sqrd_FUNC_,
//                _MUX_FUNC_,
//                pSigXQ_to,
//                ppMuXQ_to,
//                nX_Begin,
//                nAuxX_Begin);
//
//            _Update_Slice(
//                pdT,
//                pdAuxX_Begin,
//                pdAuxX_End,
//                pdX_Begin,
//                pdX_End,
//                pFuncArg,
//                _SIGMAAuxX_Sqrd_FUNC_,
//                _MUAuxX_FUNC_,
//                pSigAuxXQ_to,
//                ppMuAuxXQ_to,
//                nAuxX_Begin,
//                nX_Begin);
//
//            _ADI_SingleSlice(
//                pdX_Begin,
//                pdX_End,
//                pdAuxX_Begin,
//                pdAuxX_End,
//                ppdU,
//                nX_Begin,
//                nAuxX_Begin,
//                dInvHalfDeltaT,
//                *pSigXQ_from,
//                *pSigAuxXQ_from,
//                ppMuXQ_from,
//                ppMuAuxXQ_from,
//                *pSigXQ_to,
//                *pSigAuxXQ_to,
//                ppMuXQ_to,
//                ppMuAuxXQ_to,
//                ppd_Storage,
//                ppd_AuxXLX,
//                ppd_AuxXDInvX,
//                ppd_AuxXUX,
//                ppd_XLAuxX,
//                ppd_XDInvAuxX,
//                ppd_XUAuxX);
//
//            // SWAP information of the 2 slices
//            memswap(&pSigXQ_from, &pSigXQ_to, sizeof(double*));
//            memswap(&pSigAuxXQ_from, &pSigAuxXQ_to, sizeof(double*));
//            memswap(&ppMuXQ_from, &ppMuXQ_to, sizeof(double**));
//            memswap(&ppMuAuxXQ_from, &ppMuAuxXQ_to, sizeof(double**));
//        }
//    }
//}
//
//#ifdef _DEBUG
///////////////////////// Debug related functions  - Begin
// static double _alpha_;
// static double _gamma_;
// static double _rho_;
// static double _sigma_;
// static double _lambda_;
//
//// Sigma is a function of time only
// static double __stdcall _SIGMA1_Sqrd_TEST_(
//    const double* pdX1, const double* pdX2, const double* pdT, void* pAux)
//{
//    pdT, pAux;
//    return _sigma_ * _sigma_;
//    // unreferenced
//    pdX1, pdX2;
//}
//
// static double _LAMBDA1_TEST_()
//{
//    return _lambda_;
//}
//
// static double _ALPHA_TEST_()
//{
//    return _alpha_;
//}
//
// static double _GAMMA_TEST_()
//{
//    return _gamma_;
//}
//
// static double _RHO_TEST_()
//{
//    return _rho_;
//}
//
///// drift can be a function of states and time
// static double __stdcall _MU1_TEST_(
//    const double* pdX1, const double* pdX3, const double* pdT, void* pAux)
//{
//    const double dLambda1 = _LAMBDA1_TEST_();
//    pdX3, pdT, pAux;
//    // return -dLam1 * dX1;
//    // return 5.;
//    return -dLambda1 * (*pdX1);
//}
//
//// Sigma is a function of time only
// static double _SIGMA2_Sqrd_TEST_(
//    const double* pdX1, const double* pdX2, const double* pdT, void* pAux)
//{
//    const double dAlpha = _ALPHA_TEST_();
//    return _SIGMA1_Sqrd_TEST_(pdX1, pdX2, pdT, pAux) * dAlpha * dAlpha;
//}
//
//// Sigma is a function of time only
// static double __stdcall _SIGMA3_Sqrd_TEST_(
//    const double* pdX1, const double* pdX2, const double* pdT, void* pAux)
//{
//    return _SIGMA1_Sqrd_TEST_(pdX1, pdX2, pdT, pAux);
//}
//
///// drift can be a function of states and time
// static double __stdcall _MU3_TEST_(
//    const double* pdX3, const double* pdX1, const double* pdT, void* pAux)
//{
//    const double dGamma   = _GAMMA_TEST_();
//    const double dLambda1 = _LAMBDA1_TEST_();
//    const double dLambda2 = dLambda1 + dGamma;
//    const double dRho     = _RHO_TEST_();
//
//    pdT, pAux;
//
//    return -dLambda2 * (*pdX3) - dRho * dGamma * (*pdX1) / sqrt(1 - dRho * dRho);
//    // return dLambda3*dX3;
//}
//
// static const char* __adi_tester_(int nR_Begin, int nR_End, int nC_Begin, int nC_End, int nFlat)
//{
//    const char*  szErr = 0;
//    int          nr, nc;
//    int          nNumX1 = nC_End - nC_Begin;
//    int          nNumX2 = nR_End - nR_Begin;
//    const int    nNumT  = 100;
//    const double dTL    = 0;
//    const double dTH    = 1.;
//    const double dT     = dTH - dTL;
//
//    double** ppdU  = dmatrix(nR_Begin, nR_End - 1, nC_Begin, nC_End - 1);
//    double** ppd_U = dmatrix(nC_Begin, nC_End - 1, nR_Begin, nR_End - 1);
//    // alias
//    _Grid        grid     = {0};
//    const double dAlpha   = _ALPHA_TEST_();
//    const double dGamma   = _GAMMA_TEST_();
//    const double dRho     = _RHO_TEST_();
//    const double dLambda1 = _LAMBDA1_TEST_();
//    const double dLambda2 = dLambda1 + dGamma;
//
//    const double dTerminalVarX1 =
//        (dLambda1 == 0.)
//            ? dT * _SIGMA1_Sqrd_TEST_(0, 0, 0, 0)
//            : _SIGMA1_Sqrd_TEST_(0, 0, 0, 0) * (1 - exp(2 * dLambda1 * dT)) / (-2 * dLambda1);
//
//    const double dTerminalVarX2 =
//        (dLambda2 == 0.)
//            ? dT * _SIGMA2_Sqrd_TEST_(0, 0, 0, 0)
//            : _SIGMA2_Sqrd_TEST_(0, 0, 0, 0) * (1 - exp(2 * dLambda2 * dT)) / (-2 * dLambda2);
//
//    const double dTerminalCovar12 =
//        (dLambda1 + dLambda2) == 0.
//            ? dT * sqrt(_SIGMA1_Sqrd_TEST_(0, 0, 0, 0) * _SIGMA2_Sqrd_TEST_(0, 0, 0, 0)) * dRho
//            : sqrt(_SIGMA1_Sqrd_TEST_(0, 0, 0, 0) * _SIGMA2_Sqrd_TEST_(0, 0, 0, 0)) * dRho *
//                  (1 - exp(dLambda1 * dT + dLambda2 * dT)) / (-dLambda1 - dLambda2);
//
//    const double dTerminalVarX3 = (dTerminalVarX2 + dRho * dRho * dAlpha * dAlpha * dTerminalVarX1
//    -
//                                   2. * dRho * dAlpha * dTerminalCovar12) /
//                                  (dAlpha * dAlpha * (1. - dRho * dRho));
//
//    const double  dDeltaX1  = (2. * 7. * sqrt(dTerminalVarX1)) / (double)(nNumX1 - 1);
//    const double  dDeltaX3  = (2. * 7. * sqrt(dTerminalVarX3)) / (double)(nNumX2 - 1);
//    const double  dDeltaT   = (dTL - dTH) / ((double)nNumT);
//    double *      pdT_Begin = (double*)malloc((nNumT + 1) * sizeof(double)), *pdT = 0;
//    const double* pdT_End = pdT_Begin + nNumT + 1;
//
//    int *pnSize_X1 = (int*)malloc((nNumT + 1) * sizeof(int)), *pnX1;
//    int *pnSize_X3 = (int*)malloc((nNumT + 1) * sizeof(int)), *pnX3;
//
//    _init_Grid(nNumX1, dDeltaX1, nNumX2, dDeltaX3, nNumT, dTH - dTL, &grid);
//
//    // update pdU
//    for (nr = nR_Begin; nr < nR_End; ++nr)
//    {
//        for (nc = nC_Begin; nc < nC_End; ++nc)
//        {
//            ppdU[nr][nc] = *(grid.m_pdX3_Begin + nr - nR_Begin);
//            ppdU[nr][nc] *= ppdU[nr][nc];
//            ppdU[nr][nc] = nFlat ? 1. : ppdU[nr][nc];
//        }
//    }
//
//    for (pdT = pdT_Begin; pdT < pdT_End; ++pdT)
//        *pdT = dTH + dDeltaT * (pdT - pdT_Begin);
//
//    _ASSERTE(*(pdT_End - 1) == dTL);
//
//    for (pdT = pdT_Begin, pnX1 = pnSize_X1, pnX3 = pnSize_X3; pnX1 < pnSize_X1 + nNumT + 1;
//         ++pdT, ++pnX1, ++pnX3)
//    {
//        const double dT = *pdT;
//        const double dTerminalVarX1 =
//            (dLambda1 == 0.)
//                ? dT * _SIGMA1_Sqrd_TEST_(0, 0, 0, 0)
//                : _SIGMA1_Sqrd_TEST_(0, 0, 0, 0) * (1 - exp(2 * dLambda1 * dT)) / (-2 * dLambda1);
//
//        const double dTerminalVarX2 =
//            (dLambda2 == 0.)
//                ? dT * _SIGMA2_Sqrd_TEST_(0, 0, 0, 0)
//                : _SIGMA2_Sqrd_TEST_(0, 0, 0, 0) * (1 - exp(2 * dLambda2 * dT)) / (-2 * dLambda2);
//
//        const double dTerminalCovar12 =
//            (dLambda1 + dLambda2) == 0.
//                ? dT * sqrt(_SIGMA1_Sqrd_TEST_(0, 0, 0, 0) * _SIGMA2_Sqrd_TEST_(0, 0, 0, 0)) *
//                dRho : sqrt(_SIGMA1_Sqrd_TEST_(0, 0, 0, 0) * _SIGMA2_Sqrd_TEST_(0, 0, 0, 0)) *
//                dRho *
//                      (1 - exp(dLambda1 * dT + dLambda2 * dT)) / (-dLambda1 - dLambda2);
//
//        const double dTerminalVarX3 =
//            (dTerminalVarX2 + dRho * dRho * dAlpha * dAlpha * dTerminalVarX1 -
//             2. * dRho * dAlpha * dTerminalCovar12) /
//            (dAlpha * dAlpha * (1. - dRho * dRho));
//
//        // const double dDeltaX1 = (2.*7.*sqrt(dTerminalVarX1))/(double)(nNumX1-1);
//        // const double dDeltaX3 = (2.*7.*sqrt(dTerminalVarX3))/(double)(nNumX2-1);
//        _ASSERTE(dDeltaX1 != 0. && dDeltaX3 != 0.);
//        *pnX1 = 1 + (int)((2. * 7. * sqrt(dTerminalVarX1)) / dDeltaX1);
//        *pnX3 = 1 + (int)((2. * 7. * sqrt(dTerminalVarX3)) / dDeltaX3);
//        *pnX1 = _validate_gridpoints(*pnX1);
//        *pnX3 = _validate_gridpoints(*pnX3);
//
//        if (dT == (dTH - dTL))
//        {
//            _ASSERTE(*pnX1 == nNumX1);
//            _ASSERTE(*pnX3 == nNumX2);
//        }
//    }
//
//    szErr = ADI_MultipleSlices(
//        pdT_Begin,  // backwards
//        pdT_End,
//        grid.m_pdX1_Begin,
//        grid.m_pdX1_End,
//        grid.m_pdX3_Begin,
//        grid.m_pdX3_End,
//        0,  // FuncArg
//        _SIGMA1_Sqrd_TEST_,
//        _MU1_TEST_,
//        _SIGMA3_Sqrd_TEST_,
//        _MU3_TEST_,
//        ppdU,
//        nC_Begin,
//        nR_Begin);
//
//    if (szErr)
//        return szErr;
//
//    _print_matrix(ppdU, nR_Begin, nR_End, nC_Begin, nC_End);
//    ppdU[nR_Begin + grid.m_nX3_Zero][nC_Begin + grid.m_nX1_Zero];
//    free_dmatrix(ppdU, nR_Begin, nR_End - 1, nC_Begin, nC_End - 1);
//    return 0;
//}
//
// void _adi_tester_()
//{
//    int nNumX1 = 9;
//    int nNumX2 = 19;
//
//    const int nR_Begin = 0;
//    const int nR_End   = nR_Begin + nNumX2;
//    const int nC_Begin = 0;
//    const int nC_End   = nC_Begin + nNumX1;
//
//    //// test with no mean reversion
//    _alpha_  = 1.2;
//    _gamma_  = 0.;
//    _rho_    = -0.85;
//    _sigma_  = 2;
//    _lambda_ = 0;
//    // not flat
//    __adi_tester_(nR_Begin, nR_End, nC_Begin, nC_End, 0);
//
//    //// test with mean reversion
//    _alpha_  = 1.2;
//    _gamma_  = 0.2;
//    _rho_    = -0.85;
//    _sigma_  = 0.02;
//    _lambda_ = -0.5;
//    // flat
//    __adi_tester_(nR_Begin, nR_End, nC_Begin, nC_End, 1);
//}
///////////////////////// Debug related functions - End
//#endif  //_DEBUG
