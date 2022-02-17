//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : QuasiMQ.cpp
//
//   Description : First attempt at implementing quasi-vanilla in qlib
//                 subject to change in the future
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QuasiMQ.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IRModelVNFM.hpp"
#include "edginc/IRSmile2Q.hpp"
#include "edginc/IRSmileQuasi2Q.hpp"
#include "edginc/IRSmileMQ.hpp"

#include "edginc/MaturityPeriod.hpp"

#include <q3.h>
#include <q3_mqquasi.h>
#include <q3_mqbv.h>

DRLIB_BEGIN_NAMESPACE

///////////////////////////////
//
//  Quasi Vanilla Model
//
///////////////////////////////

/*f----------------------------------------------------------------------------
 * Q3MQQuasiInit
 *
 * Select payoff function and populate payoff parameters.
 */
//static void Q3MQQuasiInit(
//    long      optType,         /* (I) payoff type: see q3.h   */
//    double   *strike,          /* (I) strike                  */  
//    int      numPayoffParams, /* (I) number of payoff params */
//    double   *payoffParams,    /* (I) payoff coefficients     */
//    PAYOFF   *pf,              /* (O) option payoff structure */
//    FPAYOFF **payFunc          /* (O) payoff function         */
//    )
//{
//    try {
//        long  i;
//
//        /* number of payoff params */
//        if (numPayoffParams > Q3_MAX_PAY_PARAMS) {
//            throw ModelException(
//                "Number of payoff params ("+Format::toString(numPayoffParams)
//                +") exceeds max = "+Format::toString(Q3_MAX_PAY_PARAMS));
//        }
//
//        switch (optType) {
//
//            case Q3_ADJ_FWD:
//                *payFunc  = Q3Pay1D_Yield;
//                break;
//
//            case (Q3_VNL + Q3_PUT):
//            case (Q3_VNL + Q3_CALL):
//                *payFunc = Q3Pay1D_Vnl;
//                break;
//
//            case (Q3_ANN + Q3_PUT):
//            case (Q3_ANN + Q3_CALL):
//                if (numPayoffParams != 2) {
//                    throw ModelException(
//                        "Rate lock options require 2 payoff params (input: "
//                        + Format::toString(numPayoffParams)+").");
//                }
//
//                if (payoffParams[1] < TINY) {
//                    throw ModelException(
//                        "Rate maturity "+Format::toString(payoffParams[1])
//                        +" must be positive.");
//                }
//                *payFunc  = Q3Pay1D_Ann;
//                break;
//
//            case Q3_TEC:
//                if (numPayoffParams != 3) {
//                    throw ModelException("TEC options require 3 payoff params (input: "
//                        +Format::toString(numPayoffParams)+")");
//                }
//                if (*strike <= -1.) {
//                    throw ModelException("For TEC options, the strike ("
//                        +Format::toString(*strike)+") must be > -1");
//                }
//                *payFunc  = Q3Pay1D_TecFwd;
//                break;
//
//            case (Q3_TEC + Q3_PUT):
//            case (Q3_TEC + Q3_CALL):
//                if (numPayoffParams != 3) {
//                    throw ModelException("TEC options require 3 payoff params (input: "
//                        +Format::toString(numPayoffParams)+")");
//                }
//                if (*strike <= -1.) {
//                    throw ModelException("For TEC options, the strike ("
//                        +Format::toString(*strike)+") must be > -1");
//                }
//                *payFunc  = Q3Pay1D_Tec;
//                break;
//
//            case (Q3_BIN + Q3_PUT):
//            case (Q3_BIN + Q3_CALL):
//                *payFunc = Q3Pay1D_Bin;
//                break;
//
//            case (Q3_BIN + Q3_IN):
//            case (Q3_IN_BIRIB):
//                if (numPayoffParams != 4) {
//                    throw ModelException("Fixed RIB requires 4 payoff params (input: "
//                        +Format::toString(numPayoffParams)+").");
//                }
//                if (payoffParams[0] > payoffParams[1]) {
//                    throw ModelException("For fixed RIB, low barrier ("
//                        +Format::toString(payoffParams[0])+") must be less"
//                        "than high barrier ("+Format::toString(payoffParams[1])+")");
//                }
//                *payFunc = Q3Pay1D_FixedRibIn;
//                break;
//
//            case (Q3_BIN + Q3_OUT):
//            case (Q3_OUT_BIRIB):
//                if (numPayoffParams != 4) {
//                    throw ModelException("Fixed RIB requires 4 payoff params (input: "
//                        +Format::toString(numPayoffParams)+").");
//                }
//                if (payoffParams[0] > payoffParams[1]) {
//                    throw ModelException("For fixed RIB, low barrier ("
//                        +Format::toString(payoffParams[0])+") must be less"
//                        "than high barrier ("+Format::toString(payoffParams[1])+").");
//                }
//                *payFunc = Q3Pay1D_FixedRibOut;
//                break;
//
//            case Q3_POW:
//                if (numPayoffParams != 1) {
//                    throw ModelException("Pow requires 1 payoff param (input: "
//                        +Format::toString(numPayoffParams)+")");
//                }
//                *payFunc = Q3Pay1D_Pow;
//                break;
//
//            case (Q3_QUAD + Q3_PUT):
//            case (Q3_QUAD + Q3_CALL):
//                *payFunc = Q3Pay1D_Quad;
//                break;
//
//            default:
//                throw ModelException("Unsupported optType "+Format::toString((int)optType));
//        }
//
//        /* fill in payoff parameters */
//        pf->optType = optType;
//        pf->cop     = (Q3_COP_TYPE(optType) == Q3_CALL ? 1 : -1);/* cop_to_cop */
//    	pf->strike  = (strike == NULL ? 0. : *strike);
//        for (i=0; i<numPayoffParams; i++) pf->params[i] = payoffParams[i];
//    }
//    catch (exception& e) {
//        throw ModelException(e, __FUNCTION__);
//    }
//} /* Q3MQQuasiInit */
//
//
///**----------------------------------------------------------------------------
// * Q3MQQuasiPricer
// *
// * Multiq quasi-vanilla pricer: single-currency, 1d option
// */
//
//static void Q3MQQuasiPricer(
//    long              VolBaseDate,      /* 1  (I) ATM vol base date         */
//    long              SmileBaseDate,    /* 2  (I) Smile base date           */
//    long              CurveBaseDate,    /* 3  (I) Index curve base date     */
//    long              numIndxPts,       /* 4  (I) Nb of index points        */
//    long              *indxDates,       /* 5  (I) Index dates               */
//    double            *indxRates,       /* 6  (I) Act/365F index rates      */
//    long              numDiscPts,       /* 7  (I) Nb of discount points     */
//    long              *discDates,       /* 8  (I) Discount dates            */
//    double            *discRates,       /* 9 (I) Act/365F discount rates   */
//    long              numVnfmParams,    /* 10 (I) Number of VNFM params     */
//    double            *vnfmParams,      /* 11 (I) VNFM params               */
//    long              rateFreq,         /* 12 (I) rate frequency            */
//    char              *rateDCC,         /* 13 (I) rate day count conv       */
//    long              expiryDate,       /* 14 (I) rate reset date           */
//    long              startDate,        /* 15 (I) rate effective date       */
//    long              matDate,          /* 16 (I) rate maturity date        */
//    double            ScaledSigATM,     /* 17 (I) Scaled ATM volatility     */
//    double            bbATM,            /* 18 (I) ATM backbone              */
//    double            *smile,           /* 19 (I) MQ smile                  */
//    long              optType,          /* 20 (I) option type               */
//    double            *strike,          /* 21 (I) strike (NULL for Q3_ADJ)  */
//    long              numPayoffParams,  /* 22 (I) number of payoff coeffs   */
//    double            *payoffParams,    /* 23 (I) payoff coeffs             */
//    long              payDate,          /* 24 (I) option payment date       */
//    long              setlType,         /* 25 (I) cash or phys settle       */
//    char              *holidayFile,     /* 26 (I) holiday file name         */
//    char              *BusVolDCC,       /* 27 (I) bus day disc convention   */
//    double            *outputs)         /* 28 (O) fwd price and AA par rate */
//{
//    try
//    {
//        MQDATA           mq;
//        FADATA           fa;
//        PAYOFF           pf;
//        FPAYOFF          *payFunc = NULL;
//
//        double rate[11];
//        double fwdRate, expiry, expiryVolvol;
//        double fwdRateVNFM, fwdAnnVNFM, zeroRateSwap, zeroRatePay;
//        double start, swapMat, payDelay, freq, delayCapRate;
//        long   flagMethod;
//        double price;
//        
//        /* set up dates and forward rates */
//        Q3SwapRateCalc2(
//            VolBaseDate,
//            SmileBaseDate,
//            CurveBaseDate,
//            numIndxPts,
//            indxDates,
//            indxRates,
//            numDiscPts,
//            discDates,     
//            discRates,     
//            rateFreq,
//            rateDCC,
//            expiryDate,
//            startDate,
//            matDate,
//            payDate,
//            holidayFile,
//            BusVolDCC,
//            rate);
//
//        fwdRate      = rate[0];
//        fwdRateVNFM  = rate[1];
//        fwdAnnVNFM   = rate[2];
//        zeroRateSwap = rate[3];
//        zeroRatePay  = rate[4];
//        expiry       = rate[5];
//        expiryVolvol = rate[6];
//        start        = rate[7];
//        swapMat      = rate[8];
//        payDelay     = rate[9];
//        freq         = rate[10];
//
//        flagMethod = (long)(smile[0]+0.5) % 10;
//
//    /* smile: initialize */
//    if (Q3SmileInit(
//                    fwdRate,
//                    ScaledSigATM,
//                    bbATM,
//                    expiry,
//                    expiryVolvol,
//                    1.,
//                    smile,
//                    &mq) == FAILURE) throw QuasiMQ::makeException();
//
//        /* calibrate MQ */
//        if (Q3MQCalib(&mq) == FAILURE)
//            throw QuasiMQ::makeException();
//
//        /* FA measure initialization: sets alpha and power */
//        fa.mq = &mq;
//
//	    switch (flagMethod)
//	    {
//	     case Q3_VNFM_DELAY:
//	        if (Q3FASmileInit(
//                expiry,
//            	mq.sigATM,
//                start,
//                (long)freq,
//                swapMat,
//                fwdRateVNFM,
//                fwdAnnVNFM,
//                zeroRateSwap,
//                payDelay,
//                zeroRatePay,
//                numVnfmParams,
//                vnfmParams,
//                setlType,
//                Q3_PHYS_SETL,
//                &fa) == FAILURE) throw QuasiMQ::makeException();
//	        break;
//	    
//	     case Q3_VSK_DELAY:
//	        delayCapRate = smile[10];   // set up the cap rate for the delay
//	        if (Q3FASmileInitVsk(
//	            expiry,
//	            mq.sigATM,
//	            start,
//	            (long)freq,
//	            swapMat,
//	            fwdRateVNFM,
//	            fwdAnnVNFM,
//	            zeroRateSwap,
//	            payDelay,
//	            delayCapRate,
//	            zeroRatePay,
//	            numVnfmParams,
//	            vnfmParams,
//	            setlType,
//	            &fa) == FAILURE) throw QuasiMQ::makeException();
//	        break;
//	
//	    default:
//	        throw QuasiMQ::makeException();
//	    
//	    }
//     
//        /* select payoff function, populate payoff structure */
//        Q3MQQuasiInit(
//            optType,
//            strike,
//            numPayoffParams,
//            payoffParams,
//            &pf,
//            &payFunc);
//        
//        /* price option (undiscounted) */
//	    switch (flagMethod)
//	    {
//	    case Q3_VNFM_DELAY:
//	        if (Q3FAGridPricer(
//                &pf, 
//                payFunc,  
//                &fa,
//                &price) == FAILURE) throw QuasiMQ::makeException();
//        break;
//    
//	    case Q3_VSK_DELAY:
//	        if (Q3FAGridPricerVsk(
//	            &pf, 
//	            payFunc,  
//	            &fa,
//	            &price
//	            ) == FAILURE) throw QuasiMQ::makeException();
//	        break;
//	
//	    default:
//	        throw QuasiMQ::makeException();
//	    
//	    }
//            
//        /* AA par rate (swap rate) */
//        outputs[0] = price ;
//        outputs[1] = fwdRate;
//    }
//    catch (exception& e) {
//        throw ModelException(e, __FUNCTION__);
//    }
//} /* Q3MQQuasiPricer */
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////  Bivariate Quasi Model
////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//static const int Q3_BV_NCK = 1090505665; // calibration control data
//
//// function prototypes
//static void Q3MQCalibJoint(double* corr, long index, MQDATA** mq);
//static void Q3MQCalibJointFwd(double target, double* corr, long index, MQDATA**  mq);
//
///*-----------------------------------------------------------------------------
// * Q3MQBivarInit
// *
// * Select payoff function and populate payoff parameters. 
// * Routine may reorder measure data elements.
// *
// */
//static void Q3MQBivarInit(long      optType,    /* (I)   option type             */
//                          long      numPayPrms, /* (I)   number of payoff params */
//                          double*   payPrms,    /* (I)   payoff coefficients     */
//                          double*   corr,       /* (I)   index correlation       */
//                          MQDATA**  mq,         /* (I/O) measure data            */
//                          PAYOFF*   pf,         /* (O)   option payoff structure */
//                          FPAYOFF** payFunc)    /* (O)   payoff function         */
//{
//    try
//    {
//        long        calJoint    = FALSE;
//        long        calIndex = 0;
//        MQDATA     *mqTmp;
//        long        i;
//        double      RIB1, RIB2;
//            
//        /* number of payoff params */
//        if (numPayPrms > Q3_MAX_PAY_PARAMS)
//        {
//            throw ModelException("Number of payoff params exceeds maximum.");
//        }
//
//        switch (optType)
//        {   
//
//        case Q3_ADJ_FWD:
//
//            *payFunc = Q3Pay1D_YldNull;
//            break;
//
//        case Q3_JOINT_FWD:
//        case Q3_BS_JOINT_FWD:
//            
//            *payFunc = Q3Pay1D_YldYld;
//            break;
//
//        case Q3_VNL + Q3_CALL:
//
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_VnlNull;
//            pf->cop  = 1;
//            break;
//
//        case Q3_VNL + Q3_PUT:
//
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_VnlNull;
//            pf->cop  = -1;
//            break;
//
//        case Q3_CALL_SUM:           
//
//            if (numPayPrms < 2)
//            {
//                throw ModelException("Payoff requires 2 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_Sum;
//            pf->cop  = 1;
//            break;
//
//        case Q3_PUT_SUM:           
//
//            if (numPayPrms < 2)
//            {
//                throw ModelException("Payoff requires 2 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_Sum;
//            pf->cop  = -1;
//            break;
//
//        case Q3_CALL_PROD:           
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_Prod;
//            pf->cop  = 1;
//            break;
//
//        case Q3_PUT_PROD:           
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_Prod;
//            pf->cop  = -1;
//            break;
//        
//        case Q3_CALL_PERC:           
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 params.");
//            }
//            pf->params[0] = payPrms[0];
//            pf->params[1] = 0.;  /* zero spread */
//            *payFunc = Q3Pay1D_Perc;
//            pf->cop  = 1;
//            break;
//        
//        case Q3_PUT_PERC:           
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 params.");
//            }
//            pf->params[0] = payPrms[0];
//            pf->params[1] = 0.; /* zero spread */
//            *payFunc = Q3Pay1D_Perc;
//            pf->cop  = -1;
//            break;
//
//        case Q3_CALL_PERC_WGT:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_PercWgt;
//            pf->cop  = 1;
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//        
//        case Q3_PUT_PERC_WGT:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_PercWgt;
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            pf->cop  = -1;
//            break;
//
//        case Q3_FLR_W_FLR:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrFlrOrCapCap;
//            pf->cop  = 1;
//            break;
//
//        case Q3_CAP_W_CAP:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrFlrOrCapCap;
//            pf->cop  = -1;
//            break;
//
//        case Q3_FLR_W_CAP:           
//            
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrCapOrCapFlr;
//            pf->cop  = 1;
//            break;
//        
//        case Q3_CAP_W_FLR:           
//            
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrCapOrCapFlr;
//            pf->cop  = -1;
//            break;    
//        
//        case Q3_FLR_W_FLR_EMBED:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrFlrOrCapCapEmbedFlt;
//            pf->cop  = 1;
//            break;
//        
//        case Q3_CAP_W_CAP_EMBED:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrFlrOrCapCapEmbedFlt;
//            pf->cop  = -1;
//            break;
//
//        case Q3_FLR_W_CAP_EMBED:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrCapOrCapFlrEmbedFlt;
//            pf->cop  = 1;
//            break;
//        
//        case Q3_CAP_W_FLR_EMBED:           
//
//            if (numPayPrms < 3)
//            {
//                throw ModelException("Payoff requires 3 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrCapOrCapFlrEmbedFlt;
//            pf->cop  = -1;
//            break;
//
//        case Q3_FLR_CAP_SUM:
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_FlrCapSum;
//            pf->cop  = 1;
//            break;
//
//        case Q3_BS_PERC_RATE_CALL:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_Prod;
//            pf->cop  = 1;
//            calJoint = TRUE;
//            calIndex = 1;
//            break;
//
//        case Q3_BS_PERC_RATE_PUT:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_Prod;
//            pf->cop  = -1;
//            calJoint = TRUE;
//            calIndex = 1;
//            break;
//
//        case Q3_BS_PERC_CALL:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires at least 1 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            if ( numPayPrms == 1)
//            {
//                /* Backward compatible */
//                pf->params[1] = 0.;
//            }
//            *payFunc = Q3Pay1D_Perc;
//            pf->cop  = 1;
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            calIndex = 0;
//            calJoint = TRUE;
//            break;
//
//        case Q3_BS_PERC_PUT:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires at least 1 params.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            if ( numPayPrms == 1)
//            {
//                /* Backward compatible */
//                pf->params[1] = 0.;
//            }
//            *payFunc = Q3Pay1D_Perc;
//            pf->cop  = -1;
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            calIndex = 0;
//            calJoint = TRUE;
//            break;
//
//        case Q3_BS_SPRD_RATE_CALL:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            pf->params[1] = 1.;         /* payFunc expects leverage */
//            *payFunc      = Q3Pay1D_Sum;
//            pf->cop       = 1;
//            break;
//            
//        case Q3_BS_SPRD_RATE_PUT:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            pf->params[1] = 1.;         /* payFunc expects leverage */
//            *payFunc      = Q3Pay1D_Sum;
//            pf->cop       = -1;
//            break;
//
//        case Q3_BS_PERC_SPRD_CALL:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_VnlNull; /* Vanilla on 1st variable */
//            pf->cop  = 1.;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            calJoint = TRUE;
//            calIndex = 0;
//            break;
//
//        case Q3_BS_PERC_SPRD_PUT:
//
//            if (numPayPrms < 1)
//            {
//                throw ModelException("Payoff requires 1 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_VnlNull;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            pf->cop  = -1.;
//            calJoint = TRUE;
//            calIndex = 0;
//            break;
//
//        case Q3_BS_SPRD_DIG_RIBIN:
//            if (numPayPrms < 2)
//            {
//                throw ModelException("Payoff requires 2 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            /* payFunc expects leps, heps, weight1, weight2 */
//            pf->params[2] = 0.;     
//            pf->params[3] = 0.;
//            pf->params[4] = 1.;
//            pf->params[5] = 1.;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            *payFunc      = Q3Pay1D_DigSpdInRIB_EPS;
//            break;
//
//        case Q3_BS_SPRD_DIG_RIBOUT:
//            if (numPayPrms < 2)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            /* payFunc expects leps, heps, weight1, weight2 */
//            pf->params[2] = 0.;     
//            pf->params[3] = 0.;
//            pf->params[4] = 1.;
//            pf->params[5] = 1.;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            *payFunc      = Q3Pay1D_DigSpdOutRIB_EPS;
//            break;
//
//        case Q3_BS_SPRD_DIG_RIBIN_EPS:
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            /* payFunc expects weight1, weight2 */
//            pf->params[4] = 1.;
//            pf->params[5] = 1.;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            *payFunc      = Q3Pay1D_DigSpdInRIB_EPS;
//
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
//                (pf->params[3] < 0 && -pf->params[3] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            break;
//
//        case Q3_BS_SPRD_DIG_RIBOUT_EPS:
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 1 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            /* payFunc expects weight1, weight2 */
//            pf->params[4] = 1.;
//            pf->params[5] = 1.;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            *payFunc      = Q3Pay1D_DigSpdOutRIB_EPS;
//            
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
//                (pf->params[3] < 0 && -pf->params[3] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            break;
//
//        case Q3_BS_PERC_DIG_RIBIN:
//            if (numPayPrms < 2)
//            {
//                throw ModelException("Payoff requires 2 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            /* payFunc expects leps, heps */
//            pf->params[2] = 0.;     
//            pf->params[3] = 0.;
//            *payFunc      = Q3Pay1D_DigPercInRIB_EPS;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            calJoint = TRUE;
//            calIndex = 0;
//            break;
//
//        case Q3_BS_PERC_DIG_RIBOUT:
//            if (numPayPrms < 2)
//            {
//                throw ModelException("Payoff requires 2 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            /* payFunc expects leps, heps */
//            pf->params[2] = 0.;     
//            pf->params[3] = 0.;
//            *payFunc      = Q3Pay1D_DigPercOutRIB_EPS;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            calJoint = TRUE;
//            calIndex = 0;
//            break;
//
//        case Q3_BS_PERC_DIG_RIBIN_EPS:
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc      = Q3Pay1D_DigPercInRIB_EPS;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            calJoint = TRUE;
//            calIndex = 0;
//
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
//                (pf->params[3] < 0 && -pf->params[3] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            break;
//
//        case Q3_BS_PERC_DIG_RIBOUT_EPS:
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc      = Q3Pay1D_DigPercOutRIB_EPS;
//            /* payoff function expects order: rate1, rate2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            calJoint = TRUE;
//            calIndex = 0;
//
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[2] < 0 && -pf->params[2] > RIB1) ||
//                (pf->params[3] < 0 && -pf->params[3] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            break;
//
//        case Q3_IN_BIRIB:
//
//            if (numPayPrms < 6)
//            {
//                throw ModelException("Payoff requires 6 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_InRIB;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_BIRIB_EPS:
//            
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            /* Ignore small epsilons */
//            if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
//                (fabs(payPrms[7]) < Q3_MIN_EPS))
//            {
//                *payFunc = Q3Pay1D_InRIB;
//            }
//            else
//            {
//                *payFunc = Q3Pay1D_InRIB_EPS;
//            }
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//                
//        case Q3_OUT_BIRIB:
//
//            if (numPayPrms < 6)
//            {
//                throw ModelException("Payoff requires 6 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_OutRIB;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_BIRIB_EPS:
//            
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            /* Ignore small epsilons */
//            if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
//                (fabs(payPrms[7]) < Q3_MIN_EPS))
//            {
//                *payFunc = Q3Pay1D_OutRIB;
//            }
//            else
//            {
//                *payFunc = Q3Pay1D_OutRIB_EPS;
//            }
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_SPDRIB:
//
//            if (numPayPrms < 6)
//            {
//                throw ModelException("Payoff requires 6 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_InSpdRIB;
//            /* swap measures to condition on 1st index in payoff*/
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            break;
//
//        case Q3_IN_SPDRIB_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            /* Ignore small epsilons */
//            if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
//                (fabs(payPrms[7]) < Q3_MIN_EPS))
//            {
//                *payFunc = Q3Pay1D_InSpdRIB;
//            }
//            else
//            {
//                *payFunc = Q3Pay1D_InSpdRIB_EPS;
//            }
//            /* swap measures to condition on 1st index in payoff*/
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            break;
//
//        case Q3_OUT_SPDRIB:
//
//            if (numPayPrms < 6)
//            {
//                throw ModelException("Payoff requires 6 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_OutSpdRIB;
//            /* swap measures to condition on 1st index in payoff*/
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_SPDRIB_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            if(RIB1 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            if( (pf->params[6] < 0 && -pf->params[6] > RIB1) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB1))
//            {
//                throw ModelException("When leps < 0 or heps < 0, leps (or heps) < 0.5 *(HB - LB).");
//            }
//            /* Ignore small epsilons */
//            if ((fabs(payPrms[6]) < Q3_MIN_EPS) &&
//                (fabs(payPrms[7]) < Q3_MIN_EPS))
//            {
//                *payFunc = Q3Pay1D_OutSpdRIB;
//            }
//            else
//            {
//                *payFunc = Q3Pay1D_OutSpdRIB_EPS;
//            }
//            /* swap measures to condition on 1st index in payoff*/
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;  
//            break;
//
//        case Q3_IN_AND_IN:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_InAndIn;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_AND_IN_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//            *payFunc = Q3Pay1D_InAndIn;
//
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_OR_IN:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_InOrIn;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_OR_IN_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_InOrIn;
//    	    
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_AND_OUT:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_InAndOut;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_AND_OUT_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_InAndOut;
//    	    
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_OR_OUT:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_InOrOut;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_IN_OR_OUT_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_InOrOut;
//    	    
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_AND_IN:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_OutAndIn;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_AND_IN_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_OutAndIn;
//
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_OR_IN:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_OutOrIn;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_OR_IN_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_OutOrIn;
//    	    
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_AND_OUT:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_OutAndOut;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_AND_OUT_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_OutAndOut;
//    	    
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_OR_OUT:
//
//            if (numPayPrms < 4)
//            {
//                throw ModelException("Payoff requires 4 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            for (i=4; i<8; i++) pf->params[i] = 0.0;
//            *payFunc = Q3Pay1D_OutOrOut;
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_OUT_OR_OUT_EPS:
//
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 8 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            RIB1 = pf->params[1] - pf->params[0];
//            RIB2 = pf->params[3] - pf->params[2];
//            if( RIB1 < 0 || RIB2 < 0)
//            {
//                throw ModelException("Requires HB > LB.");
//            }
//            RIB1 /= 2.0;
//            RIB2 /= 2.0;
//            if( (pf->params[4] < 0 && -pf->params[4] > RIB1) ||
//                (pf->params[5] < 0 && -pf->params[5] > RIB1) ||
//                (pf->params[6] < 0 && -pf->params[6] > RIB2) ||
//                (pf->params[7] < 0 && -pf->params[7] > RIB2) )
//            {
//                throw ModelException("When leps < 0 or heps < 0, epsilon < 0.5 *(HB - LB). ");
//            }
//    	    *payFunc = Q3Pay1D_OutOrOut;
//    	    
//            /* payoff function expects order: rate1, rate 2 */
//            mqTmp = mq[0];
//            mq[0] = mq[1];
//            mq[1] = mqTmp;
//            break;
//
//        case Q3_COMPLEX_SPD:
//
//            if (numPayPrms < 19)
//            {
//                throw ModelException("Payoff requires 19 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_ComplexSPD;
//            break;
//
//        case Q3_BS_PERC_YLD_ANN:
//            
//            if (numPayPrms < 8)
//            {
//                throw ModelException("Payoff requires 1 param.");
//            }
//            for (i=0; i<numPayPrms; i++) pf->params[i] = payPrms[i];
//            *payFunc = Q3Pay1D_SpdYldAnn;
//            break;
//
//        default: 
//            
//            throw ModelException("Unrecognized option type.");
//
//        }
//        
//        /* joint calibration */
//        if (calJoint == TRUE)
//        { 
//            /* if no user calibration parameter, fwd cal. only */
//            if (numPayPrms > 2 && fabs(payPrms[2] - Q3_BS_CAL_VOL) < TINY)
//            {
//                Q3MQCalibJoint(corr, calIndex, mq);
//            }
//            else
//            {
//                Q3MQCalibJointFwd((mq[0])->fwdRate * (mq[1])->fwdRate, corr, calIndex, mq);
//            }
//        }
//
//        /* fill in common payoff parameters */
//        pf->optType = optType;
//        pf->corr    = corr[0];
//        (pf->mq)[0] = mq[0];
//        (pf->mq)[1] = mq[1];
//    }
//    catch (exception& e)
//    {
//        throw ModelException(e, __FUNCTION__);
//    }
//
//} /* Q3MQBivarInit */
//
///*-----------------------------------------------------------------------------
// * Q3MQCalibJoint
// */
//#define Q3_S_DELTA1 10*Q3_S_DELTA
//
//static void Q3MQCalibJoint(double* corr, long index, MQDATA** mq)
//{
//    try
//    {
//        double targetVolPrice, targetFwd, compositeVol;
//        double sigATM0, sigATM1, expiry;
//        double volPrice, volPriceTwk, volPriceDelta, diff;
//        long   i;
//
//        PAYOFF   pf;
//        FPAYOFF *payFunc;    
//
//        /* local variables */
//        sigATM0 = (mq[0])->sigATM; sigATM1 = (mq[1])->sigATM;
//        expiry  = (mq[0])->expiry;
//
//        /* compute forward target */
//        targetFwd = (mq[0])->fwdRate * (mq[1])->fwdRate;
//
//        /* composite volatility */
//        compositeVol = sigATM0*sigATM0 + sigATM1*sigATM1 + 2.*corr[0]*sigATM0*sigATM1;
//        compositeVol = sqrt(compositeVol);
//
//        /* compute volatility target */
//        if (BSQPricer(
//            targetFwd,    /* fwd       */
//            targetFwd,    /* strike    */
//            expiry,
//            compositeVol,
//            1,            /* lognormal */
//            Q3_CALL,
//            &targetVolPrice) == FAILURE) throw ModelException("BSQPricer failed");
//
//        /* initialization for bivariate pricing */
//        Q3MQBivarInit(Q3_CALL_PROD,
//                      1,
//                      &targetFwd,
//                      corr,
//                      mq,
//                      &pf,
//                      &payFunc);
//
//        /* Newton-Raphson to match targets */
//        for (i = 0; i < Q3_S_STEPS; ++i)
//        {
//            /* calibrate forward */
//            Q3MQCalibJointFwd(targetFwd, corr, index, mq);
//
//            /* price vol */
//            if (Q3MQGridPricer(
//                &pf,
//                payFunc,
//                mq[0],
//                &volPrice) == FAILURE) throw ModelException("Q3MQGridPricer failed");
//        
//            diff = targetVolPrice - volPrice;
//
//            if (fabs(diff) > .01 * targetVolPrice)
//            {
//                /* N-R */
//                (mq[index])->sigMQ += Q3_S_DELTA1;
//
//                Q3MQCalibJointFwd(targetFwd, corr, index, mq);
//
//                if (Q3MQGridPricer(
//                    &pf,
//                    payFunc,
//                    mq[0],
//                    &volPriceTwk) == FAILURE) throw ModelException("Q3MQGridPricer failed");
//                
//                (mq[index])->sigMQ -= Q3_S_DELTA1;
//                volPriceDelta = volPrice - volPriceTwk;
//                if (fabs(volPriceDelta) < Q3_MQ_RESN * Q3_S_DELTA1)
//                {
//                    throw ModelException("df/ds = 0 in Newton-Raphson.");
//                }
//                
//                (mq[index])->sigMQ -= diff*Q3_S_DELTA1/volPriceDelta;
//            }
//            else
//            {
//                Q3MQCalibJointFwd(targetFwd, corr, index, mq);
//                break;
//            }
//        } /* N-R loop */
//    }
//    catch (exception& e)
//    {
//        throw ModelException(e, __FUNCTION__);
//    }
//
//} /* Q3MQCalibJoint */
//
///*-----------------------------------------------------------------------------
// * Q3MQCalibJointFwd
// *
// * Calibrate joint forward to target by adjusting forward of selected index.
// *
// */
//static void Q3MQCalibJointFwd(double   target, /* (I)   calibration target     */
//                              double*  corr,   /* (I)   correlation of indices */
//                              long     index,  /* (I)   index to adjust (0,1)  */
//                              MQDATA** mq)     /* (I/O) index measure data     */
//{
//    try
//    {
//        PAYOFF   pf;
//        FPAYOFF *payFunc;    
//        double   jointFwd;
//
//        /* calculate joint forward */
//        Q3MQBivarInit(Q3_JOINT_FWD,
//                      0,              /* num pay params */
//                      NULL,
//                      corr,
//                      mq,
//                      &pf,
//                      &payFunc);
//
//        if (Q3MQGridPricer(
//            &pf,
//            payFunc,
//            mq[0],
//            &jointFwd) == FAILURE) throw ModelException("Q3MQGridPricer failed");
//
//        /* fail on zero forward */
//        if (jointFwd < TINY) throw ModelException("jointFwd < TINY");
//
//        /* adjust selected forward */
//        (mq[index])->fwdRate *= target / jointFwd;
//    }
//    catch (exception& e)
//    {
//        throw ModelException(e, __FUNCTION__);
//    }
//} /* Q3MQCalibJointFwd */
//
///*-----------------------------------------------------------------------------
// * Q3MQPSAPricer
// */
//static void Q3MQPSAPricer(
//    long    today,           /*  1 (I) Vol base date                        */
//    long    valueDate,       /*  2 (I) Value date                           */
//    long    bsFixLegFreq,    /*  3 (I) BMA and CMS fixed leg frequency      */
//    char   *bsFixLegDCC,     /*  4 (I) BMA and CMS fixed Leg DCC            */
//    long    bsFltLegFreq,    /*  5 (I) BMA and CMS float leg frequency      */
//    char   *bsFltLegDCC,     /*  6 (I) BMA and CMS float leg DCC            */
//    long    numLiborZeroPts, /*  7 (I) LIBOR zero curve num points          */
//    long   *liborZeroDates,  /*  8 (I) LIBOR zero curve dates               */
//    double *liborZeroRates,  /*  9 (I) LIBOR zero curve rates               */
//    long    numBsZeroPts,    /* 10 (I) basis zero curve num points          */
//    long   *bsZeroDates,     /* 11 (I) basis zero curve dates               */
//    double *bsZeroRates,     /* 12 (I) basis zero curve rates               */
//    long    expiryDate,      /* 13 (I) reset date                           */
//    long    startDate,       /* 14 (I) start date                           */
//    long    matDate,         /* 15 (I) end date                             */
//    double *cmsSmile,        /* 16 (I) CMS smile parameters                 */
//    double  cmsSigATM,       /* 17 (I) CMS ATM vol                          */
//    long    numVNFMParams,   /* 18 (I) 1st index num VNFM parameters        */
//    double *vnfmParams,      /* 19 (I) 1st index VNFM parameters            */
//    double *bsSpreadSmile,   /* 20 (I) spread smile parameters              */
//    double  bsSpreadSigATM,  /* 21 (I) spread ATM vol                       */
//    long    numDiscPts,      /* 22 (I) Discount curve num points            */
//    long   *discDates,       /* 23 (I) Discount curve dates                 */
//    double *discRates,       /* 24 (I) Discount curve rates                 */
//    double  corr,            /* 25 (I) Index correlation data               */
//    long    optType,         /* 26 (I) Option type                          */
//    long    payDate,         /* 27 (I) Option payment date                  */
//    long    setlType,        /* 28 (I) Cash or physical settlement          */
//    long    numPayoffParams, /* 29 (I) Num payment parameters               */
//    double *payoffParams,    /* 28 (I) Payment parameters                   */
//    char   *holidayfile,     /* 29 (I) holiday file name                    */
//    char   *BusVolDCC,       /* 30 (I) BUS/251F or BUS/BUS                  */
//    double *results          /* 31 (O) Pricing results.                     */ 
//    )
//{
//    try
//    {
//        union{
//            struct{
//                double fwdRate;
//                double fwdRate30360;
//                double fwdAnn30360;
//                double zeroRateSwap;
//                double zeroRatePay;
//                double expiry;
//                double expiryVolvol;
//                double start;
//                double swapMat;
//                double payDelay;
//                double freq;
//            } s;
//            double a[11];
//        } rate[2];
//
//        FPAYOFF  *payFunc       = NULL;
//        double   *smile         = NULL;
//
//        PAYOFF    pf;
//        FADATA    fa;       /* CMS Yield distribution in FA measure*/
//        MQDATA    mq[2];
//        MQDATA    pa;       /* re-parameterize CMS Yield distribution in FA */
//        MQDATA   *ppa[2];
//        double    sigATM;
//        double    price;
//        int       i;
//
//        double    basisParRate, basisLegFv;
//        double    payPrms[8];
//        double    AnnFv, zeroDelay;
//        double    FADensNorm, AAtoFAScaleAdj = 1.;
//
//        /* Check optType validity */ 
//        if (optType != Q3_BS_PERC_CALL &&
//            optType != Q3_BS_PERC_PUT &&
//            optType != Q3_BS_PERC_RATE_CALL &&
//            optType != Q3_BS_PERC_RATE_PUT &&
//            optType != Q3_BS_JOINT_FWD)
//        {
//            throw ModelException("Invalid option type.");
//        }
//
//        /* Fwd dates and rates for CMS. */
//        if (Q3SwapRateCalc2(
//                today,
//                today,
//                valueDate,
//                numLiborZeroPts,
//                liborZeroDates,
//                liborZeroRates,
//                numDiscPts,
//                discDates,
//                discRates,
//                bsFixLegFreq,
//                bsFixLegDCC,
//                expiryDate,
//                startDate,
//                matDate,
//                payDate,
//                holidayfile,
//                BusVolDCC,
//                (rate[0]).a) == FAILURE) throw ModelException("Q3SwapRateCalc2 failed");
//
//        /* Basis par swap yield, basis leg PV*/
//        if (Q3BasisLegCalc(
//                valueDate,
//                numBsZeroPts,
//                bsZeroDates,
//                bsZeroRates,
//                numDiscPts,
//                discDates,
//                discRates,
//                (long)rate[0].s.freq,
//                bsFixLegDCC,
//                bsFltLegFreq,
//                bsFltLegDCC,
//                expiryDate,
//                startDate,
//                matDate,
//                payDate,
//                holidayfile,
//                (rate[1].a)) == FAILURE) throw ModelException("Q3BasisLegCalc failed");
//
//
//        if (cmsSigATM * sqrt((rate[0].s.expiry)) < Q3_MIN_VOL_CALIB ||
//            bsSpreadSigATM * sqrt(rate[0].s.expiry) < Q3_MIN_VOL_CALIB) corr = 0.;
//
//        basisParRate = rate[1].a[0];
//        basisLegFv   = rate[1].a[1];
//        AnnFv        = rate[1].a[2];
//        zeroDelay      = rate[1].a[3];
//
//        /* Assign spread expiry, expiryVolvol, etc.*/
//        rate[1].s.expiry = rate[0].s.expiry;
//        rate[1].s.expiryVolvol = rate[0].s.expiryVolvol;
//
//        /* spread fwd = basisSwapYld / CMSYld*/
//        rate[1].s.fwdRate = basisParRate / rate[0].s.fwdRate;
//
//        for (i = 0; i <= 1; i++)
//        {
//            /* select index, note ordering is reversed */
//            if (i == 0)
//            {
//                sigATM          = cmsSigATM;
//                smile           = cmsSmile;
//            }
//            else
//            {
//                sigATM          = bsSpreadSigATM;
//                smile           = bsSpreadSmile;
//            }
//            
//	         /* smile: initialize */
//	         if (Q3SmileInit(
//	                    (rate[i]).s.fwdRate,                 
//	                    sigATM,
//	                    0.,
//	                    (rate[i]).s.expiry,
//	                    (rate[i]).s.expiryVolvol,
//	                    1.,
//	                    smile,
//                         &(mq[i])) == FAILURE)
//            {
//                throw ModelException("Error initialising smile");
//            }
//            
//            /* Calibrate MQ distribution to SV option prices. */
//            if (Q3MQCalib(&(mq[i])) == FAILURE) throw ModelException("Q3MQCalib failed");
//        }
//
//        /* Calibrate CMS yield distribution in FA measure    */
//        fa.mq = &(mq[0]);
//        if (Q3FASmileInit(
//                (rate[0]).s.expiry,
//            	mq[0].sigATM,
//                (rate[0]).s.start,
//                (long) (rate[0]).s.freq,
//                (rate[0]).s.swapMat,
//                (rate[0]).s.fwdRate30360,
//                (rate[0]).s.fwdAnn30360,
//                (rate[0]).s.zeroRateSwap,
//                (rate[0]).s.payDelay, 
//                (rate[0]).s.zeroRatePay,
//                numVNFMParams,
//                vnfmParams,
//                setlType,
//                setlType,
//                &fa) == FAILURE) throw ModelException("Q3FASmileInit failed");
//         
//        /* initialize PA MQ */
//        if (Q3MQCopySmileFromMQ(
//                &(mq[0]),
//                &pa) == FAILURE) throw ModelException("Q3MQCopySmileFromMQ failed");
//           
//        /* set NCK parameters */
//        if (Q3DecodeNCK(
//                Q3_BV_NCK, 
//                &pa) == FAILURE) throw ModelException("Q3DecodeNCK failed");
//
//        /* compute calibration targets */
//        if (Q3MQTargetFA( 
//                &fa,
//                &pa) == FAILURE) throw ModelException("Q3MQTargetFA failed");
//
//        /* bootstrap PA multi-q measure */
//        if (Q3MQBootstrapQ(&pa) == FAILURE) throw ModelException("Q3MQBootstrapQ failed");
//
//        /* calibrate PA */
//        if (Q3MQCalib(&pa) == FAILURE) throw ModelException("Q3MQCalib failed");
//
//        if (Q3FADensNorm(&fa, &FADensNorm) == FAILURE) throw ModelException("Q3FADensNorm failed");
//        
//        AAtoFAScaleAdj = basisParRate / ( FADensNorm * basisLegFv );
//       
//        /* Adjust spread forward:
//        * Compute E[spread * (1 - Z(T, T+M)/Z(T, Pay)] 
//        */
//        ppa[0] = &pa;
//        ppa[1] = &(mq[1]);
//        payPrms[0] = fa.alphaAnn;
//        payPrms[1] = fa.powerAnn;
//        payPrms[2] = fa.freqAnn;
//        payPrms[3] = fa.matAnn;
//        payPrms[4] = fa.alphaDel;
//        payPrms[5] = fa.powerDel;
//        payPrms[6] = fa.freqDel;
//        payPrms[7] = fa.matDel;
//
//        Q3MQBivarInit(Q3_BS_PERC_YLD_ANN,
//                      8,
//                      payPrms,
//                      &corr,
//                      &(ppa[0]),
//                      &pf,
//                      &payFunc);
//        /* integrate payoff over density. */
//        if (Q3MQGridPricer(
//            &pf,
//            payFunc,
//            ppa[0],
//            &price) == FAILURE) throw ModelException("Q3MQGridPricer failed");
//        
//        /* Adjust spread forward :
//        * The basis float leg FV is scaled first due to
//        * numerical scaling from AA to FA.
//        * See details in Changhong He's Doc.
//        */
//        (mq[1]).fwdRate *= AAtoFAScaleAdj * basisLegFv / price;
//
//        /* Price under the payment measure */
//        ppa[0] = &pa;
//        ppa[1] = &(mq[1]);      
//        Q3MQBivarInit(optType,
//                      numPayoffParams,
//                      payoffParams,
//                      &corr,
//                      &(ppa[0]),
//                      &pf,
//                      &payFunc);
//        /* integrate payoff over density. */
//        if (Q3MQGridPricer(
//            &pf,
//            payFunc,
//            ppa[0],
//            &price) == FAILURE) throw ModelException("Q3MQGridPricer failed");
//
//        results[0] = price;
//        results[1] = mq[1].fwdRate;     /* spread forward */
//        results[2] = pa.fwdRate;        /* CMS forward    */
//    }
//    catch (exception& e)
//    {
//        throw ModelException(e, __FUNCTION__);
//    }
//}/* Q3MQPSAPricer */
//
//
///*-----------------------------------------------------------------------------
// * Q3MQBivarPricer
// * 
// * Multi-q bivariate option pricer.
// *
// */
//static void Q3MQBivarPricer(
//    long    today,           /* (I) Vol base date                        */
//    long    valueDate,       /* (I) Value date                           */
//    long    numIndxPts1,     /* (I) 1st index curve num points           */
//    long   *indxDates1,      /* (I) 1st index curve dates                */
//    double *indxRates1,      /* (I) 1st index curve rates                */
//    long    rateFreq1,       /* (I) 1st index curve frequency            */
//    char   *rateDCC1,        /* (I) 1st index curve day count convention */
//    long    expiryDate1,     /* (I) 1st index reset date                 */
//    long    startDate1,      /* (I) 1st index start date                 */
//    long    matDate1,        /* (I) 1st index rate end date              */
//    double *smile1,          /* (I) 1st index smile parameters           */
//    double  sigATM1,         /* (I) 1st index ATM vol                    */
//    long    numVNFMParams1,  /* (I) 1st index num VNFM parameters        */
//    double *vnfmParams1,     /* (I) 1st index VNFM parameters            */
//    long    numIndxPts2,     /* (I) 2nd index curve num points           */
//    long   *indxDates2,      /* (I) 2nd index curve dates                */
//    double *indxRates2,      /* (I) 2nd index curve rates                */
//    long    rateFreq2,       /* (I) 2nd index curve frequency            */
//    char   *rateDCC2,        /* (I) 2nd index curve day count convention */
//    long    expiryDate2,     /* (I) 2nd index reset date                 */
//    long    startDate2,      /* (I) 2nd index start date                 */
//    long    matDate2,        /* (I) 2nd index rate end date              */
//    double *smile2,          /* (I) 2nd index smile parameters           */
//    double  sigATM2,         /* (I) 2nd index ATM vol                    */
//    long    numVNFMParams2,  /* (I) 2nd index num VNFM parameters        */
//    double *vnfmParams2,     /* (I) 2nd index VNFM parameters            */
//    long    numDiscPts,      /* (I) Discount curve num points            */
//    long   *discDates,       /* (I) Discount curve dates                 */
//    double *discRates,       /* (I) Discount curve rates                 */
//    double  corr,            /* (I) Index correlation data               */
//    double *fwdVolInput,     /* (I) Forward vol betw rate start dates    */
//    long    optType,         /* (I) Option type                          */
//    long    payDate,         /* (I) Option payment date                  */
//    long    setlType,        /* (I) Cash or physical settlement          */
//    long    numPayoffParams, /* (I) Num payment parameters               */
//    double *payoffParams,    /* (I) Payment parameters                   */
//    char   *holidayfile,     /* (I) holiday file name                    */
//    char   *BusVolDCC,       /* (I) BUS/251F or BUS/BUS                  */
//    double *results          /* (O) Pricing results.                     */
//    )
//{
//    try
//    {
//        union{
//            struct{
//                double fwdRate;
//                double fwdRate30360;
//                double fwdAnn30360;
//                double zeroRateSwap;
//                double zeroRatePay;
//                double expiry;
//                double expiryVolvol;
//                double start;
//                double swapMat;
//                double payDelay;
//                double freq;
//            } s;
//            double a[11];
//        } rate[2];
//
//        FPAYOFF  *payFunc       = NULL;
//        double   *smile         = NULL;
//        double   *indxRates     = NULL;
//        double   *vnfmParams    = NULL;
//        long     *indxDates     = NULL;
//        long      numVNFMParams = 0;
//        char     *rateDCC       = NULL;
//
//        PAYOFF    pf;
//        FADATA    fa[2];
//        MQDATA    mq[2];
//        MQDATA    pa[2];
//        MQDATA   *ppa[2];
//        double    sigATM;
//        double    price;
//    	long      numIndxPts;
//        long      rateFreq;
//        long      startDate, matDate, expiryDate;
//        int       i;
//
//        /* CMS BMA option */
//        if (optType == Q3_BS_PERC_CALL ||
//            optType == Q3_BS_PERC_PUT  ||
//            optType == Q3_BS_PERC_RATE_CALL ||
//            optType == Q3_BS_PERC_RATE_PUT ||
//            optType == Q3_BS_JOINT_FWD)
//        {
//            /* Consistency check */
//            if (expiryDate1 != expiryDate2 ||
//                startDate1  != startDate2 ||
//                matDate1    != matDate2) 
//            {
//                throw ModelException("The percentage spread and the rate require same expiry, start & maturity dates.");  
//            }
//
//            Q3MQPSAPricer(
//                today,           /*  1 (I) Vol base date                        */
//                valueDate,       /*  2 (I) Value date                           */
//                rateFreq2,       /*  3 (I) BMA and CMS fixed leg frequency      */
//                rateDCC2,        /*  4 (I) BMA and CMS fixed Leg DCC            */
//                rateFreq1,       /*  5 (I) BMA float leg frequency              */
//                rateDCC1,        /*  6 (I) BMA float leg DCC                    */
//                numIndxPts2,     /*  7 (I) LIBOR zero curve num points          */
//                indxDates2,      /*  8 (I) LIBOR zero curve dates               */
//                indxRates2,      /*  9 (I) LIBOR zero curve rates               */
//                numIndxPts1,     /* 10 (I) basis zero curve num points          */
//                indxDates1,      /* 11 (I) basis zero curve dates               */
//                indxRates1,      /* 12 (I) basis zero curve rates               */
//                expiryDate1,     /* 13 (I) reset date                           */
//                startDate1,      /* 14 (I) start date                           */
//                matDate1,        /* 15 (I) end date                             */
//                smile2,          /* 16 (I) CMS smile parameters                 */
//                sigATM2,         /* 17 (I) CMS ATM vol                          */
//                numVNFMParams2,  /* 18 (I) 1st index num VNFM parameters        */
//                vnfmParams2,     /* 19 (I) 1st index VNFM parameters            */
//                smile1,          /* 20 (I) spread smile parameters              */
//                sigATM1,         /* 21 (I) spread ATM vol                       */
//                numDiscPts,      /* 22 (I) Discount curve num points            */
//                discDates,       /* 23 (I) Discount curve dates                 */
//                discRates,       /* 24 (I) Discount curve rates                 */
//                corr,            /* 25 (I) Index correlation data               */
//                optType,         /* 26 (I) Option type                          */
//                payDate,         /* 27 (I) Option payment date                  */
//                setlType,        /* 28 (I) Cash or physical settlement          */
//                numPayoffParams, /* 29 (I) Number of payment parameters         */
//                payoffParams,    /* 28 (I) Payment parameters                   */
//                holidayfile,     /* 29 (I) holiday file name                    */
//                BusVolDCC,       /* 30 (I) BUS/251F or BUS/BUS                  */
//                results          /* 31 (O) Pricing results.                     */ 
//                );
//            return;
//        }
//
//
//        /* calibrate PA measure for each index */
//        for (i = 1; i >= 0; i--)
//        {
//            
//            /* select index, note ordering is reversed */
//            if (i == 1)
//            {
//                sigATM          = sigATM1;
//                smile           = smile1;
//                numIndxPts      = numIndxPts1;
//                indxDates       = indxDates1;
//                indxRates       = indxRates1;
//                rateFreq        = rateFreq1;
//                rateDCC         = rateDCC1;
//                expiryDate      = expiryDate1;
//                startDate       = startDate1;
//                matDate         = matDate1;
//                numVNFMParams   = numVNFMParams1;
//                vnfmParams      = vnfmParams1;
//            }
//            else
//            {
//                sigATM          = sigATM2;
//                smile           = smile2;
//                numIndxPts      = numIndxPts2;
//                indxDates       = indxDates2;
//                indxRates       = indxRates2;
//                rateFreq        = rateFreq2;
//                rateDCC         = rateDCC2;
//                expiryDate      = expiryDate2;
//                startDate       = startDate2;
//                matDate         = matDate2;
//                numVNFMParams   = numVNFMParams2;
//                vnfmParams      = vnfmParams2;
//            }
//
//            /* Fwd dates and rates for index 1. */
//            if (Q3SwapRateCalc2(
//                today,
//                today,
//                valueDate,
//                numIndxPts,
//                indxDates,
//                indxRates,
//                numDiscPts,
//                discDates,
//                discRates,
//                rateFreq,
//                rateDCC,
//                expiryDate,
//                startDate,
//                matDate,
//                payDate,
//                holidayfile,
//                BusVolDCC,
//                (rate[i]).a) == FAILURE) throw ModelException("Q3SwapRateCalc2 failed");
//            
//            if (sigATM * sqrt((rate[i].s.expiry)) < Q3_MIN_VOL_CALIB) corr = 0.;
//
//        /* smile: initialize */
//        if (Q3SmileInit(
//                    (rate[i]).s.fwdRate,                 
//                    sigATM,
//                    0.,
//                    (rate[i]).s.expiry,
//                    (rate[i]).s.expiryVolvol,
//                    1.,
//                    smile,
//                        &(mq[i])) == FAILURE) 
//            {
//                throw ModelException("Error initialising smile.");
//            }
//            
//            /* Calibrate MQ distribution to SV option prices. */
//            if (Q3MQCalib(&(mq[i])) == FAILURE) throw ModelException("Q3MQCalib failed");
//            
//            /* Initialize FA parameters. */
//            (fa[i]).mq = &(mq[i]);
//            if (Q3FASmileInit(
//                (rate[i]).s.expiry,
//                mq[i].sigATM,
//                (rate[i]).s.start,
//                (long) (rate[i]).s.freq,
//                (rate[i]).s.swapMat,
//                (rate[i]).s.fwdRate30360,
//                (rate[i]).s.fwdAnn30360,
//                (rate[i]).s.zeroRateSwap,
//                (rate[i]).s.payDelay,
//                (rate[i]).s.zeroRatePay,
//                numVNFMParams,
//                vnfmParams,
//                setlType,
//                setlType,
//                &(fa[i])) == FAILURE) throw ModelException("Q3FASmileInit failed");
//         
//            /* initialize PA MQ */
//            if (Q3MQCopySmileFromMQ(
//                &(mq[i]),
//                &(pa[i])) == FAILURE) throw ModelException("Q3MQCopySmileFromMQ failed");
//           
//            /* set NCK parameters */
//            if (Q3DecodeNCK(
//                Q3_BV_NCK, 
//                &(pa[i])) == FAILURE) throw ModelException("Q3DecodeNCK failed");
//
//            /* compute calibration targets */
//            if (Q3MQTargetFA( 
//                &(fa[i]),
//                &(pa[i])) == FAILURE) throw ModelException("Q3MQTargetFA failed");
//
//            /* bootstrap PA multi-q measure */
//            if (Q3MQBootstrapQ(&(pa[i])) == FAILURE) throw ModelException("Q3MQBootstrapQ failed");
//
//            /* calibrate PA */
//            if (Q3MQCalib(&(pa[i])) == FAILURE) throw ModelException("Q3MQCalib failed");
//        } /* i */
//
//        /* adjust correlation for rate start date difference */
//        if (startDate1 != startDate2)
//        {
//            double vol, fwdVol, fwdTenor, swapVol, swapFwdVol;
//            double zeroVol, zeroSwapCorr, adj;
//            int    iLate;
//
//            fwdTenor = (rate[1]).s.start - (rate[0]).s.start;
//            if (fwdTenor > 0) 
//            {
//                iLate = 1;
//                vol   = sigATM2;
//            }
//            else
//            {
//                iLate    = 0;
//                vol      = sigATM1;
//                fwdTenor = fabs(fwdTenor);
//            }
//
//            if (fwdVolInput == NULL)
//            {
//                /* estimate forward vol */   
//        
//                if (vol < TINY) throw ModelException("vol < TINY");
//
//                if (Q3VNFMZero2Swap(
//                    rate[iLate].s.expiry,
//                    rate[iLate].s.expiry,
//                    rate[iLate].s.start,
//                    (long) rate[iLate].s.freq,
//                    rate[iLate].s.swapMat,
//                    rate[iLate].s.fwdRate30360,
//                    rate[iLate].s.fwdAnn30360,
//                    rate[iLate].s.swapMat,
//                    rate[iLate].s.zeroRateSwap,
//                    numVNFMParams,
//                    vnfmParams,
//                    &swapVol,
//                    &zeroVol,
//                    &zeroSwapCorr) == FAILURE) throw ModelException("Q3VNFMZero2Swap failed");
//
//                if (Q3VNFMZero2Swap(
//                    rate[iLate].s.expiry,
//                    rate[iLate].s.expiry - fwdTenor,
//                    rate[iLate].s.start,
//                    (long) rate[iLate].s.freq,
//                    rate[iLate].s.swapMat,
//                    rate[iLate].s.fwdRate30360,
//                    rate[iLate].s.fwdAnn30360,
//                    rate[iLate].s.swapMat,
//                    rate[iLate].s.zeroRateSwap,
//                    numVNFMParams,
//                    vnfmParams,
//                    &swapFwdVol,
//                    &zeroVol,
//                    &zeroSwapCorr) == FAILURE) throw ModelException("Q3VNFMZero2Swap failed");
//
//                fwdVol = vol * swapFwdVol / swapVol;
//            }
//            else
//            {
//                fwdVol = fwdVolInput[0];
//            }
//
//            /* correlation adjustment */
//            if (rate[iLate].s.start < TINY) throw ModelException("rate[iLate].s.start < TINY");
//            adj = (fwdVol*fwdVol*fwdTenor)/(vol*vol*rate[iLate].s.start);
//            if (adj > 1.) throw ModelException("adj > 1.");
//            corr *= sqrt(1. - adj);
//        }
//        
//        /* populate payoff structure, select payoff function      */
//        /* note: ordering of ppa elements may be changed by call! */
//        ppa[0] = &(pa[0]);
//        ppa[1] = &(pa[1]);
//        Q3MQBivarInit(
//            optType,
//            numPayoffParams,
//            payoffParams,
//            &corr,
//            &(ppa[0]),
//            &pf,
//            &payFunc);
//        /* integrate payoff over density. */
//        if (Q3MQGridPricer(
//            &pf,
//            payFunc,
//            ppa[0],
//            &price) == FAILURE) throw ModelException("Q3MQGridPricer failed");
//
//        results[0] = price;
//        results[1] = pa[1].fwdRate; /* index 1 */
//        results[2] = pa[0].fwdRate; /* index 2 */
//    }
//    catch (exception& e)
//    {
//        throw ModelException(e, __FUNCTION__);
//    }
//} /* Q3MQBivarPricer */


/*****************************************************************************/
/********************************** QuasiMQ **********************************/
/*****************************************************************************/

class Q3VariateData
{
public:

    Q3VariateData();
    Q3VariateData(const DateTime& volBaseDate, const IndexSpecIR& rateSpec, const DateTime& expiry, bool legacy);
    
    // accessors for use to pass into Q3 pricing functions
    long numCurvePoints() const;
    long* curveDates();
    double* curveRates();
    const DateTime& expiryDate() const;
    const DateTime& startDate() const;
    const DateTime& maturityDate() const;
    long freq() const;
    char* dcc();
    double* smile();
    double sigATM() const;

    // hack to allow population of smileData
    DoubleArray& smileRef();

private:
    std::vector<long>  curveDates_;
    DoubleArray        curveRates_;
    DateTime           expiryDate_;
    DateTime           startDate_;
    DateTime           maturityDate_;
    long               freq_;
    std::string        dcc_;
    DoubleArray        smile_;
    double             atmVol_;
};

struct GetLongDate : public std::unary_function<DateTime, long>
{
    long operator()(const DateTime& dt)
    {
        return static_cast<long>(dt.getDate());
    }
};

void QuasiMQ::univarPricer(const IndexSpecIR&            rateSpec,        /*  1 (I) rate definition           */
                           const DateTime&               resetDate,       /*  2 (I) rate effective date       */
                           const DateTime&               payDate,         /*  3 (I) option payment date       */
                           double                        strike,          /*  4 (I) strike (0. for Q3_ADJ)    */
                           const DayCountConventionSP&   busDcc,          /*  5 (i) between today to reset date */
                           long                          optType,         /*  6 (I) option type               */
                           const DoubleArray&            payoffParams,    /*  7 (I) payoff coeffs             */
                           const InstrumentSettlementSP& setlType,        /*  8 (I) cash or phys settle       */
                           double&                       price,           /*  9 (O) fwd price and AA par rate */
                           double&                       fwdRate,         /* 10 (O) fwd price and AA par rate */
                           VariateDebugData&             debug) const
{
    try
    {
        DateTime volBaseDate = today;
        DateTime smileBaseDate = getValueDate();

		Q3VariateData rateData(volBaseDate, rateSpec, resetDate, legacyPricingMode);
		populateSmileData(volBaseDate, resetDate, rateSpec.tenor.get(), rateData.smileRef());

        // ??? temporary check - remove for quanto as spot days may be different
        if (rateSpec.factor->getSpotDate() != spotDate)
        {
            throw ModelException("At this stage, the spot Date (" + rateSpec.factor->getSpotDate().toString() +
                                 ") for index curve " + rateSpec.factor->getName() + " must be the same "
                                 " as the spot Date (" + spotDate.toString() + ") for the discount curve " +
                                 domesticYC->getName());
        }
		
        std::vector<long> domLDates;
        std::transform(domDates.begin(), domDates.end(), std::back_inserter(domLDates), GetLongDate());

        double bbATM = 0.;  // 0 = yield vol, 1 = bp vol

        // Call Pricer
        std::vector<double> outputs(2);
        Q3MQQuasiPricer(today.getDate(),  // volBaseDate.getDate(),
                        smileBaseDate.getDate(),
                        spotDate.getDate(),   // domesticYC->getSpotDate().getDate(),
                        rateData.numCurvePoints(),
                        rateData.curveDates(),
                        rateData.curveRates(),
                        domLDates.size(),
                        static_cast<long*>(&domLDates[0]),
                        (double*)&domRates[0],
                        vnfm.size(),
                        (double*)&vnfm[0],
                        rateData.freq(),
                        rateData.dcc(),//rateDcc,
                        rateData.expiryDate().getDate(),
                        rateData.startDate().getDate(), // to be changed
                        rateData.maturityDate().getDate(),
                        rateData.sigATM(),
                        bbATM,
                        rateData.smile(),
                        optType,
                        &strike,
                        payoffParams.size(),
                        (double*)&payoffParams[0],
                        payDate.getDate(),
                        (setlType->isPhysical() ? Q3_PHYS_SETL : Q3_CASH_SETL),
                        NULL,
                        NULL,
                        &outputs[0]);

        price = outputs[0];
        fwdRate = outputs[1];

        // debug...
        debug.expiryDates->push_back(rateData.expiryDate());
        debug.startDates->push_back(rateData.startDate());
        debug.endDates->push_back(rateData.maturityDate());
        debug.payDates->push_back(payDate);
        debug.atmVols->push_back(rateData.sigATM());
        debug.prices->push_back(price);
        debug.fwdRates->push_back(fwdRate);
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}

void QuasiMQ::bivarPricer(const IndexSpecIR&            rate1Spec,       /*  1 (I) rate1 definition           */
                          const DateTime&               rate1Expiry,     /*  2 (I) rate1 effective date       */
                          const IndexSpecIR&            rate2Spec,       /*  3 (I) rate2 definition           */
                          const DateTime&               rate2Expiry,     /*  4 (I) rate2 effective date       */
                          const DateTime&               payDate,         /*  5 (I) option payment date       */
                          const DayCountConventionSP&   busDcc,          /*  7 (i) between today to reset date */
                          long                          optType,         /*  8 (I) option type               */
                          const DoubleArray&            payoffParams,    /*  9 (I) payoff coeffs             */
                          const InstrumentSettlementSP& setlType,        /* 10 (I) cash or phys settle       */
                          double&                       price,           /* 11 (O) fwd price and AA par rate */
                          double&                       fwdRate,         /* 12 (O) fwd price and AA par rate */
                          VariateDebugData&             debug1,
                          VariateDebugData&             debug2) const
{
    try
    {
        DateTime volBaseDate = today;
        DateTime smileBaseDate = getValueDate();

        // index1
        Q3VariateData rate1Data(volBaseDate, rate1Spec, rate1Expiry, legacyPricingMode);
        populateSmileData(volBaseDate, rate1Expiry, rate1Spec.tenor.get(), rate1Data.smileRef());

        // index2
        Q3VariateData rate2Data(volBaseDate, rate2Spec, rate2Expiry, legacyPricingMode);
        populateSmileData(volBaseDate, rate2Expiry, rate2Spec.tenor.get(), rate2Data.smileRef());

        // ??? temporary check - remove for quanto as spot days may be different
        if (rate1Spec.factor->getSpotDate() != spotDate || rate2Spec.factor->getSpotDate() != spotDate)
        {
            throw ModelException("SpotDate (" + rate1Spec.factor->getSpotDate().toString() + ") for curve '" + rate1Spec.factor->getName() +
                                 "' and SpotDate (" + rate2Spec.factor->getSpotDate().toString() + ") for curve '" + rate2Spec.factor->getName() +
                                 "'must be the same as the SpotDate (" + spotDate.toString() + ") for the discount curve " + domesticYC->getName());
        }

        // domDates
        std::vector<long> domLDates;
        std::transform(domDates.begin(), domDates.end(), std::back_inserter(domLDates), GetLongDate());

        double  corr   = 0; // correlation
        double fwdVol  = 0;
       
        // Call Pricer
        std::vector<double> outputs(3);
        Q3MQBivarPricer(today.getDate(),  // volBaseDate.getDate(),
                        spotDate.getDate(),
                        rate1Data.numCurvePoints(),
                        rate1Data.curveDates(),
                        rate1Data.curveRates(),
                        rate1Data.freq(),
                        rate1Data.dcc(),
                        rate1Data.expiryDate().getDate(),
                        rate1Data.startDate().getDate(),
                        rate1Data.maturityDate().getDate(),
                        rate1Data.smile(),
                        rate1Data.sigATM(),
                        static_cast<long>(vnfm.size()),
                        const_cast<double*>(&vnfm[0]),
                        rate2Data.numCurvePoints(),
                        rate2Data.curveDates(),
                        rate2Data.curveRates(),
                        rate2Data.freq(),
                        rate2Data.dcc(),
                        rate2Data.expiryDate().getDate(),
                        rate2Data.startDate().getDate(),
                        rate2Data.maturityDate().getDate(),
                        rate2Data.smile(),
                        rate2Data.sigATM(),
                        static_cast<long>(vnfm.size()),
                        const_cast<double*>(&vnfm[0]),
                        static_cast<long>(domLDates.size()),
                        const_cast<long*>(&domLDates[0]),
                        const_cast<double*>(&domRates[0]),
                        corr,
                        &fwdVol,
                        optType,
                        payDate.getDate(),
                        (setlType->isPhysical() ? Q3_PHYS_SETL : Q3_CASH_SETL),
                        static_cast<long>(payoffParams.size()),
                        const_cast<double*>(&payoffParams[0]),
                        NULL,
                        NULL,
                        static_cast<double*>(&outputs[0]));

        price = outputs[0];
        fwdRate = outputs[1];

        // debug...
        debug1.expiryDates->push_back(rate1Data.expiryDate());
        debug1.startDates->push_back(rate1Data.startDate());
        debug1.endDates->push_back(rate1Data.maturityDate());
        debug1.payDates->push_back(payDate);
        debug1.atmVols->push_back(rate1Data.sigATM());
        debug1.prices->push_back(outputs[0]);
        debug1.fwdRates->push_back(outputs[1]);
        
        debug2.expiryDates->push_back(rate2Data.expiryDate());
        debug2.startDates->push_back(rate2Data.startDate());
        debug2.endDates->push_back(rate2Data.maturityDate());
        debug2.payDates->push_back(payDate);
        debug2.atmVols->push_back(rate2Data.sigATM());
        debug2.prices->push_back(outputs[0]);
        debug2.fwdRates->push_back(outputs[2]);
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}

void QuasiMQ::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments)
{
    static const string method = "QuasiMQ::getMarket";
    try
    {
        today = market->GetReferenceDate();
        if (!smileTable.isEmpty())
            smileTable.getData(this, market);
        if (!modelTable.isEmpty())
            modelTable.getData(this, market);
        // Phasing out of irCalib - do not use irCalib if using smileTable or modelTable.
        if (!irCalib.isEmpty() && smileTable.isEmpty() && modelTable.isEmpty())
            irCalib.getData(this, market);  // populate market wrapper

        // recurse all the potential components in the instrument and look for the type
        // IMarketFactor::TYPE, which represent marketData type fields (yield curves, 
        // index Specs etc) and store them for analysis later on
        class RetrieveFactors : public ObjectIteration::IAction
        {
            const MarketData* market;
            IMarketFactorArray& factors;

        public:
            RetrieveFactors(const MarketData* market, IMarketFactorArray& factors ) :
                market( market ), factors( factors ) { factors.clear(); }

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
            {
                IMarketFactor* factor = dynamic_cast<IMarketFactor*>(state.getObject().get());
                string name = factor->getName();
                string type = factor->getClass()->getName();

                int i;
                // if already encountered, don't store again
                for (i = 0; i < factors.size(); ++i )
                {
                    if(name == factors[i]->getName() && 
                       type == factors[i]->getClass()->getName())
                       break;
                }
                if (i >= factors.size())
                {
                    factors.push_back(IMarketFactorSP::attachToRef(factor));
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveFactors retrieveFactors(market, factors);
        ObjectIteration iteration(IMarketFactor::TYPE);
        iteration.recurse(retrieveFactors, instruments);

        // recurse all the potential components in the instrument and look for the type
        // IndexSpec::TYPE, and store them for later 
        class RetrieveIndexSpecs : public ObjectIteration::IAction
        {
            IndexSpecArray& indexSpecs;

        public:
            RetrieveIndexSpecs(IndexSpecArray& indexSpecs) : indexSpecs(indexSpecs) 
                {indexSpecs.clear();}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
            {
                int i;
                IndexSpec * indexSpec = dynamic_cast<IndexSpec*>(state.getObject().get());

                string name = indexSpec->getName();
                string type = indexSpec->getClass()->getName();

                for (i = 0; i < indexSpecs.size(); ++i )
                {
                    if (name == indexSpecs[i]->getName() && type == indexSpecs[i]->getClass()->getName())
                        break;
                }
                if (i >= indexSpecs.size())
                {
                    indexSpecs.push_back(IndexSpecSP::attachToRef(indexSpec));
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveIndexSpecs retrieveIS(indexSpecs);
        ObjectIteration retrieveISIter(IndexSpec::TYPE);
        retrieveISIter.recurse(retrieveIS, instruments);

        // get populate domesticYC  ??? do more elegant way
        domesticYC.setName((*instruments)[0]->discountYieldCurveName());
        domesticYC.getData(this, market);

        if (legacyPricingMode) {
        IrConverter::q3GetZerosData(domDates, domRates, *domesticYC.get());
            spotDate = domesticYC->getSpotDate();
        }
        else
            IrConverter::q3GetandExtendZerosData(domDates, domRates, spotDate, 
                                                 *domesticYC.get(), today);

        // retrieve 2Q smile parameters: qLeft, qRight, fwdShift
        if (smileTable.isEmpty() && modelTable.isEmpty())
        {
            // Using IRCALIB
            // retrieve 1D vnfm model parameters: oneFactorMR, oneFactorWeight
            IRCalib::ModelRequest modelRequest(modelSet);
            CVolProcessedSP volProcessed2(irCalib->getProcessedVol(&modelRequest,0));
            IRCalib::VolProcessed* volData = dynamic_cast<IRCalib::VolProcessed*>(volProcessed2.get());
            if (!volData)
                throw ModelException(method, "volProcessed should be of type IRCalib::VolProcessed");

            const DoubleArray& modelParams = volData->getParams();
            const StringArray& paramLabel = volData->getParamLabel();

            if (paramLabel.empty())
                throw ModelException(method, "Must supply optional paramLabel array field for "
                    "IRCalib::Model[1-3]FL to calculate the VNFM");

            int nbFactors = 1;
            if (nbFactors == 1)
            {
                vnfm.resize(2);
                vnfm[0] = modelParams[0];
                vnfm[1] = modelParams[1];
            }
            else 
                throw ModelException("Currently only one factor IR model parameters supported");
        }
        else
        {
            // Initialise if SmileMQ Object
            IRExoticParamTable& smileTableObj = *(smileTable.get());
            MarketObjectSP ptrSmile = smileTableObj.getMarketObject(smileSet);
            MarketObject*  pSmile = ptrSmile.get();
            if (!pSmile)
                throw ModelException(method, " - Smile object " + smileSet + 
                    " not found in Smile Table " + smileTableObj.getName() + "!");

            IRSmileMQ* pIRSmileMQ = dynamic_cast<IRSmileMQ*>(pSmile);
            if (pIRSmileMQ)
                pIRSmileMQ->initialiseData(market, this);

            // Lookup Model
            IRExoticParamTable& modelTableObj = *(modelTable.get());
            MarketObjectSP ptrModel = modelTableObj.getMarketObject(modelSet);
            MarketObject*  pModel = ptrModel.get();
            if (!pModel)
                throw ModelException(method, " - Model object " + modelSet + 
                    " not found in Model Table " + modelTableObj.getName() + "!");
            const IRModelVNFM* pIRModelVNFM = dynamic_cast<const IRModelVNFM*>(pModel);
            if (!pIRModelVNFM)
                throw ModelException(method, " - Object " + modelSet + " is not of type IRModelVNFM!");
            
            // Populate model parameters
            int nbFactors = pIRModelVNFM->numFactors();
            const DoubleMatrix& mMR = pIRModelVNFM->meanReversion();
            const DoubleMatrix& mWeight = pIRModelVNFM->weight();
            const DoubleMatrix& mCorr = pIRModelVNFM->correlation();

            // Setup vnfm array (multi-factor). Size = n MR's + n Weights + n*(n-1)/2 Correlations = n(n+3)/2
            int i;
            int j = 0;
            int nbCorrelations = nbFactors * (nbFactors - 1) / 2;
            vnfm.resize(nbFactors + nbFactors + nbCorrelations);
            for (i = 0; i < nbFactors; ++i)
                vnfm[j++] = mMR[i][0];
            for (i = 0; i < nbFactors; ++i)
                vnfm[j++] = mWeight[i][0];
            for (i = 0; i < nbCorrelations; ++i)
                vnfm[j++] = mCorr[i][0];
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

double QuasiMQ::getZero(DateTime useDate, DateTime matDate, string curveName)
{
    try
    {
        if (curveName != domesticYC->getName())
            throw ModelException("Can only calculate discount factor on the curve "
                "name specified as the model discount curve (" + domesticYC->getName() +
                ").  curveName supplied = " + curveName);

        if (useDate > matDate)
            throw ModelException("useDate " + useDate.toString() + " must be <= matDate " +
                matDate.toString());
        // ??? for now, becuase esl is in a state of flux, construct old style yyyymmdd zerodates
        // and rates and use esl to interp zero rates from curve in order to calculate df
        std::vector<IRDate> localDomDates(domDates.size());

        for (int i = 0; i < domDates.size(); ++i)
        {
            localDomDates[i] = domDates[i].toIrDate();
        }

        // ??? for now use the old rates function to    calculate discount factor
        return ZeroPrice(matDate.toIrDate(), useDate.toIrDate(), localDomDates.size(),
                                  (IRDate const*)&localDomDates[0], (double const*)&domRates[0]);
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}


void QuasiMQ::populateSmileData(const DateTime& baseDate,
                                const DateTime& expiryDate,
                                const Expiry*   pMaturityTenor,
                                DoubleArray&    smileData) const // output
{
    static const string method = "QuasiMQ::populateSmileData";
    if (smileTable.isEmpty() && modelTable.isEmpty())
    {
        // Use IRCalib
        IRCalib::SmileRequest smileRequest(smileSet);
        CVolProcessedSP volProcessed(irCalib->getProcessedVol(&smileRequest,0));
        IRCalib::VolProcessed *volData = dynamic_cast<IRCalib::VolProcessed*>(volProcessed.get());
        if (!volData)
            throw ModelException("volProcessed should be of type IRCalib::VolProcessed");
        const DoubleArray& smileParams = volData->getParams();
        if (smileParams.size() != 4)
            throw ModelException("Number of ir vol smile params supplied = " +
                                Format::toString(smileParams.size()) +
                                " must supply 4 (QLeft, QRight, FwdShift, nbCETIters)");

        smileData[0] = static_cast<double>(SMILE_2Q);
        smileData[1] = smileParams[0];
        smileData[2] = smileParams[1];
        smileData[3] = smileParams[2];
        return;
    }

    // Lookup Smile
    MarketObjectConstSP ptrSmile = smileTable->getMarketObject(smileSet);
    if (!ptrSmile.get())
        throw ModelException(method, " - Smile object " + smileSet + " not found in Smile Table " + smileTable->getName() + "!");

    // Smile2Q
    const IRSmile2Q* pIRSmile2Q = dynamic_cast<const IRSmile2Q*>(ptrSmile.get());
    if (pIRSmile2Q)
    {
        // Populate smile parameters
        smileData[0] = static_cast<double>(SMILE_2Q);
        smileData[1] = pIRSmile2Q->qLeft();
        smileData[2] = pIRSmile2Q->qRight();
        smileData[3] = pIRSmile2Q->fwdShift();
    }
    
    // SmileQuasi2Q
    const IRSmileQuasi2Q* pIRSmileQuasi2Q = dynamic_cast<const IRSmileQuasi2Q*>(ptrSmile.get());
    if (pIRSmileQuasi2Q)
    {
        // Populate smile parameters
        smileData[0] = static_cast<double>(SMILE_2Q);
        smileData[1] = pIRSmileQuasi2Q->getQLeft(baseDate, expiryDate, pMaturityTenor);
        smileData[2] = pIRSmileQuasi2Q->getQRight(baseDate, expiryDate, pMaturityTenor);
        smileData[3] = pIRSmileQuasi2Q->getFwdShift(baseDate, expiryDate, pMaturityTenor);
    }
    
    // SmileMQ
    const IRSmileMQ* pIRSmileMQ = dynamic_cast<const IRSmileMQ*>(ptrSmile.get());
    if (pIRSmileMQ)
    {
        // Populate smile parameters
        smileData[0] = static_cast<double>(SMILE_MQ);
        smileData[1]  = pIRSmileMQ->getSkew(baseDate, expiryDate, pMaturityTenor);
        smileData[2]  = pIRSmileMQ->getVolOfVol(baseDate, expiryDate, pMaturityTenor);
        smileData[3]  = pIRSmileMQ->getBBRP(baseDate, expiryDate, pMaturityTenor);
        smileData[4]  = pIRSmileMQ->getBBVP(baseDate, expiryDate, pMaturityTenor);
        smileData[5]  = pIRSmileMQ->getDeltaLeft(baseDate, expiryDate, pMaturityTenor);
        smileData[6]  = pIRSmileMQ->getTauLeft(baseDate, expiryDate, pMaturityTenor);
        smileData[7]  = pIRSmileMQ->getDeltaRight(baseDate, expiryDate, pMaturityTenor);
        smileData[8]  = pIRSmileMQ->getTauRight(baseDate, expiryDate, pMaturityTenor);
        smileData[9]  = pIRSmileMQ->getNormalCutoff();
        smileData[10] = pIRSmileMQ->getNCK();
        
        // Used by Q3SmileInit
        smileData[11] = 0.10; //-999.0;
        smileData[12] = 0.10; //-999.0;
    }
    
    // Invalid smile type    
    if (!pIRSmile2Q && !pIRSmileQuasi2Q && !pIRSmileMQ)
        throw ModelException(method, " - Object " + smileSet + " is not of type IRSmile2Q or IRSmileMQ!");
}

void QuasiMQ::Price(CInstrument* instrument,
    CControl*    control,
    CResults*    results)
{
    try
    {
        QuasiMQ::IIntoProduct* creator = dynamic_cast<QuasiMQ::IIntoProduct*>(instrument);
        if (!creator)
        {
            throw ModelException("Instrument of type " + instrument->getClass()->getName() +
                " does not support QuasiMQ::IIntoProduct");
        }

        // create product
        QuasiMQ::ProductSP product = creator->createProduct(this);

        // ??? need to think more about this
        // check what model mode we may require for pricing
        // check all instrument market factors and indexSpecs, and then can do
        // if 1 unique indexSpecIR, modelMode is univariate
        // if 2 unique indexSpecIR, modelMode is bivariate
        // if currency of indexSpecIR != domestic currency (ie. quanto) return error
 /*       int i;
        
        for (i = 0; i < factors.size(); i++) {

            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());

            if (yc) {
                string ycCcy = yc->getCcy();
                if (ycCcy != domesticYC->getCcy())
                    throw ModelException("this model does not support quanto interest "
                        "rate curves (currency = " + ycCcy + ") where the domestic/ "
                        "pricing currency is " + domesticYC->getCcy() + " for yield curve " +
                        yc->getName());
            }
            else
                throw ModelException("Only yieldCurve market factors supported - invalid type "
                    "supplied = " + mf->getClass()->getName() + ", name = " +
                    mf->getName());
        }

        size_t j;
        vector<string> indexName;
        for (i = 0; i < indexSpecs.size(); i++) {

            IndexSpecSP is = indexSpecs[i];
            const IndexSpecIR* ir = dynamic_cast<const IndexSpecIR*>(is.get());
            const IndexSpecFX* fx = dynamic_cast<const IndexSpecFX*>(is.get());

            if (ir) {
                // check number of separate index specs - these fields define
                // a unique indexSpecIR
                string tenor = ir->tenor->toString();
                string freq = ir->frequency->toString();
                string dcc = ir->dcc->toString();
                string name = tenor + freq + dcc + ir->factor->getName();
                // check if it exists before
                for (j = 0; j < indexName.size(); j++) {
                    if (name == indexName[j])
                        break;
                }
                if (j >= indexName.size())
                    indexName.push_back(name);
            }
            else if (fx) {
                // check if implicit indexSpecFX ie.  dom-dom and thus guaranteed a value = 1.0
                if (fx->sameCurve == false)
                    throw ModelException("IndexSpecFX " + fx->getName() +
                        "must be a single currency currency FX (ie. rate = 1.0, "
                        "same foreign and domestic curve) for this model");
            }
            else
                throw ModelException("Only indexSpecIR type indexSpecs are supported - "
                    "type supplied = " + is->getClass()->getName());
        }

        int nbIRIndex = indexName.size();
        if (nbIRIndex != 1)
            throw ModelException("Currently only support univariate pricing (ie. one "
                "unique interest rate spec defined in the product, found " +
                Format::toString(nbIRIndex));
*/
        EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::LINEAR ?
            ESL_INTERP_LINEAR : ESL_INTERP_FLATFWD);
        product->price(control, results);
    }
    catch (exception& e)
    {
        throw ModelException(e, __FUNCTION__);
    }
}


QuasiMQ::QuasiMQ(const CClassConstSP &type)
    : CModel(TYPE), zeroInterpStyle(ZeroInterpStyle::FLAT_FWD), legacyPricingMode(false)
{}


void QuasiMQ::registerErrorCallback()
{
    if (Q3ErrCallbackSet(errCallBack) != SUCCESS)
    {
        throw ModelException(__FUNCTION__, "Unable to register Q3 error handling callback function");
    }
}


string QuasiMQ::popErrMsg()
{
    string tmp = errMsg;
    errMsg.clear();
    size_t s;
    while ( (s = tmp.size()) != 0 && tmp[s - 1] == '\n' )
        tmp.resize(s-1);
    return tmp;
}


ModelException QuasiMQ::makeException()
{
    return ModelException(popErrMsg());
}


int QuasiMQ::errCallBack(const char *msg)
{
        errMsg += msg;
        return SUCCESS;
}


void QuasiMQ::load(CClassSP& clazz)
{
    REGISTER(QuasiMQ, clazz);
    clazz->setPublic(); // make externally visible
    clazz->setDescription("Q3 based 2Q smile pricing model");
    SUPERCLASS(CModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(modelMode, "");
    FIELD(irCalib, "Vol/model parameters market wrapper name")
    FIELD_MAKE_OPTIONAL(irCalib);
    FIELD(smileSet, "Name of 2Q smile parameters set to use from the market cache");
    FIELD(modelSet, "Name of 2Q smile parameters set to use from the market cache");
    FIELD(smileTable,"A MarketTable object containing collection of IR smile objects, "
          "accessed using the smile key given above.");
    FIELD_MAKE_OPTIONAL(smileTable);
    FIELD(modelTable,"A MarketTable object containing collection of IR model objects, "
          "accessed using the model key given above.");
    FIELD_MAKE_OPTIONAL(modelTable);
    FIELD(zeroInterpStyle, "Zero curve interpolation style.  Defaults to "
          "FLAT_FWD if not supplied");
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(legacyPricingMode, "If true set valueDate = IR spot date (as defined by curves) and "
          "do not adjust curves to start from today.  False by default, which means curves are "
          "extended to today and price is discounted to today.");
    /*In both cases time to expiry "
          "(ie. volatility) starts from today"); */
    FIELD_MAKE_OPTIONAL(legacyPricingMode);

    // transient fields
    FIELD(today,"");            FIELD_MAKE_TRANSIENT(today);
    FIELD(spotDate, "");        FIELD_MAKE_TRANSIENT(spotDate);
    FIELD(domesticYC,"");       FIELD_MAKE_TRANSIENT(domesticYC);
    FIELD(vnfm,"");             FIELD_MAKE_TRANSIENT(vnfm);
    FIELD(domDates,"");         FIELD_MAKE_TRANSIENT(domDates);
    FIELD(domRates,"");         FIELD_MAKE_TRANSIENT(domRates);

    registerErrorCallback();
}

CClassConstSP const QuasiMQ::TYPE = CClass::registerClassLoadMethod(
    "QuasiMQ", typeid(QuasiMQ), load);

string QuasiMQ::errMsg;

START_PUBLIC_ENUM_DEFINITION(QuasiMQ::ModelMode::Enum, "");
ENUM_VALUE_AND_NAME(QuasiMQ::ModelMode::Q2Q, "Q2Q", "");
ENUM_VALUE_AND_NAME(QuasiMQ::ModelMode::QMULTIQ, "QMULTIQ", "");
END_ENUM_DEFINITION(QuasiMQ::ModelMode::Enum);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// START Q3VariateData Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Q3VariateData::Q3VariateData()
{}

Q3VariateData::Q3VariateData(const DateTime& volBaseDate,
                             const IndexSpecIR& rateSpec,
                             const DateTime& expiry,
                             bool legacy)
    : expiryDate_(expiry),
      startDate_(expiry),
      maturityDate_(rateSpec.tenor->toDate(expiry)),
      freq_(rateSpec.frequency->equals(rateSpec.tenor.get()) ? 0 : rateSpec.frequency->annualFrequency()),
      dcc_(rateSpec.dcc->toString()),
      smile_(13)
{
    // Aladdin has: expiryDate = rateStart - numDaysToSpot
    // ... so in the interests of "replicating legacy models" - this is about the
    // the only way I can think of, as you can't do "expiryDate_ -= *offSet;"
    MaturityPeriod* offset = dynamic_cast<MaturityPeriod*>(rateSpec.fwdRateOffset.get());
    if (offset)
    {
        int count;
        std::string step;
        offset->decompose(count, step);
        expiryDate_ = MaturityPeriod::toDate(-count, step, expiryDate_);
    }

    if (!rateSpec.factor.get())
    {
        throw ModelException("IndexSpec Rate factor for '" + rateSpec.getName() + "' not populated");
    }

    DateTimeArray dtArray;
    if (legacy)
    {
        IrConverter::q3GetZerosData(dtArray, curveRates_, *rateSpec.factor.get());
    }
    else
    {
        DateTime spotDate;
        IrConverter::q3GetandExtendZerosData(dtArray, curveRates_, spotDate, *rateSpec.factor.get(), volBaseDate);
    }

    std::transform(dtArray.begin(), dtArray.end(), std::back_inserter(curveDates_), GetLongDate());
    
    // ATM Vol
    DateTime tenorDate   = rateSpec.tenor->toDate(volBaseDate);
    DateTime tenorDate1Y = MaturityPeriod(1).toDate(volBaseDate);
    const IRVolBase* pVol = IrConverter::getVol(rateSpec.factor.get(), tenorDate < tenorDate1Y);
    
    SwapMaturityVolRequest volReqSwapMaturity(rateSpec.tenor.get());
    CVolProcessedBSSP volProcessedBS = CVolProcessedBSSP::dynamicCast(CVolProcessedSP(pVol->getProcessedVol(&volReqSwapMaturity, 0)));
    
    DateTimeArray expiryDates;
    expiryDates.push_back(volBaseDate);
    expiryDates.push_back(expiry);
    
    DoubleArray vols(expiryDates.size() - 1);
    
    volProcessedBS->CalcVol(expiryDates, CVolProcessedBS::forward, vols);
    
    atmVol_ = vols[0];
}
    
long
Q3VariateData::numCurvePoints() const
{
    return static_cast<long>(curveDates_.size());
}

long*
Q3VariateData::curveDates()
{
    return static_cast<long*>(&curveDates_[0]);
}

double*
Q3VariateData::curveRates()
{
    return static_cast<double*>(&curveRates_[0]);
}

const DateTime&
Q3VariateData::expiryDate() const
{
    return expiryDate_;
}

const DateTime&
Q3VariateData::startDate() const
{
    return startDate_;
}

const DateTime&
Q3VariateData::maturityDate() const
{
    return maturityDate_;
}

long
Q3VariateData::freq() const
{
    return freq_;
}

char*
Q3VariateData::dcc()
{
    return const_cast<char*>(dcc_.c_str());
}

double*
Q3VariateData::smile()
{
    return static_cast<double*>(&smile_[0]);
}

DoubleArray&
Q3VariateData::smileRef()
{
    return smile_;
}

double
Q3VariateData::sigATM() const
{
    return atmVol_;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// END Q3VariateData Implementation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// type loading
bool QuasiMQLoad()
{
    return (QuasiMQ::TYPE != 0);
}


DRLIB_END_NAMESPACE
