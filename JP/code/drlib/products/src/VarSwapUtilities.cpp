//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GammaSwap.cpp
//
//   Description : Gamma Swap Instrument
//
//   Author      : Manos Venardos
//
//   Date        : 1 February 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/VarSwapUtilities.hpp"
#include "edginc/Format.hpp"
#include "edginc/Black.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/FixingType.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/ProtAsset.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/VegaMatrixLite.hpp"
#include "edginc/VegaParallel.hpp"


DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////////////////////////////////////


VanillaInfo::VanillaInfo(): CObject(TYPE) {}


VanillaInfo::VanillaInfo(const DateTime& maturity,
                         double nbContracts,
                         bool isCall,
                         double fwd,
                         double strike,
                         double vol,
                         double yearFrac):
CObject(TYPE), maturity(maturity), nbContracts(nbContracts), isCall(isCall),
fwd(fwd), strike(strike), vol(vol), yearFrac(yearFrac) {}


void VanillaInfo::scaleNbContracts(double scalingFactor) {
    nbContracts *= scalingFactor;
}


double VanillaInfo::price() const {
    double var = Maths::square(vol) * yearFrac;
    double thisPrice = Black::price(isCall, fwd, strike, 1.0, var) * nbContracts;
    return thisPrice;
}


void VanillaInfo::bracketArray(DoubleArray& xArray, double x, int& lowerBound, int& upperBound) {
    // Use locate method to bracket
    unsigned long index;
    int size = xArray.size();
    locate(&xArray[0] - 1, size, x, &index);        // index is in [0, size]

    // Bracket indices
    if(index == 0) {
        lowerBound = 0;
        upperBound = 0;
    } else if(index == size) {
        lowerBound = size-1;
        upperBound = size-1;
    } else {
        lowerBound = index-1;
        upperBound = index;
    }
}


double VanillaInfo::vega(double shift) const {
    double vega = 0.0;
    if(!Maths::isZero(shift)) {
        // Compute Vega as a Finite Difference
        double var0 = Maths::square(vol) * yearFrac;
        double var1 = Maths::square(vol + shift) * yearFrac;

        double price0 = Black::price(isCall, fwd, strike, 1.0, var0);
        double price1 = Black::price(isCall, fwd, strike, 1.0, var1);
        
        double deriv = (price1 - price0) / shift;
        vega += nbContracts * deriv;
    }
    return vega;
}


double VanillaInfo::volInterpMaturity(double t1, double t2, 
                                      double vol_t1, double vol_t2,
                                      double t) {
    if(Maths::equals(t1, t2)) {
        // Assume the user is smart enough to have
        // t1 == t2 && vol_t1 == vol_t2 && t == t1
        return vol_t1;
    }

    // Linear Var between maturities
    double w_low = (t2 - t) / (t2 - t1);
    double var = w_low * Maths::square(vol_t1) * t1 + (1.0 - w_low) * Maths::square(vol_t2) * t2;
    double vol = sqrt(var / t);
    return vol;
}


double VanillaInfo::volInterpStrike(double k1, double k2, 
                                    double vol_k1, double vol_k2,
                                    double k) {
    if(Maths::equals(k1, k2)) {
        // Assume the user is smart enough to have
        // k1 == k2 && vol_k1 == vol_k2 && k == k1
        return vol_k1;
    }
    
    // Linear vol between strikes
    double w_low = (k2 - k) / (k2 - k1);
    double vol = w_low * vol_k1 + (1.0 - w_low) * vol_k2;
    return vol;
}


double VanillaInfo::volInterp(double t1, double t2, double k1, double k2,
                              double vol_t1_k1, double vol_t1_k2, double vol_t2_k1, double vol_t2_k2,
                              double t, double k) {
    
    // Interpolate vol at benchmark strikes between benchmark dates
    double vol_t_k1 = volInterpMaturity(t1, t2, vol_t1_k1, vol_t2_k1, t);
    double vol_t_k2 = volInterpMaturity(t1, t2, vol_t1_k2, vol_t2_k2, t);
                                
    // Interpolate vol between strikes
    double vol = volInterpStrike(k1, k2, vol_t_k1, vol_t_k2, k);
    return vol;
}


double VanillaInfo::price(VanillaContractsRecorderSP recorder) {
    // Price vanilla portfolio
    VanillaInfoArraySP options = recorder->getResults();
    double price = 0.0;
    for(int i = 0; i < options->size(); i++) {
        const VanillaInfo& option = *(*options)[i];
        price += option.price();
    }
    return price;
}


double VanillaInfo::vegaParallel(VegaParallelSP sens, VanillaContractsRecorderSP recorder) {
    // Compute Vega Parallel of vanilla portfolio
    double vega = 0.0;
    double shiftSize = sens->getShiftSize();
    if(!Maths::isZero(shiftSize)) {
        VanillaInfoArraySP options = recorder->getResults();
        for(int i = 0; i < options->size(); i++) {
            const VanillaInfo& option = *(*options)[i];
            vega += option.vega(shiftSize);
        }
    }
    return vega;
}


#if 0
void VanillaInfo::storeVegaMatrix(VegaMatrixLiteSP sens,
                                  CInstrument* inst,
                                  CModelLN* model,
                                  VanillaContractsRecorderSP recorder,
                                  CAssetSP asset,
                                  const DateTime& valueDate,
                                  const DateTime& instExpiry,
                                  CResults* results) {
    static const string method = "VanillaInfo::storeVegaMatrix";
    try {
        // Make sure (again) instrument supports the appropriate interfaces
        ISensitiveStrikes* p = dynamic_cast<ISensitiveStrikes*>(inst);
        ISupportVegaMatrixLite* vmlInstrument = dynamic_cast<ISupportVegaMatrixLite*>(inst);

        QLIB_VERIFY(p!=NULL, 
            "The instrument does not support the interface ISensitiveStrikes.");
        QLIB_VERIFY(!p->avoidVegaMatrix(model), 
            "The instrument is set to avoid the vega matrix computation.");
        QLIB_VERIFY(vmlInstrument!=NULL, 
            "The instrument does not support the interface ISupportVegaMatrixLite.");

        // Create a VegaMatrix to figure out the sensitive assets & expiries
        double shiftSize = sens->getShiftSize();
        VegaMatrix vm(shiftSize);
        SensMgrConst sensMgr(asset);
        OutputNameArrayConstSP names;
        names = sensMgr.allNames(&vm);
        if (names->empty()) {
            results->storeNotApplicable(sens.get());
            return;
        }

        // Allow only single names i.e. this methodology does not work for XCBs, Funds etc.
        if (names->size() != 1) {
            throw ModelException(method, "Only single name underlyings are supported");
        }

        // Work with single name now
        OutputNameSP name = (*names)[0];

        // Extract sensitive strikes: unfortunately this is a copy past of the VegaMatrix code
        // Get them from the instument, discard close strikes
        DoubleArraySP raw = p->getSensitiveStrikes(name, model);
        // now strip duplicates
        sort(raw->begin(), raw->end());
        int xLowerIdx = 0;
        int xUpperIdx = 0;
        CDoubleArraySP sensStrikes(new DoubleArray(0));
        while (xUpperIdx < raw->size()) {
            xLowerIdx = xUpperIdx;

            /* Find the upper x idx to tweak. This will do
                nothing if adjacent strikes are
                sufficiently far apart */
            while( (xUpperIdx+1 < raw->size() ) &&
                    (((*raw)[xUpperIdx+1]/(*raw)[xUpperIdx])
                    <= MatrixShift::BETA)) {
                xUpperIdx++;
            }
            sensStrikes->push_back((*raw)[xLowerIdx]);
            xUpperIdx++;
        }
        DoubleArray& strikes = *sensStrikes;

        // Extract sensitive maturities
        TweakGroupSP group(new TweakGroup(CInstrumentSP(inst), IModelSP(model)));
        vm.setMarketDataName(name);
        ExpiryArrayConstSP expiries = vm.getExpiries(group.get());

        // Get dimensions for output
        int nbStrikes = strikes.size();
        int nbMaturities = expiries->size();
        CDoubleMatrixSP matrixOut(new CDoubleMatrix(nbStrikes, nbMaturities));
            
        // Compute if we are sensitive to something and shift size is not zero
        if(nbStrikes > 0 && nbMaturities > 0 && !Maths::isZero(shiftSize)) {
            // Compute implied vol at benchmark strikes and maturities
            LinearStrikeVolRequestSP volReq(new 
                LinearStrikeVolRequest(1.0, valueDate, instExpiry, false));
            volReq->allowNegativeFwdVar(model->negativeFwdVarAllowed());
            CVolProcessedBSSP volBS(asset->getProcessedVol(volReq.get()));
            TimeMetricSP metric = TimeMetricSP::constCast(volBS->GetTimeMetric());
            
            CDoubleMatrixSP backboneVols(new CDoubleMatrix(nbStrikes, nbMaturities));
            for(int k = 0; k < nbStrikes; k++) {
                volReq->setStrike(strikes[k]);
                volBS = CVolProcessedBSSP(asset->getProcessedVol(volReq.get()));
                for(int t = 0; t < nbMaturities; t++) {
                    DateTime maturity = (*expiries)[t]->toDate(valueDate);
                    double vol = volBS->CalcVol(valueDate, maturity);
                    (*backboneVols)[k][t] = vol;
                }
            }
            
            // Compute vega matrix of vanilla portfolio
            recorder->scaleNbContracts(0.01);
            VanillaInfoArraySP options = recorder->getResults();

            CDoubleMatrix& vegaMatrix = *matrixOut;
            const CDoubleMatrix& vols = *backboneVols;

            // Create maturities
            DoubleArray yearFracs((nbMaturities));
            unsigned long iMat, iStrike;
            for(iMat = 0; iMat < nbMaturities; iMat++) {
                DateTime thisDate = (*expiries)[iMat]->toDate(valueDate);
                yearFracs[iMat] = metric->yearFrac(valueDate, thisDate);
            }
            
            // Loop over options and bucket them
            for(int i = 0; i < options->size(); i++) {
                const VanillaInfo& option = *(*options)[i];
                double strike = option.strike;
                double yearFrac = option.yearFrac;
                double vol = option.vol;
                
                // Find buckets
                iMat = 0;
                iStrike = 0;
                locate(&yearFracs[0] - 1, nbMaturities, yearFrac, &iMat);   // iMat is in [0, nbMaturities]
                locate(&strikes[0] - 1, nbStrikes, strike, &iStrike);       // iStrike is in [0, nbStrikes]

                if(iMat == 0) {
                    if(iStrike == 0) {
                        // Lower left
                        vegaMatrix[0][0] += option.vega(shiftSize);
                    } else if(iStrike == nbStrikes) {
                        // Lower right
                        vegaMatrix[nbStrikes-1][0] += option.vega(shiftSize);
                    } else {
                        // Distribute to neighboring strikes
                        double volBase = 
                            volInterpStrike(strikes[iStrike-1], strikes[iStrike], vols[iStrike-1][0], vols[iStrike][0], strike);
                        double dvolLeft = 
                            volInterpStrike(strikes[iStrike-1], strikes[iStrike], vols[iStrike-1][0] + shiftSize, vols[iStrike][0], strike) - volBase;
                        double dvolRight = 
                            volInterpStrike(strikes[iStrike-1], strikes[iStrike], vols[iStrike-1][0], vols[iStrike][0] + shiftSize, strike) - volBase;
                        vegaMatrix[iStrike-1][0] += option.vega(dvolLeft) * dvolLeft / shiftSize;
                        vegaMatrix[iStrike][0]   += option.vega(dvolRight) * dvolRight / shiftSize;
                    }
                } else if(iMat == nbMaturities) {
                    if(iStrike == 0) {
                        // Upper left
                        vegaMatrix[0][nbMaturities-1] += option.vega(shiftSize);
                    } else if(iStrike == nbStrikes) {
                        // Upper right
                        vegaMatrix[nbStrikes-1][nbMaturities-1] += option.vega(shiftSize);
                    } else {
                        // Distribute to neighboring strikes
                        double volBase = 
                            volInterpStrike(strikes[iStrike-1], strikes[iStrike], vols[iStrike-1][nbMaturities-1], vols[iStrike][nbMaturities-1], strike);
                        double dvolLeft = 
                            volInterpStrike(strikes[iStrike-1], strikes[iStrike], vols[iStrike-1][nbMaturities-1] + shiftSize, vols[iStrike][nbMaturities-1], strike) - volBase;
                        double dvolRight = 
                            volInterpStrike(strikes[iStrike-1], strikes[iStrike], vols[iStrike-1][nbMaturities-1], vols[iStrike][nbMaturities-1] + shiftSize, strike) - volBase;
                        vegaMatrix[iStrike-1][nbMaturities-1] += option.vega(dvolLeft) * dvolLeft / shiftSize;
                        vegaMatrix[iStrike][nbMaturities-1]   += option.vega(dvolRight) * dvolRight / shiftSize;
                    }
                } else {
                    if(iStrike == 0) {
                        // Distribute to neighboring maturities
                        double volBase = volInterpMaturity(yearFracs[iMat-1], yearFracs[iMat], vols[0][iMat-1], vols[0][iMat], yearFrac);
                        double dvolFront = 
                            volInterpMaturity(yearFracs[iMat-1], yearFracs[iMat], vols[0][iMat-1] + shiftSize, vols[0][iMat], yearFrac) - volBase;
                        double dvolBack = 
                            volInterpMaturity(yearFracs[iMat-1], yearFracs[iMat], vols[0][iMat-1], vols[0][iMat] + shiftSize, yearFrac) - volBase;
                        
                        vegaMatrix[0][iMat-1] += option.vega(dvolFront) * dvolFront / shiftSize;
                        vegaMatrix[0][iMat]   += option.vega(dvolBack) * dvolBack / shiftSize;
                    } else if(iStrike == nbStrikes) {
                        // Distribute to neighboring maturities
                        double volBase = volInterpMaturity(yearFracs[iMat-1], yearFracs[iMat], vols[nbStrikes-1][iMat-1], vols[nbStrikes-1][iMat], yearFrac);
                        double dvolFront = 
                            volInterpMaturity(yearFracs[iMat-1], yearFracs[iMat], vols[nbStrikes-1][iMat-1] + shiftSize, vols[nbStrikes-1][iMat], yearFrac) - volBase;
                        double dvolBack = 
                            volInterpMaturity(yearFracs[iMat-1], yearFracs[iMat], vols[nbStrikes-1][iMat-1], vols[nbStrikes-1][iMat] + shiftSize, yearFrac) - volBase;
                        
                        vegaMatrix[nbStrikes-1][iMat-1] += option.vega(dvolFront) * dvolFront / shiftSize;
                        vegaMatrix[nbStrikes-1][iMat]   += option.vega(dvolBack) * dvolBack / shiftSize;
                    } else {
                        // Distribute to neighboring strikes & maturities
                        double volBase = 
                            volInterp(yearFracs[iMat-1], yearFracs[iMat], strikes[iStrike-1], strikes[iStrike],
                                    vols[iStrike-1][iMat-1], vols[iStrike][iMat-1], vols[iStrike-1][iMat], vols[iStrike][iMat],
                                    yearFrac, strike);
                        double dvolFrontLeft = 
                            volInterp(yearFracs[iMat-1], yearFracs[iMat], strikes[iStrike-1], strikes[iStrike],
                                    vols[iStrike-1][iMat-1] + shiftSize, vols[iStrike][iMat-1], vols[iStrike-1][iMat], vols[iStrike][iMat],
                                    yearFrac, strike) - volBase;
                        double dvolFrontRight = 
                            volInterp(yearFracs[iMat-1], yearFracs[iMat], strikes[iStrike-1], strikes[iStrike],
                                    vols[iStrike-1][iMat-1], vols[iStrike][iMat-1] + shiftSize, vols[iStrike-1][iMat], vols[iStrike][iMat],
                                    yearFrac, strike) - volBase;
                        double dvolBackLeft = 
                            volInterp(yearFracs[iMat-1], yearFracs[iMat], strikes[iStrike-1], strikes[iStrike],
                                    vols[iStrike-1][iMat-1], vols[iStrike][iMat-1], vols[iStrike-1][iMat] + shiftSize, vols[iStrike][iMat],
                                    yearFrac, strike) - volBase;
                        double dvolBackRight = 
                            volInterp(yearFracs[iMat-1], yearFracs[iMat], strikes[iStrike-1], strikes[iStrike],
                                    vols[iStrike-1][iMat-1], vols[iStrike][iMat-1], vols[iStrike-1][iMat], vols[iStrike][iMat] + shiftSize,
                                    yearFrac, strike) - volBase;

                        vegaMatrix[iStrike-1][iMat-1] += option.vega(dvolFrontLeft) * dvolFrontLeft / shiftSize;
                        vegaMatrix[iStrike][iMat-1]   += option.vega(dvolFrontRight) * dvolFrontRight / shiftSize;
                        vegaMatrix[iStrike-1][iMat]   += option.vega(dvolBackLeft) * dvolBackLeft / shiftSize;
                        vegaMatrix[iStrike][iMat]     += option.vega(dvolBackRight) * dvolBackRight / shiftSize;
                    }
                }
            }
        }
        
        // Package results and return
        MatrixResultSP matrixResult(new MatrixResult(sensStrikes, expiries, matrixOut));
            
        // Store result
        results->storeGreek(matrixResult, sens->getSensOutputName(), name);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}
#endif


void VanillaInfo::storeVegaMatrix(VegaMatrixLiteSP sens,
                                  CInstrument* inst,
                                  CModelLN* model,
                                  VanillaContractsRecorderSP recorder,
                                  CAssetSP asset,
                                  const DateTime& valueDate,
                                  const DateTime& instExpiry,
                                  CResults* results) {
    static const string method = "VanillaInfo::storeVegaMatrix";
    try {
        // Make sure (again) instrument supports the appropriate interfaces
        ISensitiveStrikes* p = dynamic_cast<ISensitiveStrikes*>(inst);
        ISupportVegaMatrixLite* vmlInstrument = dynamic_cast<ISupportVegaMatrixLite*>(inst);

        QLIB_VERIFY(p!=NULL, 
            "The instrument does not support the interface ISensitiveStrikes.");
        QLIB_VERIFY(!p->avoidVegaMatrix(model), 
            "The instrument is set to avoid the vega matrix computation.");
        QLIB_VERIFY(vmlInstrument!=NULL, 
            "The instrument does not support the interface ISupportVegaMatrixLite.");

        // Create a VegaMatrix to figure out the sensitive assets & expiries
        double shiftSize = sens->getShiftSize();
        VegaMatrix vm(shiftSize);
        SensMgrConst sensMgr(asset);
        OutputNameArrayConstSP names;
        names = sensMgr.allNames(&vm);
        if (names->empty()) {
            results->storeNotApplicable(sens.get());
            return;
        }

        // Allow only single names i.e. this methodology does not work for XCBs, Funds etc.
        if (names->size() != 1) {
            throw ModelException(method, "Only single name underlyings are supported");
        }

        // Work with single name now
        OutputNameSP name = (*names)[0];

        // Extract sensitive strikes: unfortunately this is a copy past of the VegaMatrix code
        // Get them from the instument, discard close strikes
        DoubleArraySP raw = p->getSensitiveStrikes(name, model);
        DoubleArray& rawStrikes = *raw;
        int nbRawStrikes = rawStrikes.size();
        sort(rawStrikes.begin(), rawStrikes.end());
        int xLowerIdx = 0;
        int xUpperIdx = 0;
        CDoubleArraySP sensStrikes(new DoubleArray(0));
        while (xUpperIdx < rawStrikes.size()) {
            xLowerIdx = xUpperIdx;

            /* Find the upper x idx to tweak. This will do
                nothing if adjacent strikes are
                sufficiently far apart */
            while( (xUpperIdx+1 < rawStrikes.size() ) &&
                    (rawStrikes[xUpperIdx+1] / rawStrikes[xUpperIdx]) 
                    <= MatrixShift::BETA) {
                xUpperIdx++;
            }
            sensStrikes->push_back(rawStrikes[xLowerIdx]);
            xUpperIdx++;
        }
        DoubleArray& strikes = *sensStrikes;

        // Extract sensitive maturities
        TweakGroupSP group(new TweakGroup(CInstrumentSP(inst), IModelSP(model)));
        vm.setMarketDataName(name);
        ExpiryArrayConstSP expiries = vm.getExpiries(group.get());

        // Get dimensions for output
        int nbStrikes = strikes.size();
        int nbMaturities = expiries->size();
        CDoubleMatrixSP matrixOut(new CDoubleMatrix(nbStrikes, nbMaturities));
            
        // Compute if we are sensitive to something and shift size is not zero
        if(nbStrikes > 0 && nbMaturities > 0 && !Maths::isZero(shiftSize)) {
            // Compute implied vol at benchmark strikes and maturities
            LinearStrikeVolRequestSP volReq(new 
                LinearStrikeVolRequest(1.0, valueDate, instExpiry, false));
            volReq->allowNegativeFwdVar(model->negativeFwdVarAllowed());
            CVolProcessedBSSP volBS(asset->getProcessedVol(volReq.get()));
            TimeMetricSP metric = TimeMetricSP::constCast(volBS->GetTimeMetric());
            
            CDoubleMatrixSP backboneVols(new CDoubleMatrix(nbRawStrikes, nbMaturities));
            for(int k = 0; k < nbStrikes; k++) {
                volReq->setStrike(rawStrikes[k]);
                volBS = CVolProcessedBSSP(asset->getProcessedVol(volReq.get()));
                for(int t = 0; t < nbMaturities; t++) {
                    DateTime maturity = (*expiries)[t]->toDate(valueDate);
                    double vol = volBS->CalcVol(valueDate, maturity);
                    (*backboneVols)[k][t] = vol;
                }
            }
            
            // Compute vega matrix of vanilla portfolio
            recorder->scaleNbContracts(0.01);
            VanillaInfoArraySP options = recorder->getResults();

            CDoubleMatrix& vegaMatrix = *matrixOut;
            const CDoubleMatrix& vols = *backboneVols;

            // Create maturities
            DoubleArray yearFracs((nbMaturities));
            for(int iMat = 0; iMat < nbMaturities; iMat++) {
                DateTime thisDate = (*expiries)[iMat]->toDate(valueDate);
                yearFracs[iMat] = metric->yearFrac(valueDate, thisDate);
            }
            
            // Loop over options and bucket them
            for(int i = 0; i < options->size(); i++) {
                // Copy over some info from vanilla
                const VanillaInfo& option = *(*options)[i];
                double strike = option.strike;
                double yearFrac = option.yearFrac;
                
                // Bracket maturity, sensStrike, rawStrike
                int shortMat, longMat;
                int leftSensStrike, rightSensStrike;
                int leftRawStrike, rightRawStrike;
                VanillaInfo::bracketArray(yearFracs, yearFrac, shortMat, longMat);
                VanillaInfo::bracketArray(strikes, strike, leftSensStrike, rightSensStrike);
                VanillaInfo::bracketArray(rawStrikes, strike, leftRawStrike, rightRawStrike);
                
                // Distribute to neighboring strikes & maturities
                
                // 1) Untweaked vol
                double volBase = 
                    volInterp(yearFracs[shortMat], yearFracs[longMat], rawStrikes[leftRawStrike], rawStrikes[rightRawStrike],
                              vols[leftRawStrike][shortMat], vols[rightRawStrike][shortMat], vols[leftRawStrike][longMat], vols[rightRawStrike][longMat],
                              yearFrac, strike);
                
                // 2) Shift due to short maturity, left strike
                double dvolShortLeft = 
                    volInterp(yearFracs[shortMat], yearFracs[longMat], rawStrikes[leftRawStrike], rawStrikes[rightRawStrike],
                              vols[leftRawStrike][shortMat] + shiftSize, vols[rightRawStrike][shortMat], vols[leftRawStrike][longMat], vols[rightRawStrike][longMat],
                              yearFrac, strike) - volBase;
                
                // 3) Shift due to short maturity, right strike
                double dvolShortRight = 
                    volInterp(yearFracs[shortMat], yearFracs[longMat], rawStrikes[leftRawStrike], rawStrikes[rightRawStrike],
                              vols[leftRawStrike][shortMat], vols[rightRawStrike][shortMat] + shiftSize, vols[leftRawStrike][longMat], vols[rightRawStrike][longMat],
                              yearFrac, strike) - volBase;
                
                // 4) Shift due to long maturity, left strike
                double dvolLongLeft = 
                    volInterp(yearFracs[shortMat], yearFracs[longMat], rawStrikes[leftRawStrike], rawStrikes[rightRawStrike],
                              vols[leftRawStrike][shortMat], vols[rightRawStrike][shortMat], vols[leftRawStrike][longMat] + shiftSize, vols[rightRawStrike][longMat],
                              yearFrac, strike) - volBase;
                
                // 5) Shift due to long maturity, right strike
                double dvolLongRight = 
                    volInterp(yearFracs[shortMat], yearFracs[longMat], rawStrikes[leftRawStrike], rawStrikes[rightRawStrike],
                              vols[leftRawStrike][shortMat], vols[rightRawStrike][shortMat], vols[leftRawStrike][longMat], vols[rightRawStrike][longMat] + shiftSize,
                              yearFrac, strike) - volBase;
 
                // Record the 4 vegas
                vegaMatrix[leftSensStrike][shortMat]    += option.vega(dvolShortLeft) * dvolShortLeft / shiftSize;
                vegaMatrix[rightSensStrike][shortMat]   += option.vega(dvolShortRight) * dvolShortRight / shiftSize;
                vegaMatrix[leftSensStrike][longMat]     += option.vega(dvolLongLeft) * dvolLongLeft / shiftSize;
                vegaMatrix[rightSensStrike][longMat]    += option.vega(dvolLongRight) * dvolLongRight / shiftSize;
            }
        }
        
        // Package results and store them
        MatrixResultSP matrixResult(new MatrixResult(sensStrikes, expiries, matrixOut));
        results->storeGreek(matrixResult, sens->getSensOutputName(), name);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


void VanillaInfo::load(CClassSP& clazz){
    REGISTER(VanillaInfo, clazz);
    SUPERCLASS(CObject);
    FIELD(maturity, "Maturity");
    FIELD(nbContracts, "Number of contracts");
    FIELD(isCall, "True: call, False: put");
    FIELD(fwd, "Forward price");
    FIELD(strike, "Absolute strike");
    FIELD(vol, "Implied vol");
    FIELD(yearFrac, "Year fraction");
    EMPTY_SHELL_METHOD(defaultVanillaInfo);
    clazz->setPublic();
}


IObject* VanillaInfo::defaultVanillaInfo() {
    return new VanillaInfo();
}


CClassConstSP const VanillaInfo::TYPE = CClass::registerClassLoadMethod(
    "VanillaInfo", typeid(VanillaInfo), VanillaInfo::load);


DEFINE_TEMPLATE_TYPE(VanillaInfoArray);


//////////////////////////////////////////////////////////////////////////////////////////////////////


VanillaContractsRecorder::VanillaContractsRecorder(double scalingFactor, bool record):
CObject(TYPE), scalingFactor(scalingFactor), record(record) {
    options = VanillaInfoArraySP(new VanillaInfoArray());
}


VanillaContractsRecorderSP VanillaContractsRecorder::createVanillaOptionRecorder(const Control* control) {
    static const string method = "VanillaContractsRecorder::createVanillaOptionRecorder";
    
    try {
        // Figure out if we are doing VegaMatrixLite
        bool doingVegaMatrixLite = false;
        if (control)
        {
            SensitivitySP vml = control->getCurrentSensitivity();
            doingVegaMatrixLite = !!vml;
        }        
        // Create recorder and return it
        return VanillaContractsRecorderSP(new VanillaContractsRecorder(1.0, doingVegaMatrixLite));
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


void VanillaContractsRecorder::recordContract(VanillaInfoSP optionInfo) {
    if(record) {
        optionInfo->scaleNbContracts(scalingFactor);
        options->push_back(optionInfo);
    }
}


void VanillaContractsRecorder::recordContract(
    const DateTime& maturity,
    double nbContracts,
    bool isCall,
    double fwd,
    double strike,
    double vol,
    double yearFrac) {
    if(record) {
        double nbOptions = scalingFactor * nbContracts;
        // Don't bother if it's a zero position
        if(!Maths::isZero(nbOptions)) {
            VanillaInfoSP optionInfo(new VanillaInfo(
                maturity, nbOptions, isCall, fwd, strike, vol, yearFrac));
            options->push_back(optionInfo);
        }
    }
}


double VanillaContractsRecorder::getScalingFactor() const {
    return scalingFactor;
}


void VanillaContractsRecorder::setScalingFactor(double newScalingFactor) {
    scalingFactor = newScalingFactor;
}


void VanillaContractsRecorder::scaleNbContracts(double x) {
    for(int i = 0; i < options->size(); i++) {
        (*options)[i]->scaleNbContracts(x);
    }
}


VanillaInfoArraySP VanillaContractsRecorder::getResults() const {
    return options;
}


void VanillaContractsRecorder::load(CClassSP& clazz) {
    REGISTER(VanillaContractsRecorder, clazz);
    SUPERCLASS(CObject);
    FIELD(scalingFactor, "Scaling factor");
    FIELD(record, "True: record entries, False: ignore");
    EMPTY_SHELL_METHOD(defaultVanillaContractsArray);
    clazz->setPublic();
}


VanillaContractsRecorder::VanillaContractsRecorder(): CObject(TYPE) {}


IObject* VanillaContractsRecorder::defaultVanillaContractsArray() {
    return new VanillaContractsRecorder();
}


CClassConstSP const VanillaContractsRecorder::TYPE = CClass::registerClassLoadMethod(
    "VanillaContractsRecorder", typeid(VanillaContractsRecorder), VanillaContractsRecorder::load);


//////////////////////////////////////////////////////////////////////////////////////////////////////


double BlackPrice(bool isCall, double fwd, double strike, double pv, double vol, double yearFrac,
    double nbContracts, const DateTime& maturity, VanillaContractsRecorderSP recorder) {
    // Price & record
    double price = 0.0;
    if(!Maths::isZero(nbContracts)) {
        price = Black::price(isCall, fwd, strike, pv, Maths::square(vol) * yearFrac);
    }
    recorder->recordContract(maturity, nbContracts, isCall, fwd, strike, vol, yearFrac);
    
    // Return value of position
    double value = price * nbContracts;
    return value;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////


StaticReplicationIntegrand::StaticReplicationIntegrand(
    Function1DDoubleConstSP     relativeStrikeWeight,
    VanillaContractsRecorderSP  recorder,
    CAssetConstSP               asset,
    const DateTime&             valueDate,
    const DateTime&             maturity,
    bool                        negativeVar):
    Function1DDouble(relativeStrikeWeight->getInterval()), 
    relativeStrikeWeight(relativeStrikeWeight), recorder(recorder), asset(asset), valueDate(valueDate), 
    maturity(maturity) {
    static const string routine = "StaticReplicationIntegrand::StaticReplicationIntegrand";
    try {
        // Compute some transient stuff
        fwd = asset->fwdValue(maturity);
        volRequest = LinearStrikeVolRequestSP(new 
            LinearStrikeVolRequest(1.0, valueDate, maturity, false));
        //test for negativeVar
        if(negativeVar){
            volRequest->allowNegativeFwdVar(true);
        }
        // Validate the boundaries
        const Range& interval = relativeStrikeWeight->getInterval();
        const Boundary& lowerBound = interval.getLower();
        const Boundary& upperBound = interval.getUpper();

        // Lower bound cannot be minus infinity
        if(lowerBound.isInfinite()) {
            throw ModelException("Lower bound cannot be infinite");
        }

        // Check value of lower bound
        double low = lowerBound.getValue();
        if(low<=0.0) {
            if(low<0.0) {
                throw ModelException(routine, 
                    "Lower bound " + Format::toString(low) + "cannot be negative.");
            } else if(lowerBound.isClosedBracket()) {
                throw ModelException(routine, "Lower bound cannot be closed 0.");
            }
        }

        // Finally check that upper bound is greater than lower bound
        if(!upperBound.isInfinite()) {
            double high = upperBound.getValue();
            if(high < low) {
                throw ModelException(routine, 
                    "Upper bound " + Format::toString(high) +
                    "must be greater then lower bound" + Format::toString(low));
            }
        }

        // Get trading time
        ATMVolRequest volReq;
        IVolProcessedSP vol(asset->getProcessedVol(&volReq));
        yearFrac = vol->calcTradingTime(valueDate, maturity);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


double StaticReplicationIntegrand::operator()(double relativeStrike) const {
    static const string routine = "StaticReplicationIntegrand::operator()";
    try{
        if (IMSLError::isError()) {
            // if there is already an error, dont continue, but just return value to integrator
            return 0.0;
        }
        bool isCall;
        if (relativeStrike < 1.0) {
            isCall = false;
        } else {
            isCall = true;
        }
        
        double absoluteStrike = relativeStrike * fwd;
        
        volRequest->setStrike(absoluteStrike);
        CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
        double vol = volBS->CalcVol(valueDate, maturity);
        
        double weight = (*relativeStrikeWeight)(relativeStrike) / fwd;
        double integrand = BlackPrice(isCall, fwd, absoluteStrike, 1.0, vol, yearFrac,
            weight, maturity, recorder);
        
        return integrand;
    } catch (exception& e) {
        // do not throw an exception here,
        // just store it and return arbitraty value to integrator so that there is no failure
        string message = routine + ": Failed";
        IMSLError::appendToStack(ModelException::addTextToException(e, message));
        return 0.0;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////


Trapez1DSimpleRecorder::Trapez1DSimpleRecorder(int nbSteps, VanillaContractsRecorderSP recorder): 
Trapez1DSimple(TYPE, nbSteps), recorder(recorder) {}


double Trapez1DSimpleRecorder::integrate(const Function1DDouble& func) const {
    const StaticReplicationIntegrand* ptr = dynamic_cast<const StaticReplicationIntegrand*>(&func);
    if(ptr) {
        return integrateAndRecord(*ptr);
    } else {
        return Trapez1DSimple::integrate(func);
    }
}


double Trapez1DSimpleRecorder::integrateAndRecord(const StaticReplicationIntegrand& func) const {
    static const string method = "Trapez1DSimpleRecorder::integrateAndRecord";
    try{
        const Range& interval = func.getInterval();

        Range::checkIsNonEmpty(interval);

        if (interval.isSingleton()){
            return 0.0;
        }

        const Boundary& lower = interval.getLower();
        const Boundary& upper = interval.getUpper();

        if (!lower.isClosedBracket() || lower.isInfinite()
            || !upper.isClosedBracket() || upper.isInfinite()){
            throw ModelException(method,
                                 "(Semi-)Open and / or infinite intervals are not supported; got " + interval.toString());
        }

        double a = lower.getValue();
        double b = upper.getValue();

        double del = (b - a) / nbSteps;

        // Change scaling for mid
        double scalingFactor = recorder->getScalingFactor();
        double midScalingFactor = scalingFactor * del;
        recorder->setScalingFactor(midScalingFactor);
        
        // Interior Sum
        double intSum = 0.0;
        double x = a;
        for (int j = 0; j < nbSteps - 1; j++) {
            x += del;
            intSum += func(x);
        }
        intSum *= 2.0;
        
        // Change scaling for tail
        double tailScalingFactor = scalingFactor * del * 0.5;
        recorder->setScalingFactor(tailScalingFactor);
        
        // Endpoints Contribution
        double extSum = func(a) + func(b);
        
        // Reset scaling factor
        recorder->setScalingFactor(scalingFactor);
        
        // Integral
        double integral = 0.5 * del * (intSum + extSum);
        return integral;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}


void Trapez1DSimpleRecorder::load(CClassSP& clazz) {
    REGISTER(Trapez1DSimpleRecorder, clazz);
    SUPERCLASS(Trapez1DSimple);
    EMPTY_SHELL_METHOD(defaultTrapez1DSimpleRecorder);
    FIELD(recorder, "Option recorder");
    clazz->setPublic();
}


IObject* Trapez1DSimpleRecorder::defaultTrapez1DSimpleRecorder() {
    return new Trapez1DSimpleRecorder();
}


Trapez1DSimpleRecorder::Trapez1DSimpleRecorder(): Trapez1DSimple(TYPE) {}


CClassConstSP const Trapez1DSimpleRecorder::TYPE = CClass::registerClassLoadMethod(
    "Trapez1DSimpleRecorder", typeid(Trapez1DSimpleRecorder), load);


//////////////////////////////////////////////////////////////////////////////////////////////////////


DoubleArraySP VarSwapUtilities::computeHistoricLogReturns(const CAsset*        asset,
                                                          const DateTimeArray& obsDates,
                                                          const DoubleArray&   obsSamples,
                                                          const DateTime&      valueDate,
                                                          bool                 computeCurrentReturn,
                                                          bool                 dividendAdjusted,
                                                          bool                 divAdjOnExDate) {
    static const string routine = "VarSwapUtilities::computeHistoricLogReturns";
    try {
        DoubleArraySP logReturns(new DoubleArray(obsDates.size()));
        
        // Return if no history
        if (obsDates.size() == 0 || valueDate <= obsDates[0]) {
            return logReturns;
        }
        
        VolVarShell::DivAdjuster divAdjuster;
        bool onExpiry = obsDates.back().equals(valueDate);
        for(int iStep = 1; iStep < obsDates.size(); iStep++) {
            if(valueDate < obsDates[iStep-1]) {
                // Strictly future return so nothing to do here
            } else {
                // Start samples
                DateTime startDate = obsDates[iStep-1];
                double startSample = obsSamples[iStep-1];
                DateTime endDate;
                double endSample;
                if(obsDates[iStep] < valueDate || onExpiry) {
                    // Strictly past return so compute based on samples
                    endDate = obsDates[iStep];
                    endSample = obsSamples[iStep];
                } else {
                    // Partly past & partly future: use current spot
                    endDate = valueDate;
                    endSample = asset->getSpot();
                }
                
                if(dividendAdjusted) {
                    // Compute cumulative dividend in (startDate, endDate]
                    double divAmount = getHistoricDivsBetweendates(&divAdjuster, asset,
                        valueDate, startDate, endDate);
                    
                    // Add dividend jump back to the return
                    if(divAdjOnExDate) {
                        endSample += divAmount;
                    } else {
                        startSample -= divAmount;
                    }
                }
                
                if(obsDates[iStep] < valueDate || computeCurrentReturn) {
                    // For strictly past returns or for current return when requested
                    (*logReturns)[iStep] = log(endSample / startSample);
                }
            }
        }

        return logReturns;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

// Computes realised var until today. Wraps computeHistoricLogReturns()
double VarSwapUtilities::realisedVar(const CAsset*        asset,
                                     const DateTimeArray& obsDates,
                                     const DoubleArray&   obsSamples,
                                     const DateTime&      valueDate,
                                     bool                 computeCurrentReturn,
                                     bool                 dividendAdjusted,
                                     bool                 divAdjOnExDate) {
    static const string routine = "VarSwapUtilities::realisedVar";
    try {
        // Compute log-returns for past samples
        DoubleArraySP logReturns =  VarSwapUtilities::computeHistoricLogReturns(
                asset, obsDates, obsSamples, valueDate, true, dividendAdjusted, false);

        // Compute past realized variance
        double histVar = 0.0;
        int iStep;
        for(iStep = 1; iStep < obsDates.size(); iStep++) {
            if(obsDates[iStep] < valueDate) {
                histVar += Maths::square((*logReturns)[iStep]); 
            } else {
                break;
            }
        }

        // and add current contribution
        if(valueDate <= obsDates.back()){
		    histVar += Maths::square((*logReturns)[iStep]);
	    }

        return histVar;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

/** Computes sum of dividends between two historic dates */			
double VarSwapUtilities::getHistoricDivsBetweendates(
    DividendList::IDivAdjuster* divAdjuster,
    const CAsset* asset,
    const DateTime& valueDate,
    const DateTime& startDate,
    const DateTime& endDate) {
    
    static const string routine = "VarSwapUtilities::getHistoricDivsBetweendates";
    try {
        // Compute cumulated historic divs between 2 samples
        DividendCollector collector(asset,
                                    divAdjuster,
                                    valueDate,
                                    startDate,
                                    endDate);
        
        DividendListSP divs = collector.getDividends();
        const DividendArray& divArray = divs->getArray();
        
        double divAmount = 0.0;
        for(int iDiv = 0; iDiv < divArray.size(); iDiv++) {
            divAmount += divArray[iDiv].getDivAmount();
        }  

        return divAmount;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


double VarSwapUtilities::futureDivsVariance(const CAsset*   asset,
                                            const DateTime& valueDate,
                                            const DateTime& firstDate,
                                            const DateTime& lastDate) {
    static const string routine = "VarSwapUtilities::futureDivsVariance";
    try {
                // extract all divs
        VolVarShell::DivAdjuster divAdjuster;

        // Then initialise our DividendCollector
        DividendCollector collector(asset,
                                    &divAdjuster, // optional
                                    valueDate,
                                    firstDate < valueDate ? valueDate : firstDate,
                                    lastDate);

        DividendListSP divs = collector.getDividends();

        const DividendArray& divArray = divs->getArray();

        double divVar = 0.0;
        for (int i = 0; i < divArray.size(); i++){
            // now compute the yield dividend, or the equivalent yield dividend if AMOUNT
            double amount = divArray[i].getDivAmount();
            if(divArray[i].getDivType() == Dividend::PERCENT) {
                // do nothing.
            } else if(divArray[i].getDivType() == Dividend::AMOUNT) {
                DateTime exDate = divArray[i].getExDate();
                DateTime prevDate = DateTime(exDate.getDate(), DateTime::BEFORE_EX_DIV_TIME);

                try {
                    amount = amount / (asset->fwdValue(exDate) + amount);
                } catch (exception& e) {
                    if(Asset::IQuanto::TYPE->isInstance(*asset)) {
                        const Asset::IQuanto* ptr = dynamic_cast<const Asset::IQuanto*>(asset);
                        double fwd = ptr->unadjustedFwdValue(exDate);
                        double prevFwd = ptr->unadjustedFwdValue(prevDate);
                        amount = (prevFwd - fwd) / prevFwd;
                    } else {
                        throw ModelException(e, routine);
                    }
                }
            } else {   // continuous divs get 0 contribution
                amount = 0.0;
            }
            if (!Maths::isPositive(1.0 - amount)) {
                throw ModelException(routine, "Div yield cannot be more than 100%");
            }
            double lg = log(1.0 - amount);
            divVar += Maths::square(lg);
        }

        return divVar;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}
   
    
void VarSwapUtilities::pastAndfutureExpN(const DateTimeArray&  obsDates,
                                         const HolidayWrapper& assetHols,
                                         const DateTime&       valueDate,
                                         int                   expectedN,
                                         int&                  pastExpectedN,
                                         int&                  futureExpectedN) {
    static const string routine = "VarSwapUtilities::pastAndfutureExpN";
    try {
        const DateTime& firstDate = obsDates.front();
        const DateTime& lastDate  = obsDates.back();
        
        if (lastDate <= valueDate) {
            // All past
            pastExpectedN   = expectedN;
        } else if (firstDate < valueDate) {
            // Partly past partly future
            
            // Find next sample date
            int i = 0;
            while (i < obsDates.size() && obsDates[i] < valueDate) {
                i++;
            }

            // expected number of returns from next sample date
            int numFutureTemp = assetHols->businessDaysDiff(obsDates[i], lastDate);

            // number of returns to next sample is defined as:
            // (total num returns expected on instrument creation) - (expected num returns from next sample)
            int numPastTemp = expectedN - numFutureTemp;
            pastExpectedN = Maths::max(numPastTemp, 1);
        } else {
            // Default case of everything in the future
            pastExpectedN = 0;
        }
        // NumFuture returns defined as remaining returns
        futureExpectedN = expectedN - pastExpectedN;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void VarSwapUtilities::thetaShiftCashFlows(const Theta*         shift,
                                           const CAsset*        asset,
                                           const DateTime&      valueDate,
                                           const DateTimeArray& obsDates,
                                           DoubleArray&         obsSamples) {
    static const string routine = "VarSwapUtilities::thetaShiftCashFlows";
    try {
        DateTime newDate = shift->rollDate(valueDate);
        
        // fill vol sample point if needed
        double spot = 0.0;
        bool useSpot = !shift->useAssetFwds();

        if (useSpot) {
            spot = asset->getSpot();
        }

        bool pastRollDate = false;
        int  i = 0;
        while (!pastRollDate && i < obsDates.size() ) {
            if ((obsDates[i].isGreater(valueDate) &&
                !obsDates[i].isGreater(newDate))  ||
                (obsDates[i].equals(valueDate)    &&
                Maths::isZero(obsSamples[i]))) {
                if (useSpot) {
                    obsSamples[i] = spot;
                }
                else {
                    obsSamples[i] = asset->fwdValue(obsDates[i]);
                }
            } else if (obsDates[i].isGreater(newDate)) {
                pastRollDate = true;
            }
            ++i;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void VarSwapUtilities::populateSamples(const CAsset*               asset,
                                       const ObservationSource*    assetHistorySource,
                                       const DateTime&             valueDate,
                                       const DateTimeArray&        obsDates,
                                       const ObservationTypeArray& obsTypes,
                                       DoubleArray&                obsSamples,
                                       SamplingConventionSP        sampleRule) {
    if (!sampleRule) {
        sampleRule.reset(new UnadjustedConvention());
    }

    static const string routine = "VarSwapUtilities::populateSamples";
    try {
        // Fill in historic samples
        FixingTypeSP fixType(new AssetFixType(asset->getTrueName()));
        
        // Grab strictly earlier samples
        // Will not try to grab sample if pricing at time of sample (except if pricing on expiry - see below).
        // Instead rely on products "roll to now" Theta call for populating "now" sample
        for (int i = 0; i < obsDates.size() && obsDates[i] < valueDate; i++) {
            obsSamples[i] = asset->pastValue(obsDates[i],
                                              obsTypes[i].get(),
                                              assetHistorySource,
                                              fixType.get(),
                                              0,
                                              sampleRule.get());
            // Temp validation
            if (Maths::isZero(obsSamples[i])) {
                throw ModelException(routine, "Invalid sample (0.0) found for date " +
                                              obsDates[i].toString() + " for asset " +
                                              asset->getName());
            }
        }

        // Grab today's sample if on expiry. 
        // Will fail if not in asset history
        if (valueDate == obsDates.back()) {
            obsSamples.back() = asset->pastValue(obsDates.back(),
                                              obsTypes.back().get(),
                                              assetHistorySource,
                                              fixType.get(),
                                              0,
                                              sampleRule.get());
            // Temp validation
            if (Maths::isZero(obsSamples.back())) {
                throw ModelException(routine, "Invalid sample (0.0) found for date " +
                                              obsDates.back().toString() + " for asset " +
                                              asset->getName());
            }
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

double VarSwapUtilities::futureDiscreteDivsAdjustment(const CAsset* asset,
                                                      const DateTime& valueDate,
                                                      const DateTime& startDate,
                                                      const DateTime& endDate,
                                                      int observationsPerYear,
                                                      int numSamplesInSwap)
{
    static const string method = "VarSwapUtilities::futureDiscreteDivsAdjustment";

    try {
        // extract all divs
        VolVarShell::DivAdjuster divAdjuster;

        // Then initialise our DividendCollector
        DividendCollector collector
            (asset,
            &divAdjuster, // optional
            valueDate,
            startDate < valueDate ? valueDate : startDate,
            endDate);

        DividendListSP divs = collector.getDividends();

        const DividendArray& divArray = divs->getArray();

        double sum = 0.0;
        for (int i = 0; i < divArray.size(); i++){
            // now compute the yield dividend, or the equivalent yield dividend if AMOUNT
            double amount = divArray[i].getDivAmount();
            if(divArray[i].getDivType() == Dividend::PERCENT)
            {
                // do nothing.
            }
            else if(divArray[i].getDivType() == Dividend::AMOUNT)
            {
                DateTime    exDate = divArray[i].getExDate();
                DateTime    prevDate = DateTime(exDate.getDate(), DateTime::BEFORE_EX_DIV_TIME);

                try {
                    amount = amount / (asset->fwdValue(exDate) + amount);
                } catch (exception& e) {
                    if(Asset::IQuanto::TYPE->isInstance(*asset)) {
                        const Asset::IQuanto* ptr = dynamic_cast<const Asset::IQuanto*>(asset);
                        double fwd = ptr->unadjustedFwdValue(exDate);
                        double prevFwd = ptr->unadjustedFwdValue(prevDate);
                        amount = (prevFwd - fwd)/prevFwd;
                    } else {
                        throw ModelException(e, method);
                    }
                }
            }
            else
            {   // continuous divs get 0 contribution
                amount = 0.0;
            }
            if (!Maths::isPositive(1.0 - amount))
            {
                throw ModelException(method, "Div yield cannot be more than 100%");
            }
            double lg = log(1.0 - amount);
            sum += lg * lg;
        }

        return sum * (double)observationsPerYear/(double)(numSamplesInSwap - 1);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

                                          
DRLIB_END_NAMESPACE
