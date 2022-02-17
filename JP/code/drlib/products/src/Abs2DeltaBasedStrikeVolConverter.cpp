//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Abs2DeltaBasedStrikeVolConverter.cpp
//
//   Description : Class that converts an absolute strike volbase
//                 into a delta-based volatility representation
//                 
//
//   Author      : Regis Guichard
//
//   Date        : 29 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/Enum2StringListHelper.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/Weightmatrix.hpp"
#include "edginc/BenchmarkDate.hpp"

DRLIB_BEGIN_NAMESPACE

class Abs2DeltaBasedStrikeVolConverter: public CObject{
public:
    static CClassConstSP const TYPE;

    class Output{
    public:
        enum StrikesIndices{
            TAIL_DWN = 0,
            MID_DWN,
            ATM,
            MID_UP,
            TAIL_UP,
            NB_STRIKES
        };

        // types of outputs
        struct Type{
            enum {
                VOL_MATRIX = 0,
                VHT,
                VLLP,
                NB_ENUMS
            };
            enum {
                defaultIndex = VLLP
            };
        };
        typedef Enum2StringListHelper<Type> TypeHelper;

        static string getDefaultType();

        static string getTypeList();

        static IObjectSP create(string                  type,
                                const DoubleArrayArray& strikes,
                                const DoubleArrayArray& vols);

        static bool linkIn(){
            return (VolMatrix::TYPE != 0
                    && Vht::TYPE != 0
                    && Vllp::TYPE != 0);
        }

    private:
        class VolMatrix: public CObject{
        public:
            static CClassConstSP const TYPE;
            friend class Output;

        private:
            VolMatrix(const DoubleArrayArray& strikes,
                      const DoubleArrayArray& vols):
            CObject(TYPE),
            strikes(strikes),
            vols(vols){}

            // for reflection
            VolMatrix():
            CObject(TYPE){}

            /** Invoked when Class is 'loaded' */
            static void load(CClassSP& clazz){
                clazz->setPublic(); // make visible to EAS/spreadsheet
                REGISTER(VolMatrix, clazz);
                SUPERCLASS(CObject);
                EMPTY_SHELL_METHOD(defaultCtor);
       
                FIELD(strikes, "strike matrix");
                FIELD(vols, "volatility matrix");
            }

            static IObject* defaultCtor(){
                return new VolMatrix();
            }

            // registered fields
            DoubleMatrix strikes;
            DoubleMatrix vols; 
        };

        class Vht: public CObject{
        public:
            static CClassConstSP const TYPE;
            friend class Output;

        private:
            Vht(const DoubleArrayArray& strikes,
                const DoubleArrayArray& vols):
            CObject(TYPE){
                static const double ln100To90 = -log(0.9);
                int nbDates = strikes.size();
                strikeRefs.resize(nbDates);
                atmVols.resize(nbDates);
                skews.resize(nbDates);
                convexities.resize(nbDates);
                tailSkews.resize(nbDates);
                tailConvexities.resize(nbDates);
                for(int iDate = 0; iDate < nbDates; ++iDate){
                    double strikeRef = strikes[iDate][ATM];
                    strikeRefs[iDate] = strikeRef;
                    double atmVol = vols[iDate][ATM];
                    atmVols[iDate] = atmVol;
                    double skewDwn = (atmVol - vols[iDate][MID_DWN]) 
                                     / log(strikeRef / strikes[iDate][MID_DWN])
                                     * ln100To90;
                    double skewUp = (vols[iDate][MID_UP] - atmVol) 
                                     / log(strikes[iDate][MID_UP] / strikeRef)
                                     * ln100To90;
                    skews[iDate] = 0.5 * (skewUp + skewDwn);
                    convexities[iDate] = skewUp - skewDwn;
                    double tailSkewDwn = (atmVol - vols[iDate][TAIL_DWN]) 
                                         / log(strikeRef / strikes[iDate][TAIL_DWN])
                                         * ln100To90;
                    double tailSkewUp = (vols[iDate][TAIL_UP] - atmVol) 
                                        / log(strikes[iDate][TAIL_UP] / strikeRef)
                                        * ln100To90;
                    tailSkews[iDate] = 0.5 * (tailSkewUp + tailSkewDwn);
                    tailConvexities[iDate] = tailSkewUp - tailSkewDwn;
                }
            }

            // for reflection
            Vht():
            CObject(TYPE){}

            /** Invoked when Class is 'loaded' */
            static void load(CClassSP& clazz){
                clazz->setPublic(); // make visible to EAS/spreadsheet
                REGISTER(Vht, clazz);
                SUPERCLASS(CObject);
                EMPTY_SHELL_METHOD(defaultCtor);
       
                FIELD(strikeRefs, "reference strikes");
                FIELD(atmVols, "ATM volatilities");
                FIELD(skews, "skews");
                FIELD(convexities, "convexities");
                FIELD(tailSkews, "tailSkews");
                FIELD(tailConvexities, "tailConvexities");
            }

            static IObject* defaultCtor(){
                return new Vht();
            }

            // registered fields
            DoubleArray strikeRefs;
            DoubleArray atmVols;
            DoubleArray skews;
            DoubleArray convexities;
            DoubleArray tailSkews;
            DoubleArray tailConvexities;
        };

        class Vllp: public CObject{
        public:
            static CClassConstSP const TYPE;
            friend class Output;

        private:
            Vllp(const DoubleArrayArray& strikes,
                 const DoubleArrayArray& vols):
            CObject(TYPE){
                static const double ln100To90 = -log(0.9);
                int nbDates = strikes.size();
                strikeRefs.resize(nbDates);
                atmVols.resize(nbDates);
                skews.resize(nbDates);          
                volSpreads = DoubleMatrix(nbDates, 3);
                for(int iDate = 0; iDate < nbDates; ++iDate){
                    double strikeRef = strikes[iDate][ATM];
                    strikeRefs[iDate] = strikeRef;
                    double atmVol = vols[iDate][ATM];
                    atmVols[iDate] = atmVol;
                    double beta =  (atmVol - vols[iDate][MID_DWN]) 
                                   / log(strikeRef / strikes[iDate][MID_DWN]);
                    skews[iDate] = beta * ln100To90;
                    for (int iVolSpread = 0, iStrike = 0; iVolSpread < 3; 
                         ++iVolSpread, 
                         iStrike == TAIL_DWN ? iStrike = MID_UP : ++iStrike){
                        double strike = strikes[iDate][iStrike];
                        double vol = atmVol + beta * log(strike / strikeRef);
                        volSpreads[iDate][iVolSpread] 
                            = vols[iDate][iStrike] - vol;
                    }
                }
            }

            // for reflection
            Vllp():
            CObject(TYPE){}

            /** Invoked when Class is 'loaded' */
            static void load(CClassSP& clazz){
                clazz->setPublic(); // make visible to EAS/spreadsheet
                REGISTER(Vllp, clazz);
                SUPERCLASS(CObject);
                EMPTY_SHELL_METHOD(defaultCtor);
       
                FIELD(strikeRefs, "reference strikes");
                FIELD(atmVols, "ATM volatilities");
                FIELD(skews, "skews");
                FIELD(volSpreads, "volatility spreads");
            }

            static IObject* defaultCtor(){
                return new Vllp();
            }

            // registered fields
            DoubleArray  strikeRefs;
            DoubleArray  atmVols;
            DoubleArray  skews;
            DoubleMatrix volSpreads;
        };
    };

private:
    enum StrikesIndices{
        TAIL_DWN = 0,
        MID_DWN,
        MID_UP,
        TAIL_UP,
        NB_STRIKES
    };

    CMarketDataSP           market;
    CAssetWrapper           asset;
    YieldCurveWrapper       discount;
    InstrumentSettlementSP  instSettle;
    IModelSP                model;  
    string                  volType;
    double                  deltaShiftSize;
    StringArray             BMs;
    string                  time;
    DoubleArray             deltas;
    StringArray             outputRequests;

    // transients
    int                     time2;
    DateTime                baseDate;
    DateTimeArray           dates;

    /** for reflection */
    Abs2DeltaBasedStrikeVolConverter():
    CObject(TYPE),
    instSettle(new CashSettlePeriod(0)),
    deltaShiftSize(0.5 / 100.0),    // 0.5% 
    time("EOD"),
    time2(0){
        outputRequests.resize(1);
        outputRequests[0] = Output::getDefaultType();
        deltas.resize(NB_STRIKES);
        deltas[TAIL_DWN] = 0.05; 
        deltas[MID_DWN] = 0.25; 
        deltas[MID_UP] = 0.25; 
        deltas[TAIL_UP] = 0.05; 
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Abs2DeltaBasedStrikeVolConverter, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
       
        FIELD(market, "market");        
        FIELD(asset, "asset wrapper");
        FIELD(discount, "discount");
        FIELD(instSettle, "instrument settlement");
        FIELD_MAKE_OPTIONAL(instSettle);
        FIELD(model, "model");
        FIELD_MAKE_OPTIONAL(model);
        FIELD(volType, "vol type");
        FIELD(deltaShiftSize, "delta shift size (defaulted to .05%)");
        FIELD_MAKE_OPTIONAL(deltaShiftSize);
        FIELD(BMs, "benchmarks");
        FIELD(time, "time of day");
        FIELD_MAKE_OPTIONAL(time);
        FIELD_MAKE_OPTIONAL(time);
        FIELD(deltas, "deltas");
        FIELD_MAKE_OPTIONAL(deltas);
        FIELD(outputRequests, "array of output requests ("
            + Output::getTypeList()
            + ")");
        FIELD_MAKE_OPTIONAL(outputRequests);
        // transients
        FIELD(time2, "");
        FIELD_MAKE_TRANSIENT(time2);
        FIELD(baseDate, "");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(dates, "");
        FIELD_MAKE_TRANSIENT(dates);

        Addin::registerClassObjectMethod("ABS2DELTA_BASED_STRIKE_VOL_CONVERT",
                                          Addin::MARKET,
                                          "Converts a volbase (absolute strike) into a delta-based vol representation",
                                          TYPE,
                                          true,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)output);
    }

    static IObject* defaultCtor(){
        return new Abs2DeltaBasedStrikeVolConverter();
    }

    virtual void validatePop2Object(){
        static const string method = "Abs2DeltaBasedStrikeVolConverter::validatePop2Object";
        try{
            baseDate = market->GetReferenceDate();
            time2 = DateTime::timeConvert(time);
            // check BMs and create dates
            int nbBMs = BMs.size();
            if (nbBMs == 0){
                throw ModelException(method,
                                     "The number of benchmarks should be greater than, or equal to, 1");
            }
            dates.resize(nbBMs);
            int iDate = 0;
            for (; iDate < nbBMs; ++iDate){
                MaturityPeriod period(BMs[iDate]);
                dates[iDate] = DateTime(period.toDate(baseDate).getDate(), time2);
                if (iDate > 0
                    && dates[iDate - 1] >= dates[iDate]){
                    throw ModelException(method,
                                         "The benchmarks must be increasing.\nCompare the "
                                         + Format::toString(iDate)
                                         + "-th benchmark date ("
                                         + dates[iDate - 1].toString()
                                         + ") with the "
                                         + Format::toString(iDate + 1)
                                         + "-th benchmark date ("
                                         + dates[iDate].toString()
                                         + ")");
                }
            }
            // check deltas
            if (deltas.size() != NB_STRIKES){
                throw ModelException(method,
                                     "the number of deltas should be equal to "
                                     + Format::toString(NB_STRIKES));
            }
            if (!Maths::isPositive(deltas[TAIL_DWN])
                || !Maths::isPositive(deltas[MID_DWN] - deltas[TAIL_DWN])
                || !Maths::isPositive(0.5 - deltas[MID_DWN])){
                throw ModelException(method,
                                     "lower tail delta ('tailDeltaDown'),"
                                     " lower middle delta ('midDeltaDown') should be such that "
                                     "0.0 < tailDeltaDown < midDeltaDown < 0.5;"
                                     "\ngot "
                                     + Format::toString(deltas[TAIL_DWN])
                                     + " and "
                                     + Format::toString(deltas[MID_DWN])
                                     + " respectively");
            }
            if (!Maths::isPositive(deltas[TAIL_UP])
                || !Maths::isPositive(deltas[MID_UP] - deltas[TAIL_UP])
                || !Maths::isPositive(0.5 - deltas[MID_UP])){
                throw ModelException(method,
                                     "upper tail delta ('tailDeltaUp'),"
                                     " upper middle delta ('midDeltaUp') should be such that "
                                     "0.0 < tailDeltaUp < midDeltaUp < 0.5;"
                                     "\ngot "
                                     + Format::toString(deltas[TAIL_UP])
                                     + " and "
                                     + Format::toString(deltas[MID_UP])
                                     + " respectively");
            }
            // populate model if needed
            if (!model){
                if (volType.empty()){
                    throw ModelException(method,
                                         "Need to specify a volatility type when the model is not provided");
                }
                model = IModelSP(new CClosedFormLN(volType, true /* allow neg fwd vars */));
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // populate from market, further validate and initialize
    void validate(){
        static const string method = "Abs2DeltaBasedStrikeVolConverter::validate";
        try{
            // get market
            Vanilla::DeltaImpliedStrike::Helper::getMarket(
                market.get(),
                model.get(),
                asset,
                discount,
                instSettle.get());
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    DoubleArrayArraySP calcStrikes(){
        static const string routine = "Abs2DeltaBasedStrikeVolConverter::calcStrike";
        try {
            double spot = asset->getSpot();
            double strikeAbsAcc = spot * 0.01 / 100.0;   // 0.01% of current spot
            int nbDates = dates.size();
            DoubleArrayArraySP strikes(new DoubleArrayArray(nbDates));
            for (int iDate = 0; iDate < nbDates; ++iDate){
                DoubleArray& mystrikes = (*strikes)[iDate];
                mystrikes.resize(NB_STRIKES);
                double upperStrike = spot;
                double lowerStrike = 0.9 * spot;    // 90% of current spot
                for(int iStrike = 0; iStrike < NB_STRIKES; ++iStrike){
                    bool isCall = (iStrike <= MID_DWN ? false : true);
                    double tgtDelta = (isCall ? deltas[iStrike] : -deltas[iStrike]);
                    Vanilla::DeltaImpliedStrike implier(
                        *model,
                        baseDate,
                        dates[iDate],
                        isCall,
                        asset.get(),
                        discount.get(),
                        instSettle.get(),
                        deltaShiftSize,
                        tgtDelta);
                    mystrikes[iStrike] = implier.calcStrike(
                        lowerStrike,
                        upperStrike,
                        strikeAbsAcc);
                }
            }
            return strikes;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    IObjectSP calcOutput() {
        static const string method = "Abs2DeltaBasedStrikeVolConverter::calcOutput";
        try{
            static const int nbAllStrikes = NB_STRIKES + 1;
            // get market 
            validate();
            // calc strikes
            DoubleArrayArraySP strikes(calcStrikes());
            // calc forwards
            int nbDates = dates.size();
            DoubleArray fwds(nbDates);
            asset->fwdValue(dates, fwds);
            // will need an indicative vol request to compute the vols
            SensitivityArraySP sens(new SensitivityArray(0));
            OutputRequestArraySP outReqs(new OutputRequestArray(1));
            (*outReqs)[0] = OutputRequestSP(new OutputRequest(OutputRequest::IND_VOL));
            Control ctrl(sens,
                         outReqs,
                         false,
                         "");
            // check consistency of strikes vs fwd, and
            // calculate vols
            DoubleArrayArray vols(nbDates);
            DoubleArrayArray allStrikes(nbDates);
            for (int iDate = 0; iDate < nbDates; ++iDate){
                double fwd = fwds[iDate];
                if (!Maths::isNegative((*strikes)[iDate][MID_DWN] - fwd)
                    || !Maths::isPositive((*strikes)[iDate][MID_UP] - fwd)){
                    throw ModelException(method,
                                         "the lower middle strike ("
                                         + Format::toString((*strikes)[iDate][MID_DWN])
                                         + ") at the "
                                         + Format::toString(iDate + 1)
                                         + "-th benchmark date ("
                                         + dates[iDate].toString()
                                         + ") should be smaller than the forward ("
                                         + Format::toString(fwd)
                                         + "),\n"
                                         + "which itself should be smaller than the upper middle strike ("
                                         + Format::toString((*strikes)[iDate][MID_UP])
                                         + ")");
                }
                // insert fwd inside strike array
                allStrikes[iDate].resize(nbAllStrikes);
                int iStrike = 0;
                for (; iStrike <= MID_DWN; ++iStrike){
                    allStrikes[iDate][iStrike] = (*strikes)[iDate][iStrike];
                }
                int iAllStrike = iStrike;
                allStrikes[iDate][iAllStrike] = fwd;
                for (++iAllStrike; iStrike < NB_STRIKES; ++iStrike, ++iAllStrike){
                    allStrikes[iDate][iAllStrike] = (*strikes)[iDate][iStrike];
                }
                // calc the vols (go via the VanillaGrid route because 
                // not all volbases support IVolatilityBS)
               
                // create a weightmatrix for the vanillagrid
                DateTimeArraySP maturitiesSP(new DateTimeArray(1, dates[iDate]));
                DoubleArraySP strikeSP(copy(&(allStrikes[iDate])));
                CDoubleMatrixSP weightsSP(new CDoubleMatrix(DoubleArray(nbAllStrikes, 1.0)));
                WeightMatrixSP weightMatrix(new WeightMatrix("", maturitiesSP, strikeSP, weightsSP, 
                                WeightMatrix::ABSOLUTE, WeightMatrix::CALL,baseDate));
                WeightMatrixWrapper weightMatrixWrapper(weightMatrix);

                VanillaGrid vanillaGrid(baseDate,
                                        asset,
                                        discount,
                                        weightMatrixWrapper,
                                        instSettle);
                CResultsSP results(model->Run(&vanillaGrid, &ctrl));
                IObjectConstSP resobj(results->retrieveRequestResult(OutputRequest::IND_VOL));
                VanillaGrid::OutputConstSP res(VanillaGrid::OutputConstSP::dynamicCast(resobj));                
                vols[iDate].resize(nbAllStrikes);
                for (iAllStrike = 0; iAllStrike < nbAllStrikes; ++iAllStrike){                    
                    vols[iDate][iAllStrike] = res->getValue(0, iAllStrike);
                }
            }
            // reformat the strikes and vols as needed
            ObjectArraySP res(new ObjectArray(0));
            int nbOutReqs = outputRequests.size();
            res->reserve(nbOutReqs);
            for (int iOutReq = 0; iOutReq < nbOutReqs; ++ iOutReq){
                res->push_back(
                    Output::create(outputRequests[iOutReq],
                                   allStrikes,
                                   vols));
            }
            return res;
        } 
        catch (exception& e){
            throw ModelException(e, method);
        }        
    }

    static IObjectSP output(Abs2DeltaBasedStrikeVolConverter* params) {
        return params->calcOutput();
    }
};

CClassConstSP const Abs2DeltaBasedStrikeVolConverter::TYPE =
CClass::registerClassLoadMethod("Abs2DeltaBasedStrikeVolConverter", typeid(Abs2DeltaBasedStrikeVolConverter), load);


template<> string nameForType<Abs2DeltaBasedStrikeVolConverter::Output::Type>(
        Abs2DeltaBasedStrikeVolConverter::Output::Type*){
    return "Abs2DeltaBasedStrikeVolConverter::Output::Type";
}
template<> string Abs2DeltaBasedStrikeVolConverter::Output::TypeHelper::names[
        Abs2DeltaBasedStrikeVolConverter::Output::TypeHelper::EnumList::NB_ENUMS] = {
    "VOL_MATRIX",
    "VHT",
    "VLLP"
};

string Abs2DeltaBasedStrikeVolConverter::Output::getDefaultType(){
    return TypeHelper::getDefaultName();
}

string Abs2DeltaBasedStrikeVolConverter::Output::getTypeList(){
    return TypeHelper::getNameList();
}

IObjectSP Abs2DeltaBasedStrikeVolConverter::Output::create(
        string                  type,
        const DoubleArrayArray& strikes,
        const DoubleArrayArray& vols){
    static const string method = "Abs2DeltaBasedStrikeVolConverter::Output::create";
    try{
        int typeIdx = TypeHelper::getIndex(type);
        if (typeIdx == TypeHelper::EnumList::VOL_MATRIX){
            return IObjectSP(new VolMatrix(strikes, vols));
        }
        else if (typeIdx == TypeHelper::EnumList::VHT){
            return IObjectSP(new Vht(strikes, vols));
        }
        else if (typeIdx == TypeHelper::EnumList::VLLP){
            return IObjectSP(new Vllp(strikes, vols));
        }
        else{
            throw ModelException(method,
                                 "internal error: unexpected output request of type "
                                 + type);
        }
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }        
}

CClassConstSP const Abs2DeltaBasedStrikeVolConverter::Output::VolMatrix::TYPE =
CClass::registerClassLoadMethod("Abs2DeltaBasedStrikeVolConverter_Output_VolMatrix", 
typeid(Abs2DeltaBasedStrikeVolConverter::Output::VolMatrix), load);

CClassConstSP const Abs2DeltaBasedStrikeVolConverter::Output::Vht::TYPE =
CClass::registerClassLoadMethod("Abs2DeltaBasedStrikeVolConverter_Output_Vht", 
typeid(Abs2DeltaBasedStrikeVolConverter::Output::Vht), load);

CClassConstSP const Abs2DeltaBasedStrikeVolConverter::Output::Vllp::TYPE =
CClass::registerClassLoadMethod("Abs2DeltaBasedStrikeVolConverter_Output_Vllp", 
typeid(Abs2DeltaBasedStrikeVolConverter::Output::Vllp), load);

// external symbol to allow class to be forced to be linked in
bool Abs2DeltaBasedStrikeVolConverterLoad(){
    return (Abs2DeltaBasedStrikeVolConverter::TYPE != 0
            && Abs2DeltaBasedStrikeVolConverter::Output::linkIn());
}

DRLIB_END_NAMESPACE
