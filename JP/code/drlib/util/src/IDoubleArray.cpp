//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IDoubleArray.cpp
//
//   Description : Provides a subset of array of doubles behaviour, plus embellishments
//
//   Date        : Nov 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_IDOUBLEARRAY_CPP
#include "edginc/Maths.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

/* XXX Not sure how to provide central definition of these */
#define PERF_TYPE_FORWARD   'F'
#define PERF_TYPE_CALL      'C'
#define PERF_TYPE_PUT       'P'
#define PERF_TYPE_STRADDLE  'S'
#define PERF_TYPE_DIGITAL   'D'
#define PERF_TYPE_NOTHING   'N'
#define PERF_TYPE_FORWARD_STR   "F"
#define PERF_TYPE_CALL_STR      "C"
#define PERF_TYPE_PUT_STR       "P"
#define PERF_TYPE_STRADDLE_STR  "S"
#define PERF_TYPE_DIGITAL_STR   "D"
#define PERF_TYPE_NOTHING_STR   "N"

void IDoubleArrayModifierMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IDoubleArrayModifierMaker, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IDoubleArrayModifierMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "IDoubleArrayModifierMaker", typeid(IDoubleArrayModifierMaker), load);

/*********************************************************************/

TrivialDoubleArray::TrivialDoubleArray() : CObject(TYPE), myArray(0.0) {};
TrivialDoubleArray::TrivialDoubleArray(double init): CObject(TYPE), myArray(init) {}
double& TrivialDoubleArray::operator[](const int index) {
    return myArray;
}
double& TrivialDoubleArray::operator() () {
    return myArray;
}
int TrivialDoubleArray::size() const {
    return 1;
}
int TrivialDoubleArray::begin() const {
    return 0;
}
int TrivialDoubleArray::end() const {
    return 1;
}

SimpleDoubleArray::SimpleDoubleArray(int size, double init): CObject(TYPE), myArray(size, init) {};
double& SimpleDoubleArray::operator[](const int index) {
    return myArray[index];
}
int SimpleDoubleArray::size() const {
    return myArray.size();
}
int SimpleDoubleArray::begin() const {
    return 0;
}
int SimpleDoubleArray::end() const {
    return myArray.size();
}

/*********************************************************************/

class PerfTypeNothingMaker::PerfTypeNothing: virtual public IDoubleArrayModifier {
public:
    PerfTypeNothing(PerfTypeNothingMaker* maker,
                    IDoubleArray*         p): 
    p(p) {};
    
    virtual void apply() {
        // Do nothing
    }

    virtual void apply(int index) {
        // Do nothing
    }
    
private:
    IDoubleArray*         p;
};

IDoubleArrayModifier* PerfTypeNothingMaker::getModifier(IDoubleArray* p){
    return new PerfTypeNothingMaker::PerfTypeNothing(this, p);
}

double PerfTypeNothingMaker::getInterpLevel(const int index) const {
    // Use ATM vol if LN is requested
    return 1.0;
}

PerfTypeNothingMaker::PerfTypeNothingMaker(): CObject(TYPE) {}

class PerfTypeNothingMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeNothingMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeNothingMaker);
        //clazz->setPublic(); // will always be built through PerfTypeSimpleMakerWrapper
    }
     
    static IObject* defaultPerfTypeNothingMaker(){
        return new PerfTypeNothingMaker();
    }
};

CClassConstSP const PerfTypeNothingMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeNothingMaker", 
                                typeid(PerfTypeNothingMaker), PerfTypeNothingMakerHelper::load);

/*********************************************************************/
/*********************************************************************/

class PerfTypeSimpleMaker::PerfTypeSimple : virtual public IDoubleArrayModifier {
public:
    PerfTypeSimple(PerfTypeSimpleMaker* maker,
                   IDoubleArray*         p):
        strike(maker->strike),
        participation(maker->participation),
        perfType(maker->perfType[0]),
        p(p) {};
    
    // modify the array elements in place
    virtual void apply() {
        for(int i=p->begin(); i<p->end(); i++) {
            double& a = (*p)[i];
            double fwd = a - strike;
            double perf = 0.0; // init to keep compiler quiet
            switch (perfType) {
            case PERF_TYPE_FORWARD:
                perf = fwd;
                break;
            case PERF_TYPE_CALL:
                perf = Maths::max(fwd, 0.0);
                break;
            case PERF_TYPE_PUT:
                perf = Maths::max(-fwd, 0.0);
                break;
            case PERF_TYPE_STRADDLE:
                perf = fabs(fwd);
                break;
            case PERF_TYPE_DIGITAL:
                perf = fwd>0.0?1.0:0.0;
                break;
            }
            a = participation*perf;
        }
    }

    virtual void apply(int index) {
        double& a = (*p)[index];
        double fwd = a - strike;
        double perf = 0.0; // init to keep compiler quiet
        switch (perfType) {
        case PERF_TYPE_FORWARD:
            perf = fwd;
            break;
        case PERF_TYPE_CALL:
            perf = Maths::max(fwd, 0.0);
            break;
        case PERF_TYPE_PUT:
            perf = Maths::max(-fwd, 0.0);
            break;
        case PERF_TYPE_STRADDLE:
            perf = fabs(fwd);
            break;
        case PERF_TYPE_DIGITAL:
            perf = fwd>0.0?1.0:0.0;
            break;
        }
        a = participation*perf;
    }
    
private:
    // try caching values from maker, instead of maker itself
    double                strike;
    double                participation;
    char                  perfType;
    IDoubleArray*         p;
};

IDoubleArrayModifier* PerfTypeSimpleMaker::getModifier(IDoubleArray* p){
    // XXX Could do the switch on perfType here, building one of 
    // a number of small classes such as PerfTypeSimpleForward, or
    // PerfTypeSimpleCall
    return new PerfTypeSimpleMaker::PerfTypeSimple(this, p);
}

double PerfTypeSimpleMaker::getInterpLevel(const int index) const {
    return strike;
}

// validation
void PerfTypeSimpleMaker::validatePop2Object(){
    static const string routine = "PerfTypeSimpleMaker::validatePop2Object";
    try{
        if (perfType.empty()){
            throw ModelException(routine, "Blank performance type specified -"
                                 " must supply one of F, C, P, S");
        }
        if (perfType.length()>1) {
            throw ModelException(routine,
                                 "Performance type must be a single character : "+perfType);
        }
        switch (perfType[0]){
        case PERF_TYPE_FORWARD:
        case PERF_TYPE_CALL:
        case PERF_TYPE_PUT:
        case PERF_TYPE_STRADDLE:
        case PERF_TYPE_DIGITAL:
            break;
        default:
            throw ModelException(routine,
                                 "Unrecognised performance type "+perfType+
                                 ". Must be F, C, P, D or S");
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

PerfTypeSimpleMaker::PerfTypeSimpleMaker(): CObject(TYPE), strike(0.0), participation(0.0){}

PerfTypeSimpleMaker::PerfTypeSimpleMaker(string perfType,
                                         double strike,
                                         double participation): 
    CObject(TYPE), perfType(perfType), strike(strike), participation(participation){
    
    validatePop2Object();
}

class PerfTypeSimpleMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeSimpleMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeSimpleMaker);
        FIELD(perfType, "Performance Type");
        FIELD(strike,   "Strike");
        FIELD(participation, "Participation");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     

    static IObject* defaultPerfTypeSimpleMaker(){
        return new PerfTypeSimpleMaker();
    }
};

CClassConstSP const PerfTypeSimpleMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeSimpleMaker", 
                                typeid(PerfTypeSimpleMaker), PerfTypeSimpleMakerHelper::load);

/*********************************************************************/

class PerfTypeSimpleSpreadMaker::PerfTypeSimpleSpread : virtual public IDoubleArrayModifier {
public:
    PerfTypeSimpleSpread(PerfTypeSimpleSpreadMaker* maker,
                         IDoubleArray*         p):
        maker(PerfTypeSimpleSpreadMakerSP(copy(maker))), p(p) {};
    
    // modify the array elements in place
    virtual void apply() {
        for(int i=p->begin(); i<p->end(); i++) {
            double& a = (*p)[i];
            double perf = 0.0; // init to keep compiler quiet
            switch (maker->perfType[0]) {
            case PERF_TYPE_CALL:
                // Phrased this way avoids some rounding errors
                if (a > maker->hiStrike) {
                    perf = maker->hiStrike - maker->loStrike;
                } else if (a < maker->loStrike) {
                    perf = 0.0;
                } else {
                    perf = a - maker->loStrike;
                }
                break;
            case PERF_TYPE_PUT:
                if (a < maker->loStrike) {
                    perf = maker->hiStrike - maker->loStrike;
                } else if (a > maker->hiStrike) {
                    perf = 0.0;
                } else {
                    perf = maker->hiStrike - a;
                }
                break;
            }
            a = maker->participation*perf;
        }
    }
    
    virtual void apply(int index) {
        double& a = (*p)[index];
        double perf = 0.0; // init to keep compiler quiet
        switch (maker->perfType[0]) {
        case PERF_TYPE_CALL:
            if (a > maker->hiStrike) {
                perf = maker->hiStrike - maker->loStrike;
            } else if (a < maker->loStrike) {
                perf = 0.0;
            } else {
                perf = a - maker->loStrike;
            }
            break;
        case PERF_TYPE_PUT:
            if (a < maker->loStrike) {
                perf = maker->hiStrike - maker->loStrike;
            } else if (a > maker->hiStrike) {
                perf = 0.0;
            } else {
                perf = maker->hiStrike - a;
            }
            break;
        }
        a = maker->participation*perf;
    }

private:
    PerfTypeSimpleSpreadMakerSP maker;
    IDoubleArray*               p;
};

IDoubleArrayModifier* PerfTypeSimpleSpreadMaker::getModifier(IDoubleArray* p){
    return new PerfTypeSimpleSpread(this, p);
}

// This is very weak and LN should be avoided 
double PerfTypeSimpleSpreadMaker::getInterpLevel(const int index) const {
    return loStrike;
}

// validation
void PerfTypeSimpleSpreadMaker::validatePop2Object(){
    static const string routine =
        "PerfTypeSimpleSpreadMaker::validatePop2Object";
    try{
        if (perfType.empty()){
            throw ModelException(routine, "Blank performance type specified - "
                                 "must supply one of C or P");
        }
        if (perfType.length()>1) {
            throw ModelException(routine,
                                 "Performance type must be a single character : "+perfType);
        }
        switch (perfType[0]){
        case PERF_TYPE_CALL:
        case PERF_TYPE_PUT:
            break;
        default:
            throw ModelException(routine,
                                 "Unrecognised performance type "+perfType+
                                 ". Must be C or P");
        }
        if (loStrike>hiStrike) {
            throw ModelException(routine,
                                 "loStrike (" + Format::toString(loStrike) + 
                                 ") must not be > hiStrike (" 
                                 + Format::toString(hiStrike) + ")");
            
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


// for reflection
PerfTypeSimpleSpreadMaker::PerfTypeSimpleSpreadMaker(): CObject(TYPE){
    loStrike = 0.0;
    hiStrike = 0.0;
    participation= 0.0;
}

PerfTypeSimpleSpreadMaker::PerfTypeSimpleSpreadMaker(string perfType,
                                                     double loStrike,
                                                     double hiStrike,
                                                     double participation): 
    CObject(TYPE), perfType(perfType), loStrike(loStrike), 
    hiStrike(hiStrike), participation(participation){
    try{
        validatePop2Object();
    } catch (exception&){
        throw; // gcc solaris opt fix
    }
}


class PerfTypeSimpleSpreadMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeSimpleSpreadMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeSimpleSpreadMaker);
        FIELD(perfType, "Performance Type");
        FIELD(loStrike,   "Low Strike");
        FIELD(hiStrike,   "High Strike");
        FIELD(participation, "Participation");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    static IObject* defaultPerfTypeSimpleSpreadMaker(){
        return new PerfTypeSimpleSpreadMaker();
    }
};

typedef smartPtr<PerfTypeSimpleSpreadMaker> PerfTypeSimpleSpreadMakerSP;

CClassConstSP const PerfTypeSimpleSpreadMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeSimpleSpreadMaker", 
                                typeid(PerfTypeSimpleSpreadMaker), PerfTypeSimpleSpreadMakerHelper::load);

/*********************************************************************/

class PerfTypeSimpleBandedMaker::PerfTypeSimpleBanded : virtual public IDoubleArrayModifier {
public:
    PerfTypeSimpleBanded(PerfTypeSimpleBandedMaker* maker,
                         IDoubleArray*         p):
        maker(PerfTypeSimpleBandedMakerSP(copy(maker))), p(p) {};
    
    // modify the array elements in place
    virtual void apply() {
        for(int i=p->begin(); i<p->end(); i++) {
            double& a = (*p)[i];
            double fwd = Maths::max(maker->floor, Maths::min(maker->cap, a - maker->strike));
            a = maker->participation*fwd;
        }
    }
    
    virtual void apply(int index) {
        double& a = (*p)[index];
        double fwd = Maths::max(maker->floor, Maths::min(maker->cap, a - maker->strike));
        a = maker->participation*fwd;
    }

private:
    PerfTypeSimpleBandedMakerSP maker;
    IDoubleArray*               p;
};


IDoubleArrayModifier* PerfTypeSimpleBandedMaker::getModifier(IDoubleArray* p){
    return new PerfTypeSimpleBanded(this, p);
}

// This is very weak and LN should be avoided 
double PerfTypeSimpleBandedMaker::getInterpLevel(const int index) const {
    return floor;
}

// validation
void PerfTypeSimpleBandedMaker::validatePop2Object(){
    static const string routine = 
        "PerfTypeSimpleBandedMaker::validatePop2Object";
    try{
        if (Maths::isPositive(floor-cap)) {
            throw ModelException(routine,
                                 "floor (" + Format::toString(floor) + 
                                 ") must not be greater than cap (" +
                                 Format::toString(cap) + ")" );
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

// for reflection
PerfTypeSimpleBandedMaker::PerfTypeSimpleBandedMaker(): CObject(TYPE), floor(0.0), 
    strike(0.0), cap(0.0), participation(0.0){}

PerfTypeSimpleBandedMaker::PerfTypeSimpleBandedMaker(double floor,
                                                     double strike,
                                                     double cap,
                                                     double participation): 
    CObject(TYPE), floor(floor), strike(strike), 
    cap(cap), participation(participation){
    try{
        validatePop2Object();
    } catch (exception&){
        throw; // gcc solaris opt fix
    }
}

class PerfTypeSimpleBandedMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeSimpleBandedMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeSimpleBandedMaker);
        FIELD(floor,    "Floor");
        FIELD(strike,   "Strike");
        FIELD(cap,      "Cap");
        FIELD(participation, "Participation");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    static IObject* defaultPerfTypeSimpleBandedMaker(){
        return new PerfTypeSimpleBandedMaker();
    }
};

CClassConstSP const PerfTypeSimpleBandedMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeSimpleBandedMaker", 
                                typeid(PerfTypeSimpleBandedMaker), PerfTypeSimpleBandedMakerHelper::load);

/*********************************************************************/

class PerfTypeSimpleDigitalHiLoMaker::PerfTypeSimpleDigitalHiLo : virtual public IDoubleArrayModifier {
public:
    PerfTypeSimpleDigitalHiLo(PerfTypeSimpleDigitalHiLoMaker* maker,
                         IDoubleArray*         p):
        maker(PerfTypeSimpleDigitalHiLoMakerSP(copy(maker))), p(p) {};
    
    // modify the array elements in place
    virtual void apply() {
        double dK = maker->highStrike - maker->lowStrike;
        if (dK > 0.0){
            for(int i=p->begin(); i<p->end(); i++) {
                double& a = (*p)[i];
                double cpn = maker->lowSideCoupon 
                           + (maker->highSideCoupon - maker->lowSideCoupon) 
                            *(Maths::max(a - maker->lowStrike, 0.0) - Maths::max(a - maker->highStrike,0.0) )/dK;
                a = cpn;
            }
        }else{
            for(int i=p->begin(); i<p->end(); i++) {
                double& a = (*p)[i];
                double cpn = a > maker->lowStrike ? maker->highSideCoupon : maker->lowSideCoupon;
                a = cpn;
            }
        }
    }
    
    virtual void apply(int index) {
        double& a = (*p)[index];
        double dK = maker->highStrike - maker->lowStrike;
        if (dK > 0.0){
            double cpn = maker->lowSideCoupon 
                       + (maker->highSideCoupon - maker->lowSideCoupon) 
                        *(Maths::max(a - maker->lowStrike, 0.0) - Maths::max(a - maker->highStrike,0.0) )/dK;
            a = cpn;
        }else{
            double cpn = a > maker->lowStrike ? maker->highSideCoupon : maker->lowSideCoupon;
            a = cpn;
        }
    }

private:
    PerfTypeSimpleDigitalHiLoMakerSP maker;
    IDoubleArray*               p;
};


IDoubleArrayModifier* PerfTypeSimpleDigitalHiLoMaker::getModifier(IDoubleArray* p){
    return new PerfTypeSimpleDigitalHiLo(this, p);
}

// This is very weak and LN should be avoided 
double PerfTypeSimpleDigitalHiLoMaker::getInterpLevel(const int index) const {
    return lowStrike;
}

// validation
void PerfTypeSimpleDigitalHiLoMaker::validatePop2Object(){
    static const string routine = 
        "PerfTypeSimpleDigitalHiLoMaker::validatePop2Object";
    try{
        throw ModelException(routine,"DigitalHiLo is not available as Simple Performance type." );
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

// for reflection
PerfTypeSimpleDigitalHiLoMaker::PerfTypeSimpleDigitalHiLoMaker(): CObject(TYPE), lowStrike(0.0), 
    highStrike(0.0), lowSideCoupon(0.0), highSideCoupon(0.0){}

PerfTypeSimpleDigitalHiLoMaker::PerfTypeSimpleDigitalHiLoMaker(double lowSideCoupon,
                                                     double highSideCoupon,
                                                     double lowStrike,
                                                     double highStrike): 
    CObject(TYPE), lowSideCoupon(lowSideCoupon), highSideCoupon(highSideCoupon), 
    lowStrike(lowStrike), highStrike(highStrike){
    try{
        validatePop2Object();
    } catch (exception&){
        throw; // gcc solaris opt fix
    }
}

class PerfTypeSimpleDigitalHiLoMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeSimpleDigitalHiLoMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeSimpleDigitalHiLoMaker);
        FIELD(lowSideCoupon,    "Lower Side Coupon");
        FIELD(highSideCoupon,   "higher Side Coupon");
        FIELD(lowStrike,      "lowStrike");
        FIELD(highStrike, "highStrike");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    static IObject* defaultPerfTypeSimpleDigitalHiLoMaker(){
        return new PerfTypeSimpleDigitalHiLoMaker();
    }
};

CClassConstSP const PerfTypeSimpleDigitalHiLoMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeSimpleDigitalHiLoMaker", 
                                typeid(PerfTypeSimpleDigitalHiLoMaker), PerfTypeSimpleDigitalHiLoMakerHelper::load);

/*********************************************************************/
/* This is more a "public interface" class which is capable of representing any 
   of the "Simple" ones above.
   Intended for use in Pyramid and probably spreadsheets too.
   Recognise perfTypes (plus aliases below):
   C    - PerfTypeSimple,       Call,     1 strike
   P    - PerfTypeSimple,       Put,      1 strike
   F    - PerfTypeSimple,       Forward,  1 strike
   S    - PerfTypeSimple,       Straddle, 1 strike
   CS   - PerfTypeSimpleSpread, Call,     2 strikes
   PS   - PerfTypeSimpleSpread, Put,      2 strikes
   BDF  - PerfTypeSimpleBanded, Forward,  3 strikes
   DHL  - PerfTypeSimpleDigitalHiLo, no participatoin.  2 strike, 2 coupon. (n/a as Simple)
   NONE - PerfTypeNothing.
*/
class PerfTypeSimpleMakerWrapper : public CObject,
                                   virtual public IDoubleArrayModifierMaker {
public: // how can I have this protected or private?
    string      genPerfType;      // C/P/F/S, or CS/PS, or BDF, or NONE
    DoubleArray strikesPct;    // 1, 2 or 3 strikes
    double      participation;

private:
    IDoubleArrayModifierMakerSP realMaker;
    bool                        isValid;
    string                      err; // if !isValid

    enum perfCategory {
        Unknown,
        Simple,
        Spread,
        Banded,
        Nothing
    };
    
    static void findPerf(string             genPerfType,
                         string&            perfType,
                         enum perfCategory& id) {
        if (genPerfType=="F" || genPerfType=="Forward") {
            perfType = PERF_TYPE_FORWARD_STR;
            id = Simple;
        } else if (genPerfType=="C" || genPerfType=="Call") {
            perfType = PERF_TYPE_CALL_STR;
            id = Simple;
        } else if (genPerfType=="P" || genPerfType=="Put") {
            perfType = PERF_TYPE_PUT_STR;
            id = Simple;
        } else if (genPerfType=="S" || genPerfType=="Straddle") {
            perfType = PERF_TYPE_STRADDLE_STR;
            id = Simple;
        } else if (genPerfType=="D" || genPerfType=="Digital") {
            perfType = PERF_TYPE_DIGITAL_STR;
            id = Simple;
        } else if (genPerfType=="CS" || genPerfType=="Call Spread") {
            perfType = PERF_TYPE_CALL_STR;
            id = Spread;
        } else if (genPerfType=="PS" || genPerfType=="Put Spread") {
            perfType = PERF_TYPE_PUT_STR;
            id = Spread;
        } else if (genPerfType=="BDF" || genPerfType=="Banded Forward") {
            perfType = PERF_TYPE_FORWARD_STR;
            id = Banded;
        } else if (genPerfType=="NONE" || genPerfType=="Nothing" || genPerfType=="") {
            perfType = PERF_TYPE_NOTHING_STR;
            id = Nothing;
        } else {
            throw ModelException("PerfTypeSimpleMakerWrapper::findPerf",
                                 "Unrecognised Gen Perf Type : " + genPerfType);
        }
    };

public:
    static CClassConstSP const TYPE;

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p){
        if (!isValid) {
            throw ModelException("Attempt to use invalid PerfTypeSimpleMakerWrapper!\n" + err);
        }
        return realMaker->getModifier(p);
    }

    virtual double getInterpLevel(const int index) const {
        if (!isValid) {
            throw ModelException("Attempt to use invalid PerfTypeSimpleMakerWrapper!\n" + err);
        }
        return strikesPct.size()>0?strikesPct[0]:0.0;
    }

    // validation - deferred since construction via IMS requires it
    void validatePop2Object(){
        static const string routine = 
            "PerfTypeSimpleMakerWrapper::validatePop2Object";
        // put maker outside of try/catch block - compiler was having problems
        IDoubleArrayModifierMaker* maker = 0;
        try{

            string            perfType;
            enum perfCategory id;
            findPerf(genPerfType, perfType, id);
            if (id==Nothing) {
                if (strikesPct.size() != 0) {
                    err = "genPerfType " + genPerfType + 
                        " requires no strikes, but given " + 
                        Format::toString(strikesPct.size());
                    return;
                }
                maker = new PerfTypeNothingMaker();
            } else if (id==Simple) {
                if (strikesPct.size() != 1) {
                    err = "genPerfType " + genPerfType + 
                        " requires single strike, but given " + 
                        Format::toString(strikesPct.size());
                    return;
                }
                maker = new PerfTypeSimpleMaker(perfType,
                                                strikesPct[0],
                                                participation);
            } else if (id==Spread) {
                if (strikesPct.size() != 2) {
                    err = "PerfType " + genPerfType + 
                        " requires 2 strikes, but given " +
                        Format::toString(strikesPct.size());
                    return;
                }
                maker = new PerfTypeSimpleSpreadMaker(perfType,
                                                      strikesPct[0],
                                                      strikesPct[1],
                                                      participation);
            } else if (id==Banded) {
                if (strikesPct.size() != 3) {
                    err = "PerfType " + genPerfType + 
                        " requires 3 strikes, but given " +
                        Format::toString(strikesPct.size());
                    return;
                }
                maker = new PerfTypeSimpleBandedMaker(strikesPct[0],
                                                      strikesPct[1],
                                                      strikesPct[2],
                                                      participation);
            } 
            isValid = true;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
        realMaker = IDoubleArrayModifierMakerSP(maker);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeSimpleMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeSimpleMakerWrapper);
        FIELD(genPerfType, "Generalised Performance Type");
        FIELD_MAKE_OPTIONAL(genPerfType);
        FIELD(strikesPct,  "Strikes");
        FIELD_MAKE_OPTIONAL(strikesPct);
        FIELD(participation, "Participation");
        FIELD_MAKE_OPTIONAL(participation);
        FIELD(realMaker, "realMaker");
        FIELD_MAKE_TRANSIENT(realMaker);
        FIELD(isValid, "isValid");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(err, "deferred error message");
        FIELD_MAKE_TRANSIENT(err);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    PerfTypeSimpleMakerWrapper(): CObject(TYPE), strikesPct(0), 
        participation(0.0), isValid(false), err("") {}
    
    static IObject* defaultPerfTypeSimpleMakerWrapper(){
        return new PerfTypeSimpleMakerWrapper();
    }
};

typedef smartPtr<PerfTypeSimpleMakerWrapper> PerfTypeSimpleMakerWrapperSP;

CClassConstSP const PerfTypeSimpleMakerWrapper::TYPE =
CClass::registerClassLoadMethod("PerfTypeSimpleMakerWrapper", 
                                typeid(PerfTypeSimpleMakerWrapper), load);


/*********************************************************************/

class PerfTypePerElementMaker::PerfTypePerElement : virtual public IDoubleArrayModifier {
public:
    PerfTypePerElement(PerfTypePerElementMaker* maker,
                       IDoubleArray*         p):
        maker(PerfTypePerElementMakerSP(copy(maker))), p(p) {};
    
    // modify the array elements in place
    virtual void apply() {
        for(int i=p->begin(); i<p->end(); i++) {
            double& a = (*p)[i];
            double fwd = a - maker->strikesPct[i];
            double perf = 0.0; // init to keep compiler quiet
            switch (maker->perfTypes[i][0]) {
            case PERF_TYPE_FORWARD:
                perf = fwd;
                break;
            case PERF_TYPE_CALL:
                perf = Maths::max(fwd, 0.0);
                break;
            case PERF_TYPE_PUT:
                perf = Maths::max(-fwd, 0.0);
                break;
            case PERF_TYPE_STRADDLE:
                perf = fabs(fwd);
                break;
            }
            a = maker->participations[i]*perf;
        }
    }

    virtual void apply(int index) {
        double& a = (*p)[index];
        double fwd = a - maker->strikesPct[index];
        double perf = 0.0; // init to keep compiler quiet
        switch (maker->perfTypes[index][0]) {
        case PERF_TYPE_FORWARD:
            perf = fwd;
            break;
        case PERF_TYPE_CALL:
            perf = Maths::max(fwd, 0.0);
            break;
        case PERF_TYPE_PUT:
            perf = Maths::max(-fwd, 0.0);
            break;
        case PERF_TYPE_STRADDLE:
            perf = fabs(fwd);
            break;
        }
        a = maker->participations[index]*perf;
    }
    
private:
    PerfTypePerElementMakerSP maker;
    IDoubleArray*             p;
};

IDoubleArrayModifier* PerfTypePerElementMaker::getModifier(IDoubleArray* p){
    if (!isValid) {
        throw ModelException("Attempt to use invalid PerfTypePerElementMaker!\n" + err);
    }
    // check dimensions 
    if (p->size() != perfTypes.size()) {
        throw ModelException("This PerfTypePerElementMaker is for arrays length " + 
                             Format::toString(perfTypes.size()) +
                             " but is being asked to act on array of length " +
                             Format::toString(p->size()));
    }
    return new PerfTypePerElement(this, p);
}

double PerfTypePerElementMaker::getInterpLevel(const int index) const {
    if (!isValid) {
        throw ModelException("Attempt to use invalid PerfTypePerElementMaker!\n" + err);
    }
    if (index < 0 || index >= strikesPct.size()) {
        throw ModelException("Not enough strikes (given " + 
                             Format::toString(strikesPct.size()) + 
                             " but need at least " +
                             Format::toString(index+1) + ")!\n");
    }
    return strikesPct[index];
}

// validation - deferred since construction via IMS requires it
void PerfTypePerElementMaker::validatePop2Object(){
    static const string routine = 
        "PerfTypePerElementMaker::validatePop2Object";
    try{
        int num = perfTypes.size();
        
        if (num<1) {
            err = "Require values for at least 1 element - none provided!";
            return;
        }
        if (num!=strikesPct.size()) {
            err = "Number of perfTypes (" + Format::toString(num) + 
                ") must equal number of strikes (" + 
                Format::toString(strikesPct.size()) + ")";
            return;
        }
        if (num!=participations.size()) {
            err = "Number of perfTypes (" + Format::toString(num) + 
                ") must equal number of participations (" + 
                Format::toString(participations.size()) + ")";
            return;
        }
        for(int i=0; i<num; i++) {
            if (perfTypes[i].empty()){
                err = Format::toString(i+1) + 
                    "'th performance type is empty - "
                    "must supply one of F, C, P, S";
                return;
            }
            if (perfTypes[i].length()>1) {
                throw ModelException(routine,
                                     "Performance type must be a single character but perf type [" +
                                     Format::toString(i+1) + "] is "+perfTypes[i]);
            }
            switch (perfTypes[i][0]){
            case PERF_TYPE_FORWARD:
            case PERF_TYPE_CALL:
            case PERF_TYPE_PUT:
            case PERF_TYPE_STRADDLE:
                break;
            default:
                err = "Unrecognised performance type [" + 
                    Format::toString(i+1) + "] = " 
                    + perfTypes[i] + ". Must be F, C, P or S";
                return;
            }
        }
        isValid = true;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

// for reflection
PerfTypePerElementMaker::PerfTypePerElementMaker(): CObject(TYPE), 
    perfTypes(0), isValid(false), err("") {}

PerfTypePerElementMaker::PerfTypePerElementMaker(StringArray perfTypes,
                                                 DoubleArray strikes,
                                                 DoubleArray participations) :
    CObject(TYPE), perfTypes(perfTypes), strikesPct(strikes), 
    participations(participations), isValid(false), err("") {
    validatePop2Object();
}

class PerfTypePerElementMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypePerElementMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypePerElementMaker);
        FIELD(perfTypes, "Performance Type");
        FIELD(strikesPct,   "Strike");
        FIELD(participations, "Participation");
        FIELD(isValid, "isValid");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(err, "deferred error message");
        FIELD_MAKE_TRANSIENT(err);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    static IObject* defaultPerfTypePerElementMaker(){
        return new PerfTypePerElementMaker();
    }
};

typedef smartPtr<PerfTypePerElementMaker> PerfTypePerElementMakerSP;

CClassConstSP const PerfTypePerElementMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypePerElementMaker", 
                                typeid(PerfTypePerElementMaker), PerfTypePerElementMakerHelper::load);


/*********************************************************************/

class PerfTypeBandedPerElementMaker::PerfTypeBandedPerElement : virtual public IDoubleArrayModifier {
public:
    PerfTypeBandedPerElement(PerfTypeBandedPerElementMaker* maker,
                             IDoubleArray*                  p):
        maker(PerfTypeBandedPerElementMakerSP(copy(maker))), p(p) {};
    
    // modify the array elements in place
    virtual void apply() {
        for(int i=p->begin(); i<p->end(); i++) {
            double& a = (*p)[i];
            double fwd = Maths::max(maker->floors[i], 
                                    Maths::min(maker->caps[i], 
                                               a - maker->strikesPct[i]));
            a = maker->participations[i]*fwd;
        }
    }

    virtual void apply(int index) {
        double& a = (*p)[index];
        double fwd = Maths::max(maker->floors[index], 
                                Maths::min(maker->caps[index], 
                                           a - maker->strikesPct[index]));
        a = maker->participations[index]*fwd;
    }
    
private:
    PerfTypeBandedPerElementMakerSP   maker;
    IDoubleArray*                     p;
};

IDoubleArrayModifier* PerfTypeBandedPerElementMaker::getModifier(IDoubleArray* p){
    if (!isValid) {
        throw ModelException("Attempt to use invalid PerfTypeBandedPerElementMaker!\n" + err);
    }
    // check dimensions 
    if (p->size() != strikesPct.size()) {
        throw ModelException("This PerfTypePerElementMaker is for arrays length " + 
                             Format::toString(strikesPct.size()) +
                             " but is being asked to act on array of length " +
                             Format::toString(p->size()));
    }
    return new PerfTypeBandedPerElement(this, p);
}

double PerfTypeBandedPerElementMaker::getInterpLevel(const int index) const {
    if (!isValid) {
        throw ModelException("Attempt to use invalid PerfTypeBandedPerElementMaker!\n" + err);
    }
    if (index < 0 || index >= floors.size()) {
        throw ModelException("Not enough floor levels (given " + 
                             Format::toString(floors.size()) + 
                             " but need at least " +
                             Format::toString(index+1) + ")!\n");
    }
    return floors[index];
}

// validation - deferred since construction via IMS requires it
void PerfTypeBandedPerElementMaker::validatePop2Object(){
    static const string routine = 
        "PerfTypeBandedPerElementMaker::validatePop2Object";
    try{
        int num = floors.size();
        if (num<1) {
            err = "Require values for at least 1 element - none provided!";
            return;
        }
        if (num!=strikesPct.size()) {
            err = "Number of floors (" + Format::toString(num) + 
                ") must equal number of strikes (" + 
                Format::toString(strikesPct.size()) + ")";
            return;
        }
        if (num!=caps.size()) {
            err = "Number of floors (" + Format::toString(num) + 
                ") must equal number of caps (" + 
                Format::toString(caps.size()) + ")";
            return;
        }
        if (num!=participations.size()) {
            err = "Number of floors (" + Format::toString(num) + 
                ") must equal number of participations (" + 
                Format::toString(participations.size()) + ")";
            return;
        }
        for(int i=0; i<num; i++) {
            if (Maths::isPositive(floors[i]-caps[i])) {
                err = "floor #" + Format::toString(i+1) + " = " + 
                    Format::toString(floors[i]) + 
                    " must not be greater than corresponding cap (" +
                    Format::toString(caps[i]) + ")";
                return;
            }
        }
    }catch (exception& e){
        throw ModelException(e, routine);
    }
    
    isValid = true;
}

// for reflection
PerfTypeBandedPerElementMaker::PerfTypeBandedPerElementMaker(): CObject(TYPE), 
    isValid(false), err("") {}


PerfTypeBandedPerElementMaker::PerfTypeBandedPerElementMaker(DoubleArray  floors,
                                                             DoubleArray  strikes,
                                                             DoubleArray  caps,
                                                             DoubleArray  participations):
    CObject(TYPE), floors(floors), strikesPct(strikes),
    caps(caps), participations(participations),
    isValid(false), err("") {
    validatePop2Object();
}

class PerfTypeBandedPerElementMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeBandedPerElementMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeBandedPerElementMaker);
        FIELD(floors,       "Floors");
        FIELD_MAKE_OPTIONAL(floors);
        FIELD(strikesPct,   "Strike");
        FIELD_MAKE_OPTIONAL(strikesPct);
        FIELD(caps,         "Caps");
        FIELD_MAKE_OPTIONAL(caps);
        FIELD(participations, "Participation");
        FIELD_MAKE_OPTIONAL(participations);
        FIELD(isValid, "isValid");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(err, "deferred error message");
        FIELD_MAKE_TRANSIENT(err);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    static IObject* defaultPerfTypeBandedPerElementMaker(){
        return new PerfTypeBandedPerElementMaker();
    }
};

CClassConstSP const PerfTypeBandedPerElementMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeBandedPerElementMaker", 
                                typeid(PerfTypeBandedPerElementMaker), PerfTypeBandedPerElementMakerHelper::load);

/*********************************************************************/

class PerfTypeDigHiLoPerElementMaker::PerfTypeDigitalHiLoPerElement : virtual public IDoubleArrayModifier {
public:
    PerfTypeDigitalHiLoPerElement(PerfTypeDigHiLoPerElementMaker* maker,
                             IDoubleArray*                  p):
        maker(PerfTypeDigHiLoPerElementMakerSP(copy(maker))), p(p) {};
    
    // modify the array elements in place
    virtual void apply() {

        for(int i=p->begin(); i<p->end(); i++) {
            double dK = maker->highStrikes[i] - maker->lowStrikes[i];
            if (dK > 0.0){
                double& a = (*p)[i];
                double cpn = maker->lowSideCoupons[i] 
                           + (maker->highSideCoupons[i] - maker->lowSideCoupons[i]) 
                            *(Maths::max(a - maker->lowStrikes[i], 0.0) - Maths::max(a - maker->highStrikes[i],0.0) )/dK;
                a = cpn;
            }else{
                double& a = (*p)[i];
                double cpn = a > maker->lowStrikes[i] ? maker->highSideCoupons[i] : maker->lowSideCoupons[i];
                a = cpn;
            }
        }
    }

    virtual void apply(int index) {
        double& a = (*p)[index];
        double dK = maker->highStrikes[index] - maker->lowStrikes[index];
        double cpn;
        if (dK > 0.0){
            cpn = maker->lowSideCoupons[index] 
                + (maker->highSideCoupons[index] - maker->lowSideCoupons[index]) 
                 *(Maths::max(a - maker->lowStrikes[index], 0.0) - Maths::max(a - maker->highStrikes[index],0.0) )/dK;
        }else{
            cpn = a > maker->lowStrikes[index] ? maker->highSideCoupons[index] : maker->lowSideCoupons[index];
        }
        a = cpn;
    }
    
private:
    PerfTypeDigHiLoPerElementMakerSP   maker;
    IDoubleArray*                     p;
};

IDoubleArrayModifier* PerfTypeDigHiLoPerElementMaker::getModifier(IDoubleArray* p){
    if (!isValid) {
        throw ModelException("Attempt to use invalid PerfTypeDigHiLoPerElementMaker!\n" + err);
    }
    // check dimensions 
    if (p->size() != lowStrikes.size()) {
        throw ModelException("This PerfTypePerElementMaker is for arrays length " + 
                             Format::toString(lowStrikes.size()) +
                             " but is being asked to act on array of length " +
                             Format::toString(p->size()));
    }
    return new PerfTypeDigitalHiLoPerElement(this, p);
}

double PerfTypeDigHiLoPerElementMaker::getInterpLevel(const int index) const {
    if (!isValid) {
        throw ModelException("Attempt to use invalid PerfTypeDigHiLoPerElementMaker!\n" + err);
    }
    if (index < 0 || index >= lowStrikes.size()) {
        throw ModelException("Not enough lowStrike levels (given " + 
                             Format::toString(lowStrikes.size()) + 
                             " but need at least " +
                             Format::toString(index+1) + ")!\n");
    }
    return lowStrikes[index];
}

// validation - deferred since construction via IMS requires it
void PerfTypeDigHiLoPerElementMaker::validatePop2Object(){
    static const string routine = 
        "PerfTypeDigHiLoPerElementMaker::validatePop2Object";
    try{
        int num = lowStrikes.size();
        if (num<1) {
            err = "Require values for at least 1 element - none provided!";
            return;
        }
        if (num!=lowSideCoupons.size()) {
            err = "Number of lowStrikes (" + Format::toString(num) + 
                ") must equal number of lowSideCoupons (" + 
                Format::toString(lowSideCoupons.size()) + ")";
            return;
        }
        if (num!=highSideCoupons.size()) {
            err = "Number of lowStrikes (" + Format::toString(num) + 
                ") must equal number of highSideCoupons (" + 
                Format::toString(highSideCoupons.size()) + ")";
            return;
        }
        if (num!=highStrikes.size()) {
            err = "Number of lowStrikes (" + Format::toString(num) + 
                ") must equal number of highStrikes (" + 
                Format::toString(highStrikes.size()) + ")";
            return;
        }
        for(int i=0; i<num; i++) {
            if (Maths::isPositive(lowStrikes[i]-highStrikes[i])) {
                err = "lowStrikes #" + Format::toString(i+1) + " = " + 
                    Format::toString(lowStrikes[i]) + 
                    " must not be greater than corresponding highStrikes (" +
                    Format::toString(highStrikes[i]) + ")";
                return;
            }
        }
    }catch (exception& e){
        throw ModelException(e, routine);
    }
    
    isValid = true;
}

// for reflection
PerfTypeDigHiLoPerElementMaker::PerfTypeDigHiLoPerElementMaker(): CObject(TYPE), 
    isValid(false), err("") {}


PerfTypeDigHiLoPerElementMaker::PerfTypeDigHiLoPerElementMaker(DoubleArray  lowSideCoupons,
                                                                          DoubleArray  highSideCoupons,
                                                                          DoubleArray  lowStrikes,
                                                                          DoubleArray  highStrikes):
    CObject(TYPE), lowSideCoupons(lowSideCoupons), highSideCoupons(highSideCoupons),
    lowStrikes(lowStrikes), highStrikes(highStrikes),
    isValid(false), err("") {
    validatePop2Object();
}

class PerfTypeDigHiLoPerElementMakerHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeDigHiLoPerElementMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeDigHiLoPerElementMaker);
        FIELD(lowSideCoupons,       "lowSideCoupons");
        FIELD_MAKE_OPTIONAL(lowSideCoupons);
        FIELD(highSideCoupons,   "highSideCoupons");
        FIELD_MAKE_OPTIONAL(highSideCoupons);
        FIELD(lowStrikes,         "lowStrikes");
        FIELD_MAKE_OPTIONAL(lowStrikes);
        FIELD(highStrikes, "highStrikes");
        FIELD_MAKE_OPTIONAL(highStrikes);
        FIELD(isValid, "isValid");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(err, "deferred error message");
        FIELD_MAKE_TRANSIENT(err);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    static IObject* defaultPerfTypeDigHiLoPerElementMaker(){
        return new PerfTypeDigHiLoPerElementMaker();
    }
};

CClassConstSP const PerfTypeDigHiLoPerElementMaker::TYPE =
CClass::registerClassLoadMethod("PerfTypeDigHiLoPerElementMaker", 
                                typeid(PerfTypeDigHiLoPerElementMaker), PerfTypeDigHiLoPerElementMakerHelper::load);

/*********************************************************************/

/* The ultimate wrapping of IDoubleArrayModifierMaker's mainly for use in Pyramid 
 */
#define MAKER_TYPE_SIMPLE          "Simple"
#define MAKER_TYPE_PER_ELT         "PerElement"
#define MAKER_TYPE_BANDED_PER_ELT  "BandedPerElement"
#define MAKER_TYPE_DHL_PER_ELT     "DigitalHiLoElement"
class PerfTypeMakerWrapper : public CObject,
                                   virtual public IDoubleArrayModifierMaker {
public: // how can I have this protected or private?
    string                                  perfMakerType;     // Simple or PerElement
    PerfTypeSimpleMakerWrapperSP            simpleMaker;
    PerfTypePerElementMakerSP               perElementMaker;
    PerfTypeBandedPerElementMakerSP         bandedPerElementMaker;
    PerfTypeDigHiLoPerElementMakerSP    digitalHiLoPerElementMaker;

private:
    IDoubleArrayModifierMakerSP realMaker;

public:
    static CClassConstSP const TYPE;

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p){
        return realMaker->getModifier(p);
    }

    virtual double getInterpLevel(const int index) const {
        return realMaker->getInterpLevel(index);
    }

    // validation
    void validatePop2Object(){
        static const string routine = 
            "PerfTypeMakerWrapper::validatePop2Object";
        try{
            if (perfMakerType.empty()){
                throw ModelException(routine, 
                                     "Blank Performance Maker specified!");
            }
            if (perfMakerType==MAKER_TYPE_SIMPLE) {
                if (simpleMaker.get()) {
                    realMaker = simpleMaker;
                } else {
                    throw ModelException(routine, "Expected simpleMaker "
                                         "but none supplied!");
                }
            } else if (perfMakerType==MAKER_TYPE_PER_ELT) {
                if (perElementMaker.get()) {
                    realMaker = perElementMaker;
                } else {
                    throw ModelException(routine, "Expected perElementMaker "
                                         "but none supplied!");
                }
            } else if (perfMakerType==MAKER_TYPE_BANDED_PER_ELT) {
                if (bandedPerElementMaker.get()) {
                    realMaker = bandedPerElementMaker;
                } else {
                    throw ModelException(routine, "Expected bandedPer"
                                         "ElementMaker but none supplied!");
                }
            } else if (perfMakerType==MAKER_TYPE_DHL_PER_ELT) {
                if (digitalHiLoPerElementMaker.get()) {
                    realMaker = digitalHiLoPerElementMaker;
                } else {
                    throw ModelException(routine, "Expected digitalHiLo"
                                         "ElementMaker but none supplied!");
                }
            } else {
                throw ModelException(routine, "Unrecognised Performance Maker "
                                     + perfMakerType + 
                                     ". Expected " + MAKER_TYPE_SIMPLE + ", " +
                                     MAKER_TYPE_PER_ELT + " or " +
                                     MAKER_TYPE_BANDED_PER_ELT);
            }
        }  catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PerfTypeMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDoubleArrayModifierMaker);
        EMPTY_SHELL_METHOD(defaultPerfTypeMakerWrapper);
        FIELD(perfMakerType, "Simple or PerElement");
        FIELD(simpleMaker,  "simpleMaker");
        FIELD_MAKE_OPTIONAL(simpleMaker);
        FIELD(perElementMaker,  "perElementMaker");
        FIELD_MAKE_OPTIONAL(perElementMaker);
        FIELD(bandedPerElementMaker,  "bandedPerElementMaker");
        FIELD_MAKE_OPTIONAL(bandedPerElementMaker);
        FIELD(digitalHiLoPerElementMaker,  "digitalHiLoPerElementMaker");
        FIELD_MAKE_OPTIONAL(digitalHiLoPerElementMaker);
        FIELD(realMaker, "realMaker");
        FIELD_MAKE_TRANSIENT(realMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    PerfTypeMakerWrapper(): CObject(TYPE){}

    static IObject* defaultPerfTypeMakerWrapper(){
        return new PerfTypeMakerWrapper();
    }
};

typedef smartPtr<PerfTypeMakerWrapper> PerfTypeMakerWrapperSP;

CClassConstSP const PerfTypeMakerWrapper::TYPE =
CClass::registerClassLoadMethod("PerfTypeMakerWrapper", 
                                typeid(PerfTypeMakerWrapper), load);

/*********************************************************************/

/* Tester class */

class IDoubleArrayTester: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    IDoubleArrayModifierMakerSP  modifierMaker;
    DoubleArray                  data;

    /** set an object in a data dictionary */
    static IObjectSP apply(IDoubleArrayTester* params){

        IDoubleArrayModifierSP   modifier;
        SimpleDoubleArray        myData(params->data.size(), 0.0);
        DoubleArraySP            result;
        int                      i;

        for(i=0; i<params->data.size(); i++) {
            myData[i] = params->data[i];
        }
        modifier = IDoubleArrayModifierSP(params->modifierMaker->getModifier(&myData));
        modifier->apply();
        result = DoubleArraySP(new DoubleArray(myData.size()));
        for(i=0; i<result->size(); i++) {
            (*result)[i] = myData[i];
        }
        return result;
    }

    IDoubleArrayTester(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IDoubleArrayTester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIDoubleArrayTester);
        FIELD(modifierMaker,   "modifierMaker");
        FIELD(data, "Array of values to which to apply modifier");
        Addin::registerInstanceObjectMethod(
            "IDOUBLE_ARRAY_TESTER",
            Addin::XL_TESTS,
            "Applies a double array modifier",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)apply);
    }

    static IObject* defaultIDoubleArrayTester(){
        return new IDoubleArrayTester();
    }
};

CClassConstSP const IDoubleArrayTester::TYPE = CClass::registerClassLoadMethod(
    "IDoubleArrayTester", typeid(IDoubleArrayTester), load);


DRLIB_END_NAMESPACE

    
