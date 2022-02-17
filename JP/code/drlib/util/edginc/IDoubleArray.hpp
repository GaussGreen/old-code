//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IDoubleArray.hpp
//
//   Description : Provides a subset of array of doubles behaviour, plus embellishments
//
//   Date        : Nov 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IDOUBLEARRAY_HPP
#define EDR_IDOUBLEARRAY_HPP

#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL IDoubleArray {
public:
    virtual double& operator[] (const int index) = 0;
    virtual int size() const = 0;
    virtual int begin() const = 0;
    virtual int end() const = 0;

    virtual ~IDoubleArray(){}
};

typedef refCountPtr<IDoubleArray> IDoubleArraySP;

// General operations to be performed on IDoubleArrays look like this
class UTIL_DLL IDoubleArrayModifier {
public:
    virtual void apply() = 0;
    virtual void apply(int index) = 0;

    virtual ~IDoubleArrayModifier() {}
};

typedef refCountPtr<IDoubleArrayModifier> IDoubleArrayModifierSP;

// Build IDoubleArrayModifier via a Maker class so the IDoubleArray can
// be identified once and apply() needs no parameters. It's a bit like
// building the thing we need in 2 stages.
class UTIL_DLL IDoubleArrayModifierMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p) = 0;

    // For lognormal require vol interp levels
    virtual double getInterpLevel(const int index) const = 0;

    virtual ~IDoubleArrayModifierMaker() {}

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IDoubleArrayModifierMaker> IDoubleArrayModifierMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<IDoubleArrayModifierMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<IDoubleArrayModifierMaker>);
#endif

class UTIL_DLL TrivialDoubleArray : public CObject,
                           virtual public IDoubleArray {
public:
    TrivialDoubleArray();
    TrivialDoubleArray(double init);

    double& operator[](const int index);
    int size() const;
    int begin() const;
    int end() const;

    // special one for trivial case
    double& operator() ();

private:
    double myArray;
};
typedef smartPtr<TrivialDoubleArray> TrivialDoubleArraySP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<TrivialDoubleArray>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<TrivialDoubleArray>);
#endif

class UTIL_DLL SimpleDoubleArray : public CObject,
                          virtual public IDoubleArray {
public:
    SimpleDoubleArray(int size, double init);
    double& operator[](const int index);
    int size() const;
    int begin() const;
    int end() const;
private:
    DoubleArray myArray;
};

typedef smartPtr<SimpleDoubleArray> SimpleDoubleArraySP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<SimpleDoubleArray>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<SimpleDoubleArray>);
#endif

/*********************************************************************/

class UTIL_DLL PerfTypeNothingMaker : public CObject,
                            virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeNothingMakerHelper;
    
    PerfTypeNothingMaker();

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

protected:

private:
    class PerfTypeNothing;
    friend class PerfTypeNothing;
};

typedef smartPtr<PerfTypeNothingMaker> PerfTypeNothingMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeNothingMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeNothingMaker>);
#endif

/*********************************************************************/

class UTIL_DLL PerfTypeSimpleMaker : public CObject,
                            virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeSimpleMakerHelper;
    
    PerfTypeSimpleMaker(string perfType,
                        double strike,
                        double participation);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    string perfType;
    double strike;
    double participation;

private:
    // for reflection
    PerfTypeSimpleMaker();
    
    class PerfTypeSimple;
    friend class PerfTypeSimple;

};

typedef smartPtr<PerfTypeSimpleMaker> PerfTypeSimpleMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeSimpleMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeSimpleMaker>);
#endif

/*********************************************************************/

class UTIL_DLL PerfTypeSimpleSpreadMaker : public CObject,
                                  virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeSimpleSpreadMakerHelper;

    PerfTypeSimpleSpreadMaker(string perfType,
                              double loStrike,
                              double hiStrike,
                              double participation);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    string perfType;
    double loStrike;
    double hiStrike;
    double participation;

private:
    // for reflection
    PerfTypeSimpleSpreadMaker();

    class PerfTypeSimpleSpread;
    friend class PerfTypeSimpleSpread;
};

typedef smartPtr<PerfTypeSimpleSpreadMaker> PerfTypeSimpleSpreadMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeSimpleSpreadMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeSimpleSpreadMaker>);
#endif

/*********************************************************************/

// Only "F" makes sense
class UTIL_DLL PerfTypeSimpleBandedMaker : public CObject,
                                  virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeSimpleBandedMakerHelper;

    PerfTypeSimpleBandedMaker(double floor,
                              double strike,
                              double cap,
                              double participation);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    double floor;
    double strike;
    double cap;
    double participation;

private:
    // for reflection
    PerfTypeSimpleBandedMaker();

    class PerfTypeSimpleBanded;
    friend class PerfTypeSimpleBanded;
};

typedef smartPtr<PerfTypeSimpleBandedMaker> PerfTypeSimpleBandedMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeSimpleBandedMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeSimpleBandedMaker>);
#endif

/*********************************************************************/
//  Not alloed for time being, as the input variables are not same 
//  other performance.  No participation But 2 coupon levels in this payoff.
class UTIL_DLL PerfTypeSimpleDigitalHiLoMaker : public CObject,
                                  virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeSimpleDigitalHiLoMakerHelper;

    PerfTypeSimpleDigitalHiLoMaker(double lowSideCoupon,
                              double highSideCoupon,
                              double lowStrike,
                              double highStrike);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    double lowSideCoupon;
    double highSideCoupon;
    double lowStrike;
    double highStrike;

private:
    // for reflection
    PerfTypeSimpleDigitalHiLoMaker();

    class PerfTypeSimpleDigitalHiLo;
    friend class PerfTypeSimpleDigitalHiLo;
};

typedef smartPtr<PerfTypeSimpleDigitalHiLoMaker> PerfTypeSimpleDigitalHiLoMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeSimpleDigitalHiLoMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeSimpleDigitalHiLoMaker>);
#endif

/*********************************************************************/

class UTIL_DLL PerfTypePerElementMaker : public CObject,
                                virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypePerElementMakerHelper;

    PerfTypePerElementMaker(StringArray perfTypes,
                            DoubleArray strikes,
                            DoubleArray participations);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    StringArray perfTypes;
    DoubleArray strikesPct;
    DoubleArray participations;

private:
    // for reflection
    PerfTypePerElementMaker();

    bool        isValid;
    string      err;

    class PerfTypePerElement;
    friend class PerfTypePerElement;
};

typedef smartPtr<PerfTypePerElementMaker> PerfTypePerElementMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypePerElementMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypePerElementMaker>);
#endif

/*********************************************************************/

// validation - deferred since construction via IMS requires it
class UTIL_DLL PerfTypeBandedPerElementMaker : public CObject,
                                      virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeBandedPerElementMakerHelper;

    PerfTypeBandedPerElementMaker(DoubleArray  floors,
                                  DoubleArray  strikes,
                                  DoubleArray  caps,
                                  DoubleArray  participations);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    DoubleArray  floors;
    DoubleArray  strikesPct;
    DoubleArray  caps;
    DoubleArray  participations;

private:
    // for reflection
    PerfTypeBandedPerElementMaker();

    bool        isValid;
    string      err;

    class PerfTypeBandedPerElement;
    friend class PerfTypeBandedPerElement;
};

typedef smartPtr<PerfTypeBandedPerElementMaker> PerfTypeBandedPerElementMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeBandedPerElementMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeBandedPerElementMaker>);
#endif

/*********************************************************************/

// validation - deferred since construction via IMS requires it
class UTIL_DLL PerfTypeDigHiLoPerElementMaker : public CObject,
                                      virtual public IDoubleArrayModifierMaker {
public:
    static CClassConstSP const TYPE;
    friend class PerfTypeDigHiLoPerElementMakerHelper;

    PerfTypeDigHiLoPerElementMaker(DoubleArray  lowSideCoupons,
                                  DoubleArray  highSideCoupons,
                                  DoubleArray  lowStrikes,
                                  DoubleArray  highStrikes);

    virtual IDoubleArrayModifier* getModifier(IDoubleArray* p);

    virtual double getInterpLevel(const int index) const;

    virtual void validatePop2Object();

protected:
    DoubleArray  lowSideCoupons;
    DoubleArray  highSideCoupons;
    DoubleArray  lowStrikes;
    DoubleArray  highStrikes;

private:
    // for reflection
    PerfTypeDigHiLoPerElementMaker();

    bool        isValid;
    string      err;

    class PerfTypeDigitalHiLoPerElement;
    friend class PerfTypeDigitalHiLoPerElement;
};

typedef smartPtr<PerfTypeDigHiLoPerElementMaker> PerfTypeDigHiLoPerElementMakerSP;
#ifndef QLIB_IDOUBLEARRAY_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PerfTypeDigHiLoPerElementMaker>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PerfTypeDigHiLoPerElementMaker>);
#endif

DRLIB_END_NAMESPACE

#endif

