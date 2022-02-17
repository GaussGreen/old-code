//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : WeightMatrix.hpp
//
//   Description : Weight Matrix
//
//   Author      : Regis Guichard
//
//   Date        : 03 Jan 03
//
//
//----------------------------------------------------------------------------

#ifndef EDR_WEIGHTMATRIX_HPP
#define EDR_WEIGHTMATRIX_HPP
#include "edginc/Class.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

/*  DefaultWeightMatrix / Addin class */
class PRODUCTS_DLL DefaultWeightMatrix: public CObject {
public:
    static CClassConstSP const TYPE;

    class PRODUCTS_DLL DefaultRefExpiries: public CObject {
    public:
        static CClassConstSP const TYPE;

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);

        static ExpiryArray create();

    private:
        DefaultRefExpiries(const DefaultRefExpiries& rhs);
        DefaultRefExpiries& operator=(const DefaultRefExpiries& rhs);

        static IObjectSP run(DefaultRefExpiries* params);

        DefaultRefExpiries();

        static IObject* defaultCtor();
    };

    class PRODUCTS_DLL DefaultRefLowerStrikes: public CObject {
    public:
        static CClassConstSP const TYPE;

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);

        static DoubleArray create(const string& strikeUnits,
                                  const string& instType,
                                  const string& userMode);

    private:
        DefaultRefLowerStrikes(const DefaultRefLowerStrikes& rhs);
        DefaultRefLowerStrikes& operator=(const DefaultRefLowerStrikes& rhs);

        static IObjectSP run(DefaultRefLowerStrikes* params);

        DefaultRefLowerStrikes();

        static IObject* defaultCtor();

        string strikeUnits;
        string instType;
        string userMode;
    };

    class PRODUCTS_DLL DefaultRefUpperStrikes: public CObject {
    public:
        static CClassConstSP const TYPE;

        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);

        static DoubleArray create(const string& strikeUnits,
                                  const string& instType,
                                  const string& userMode);

    private:
        DefaultRefUpperStrikes(const DefaultRefUpperStrikes& rhs);
        DefaultRefUpperStrikes& operator=(const DefaultRefUpperStrikes& rhs);

        static IObjectSP run(DefaultRefUpperStrikes* params);

        DefaultRefUpperStrikes();

        static IObject* defaultCtor();

        string strikeUnits;
        string instType;
        string userMode;
    };

    static CDoubleMatrixSP create(const DateTime&      baseDate,
                                  double               spot,
                                  const DateTimeArray& maturities,
                                  const DoubleArray&   strikes,
                                  const string&        strikeUnits,
                                  const string&        instType,
                                  const string&        userMode);

    static CDoubleMatrixSP create(const DateTime&      baseDate,
                                  double               spot,
                                  const DateTimeArray& maturities,
                                  const DoubleArray&   strikes,
                                  const ExpiryArray&   refExpiries,
                                  const DoubleArray&   refLowerStrikes,
                                  const DoubleArray&   refUpperStrikes,
                                  bool                 useFwd,
                                  const DoubleArray    fwds,
                                  const string&        strikeUnits,
                                  const string&        instType,
                                  const string&        userMode);

    virtual void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static const string OLDMODE;
    static const string NEWMODE;

private:
    DefaultWeightMatrix(const DefaultWeightMatrix& rhs);
    DefaultWeightMatrix& operator=(const DefaultWeightMatrix& rhs);

    static IObjectSP run(DefaultWeightMatrix* params);

    DefaultWeightMatrix();

    static IObject* defaultCtor();
    
    CAssetWrapper asset;
    CMarketDataSP market;
    DateTimeArray maturities;
    DoubleArray   strikes;
    ExpiryArray   refExpiries;
    DoubleArray   refLowerStrikes;
    DoubleArray   refUpperStrikes;
    bool          useFwd;
    bool          isVegaWeighted;
    YieldCurveWrapper yield;

    string        userMode;
    string        strikeUnits;
    string        instType;

    //transient
    DateTime      baseDate;
    double        spot;
    DoubleArray   fwds;
};

/** Allows specification of weights independently of calibration class
 **/
class PRODUCTS_DLL WeightMatrix : public MarketObject {
public:
    static CClassConstSP const TYPE;
    friend class WeightMatrixHelper;

    WeightMatrix(const string&          name, 
                 DateTimeArraySP        dateTimeMaturities,
                 DoubleArraySP          strikes,
                 CDoubleMatrixSP        weights,
                 const string&          strikeUnits,
                 const string&          instType,
                 const DateTime&        refDate);

    ~WeightMatrix();

    /** return copy of WeightMatrix name */
    virtual string getName() const;

    /** Hardly strong motivation for class existence, but
        this is all it's wanted for... */
    DateTimeArraySP getDateTimeMaturities() const;
    DoubleArraySP getStrikes() const;
    CDoubleMatrixSP getWeights() const;
    const string& getStrikeUnits() const;
    const string& getInstType() const; 
    CDoubleMatrixSP getInstStrikesUsed() const;
    IntArrayArraySP getInstTypesUsed() const;
    void setWeights(CDoubleMatrixSP weightsSP);

    //Compute absolute strikes and types
    void computeStrikesAndTypes(const CAssetWrapper&            asset,
                                const YieldCurveWrapper&        discount,
                                InstrumentSettlementSP          instSettle);

    void validatePop2Object();

    virtual void getMarket(const IModel* model, const MarketData* market);

    static const string DEFAULT;
    static const string EMPTY;
    static const string CALL;
    static const string PUT;
    static const string OTM;
    static const string ITM;
    static const string ABSOLUTE;
    static const string SPOTMONEYNESS;
    static const string FWDMONEYNESS;
    static const string DELTA;
	static const string ATMVARIANCE;

private:
    
    WeightMatrix();
    WeightMatrix(const WeightMatrix &rhs);
    WeightMatrix& operator=(const WeightMatrix& rhs);

    // Mandatory
    string                  name;
    ExpiryArraySP           maturities;
    DoubleArraySP           strikes;
    CDoubleMatrixSP         weights;
    string                  strikeUnits;
    string                  instType;
    
    // Transient
    DateTimeArraySP         dateTimeMaturities;
    CDoubleMatrixSP         instStrikesUsed;
    IntArrayArraySP         instTypesUsed;
    DateTime                refDate;

};

typedef smartConstPtr<WeightMatrix> WeightMatrixConstSP;
typedef smartPtr<WeightMatrix> WeightMatrixSP;

// support for wrapper class
typedef MarketWrapper<WeightMatrix> WeightMatrixWrapper;

DRLIB_END_NAMESPACE

#endif
