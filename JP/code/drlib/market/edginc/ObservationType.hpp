//----------------------------------------------------------------------------
//
//   Filename    : ObservationType.hpp
//
//   Description : Classes for observation types - for use with AssetHistory
//
//   Author      : Ian Stares   
//
//   Date        : January 19, 2006
//

//
//----------------------------------------------------------------------------

#ifndef OBSTYPE_HPP
#define OBSTYPE_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/Array.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class ObservationType;
typedef smartPtr<ObservationType> ObservationTypeSP;

class MARKET_DLL ObservationType :  public CObject {
public:
    static CClassConstSP const TYPE;
    
    virtual bool equals(const ObservationType& other) const;

    virtual string toString() const = 0;

    // Indicative observation time
    virtual string indicativeTime() const = 0;

    /** Factory method based on string. 
        QLib class Name is "Observation" + input string */
    static ObservationTypeSP make(const string& suffix);

protected:
    ObservationType(CClassConstSP clazz);

private:
    static void load(CClassSP& clazz);
};


// support for arrays
typedef array<ObservationTypeSP, ObservationType> ObservationTypeArray;
typedef smartPtr<ObservationTypeArray> ObservationTypeArraySP;

/*****************************************************************************/
// Concrete implementations

 /****************************************************************************/

class MARKET_DLL ObservationOpen : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationOpen();

    static IObject* defaultOpenObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationClose : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationClose();

    static IObject* defaultCloseObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationHigh : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationHigh();

    static IObject* defaultHighObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationLow : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationLow();

    static IObject* defaultLowObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationMid : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationMid();

    static IObject* defaultMidObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationExact : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

    // the exact rule overrides the equals method as it needs to check the
    // id and time match as well
    virtual bool equals(const ObservationType& other) const;

    // temporary for rates code which has some empty ObservationExact objects
    bool hasTimeAndId() const;

    // only public for temporary backwards compatibility with rates code
    ObservationExact();
private:
    static IObject* defaultExactObsType();

    string time;
    string id;
};

 /****************************************************************************/

class MARKET_DLL ObservationVWAP : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationVWAP();

    static IObject* defaultVWAPObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationOSPOption : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationOSPOption();

    static IObject* defaultOSPOptionObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationOSPOptionOpen : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationOSPOptionOpen();

    static IObject* defaultOSPOptionOpenObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationOSPFutureOpen : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationOSPFutureOpen();

    static IObject* defaultOSPFutureOpenObsType();
};

/****************************************************************************/

class MARKET_DLL ObservationFuture : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationFuture();

    static IObject* defaultFutureObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationAverage : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationAverage();

    static IObject* defaultAverageObsType();
};

/****************************************************************************/

class MARKET_DLL ObservationOther : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationOther();

    static IObject* defaultOtherObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationNotUsed : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationNotUsed();

    static IObject* defaultNotUsedObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationMIBOpen : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationMIBOpen();

    static IObject* defaultMIBOpenObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationNikkei : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationNikkei();

    static IObject* defaultNikkeiObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationNikkei300 : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationNikkei300();

    static IObject* defaultNikkei300ObsType();
};

 /****************************************************************************/

class MARKET_DLL ObservationTime : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationTime();

    static IObject* defaultTimeObsType();
};

  /****************************************************************************/

class MARKET_DLL ObservationPrecio : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationPrecio();

    static IObject* defaultPrecioObsType();
};

  /****************************************************************************/

class MARKET_DLL ObservationPremi : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationPremi();

    static IObject* defaultPremiObsType();
};

  /****************************************************************************/

class MARKET_DLL ObservationPrezzo : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationPrezzo();

    static IObject* defaultPrezzoObsType();
};

  /****************************************************************************/

class MARKET_DLL ObservationRifer : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationRifer();

    static IObject* defaultRiferObsType();
};

  /****************************************************************************/

class MARKET_DLL ObservationSQSP500 : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationSQSP500();

    static IObject* defaultSQSP500ObsType();
};

/****************************************************************************/

class MARKET_DLL ObservationOptionExpiry : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationOptionExpiry();

    static IObject* defaultOptionExpiryObsType();
};

/****************************************************************************/

class MARKET_DLL ObservationTopix : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

private:
    ObservationTopix();

    static IObject* defaultTopixObsType();
};

/****************************************************************************/
// NOTE THIS ONE IS CURRENTLY EXPOSED PUBLICLY
class MARKET_DLL StruckObservationType : public ObservationType {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    virtual string toString() const;

    virtual string indicativeTime() const;

    ObservationType* getPrimaryObsType() const;

    ObservationType* getFXObsType() const;

    StruckObservationType(ObservationType* primary,
                          ObservationType* fx);
private:
    StruckObservationType();

    ObservationTypeSP primaryObsType;
    ObservationTypeSP fxObsType;

    static IObject* defaultStruckObsType();
};

DRLIB_END_NAMESPACE

#endif

