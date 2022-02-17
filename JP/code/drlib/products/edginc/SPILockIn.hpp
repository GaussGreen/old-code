//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPILockIn.hpp
//
//   Description : Lock In interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_LOCKIN_HPP
#define EDR_SPI_LOCKIN_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/SPIUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// this is the actual interface of what a lock in object 
// needs to do internally
class PRODUCTS_DLL ILockInSPI {
public:

    virtual double getInitialLockIn() const = 0;

    virtual void init(const DateTimeArray& sampleDates,
                      const DateTime&      lastRebalDate) = 0; // allows mapping to lock-in dates and efficient use of "iStep"

    virtual void apply(double&  BL,
                       double   B,
                       double   BF,
                       int      iStep) const = 0;

    virtual DateTimeArray getEssentialDates() const = 0;

	virtual bool doesNothing() const = 0;
};

// this is the external interface for abstraction so that users can bolt in 
// any LockIn type they want - note this includes the SPILockinWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ILockInSPI as soon as possible
class PRODUCTS_DLL ILockInSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    // Painful - the optional strike field in SPI instrument needs to 
    // be able to get lockin info early to support backwards compatability
    // so this is needed publicly. That's the only reason - else would be private.
    // aslos used at object creation to check compatibility with coupon type
    virtual const ILockInSPI* getLockInSPINoInit() const = 0;

    // does the init as well
    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate) = 0;
private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<ILockInSPIInterface> ILockInSPIInterfaceSP;

// XXX Should have a base class for use with Std and Excess. Also the Schedule is a schedule of Std - why not of excess?
class PRODUCTS_DLL SPILockInStd : public SPIInterfaceIMS,
                     virtual public ILockInSPI,
                     virtual public ILockInSPIInterface {
public:
    static CClassConstSP const TYPE;

    double          initialLockIn;      
    double          lockInPercentage;   // lock-in level
    DateTimeArray   lockInDates;
    BoolArray       lockInFlag;         

    virtual const ILockInSPI* getLockInSPINoInit() const;

    // does the init as well
    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate);
    
    double getInitialLockIn() const;

    void init(const DateTimeArray& sampleDates,
              const DateTime&      lastRebalDate);

    void apply(double&  BL,
               double   B,
               double   BF,
               int      iStep) const;

    SPILockInStd();// for reflection

    virtual void validatePop2Object();

    virtual DateTimeArray getEssentialDates() const;

	virtual bool doesNothing() const;

private:
    SPILockInStd(const SPILockInStd& rhs); // not implemented
    SPILockInStd& operator=(const SPILockInStd& rhs); // not implemented

    static IObject* defaultSPILockInStd();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPILockInStd> SPILockInStdSP;

class PRODUCTS_DLL SPILockInExcess : public SPIInterfaceIMS,
                        virtual public ILockInSPI,
                        virtual public ILockInSPIInterface  {
public:
    static CClassConstSP const TYPE;

    double          initialLockIn;      
    double          lockInPercentage;   // lock-in level
    DateTimeArray   lockInDates;
    BoolArray       lockInFlag;         

    virtual const ILockInSPI* getLockInSPINoInit() const;

    // does the init as well
    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate);
    
    double getInitialLockIn() const;

    void init(const DateTimeArray& sampleDates,
              const DateTime&      lastRebalDate);

    void apply(double&  BL,
               double   B,
               double   BF,
               int      iStep) const;

    SPILockInExcess();// for reflection

    virtual DateTimeArray getEssentialDates() const;

    virtual bool doesNothing() const;

    virtual void validatePop2Object();

private:
    SPILockInExcess(const SPILockInExcess& rhs); // not implemented
    SPILockInExcess& operator=(const SPILockInExcess& rhs); // not implemented

    static IObject* defaultSPILockInExcess();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPILockInExcess> SPILockInExcessSP;

class PRODUCTS_DLL SPILockInSchedule : public SPIInterfaceIMS,
                          virtual public ILockInSPI,
                          virtual public ILockInSPIInterface  {
public:
    static CClassConstSP const TYPE;

    double          initialLockIn;      
    CashFlowArray   lockIns;            // lock-in levels, one per date
    DoubleArray     buffers;            // [nbLockIn], safety interval between BF & Basket which never locks in
    bool            doBestLockIn;       // false=>if condition not met then no lock-in done, else up to buffer is locked in
    IntArray        lockInIdx;          // If lockInIdx[iStep]>0 there's a lock-in at iStep == lockIns[lockInIdx[iStep]]

    virtual const ILockInSPI* getLockInSPINoInit() const;

    // does the init as well
    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate);
    
    double getInitialLockIn() const;

    void init(const DateTimeArray& sampleDates,
              const DateTime&      lastRebalDate);

    void apply(double&  BL,
               double   B,
               double   BF,
               int      iStep) const;

    SPILockInSchedule(); // for reflection
    
    virtual DateTimeArray getEssentialDates() const;

    virtual bool doesNothing() const;

    virtual void validatePop2Object();

private:
    SPILockInSchedule(const SPILockInSchedule& rhs); // not implemented
    SPILockInSchedule& operator=(const SPILockInSchedule& rhs); // not implemented

    static IObject* defaultSPILockInSchedule();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPILockInSchedule> SPILockInScheduleSP;

class PRODUCTS_DLL SPILockInGrowth : public SPIInterfaceIMS,
                     virtual public ILockInSPI,
                     virtual public ILockInSPIInterface  {
public:
    static CClassConstSP const TYPE;

    double          initialLockIn;      
    double          lockInPercentage;   // lock-in level
    DateTimeArray   lockInDates;
    BoolArray       lockInFlag;         

    virtual const ILockInSPI* getLockInSPINoInit() const;

    // does the init as well
    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate);
    
    double getInitialLockIn() const;

    void init(const DateTimeArray& sampleDates,
              const DateTime&      lastRebalDate);

    void apply(double&  BL,
               double   B,
               double   BF,
               int      iStep) const;

    SPILockInGrowth();// for reflection

    virtual void validatePop2Object();
    
    virtual DateTimeArray getEssentialDates() const;

	virtual bool doesNothing() const;

private:
    SPILockInGrowth(const SPILockInGrowth& rhs); // not implemented
    SPILockInGrowth& operator=(const SPILockInGrowth& rhs); // not implemented

    static IObject* defaultSPILockInGrowth();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPILockInGrowth> SPILockInGrowthSP;

class PRODUCTS_DLL SPILockInCappedReserve : public SPIInterfaceIMS,
                               virtual public ILockInSPI,
                     virtual public ILockInSPIInterface  {
public:
    static CClassConstSP const TYPE;

    double          initialLockIn;      
    double          lockInCap;
    DateTimeArray   lockInDates;
    BoolArray       lockInFlag;         

    virtual const ILockInSPI* getLockInSPINoInit() const;

    // does the init as well
    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate);
    
    double getInitialLockIn() const;

    void init(const DateTimeArray& sampleDates,
              const DateTime&      lastRebalDate);

    void apply(double&  BL,
               double   B,
               double   BF,
               int      iStep) const;

    SPILockInCappedReserve();// for reflection

    virtual void validatePop2Object();
    
    virtual DateTimeArray getEssentialDates() const ;

    virtual bool doesNothing() const;

private:
    SPILockInCappedReserve(const SPILockInCappedReserve& rhs); // not implemented
    SPILockInCappedReserve& operator=(const SPILockInCappedReserve& rhs); // not implemented

    static IObject* defaultSPILockInCappedReserve();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

};
typedef smartPtr<SPILockInCappedReserve> SPILockInCappedReserveSP;

/*****************************************************************************/
#define SPI_LOCK_IN_TYPE_STD                "Standard"
#define SPI_LOCK_IN_TYPE_EXCESS             "Excess"
#define SPI_LOCK_IN_TYPE_SCHEDULE           "Schedule"
#define SPI_LOCK_IN_TYPE_GROWTH             "Growth"
#define SPI_LOCK_IN_TYPE_CAPPED_RESERVE     "CappedReserve"

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPILockInWrapper : public CObject,
                     virtual public ILockInSPIInterface  {
public:
    static CClassConstSP const TYPE;

    string                       SPILockInType;
    SPILockInStdSP               lockInStd;
    SPILockInExcessSP            lockInExcess;
    SPILockInScheduleSP          lockInSchedule;
    SPILockInGrowthSP            lockInGrowth;
    SPILockInCappedReserveSP     lockInCappedReserve;

private:
    ILockInSPI*        theLockIn; // $unregistered

public:

    virtual const ILockInSPI* getLockInSPINoInit() const;

    virtual ILockInSPI* getLockInSPI(const DateTimeArray& rebalanceDates,
                                     const DateTime&      lastRebalDate);

    // validation
    void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

private:
    
    // for reflection
    SPILockInWrapper();

    static IObject* defaultSPILockInWrapper();
};
typedef smartPtr<SPILockInWrapper> SPILockInWrapperSP;

DRLIB_END_NAMESPACE

#endif
