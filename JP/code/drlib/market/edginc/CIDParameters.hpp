//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CIDParameters.hpp
//
//   Description : Parameters needed for CID model (for all names, not per name)
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CIDPARAMETERS_HPP
#define QLIB_CIDPARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ICDSParSpreads.hpp"

DRLIB_BEGIN_NAMESPACE

class CIDParameters;
class IModel;
class MarketData;

/************************************************************************/
/* An attempt to restructure the CID parameters while using the         */
/* proper input objects and separating the input data from the          */
/* model as far as possible                                             */
/************************************************************************/

class MARKET_DLL CIDCommonRecord :  public CObject
{
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const { return name; }

    virtual void getMarket(const IModel* model, const MarketData* market);

    // calculation is based purely on the CDS par spread, non-virtual(!)
    double calcSpreadSurvProbability(const DateTime& basedate, const DateTime& timepoint)  const; 

    // calculation of context-dependent survival probability
    // for the common market it is the same as one obtained from spread
    virtual double calcProperSurvProbability(const DateTime& basedate, const DateTime& timepoint) const
        { return calcSpreadSurvProbability(basedate, timepoint); }

    void validateCommon();

public: // at the moment

    /** Name this object will be stored under in the market cache ************/
    string           name;
    /** CDS wrappers *********************** *********************************/
    ICDSParSpreadsWrapper pCDSwrapper; 


    // Volatility
    double           sigma;

    // Mean reversion speed
    double           kappa;

    // If we are to set the desired discretization of thetas
    ExpiryArrayConstSP endDatesForTheta;

protected:
    CIDCommonRecord(const CClassConstSP& clazz);

private:
    CIDCommonRecord();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};
DECLARE(CIDCommonRecord);


class MARKET_DLL CIDNameRecord : public CIDCommonRecord
{
public:
    static CClassConstSP const TYPE;
    void validate(int numJSources, CIDParameters* link);

    // calculation is based on the CDS par spread, with jumps and common markets deducted (non-virtual!)
    double calcIdioSurvProbability(const DateTime& basedate, const DateTime& timepoint) const; 

    // calculation of SurvProb due to jumps from a single source (non-virtual!)
    double calcJumpsSurvProbability(const DateTime& basedate, const DateTime& timepoint, size_t sourceIdx) const;
    double calcJumpsSurvProbability(const DateTime& timepoint, size_t sourceIdx) const; 

    // calculation of context-dependent survival probability
    // for the name-specific portion it is the same as calculated for Idio
    virtual double calcProperSurvProbability(const DateTime& basedate, const DateTime& timepoint) const
        { return calcIdioSurvProbability(basedate, timepoint); }

public: // at the moment
    double                 A;  // A is the weight for common diffusion for this name
    DoubleArrayConstSP     jInitSpreads;
    DoubleArrayConstSP     jImpacts;
    DoubleArrayConstSP     jSizes;
    // recovery rates ---------------------------------------------------------
    double                 cataR;
    // mean recovery rate always come from the CDSs

private:

    // Internally generated
    CIDParameters*  pCIDParameters;

    CIDNameRecord();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};
DECLARE(CIDNameRecord);


class MARKET_DLL CIDJumpFactor :  public CObject
{
public:
    static CClassConstSP const TYPE;

    /** Returns name of vol */
    virtual string getName() const { return name; }

    void validate();

public:  // at the moment
    string name;
    double jFreqs;
    double jDecays;

private:

    CIDJumpFactor();
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};
DECLARE(CIDJumpFactor);



/******************************************************************************
*******************************************************************************
** This class is an interface of CIDParameters used by MCPathGenCID          **
** All its methods should be constant.                                       **
*******************************************************************************
******************************************************************************/
class MARKET_DLL ICIDParameters 
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

///////////////////////////////////////////////////////////////////////////////
//  A C C E S S O R S                                                        //
///////////////////////////////////////////////////////////////////////////////

    /** Returns common diffusion **********************************************/
    virtual CIDCommonRecordConstSP getCommonFactor()  const = 0;

    /** Returns idiosyncratic diffusion factor for a given name ***************/
    virtual CIDNameRecordConstSP getSingleNameRecord(const string & gName) const  = 0;

    /** Returns number of jump factors and jump factors themselves ************/
    virtual int getNumOfJumpFactors() const = 0;
    virtual CIDJumpFactorConstSP getJumpFactor(int gIdxOfFactor) const = 0;
};


class MARKET_DLL CIDParameters : public RationalisedCreditEngineParameters
                               , public ICIDParameters
{
public:
    static const string DEFAULTNAME;
    static CClassConstSP const TYPE;
    DateTime                    today;


    string getName() const { return name; }
    /** getMarket finalizes construction of 'this' given the rest of the market **/
    virtual void getMarket(const IModel* model, const MarketData* market);

    /**
      * Returns whether the name is modeling stochastic recovery
      * - true : stochastic recovery
      * - false : fixed recovery
      * in CID case returns TRUE **********************************************
    **/
    virtual bool hasStochasticRecovery() const 
        {   return true; }

    /** Returns common diffusion **********************************************/
    virtual CIDCommonRecordConstSP getCommonFactor() const
        {   return commonFactor;    }

    /** Returns idiosyncratic diffusion factor for a given name ***************/
    virtual CIDNameRecordConstSP getSingleNameRecord(const string & gName)  const
        {   return mapNameToRecord[gName]; }
    virtual int getNumOfJumpFactors() const
        {   return jumpSources->size(); }
    virtual CIDJumpFactorConstSP getJumpFactor(int gIdxOfFactor)  const
        {   return (*jumpSources)[gIdxOfFactor];   }

    bool VerifyIfCalibrationIsValid() const 
        {   return verifyIfCalibrationIsValid; }

public:
    void     validatePop2Object();
    IObject* clone() const;

private:
    CIDParameters();

    void operator=(const CIDParameters & g2Clone); // not defined
    CIDParameters( const CIDParameters & g2Clone); // not defined

    static IObject*  defaultConstructor();
    static void      load(CClassSP& clazz);


private:
    string                      name;


    // Per name storage of all idiosyncratic information
    CIDNameRecordArraySP        singleNameRecords;

    // Storage for common factor
    CIDCommonRecordSP           commonFactor;

    // jump Sources ------------------------------------------------------------
    CIDJumpFactorArraySP        jumpSources;

    // switches for diffusions ------------------------------------------------
    string                      floorTypeForCIR;
    bool                        verifyIfCalibrationIsValid;

    // internal map for easy indexing. Might go away at some point
    // not to be initialized or seen outside of this class. Used ONLY by accessors that use names

    mutable std::map< string, CIDNameRecordSP> mapNameToRecord;
};


///////////////////////////////////////////////////////////////////////////////
// Support for smart pointers
#ifndef QLIB_CIDPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CIDParameters>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CIDParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CIDParameters>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CIDParameters>);
#endif

// Support for wrapper
typedef MarketWrapper<CIDParameters> CIDParametersWrapper;
#ifndef QLIB_CIDPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CIDParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CIDParameters>);
#endif
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif //QLIB_CIDPARAMETERS_HPP
