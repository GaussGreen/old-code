//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPICutoff.hpp
//
//   Description : Cutoff interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_CUTOFF_HPP
#define EDR_SPI_CUTOFF_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/SPIUtil.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


/*****************************************************************************/
// run time post cutoff exposure behaviour
class PRODUCTS_DLL ISPIPostCutoffRT {
public:
	virtual double getPostCutoffExposure(double B, double BF, int iStep) const = 0;
	virtual double getPostCutoffCrash() const = 0;
    virtual ~ISPIPostCutoffRT();
};
typedef refCountPtr<ISPIPostCutoffRT> SPIPostCutoffRTSP;

class PRODUCTS_DLL SPIZeroCutoff : public ISPIPostCutoffRT {
public:
	SPIZeroCutoff ();
	virtual double getPostCutoffExposure(double B, double BF, int iStep) const;
	virtual double getPostCutoffCrash() const;
};

class PRODUCTS_DLL SPIBMinusBFCutoff : public ISPIPostCutoffRT {
public:
	SPIBMinusBFCutoff (double crash, int finalRebal);
	virtual double getPostCutoffExposure(double B, double BF, int iStep) const;
	virtual double getPostCutoffCrash() const;

	double crash;
	int finalRebal;
};

/*****************************************************************************/
// this is the actual interface of what a cutoff object 
// needs to do internally
class PRODUCTS_DLL ICutoffSPI {
public:
    virtual bool isCutoff(double  UE,
                          double  TE,
                          double  B,
                          double  BF) const = 0;

    virtual void crossValidate(double equityExposureMin,
							   DoubleArray crashSizes) = 0;

	virtual SPIPostCutoffRTSP getPostCutoff(int finalRebal) = 0;

    virtual ~ICutoffSPI();
};
DECLARE_REF_COUNT(ICutoffSPI);

// this is the external interface for abstraction so that users can bolt in 
// any Cutoff type they want - note this includes the SPICutoffWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ICutoffSPI as soon as possible
class PRODUCTS_DLL ICutoffSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    virtual ICutoffSPISP getCutoffSPI() = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<ICutoffSPIInterface> ICutoffSPIInterfaceSP;

// base class to incorporate post cutoff exposure calcs
class PRODUCTS_DLL SPIPostCutoff : public SPIInterfaceIMS,
                      virtual public ICutoffSPI {
public:
    static CClassConstSP const TYPE;

	static const string ZERO;
	static const string B_MINUS_BF;

	SPIPostCutoff();

	void validatePostCutoff(DoubleArray crashSizes);

	virtual SPIPostCutoffRTSP getPostCutoff(int finalRebal);

	virtual bool isCutoff(double  UE,
                          double  TE,
                          double  B,
                          double  BF) const;

    virtual void crossValidate(double equityExposureMin,
							   DoubleArray crashSizes);


protected:
	StringArray		postCutoffExposure;// after cutoff what exposure remains
	int				numNonZero; // number of assets with non-zero exposure post cutoff
	double			crash;// to fix weights of each risky asset post cutoff

	SPIPostCutoff(CClassConstSP clazz);
private:
    SPIPostCutoff (const SPIPostCutoff & rhs); // not implemented
    SPIPostCutoff & operator=(const SPIPostCutoff & rhs); // not implemented

    static IObject* defaultSPIPostCutoff();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

class PRODUCTS_DLL SPICutoffStd : public SPIPostCutoff,
                                virtual public ICutoffSPIInterface {
public:
    static CClassConstSP const TYPE;
    SPICutoffStd (); // for reflection

    bool isCutoff(double  UE,
                  double  TE,
                  double  B,
                  double  BF) const;

    void crossValidate(double equityExposureMin,
					   DoubleArray crashSizes);

    virtual ICutoffSPISP getCutoffSPI();

private:
    SPICutoffStd (const SPICutoffStd & rhs); // not implemented
    SPICutoffStd & operator=(const SPICutoffStd & rhs); // not implemented

    double                  equityExposureCutoff;        // equity exposure cutoff

    static IObject* defaultSPICutoffStd ();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPICutoffStd> SPICutoffStdSP;

////////////////////////////////////////

class PRODUCTS_DLL SPICutoffNearBondFloor : public SPIPostCutoff,
                                virtual public ICutoffSPIInterface {
public:
    static CClassConstSP const TYPE;
    SPICutoffNearBondFloor ();// for reflection

    bool isCutoff(double  UE,
                  double  TE,
                  double  B,
                  double  BF) const;

    void crossValidate(double equityExposureMin,
					   DoubleArray crashSizes);

    virtual ICutoffSPISP getCutoffSPI();

private:
    SPICutoffNearBondFloor (const SPICutoffNearBondFloor & rhs); // not implemented
    SPICutoffNearBondFloor & operator=(const SPICutoffNearBondFloor & rhs); // not implemented

    double                  bondFloorBuffer;        // equity exposure cutoff

	static IObject* defaultSPICutoffNearBondFloor ();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPICutoffNearBondFloor> SPICutoffNearBondFloorSP;

////////////////////////////////////////

#define SPI_CUTOFF_TYPE_STD   "Standard"
#define SPI_CUTOFF_TYPE_BF    "NearBondFloor"

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPICutoffWrapper : public CObject,
                                virtual public ICutoffSPIInterface {
public:
    static CClassConstSP const TYPE;

    string                    SPICutoffType;
    SPICutoffStdSP            cutoffStd;
    SPICutoffNearBondFloorSP  cutoffNearBF;

public:

    virtual ICutoffSPISP getCutoffSPI();

    // validation
    void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    // for reflection
    SPICutoffWrapper();

    static IObject* defaultSPICutoffWrapper();
};
typedef smartPtr<SPICutoffWrapper> SPICutoffWrapperSP;

DRLIB_END_NAMESPACE

#endif

