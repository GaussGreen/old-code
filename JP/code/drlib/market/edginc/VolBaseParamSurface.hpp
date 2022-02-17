//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBaseParamSurface.hpp
//
//   Description : Implements the tweak dispatch for Parametrized Vol Surfaces
//                 sensitivities
//
//   Author      : Jean-Noël Juston
//
//   Date        : 01 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_VOL_BASE_PARAM_SFC_HPP
#define EDG_VOL_BASE_PARAM_SFC_HPP

#include "edginc/VolBase.hpp"
#include "edginc/VolParamTweak.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/Theta.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE
class VolRequestTime;
class CVolRequestLN;
class CVolRequestDVF;

/** A parameterised CVolBase that describes an implied vol surface
    This class implements the tweak dispatch for such VolBases. To
    create a new parameterised vol you should derive from this class
    and implement the two pure virtual methods (ComputeImpVol and
    spotVolSurfaceFromStrikes). You also probably want to override the
    buildCache() method. The only method
    you should need to call is getBackboneSurface(). See
    VolNormalLog.cpp for an example */
class MARKET_DLL CVolBaseParamSurface: public CVolBase,
                            virtual public IPDFCalculator,
                            virtual public IRestorableWithRespectTo<VolParallel>,
                            virtual public IRestorableWithRespectTo<VolPointwise>,
                            virtual public VegaMatrix::IRestorableShift,
                            virtual public RootTimeVega::IRestorableShift,
                            virtual public VegaSkewParallel::IRestorableShift,
                            virtual public VegaSkewPointwise::IRestorableShift,
                            virtual public VolLevel::Shift,
                            virtual public VolParallelShift::Shift,
                            virtual public VolBenchmarkShift::Shift,
                            virtual public PowerVega::Shift,
                            virtual public Theta::IShift,
                            virtual public VolRelativeShift::IShift,
                            virtual public VolAbsoluteShift::IShift {
public:
    static CClassConstSP const TYPE;

    ~CVolBaseParamSurface();

    //// copies over the fields as appropriate
    virtual IObject* clone() const;

    /** Returns name of vol */
    virtual string getName() const;

    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Creates processed vol - CVolProcessedBS or DVF */
    CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      asset) const;

    /** Creates Struck processed vol - CVolProcessedBS or DVF */
    CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    /** VegaParallel */
    string sensName(const VolParallel*) const;
    TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak);
    void sensRestore(const PropertyTweak<VolParallel>& tweak);
    
    /** VegaMatrix Interface */
    string sensName(VegaMatrix* shift) const;
    ExpiryArrayConstSP sensExpiries(VegaMatrix* shift) const;
    bool sensShift(VegaMatrix* shift);
    void sensRestore(VegaMatrix* shift);

    /** VegaPointwise Interface */
    string sensName(const VolPointwise*) const;
    ExpiryWindowArrayConstSP sensQualifiers(const VolPointwise*) const;
    TweakOutcome sensShift(const PropertyTweak<VolPointwise>&);
    void sensRestore(const PropertyTweak<VolPointwise>&);

    /** RootTimeVega Interface */
    string sensName(RootTimeVega* shift) const;
    bool sensShift(RootTimeVega* shift);
    void sensRestore(RootTimeVega* shift);

    /** VegaSkewParallel Interface */
    string sensName(VegaSkewParallel* shift) const;
    bool sensShift(VegaSkewParallel* shift);
    void sensRestore(VegaSkewParallel* shift);

    /** VegaSkewPointwise Interface */
    string sensName(VegaSkewPointwise* shift) const;
    ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    bool sensShift(VegaSkewPointwise* shift);
    void sensRestore(VegaSkewPointwise* shift);

    /** VolLevel */
    string sensName(VolLevel* shift) const;
    bool sensShift(VolLevel* shift);

    /** VolParallelShift */
    string sensName(VolParallelShift* shift) const;
    bool sensShift(VolParallelShift* shift);

    /** VolBenchmarkShift Interface */
    string sensName(VolBenchmarkShift* shift) const;
    bool sensShift(VolBenchmarkShift* shift);

    /** PowerVega Shift Interface */
    string sensName(PowerVega* shift) const;
    bool sensShift(PowerVega* shift);

    /** Theta Shift Interface. */
    bool sensShift(Theta* shift);
    
    /** VolRelativeShift scenario */
    virtual string sensName(VolRelativeShift* shift) const;
    virtual bool sensShift(VolRelativeShift* shift);

    /** Implements VolAbsoluteShift scenario */
    virtual string sensName(VolAbsoluteShift* shift) const;
    virtual bool sensShift(VolAbsoluteShift* shift);

    virtual PDFCalculator* getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const;

    /** Returns the parameterised vol - needed to calls can be made directly
        to ComputeImpVol (which we want to keep out of this interface) by,
        for example, the pdfCalculator. 
        Note that either getMarket() must have been called or this object needs
        to have been constructed with a supplied vol surface for this method
        to work */
    CVolParamConstSP getVolParam();
protected:
    //// constructor
    CVolBaseParamSurface(const CClassConstSP& clazz);

    CVolBaseParamSurface(const CClassConstSP& clazz, const string& name);
    
    //// constructor (eg see VolSpline on how to calculate volParam).
    //// Clients must call buildCache() themselves (if needed) since it
    //// can't be called here (object not fully created)
    CVolBaseParamSurface(const CClassConstSP&    clazz,
                         const VolSurface&       volSurface,
                         const CVolParamSP&      volParam);

    /** This method is invoked every time the backbone is altered. The
     default implementation does nothing but most derived classes should
     override this method. */
    virtual void buildCache();

	virtual void buildCache(VolSurfaceSP volSurfaceForBackbone);
    
    /** method that builds a CVolParam. This method is invoked after
        getMarket. Any overriding implementations of this method MUST
        not take a reference to 'this' instead they should take a
        clone of 'this'. This, however, should not be necessary since
        'this' is passed as a parameter in the ComputeImpVol and
        spotVolSurfaceFromStrikes methods */
    virtual CVolParam* createVolParam() const = 0;

    /** returns a constant reference to surface to be used for the backbone */
    virtual const VolSurface* getBackboneSurface() const;

    /** Helper method - may be used by derived classes */    
    void getMarket(const IModel*     model, 
                   const MarketData* market,
                   const string&     name);


private:
    string             name;    // name of the vol
protected:
    // used to get the vol @ strikeRef - backbone changes under some scenarios
    VolSurfaceSP        volSurfaceForBackbone;
private:
    IVolProcessed* getProcessedVolTime(const VolRequestTime* volRequest,
                                       const CAsset*         asset) const;
    
    IVolProcessed* getProcessedVolLN(const CVolRequestLN* volRequest,
                                     const CAsset*        asset) const;

    IVolProcessed* getProcessedVolDVF(const CVolRequestDVF* volRequest,
                                      const CAsset*         asset) const;

    friend class ImpliedVolAddin;
    // none below are registered 
    CVolParamSP         myVolParam; // the current parameterised vol $unregistered
    /* holds all the previous vols ie before a shift, myVolParam is
       added to this array */
    vector<CVolParamSP> previousVols; // $unregistered

    CVolParamSP getMyVolParam(); // returns the current parameterised vol
    bool shiftBackbone(SensControl* shift);
    bool genericRestorableShift(ITweakID* shift);
    bool genericNonRestorableShift(ITweakID* shift);
    void genericRestore();
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif //EDG_VOL_BASE_PARAM_SFC_HPP
