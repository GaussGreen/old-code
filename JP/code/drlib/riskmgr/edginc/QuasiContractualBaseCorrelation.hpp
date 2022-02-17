//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : QuasiContractualBaseCorrelation.hpp
//
//   Description : Composite greek: Applies a multiple scenario as defined in
//                 this class, and then computes the sensitivity returned by
//                 derived classes through the "sensitivityToComputeAfterShift"
//                 virtual method (typically BetaSkewParallel/Matrix, used
//                 to compute the base correlation sensitivities of the 
//                 "quasi-contractual" capital structure representation of a
//                 trade).
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_QUASICONTRACTUALBASECORRELATION_HPP
#define QLIB_QUASICONTRACTUALBASECORRELATION_HPP

#include "edginc/SensControlAllNames.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL QuasiContractualBaseCorrelation : 
    public SensControlAllNames,
    public virtual Additive 
{
public:
    const static CClassConstSP TYPE;

    /** Returns the sensitivity to be computed in this sensitivity after 
     * applying the initial shift */
    virtual SensitivitySP sensitivityToComputeAfterShift() = 0;

    /** What an object must implement to be tweakable for this perturbation */
    class RISKMGR_DLL IShift {
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(QuasiContractualBaseCorrelation* shift) = 0;
    };
    
    /** What an object must implement to be able to perform a restorable
        tweak for THETA. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift : public virtual IShift {
    public:
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(QuasiContractualBaseCorrelation* shift) = 0;
    };

    /** How to restore the object after a tweak */
    virtual void restore(IObjectSP obj);

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    virtual CClassConstSP shiftInterface() const;

    /** Shifts the object (which supports being tweaked
        by this type of sens control) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */
    virtual bool shift(IObjectSP obj);

    /** Once used to make a shift, this reports the appropriate divisor
     * for this sensitivity */
    double divisor() const;
    
    /** Is this sensitivity made using a discrete shift (ie a jump) or a
     * an approximately continuous one. Returns true */
    virtual bool discreteShift() const;
    
    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;

    /* Computes this sensitivity */
    virtual void calculate(TweakGroup* tweakGroup, CResults* results);

    /** Returns the level to set the credit spreads to */
    double getCreditSpreadsLevel() const;

    /** Returns the compression ratio to use when rescaling the historical betas
     * dispersion. */
    double getCompressionRatio() const;

    /** Returns the value to set the yPoints of the PiecewiseMappingFunction to */
    double getPiecewiseMappingFunctionLevel() const;

    /** Note betaBasis is not configurable - its value will be hardcoded to null */

    /** Returns the value to set the skews in the skewSurface to */
    double getSkewLevel() const;
    
protected:
    // For Reflection
    QuasiContractualBaseCorrelation(CClassConstSP clazz, 
                                    const string& name);
    
private:
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);                               
                           
    // FIELDS                           
    /** Level to set the credit spreads to */
    double creditSpreadsLevel;

    /** Compression ratio to use when rescaling historical betas dispersion */
    double compressionRatio;

    /** Value to set the yPoints of the PiecewiseMappingFunction to */
    double piecewiseMappingFunctionLevel;

    /** Note index basis  is not configurable - its value is hardcoded so that
        the additive and multiplicative coefficients are set to fix values. 
        See CDSParSpreadsLegalBasis.hpp for details. */

    /** Value to set the skews in the skewSurface to */
    double skewLevel;
};

DRLIB_END_NAMESPACE

#endif
