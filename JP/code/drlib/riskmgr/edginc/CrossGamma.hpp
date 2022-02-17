//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CrossGamma.hpp
//
//   Description : CrossGamma sensitivity
//
//   Author      : Mark A Robson
//
//   Date        : 16 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CROSS_GAMMA_HPP
#define EDG_CROSS_GAMMA_HPP
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/AssetSpotGreek.hpp"
#include "edginc/SpotCrossDerivative.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for CROSS_GAMMA - a second order cross
    derivative. The Delta::Shift interface defines what needs to be
    implemented for CROSS_GAMMA to be calculated */
class RISKMGR_DLL CrossGamma: public Sensitivity,
                  virtual public Additive,
                  virtual public IAssetSpotGreek,
                  virtual public ITwoSidedDeriv,
                  virtual public ISpotCrossDerivative{
public:
    friend class CrossGammaHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    CrossGamma(double     shiftSize);

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    /** calculates crossGamma - override delta calculation  */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results);

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const;

    /** When calculating cross gamma, several pricings have to be done. These
        either involve a single shift or a double shift. This routine, called
        upon the Sensitivity returned from Control::getSens(), returns the
        shifts which have been made for the current pricing call */
    virtual ScalarShiftArray getComponentShifts() const;

    /** Returns the shifts which define the 'area of operation of the
        shifts'. This means the names of what is going to be shifted,
        the largest shift size for each shift and the nature of the
        shift (which is captured by the type of the SensControl). The
        names should be stored so that they can be retrieved using the
        names() method. */
    virtual ScalarShiftArray definingShifts(CInstrument*   instrument,
                                            IModel*        model,
                                            Control*       control) const;
    
private:
    /** for reflection */
    CrossGamma();
    CrossGamma(const CrossGamma &rhs);
    CrossGamma& operator=(const CrossGamma& rhs);

    // **** registered fields ****
    double shiftSize;
    // not registered
    ScalarShiftArray   shifts; // $unregistered
};

typedef smartConstPtr<CrossGamma> CrossGammaConstSP;
typedef smartPtr<CrossGamma> CrossGammaSP;

DRLIB_END_NAMESPACE

#endif
