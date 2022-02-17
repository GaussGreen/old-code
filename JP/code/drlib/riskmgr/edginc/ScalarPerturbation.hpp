//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : ScalarPerturbation.hpp
//
//   Description : Defines the ability to shift a named double
//
//   Author      : Mark A Robson
//
//   Date        : 29 June 2005
//
//
//----------------------------------------------------------------------------


#ifndef QRD_SCALAR_PERTURBATION_HPP
#define QRD_SCALAR_PERTURBATION_HPP
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/TweakNameListID.hpp"
#include "edginc/TweakID.hpp"
#include "edginc/ScalarPerNameShift.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ScalarPerturbation)

/** This class offers minimal functionality - it provides the two additional
    methods below (as well as declaring it implements IScalarPerNameShift).
    The cost is, though, being forced to use 'shiftSize' as the field name */
class RISKMGR_DLL ScalarPerturbation: public Perturbation,
                                      public virtual ITweakNameListID,
                                      public virtual ITweakID,
                                      public virtual IScalarPerNameShift {
public:
    static CClassConstSP const TYPE;
    ~ScalarPerturbation();

    /** Returns the scalar shift size */
    virtual double getShiftSize() const;

    /** Sets the scalar shift size */
    virtual void setShiftSize(double shiftSize);

    /** IScalarPerNameShift implementation, for ImpliedScalarShift */
    virtual IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                                    double shiftSize,
                                                    IObjectSP world);

    /** IScalarPerNameShift implementation, for ImpliedScalarShift */
    virtual OutputNameArrayConstSP allNames(const IObject*) const;

/** Returns array of output names which need to be tweaked for this
    sensitivity. In particular, if there are toTweak names, then returns
    these otherwise generates a list based on object supplied (which
    is typically, but needn't be, a TweakGroup */
    virtual OutputNameArrayConstSP names(const IObject* tweakGroup) const;

    virtual IScalarRiskPropertyConstSP asRiskProperty() const;

protected:
    ScalarPerturbation(CClassConstSP clazz);
    ScalarPerturbation(CClassConstSP clazz, double shiftSize);

    double shiftSize;
private:
    OutputNameArraySP toTweak; /* to be removed - ignored. Legacy implementation
                                  reasons */
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
