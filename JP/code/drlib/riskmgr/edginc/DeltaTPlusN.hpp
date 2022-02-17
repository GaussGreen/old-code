//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaTPlusN.hpp
//
//   Description : T+N delta shift where N is calculated by the instrument and 
//                 is NOT a user input. For instance SPI will calculate delta, 
//                 lag many rebalance dates after today 
//
//   Author      : Ian Stares
//
//   Date        : 23 Nov 2004
//
//
//----------------------------------------------------------------------------


#ifndef DeltaTPlusN_HPP
#define DeltaTPlusN_HPP

#include "edginc/T1Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** T+N delta shift */
class RISKMGR_DLL DeltaTPlusN: public T1Delta {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;
    const static int    DEFAULT_OFFSET;

    /** What an instrument must implement to allow DELTA_T_PLUS_N */
    class RISKMGR_DLL IDeltaTPlusNImnt : virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        virtual int deltaOffset(const Holiday* hols) const = 0;
    };

    /** constructor */
    DeltaTPlusN(double shift, int offset, HolidaySP hols);

    /** identifies the name used storing associated results in the output 
        (Implements pure virtual function in Sensitivity class) */
    const string& getSensOutputName() const;

    /** From INextDaySensitivity interface */
    virtual SensitivitySP requiredSensitivity(TweakGroup* tweakGroup) const;

    // override the base class offsetRequired so we can calculate offset
    // before then delegating back to base class offsetRequired!
    virtual ThetaSP offsetRequired(const CInstrument* inst) const;

protected:
    DeltaTPlusN(CClassConstSP clazz);

private:
    /** for reflection */
    DeltaTPlusN();
    DeltaTPlusN(const DeltaTPlusN &rhs);
    DeltaTPlusN& operator=(const DeltaTPlusN& rhs);
    friend class DeltaTPlusNHelper;
};

typedef smartConstPtr<DeltaTPlusN> DeltaTPlusNConstSP;
typedef smartPtr<DeltaTPlusN> DeltaTPlusNSP;

DRLIB_END_NAMESPACE

#endif
