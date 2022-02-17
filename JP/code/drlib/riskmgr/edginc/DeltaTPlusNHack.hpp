///////////////////////////////////////////////////////////////////////////////
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaTPlusNHack.hpp
//
//   Description : Delta T+n but implemented to conform to current database "shape" 
//                 since otherwise won't get used for ages
//
//   Date        : Mar 2003
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EDG_DELTA_T_PLUS_N_HACK_H
#define EDG_DELTA_T_PLUS_N_HACK_H
#include "edginc/DeltaTPlusN.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for Delta T+n */
class RISKMGR_DLL DeltaTPlusNHack: public DeltaTPlusN {
public:
    static CClassConstSP const TYPE;

    /** constructor */
    DeltaTPlusNHack(double shift, int offset, HolidaySP hols);

protected:

    DeltaTPlusNHack(const CClassConstSP& clazz,
                    const string&        sensName,
                    double               shiftSize,
                    int                  offset);

    DeltaTPlusNHack(const CClassConstSP& clazz,
                    const string&        sensName,
                    int                  offset);

private:
    /** for reflection */
    DeltaTPlusNHack();
    DeltaTPlusNHack(const DeltaTPlusNHack &rhs);
    DeltaTPlusNHack& operator=(const DeltaTPlusNHack& rhs);
    friend class DeltaTPlusNHackHelper;
};

typedef smartConstPtr<DeltaTPlusNHack> DeltaTPlusNHackConstSP;
typedef smartPtr<DeltaTPlusNHack> DeltaTPlusNHackSP;

DRLIB_END_NAMESPACE

#endif
