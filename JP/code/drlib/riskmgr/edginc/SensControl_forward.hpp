
#ifndef EDG_SENS_CONTROL_FORWARD_H
#define EDG_SENS_CONTROL_FORWARD_H

#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class SensControl;

typedef SensControl                CSensControl;
typedef smartPtr<SensControl>      SensControlSP;
typedef smartConstPtr<SensControl> SensControlConstSP;
typedef array<SensControlSP, SensControl> SensControlArray;
typedef smartPtr<SensControlArray> SensControlArraySP;
typedef smartConstPtr<SensControlArray> SensControlArrayConstSP;

DRLIB_END_NAMESPACE

#endif

