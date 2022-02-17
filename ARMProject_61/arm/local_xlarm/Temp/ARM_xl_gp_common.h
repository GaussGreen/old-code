

#ifndef ARM_XL_GP_COMMON_H
#define ARM_XL_GP_COMMON_H

#include <ARM\libarm_local\ARM_local_glob.h>
#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

struct ARM_GramFctorArg;

LPXLOPER ARM_GramFunctorToXLOPER(const ARM::ARM_GramFctorArg& gramFunctorArg, XLOPER& XL_result, ARM_result C_result, bool fromExcel = true, int index = 0);

#endif