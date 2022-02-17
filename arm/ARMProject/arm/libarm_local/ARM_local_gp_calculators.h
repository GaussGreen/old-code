/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculator.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_CALCULATORS_H
#define ARMLOCAL_GP_CALCULATORS_H

#include "ARM_local_gp_calculators_tools.h"
#include "ARM_local_gp_FXVanilla_calculators.h"
#include "ARM_local_gp_IRCallable_calculators.h"
#include "ARM_local_gp_IRPathDep_calculators.h"
#include "ARM_local_gp_PRDCandCo_calculators.h"
#include "ARM_local_gp_TARNandCo_calculators.h"


extern long ARMLOCAL_Calculator_Set(
        const long& csoId,
        const long& dataId,
        const string& setType,
		const vector< string >& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId );

#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
