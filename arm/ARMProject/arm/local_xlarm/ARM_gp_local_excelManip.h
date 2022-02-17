/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_gplocal_excelManip.h,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */

////////////////////////////////////////////////
/// This file is for manipulation of excel type
/// pricer... Everything specific to the generic
/// pricer should be included here and not in ARM_local_interglob.cpp
/// to avoid exporting the generic pricer in the 
/// local dllarm project (activex project)
/// since the interglob.h file is included in the 
/// local dllarm project (activex project)
////////////////////////////////////////////////

#ifndef _ARM_GP_LOCAL_EXCELMANIP_H
#define _ARM_GP_LOCAL_EXCELMANIP_H

#include <string>			/// for string
using std::string;

#include <GP_Base\gpbase\valuetype.h>
using ARM::ARM_GP_VALUE_TYPE;

/// function to convert a XLtype to an ARM_GP_VALUE_TYPE and vice versa
ARM_GP_VALUE_TYPE ConvertXLtoARMType( long type );
long ConvertARMToXLType( ARM_GP_VALUE_TYPE type );


/// function to convert a ARM_GP_VALUE_TYPE to the appropriate ARM_GP_VALUE_TYPE 
/// according to the corresponding string... and vice versa
ARM_GP_VALUE_TYPE ConvertToGenPricerType( ARM_GP_VALUE_TYPE type, string& val );
ARM_GP_VALUE_TYPE ConvertFromGenPricerType( ARM_GP_VALUE_TYPE type);

#endif	
