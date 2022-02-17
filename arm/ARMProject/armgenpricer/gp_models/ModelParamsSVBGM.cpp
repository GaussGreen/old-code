/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
*/

/// gpmodels
#include "gpmodels/ModelParamsSVBGM.h"

/// gpbase

/// gpinfra



CC_BEGIN_NAMESPACE( ARM )


string ARM_ModelParamsSVBGM::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	
    os << "ARM_ModelParamsSVBGM\n";
    os << "----------------------\n\n";
    os << ARM_ModelParams::toString();
	os << "\n\n";
	os << indent << "Effective number of factors : \n";
	os << itsNbFactors << "\n";

	return os.str();
}



CC_END_NAMESPACE()