/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  General Class For ODE Functions.  
 *	Abstract Class
 *  \brief
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */
 
#ifndef _INGPNUMLIB_ODEFUNCTIONS_H
#define _INGPNUMLIB_ODEFUNCTIONS_H


#include "firsttoinc.h"

#include "gpbase/countedptr.h"
#include "gpbase/functor.h"
#include "gpbase/typedef.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/curve.h"


#include <vector>
CC_USING_NS(std,vector)

#include <map>
CC_USING_NS(std,map)



CC_BEGIN_NAMESPACE(ARM)

class ARM_ODEFunc : public ARM_RootObject
{
	public:
    enum SolverType
    {
        RK4Constant=0,
        RK5Adaptative
	};

	virtual void derivs(double x, std::vector<double>* yt, std::vector<double>* dyt) const = 0;
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/