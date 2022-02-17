/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  function solution of the Riccati Equation with Piece-Wise functions A(t), B(t) and C(t)
 *	df/dt = A(t) * f(t)² + B(t) * f(t) + C(t)
 *  
 *	\file riccati.h
 *
 *  \brief
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */
 
#ifndef _GP_CF_RICCATI_H
#define _GP_CF_RICCATI_H

#include <glob/firsttoinc.h>

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

class ARM_ODEFunc:public ARM_RootObject
{

};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/