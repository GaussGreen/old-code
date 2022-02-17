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
 
#include <glob/firsttoinc.h>
#include "gpbase/rootobject.h"
#include <cmath>
#include "gpnumlib/odefunctions.h"
#include <glob/expt.h>   // for the exceptions
#include <algorithm>

using namespace std;


CC_BEGIN_NAMESPACE(ARM)

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/