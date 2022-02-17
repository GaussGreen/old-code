///////////////////////////////////////////////////////////////////////////////
// Filename: ModelLib.h 
// Author: Adam Gladstone 
// Date: Friday, August 08, 2003 
// 
// Description
//    Public #include file for ModelLib.lib
//
///////////////////////////////////////////////////////////////////////////////

#if !defined(CLSGEN_CMODELLIB_H__90572221_C9AC_11D7_BC8A_00805F9B19DC__INCLUDED_)
#define CLSGEN_CMODELLIB_H__90572221_C9AC_11D7_BC8A_00805F9B19DC__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


/////////////////////////////////////////////////////////////////////////////////
// Includes

// ----------------------------------------------------------------------------
// C Runtime includes
// ----------------------------------------------------------------------------

#include <assert.h>
#include <math.h>

// ----------------------------------------------------------------------------
// ModelLib includes
// ----------------------------------------------------------------------------

#include "TContainer.h"					// Generic container template class
#include "CMatrix.h"					// CVector and CMatrix classes
#include "smart.h"						// Smart reference counted object pointer
#include "MlEqInterpolator.h"			// Interpolator
//#include "utility.h"					// Miscellaneous utility classes
#include "lsinfo.h"						//   This header file contains 'global' structure definitions
										//   prototype declarations for lsgrgc v3.0 ANSI C version
#include "solve.h"						// LSGRGSolver class wrapper
#include "MLEqMaths.h"					// Generic math functions

#include "genericVol.h"					// Generic Strike Class
#include "edcubicspline.h"				// ummm...
#include "MonteCarlo.h"					// Forward skew Monte Carlo class


#endif		//CLSGEN_CMODELLIB_H__90572221_C9AC_11D7_BC8A_00805F9B19DC__INCLUDED_
