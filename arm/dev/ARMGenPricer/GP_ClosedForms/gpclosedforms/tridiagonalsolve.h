/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2005
 *
 *  solve the tridiagonal matrix problem 
 *
 *	\file tridiagonalsolve.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Dec 2005
 */
 
#ifndef _GP_CF_TRIDIAGONALSOLVE_H
#define _GP_CF_TRIDIAGONALSOLVE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"



#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpvector.h"
#include <vector>



CC_BEGIN_NAMESPACE(ARM)

using std::vector;

void tridiagonalsolve(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& r, vector<double>& u,vector<double>& gam);

vector<double>* Export_Util_Trigonal_Solve(ARM_GP_Vector* A_Vec,ARM_GP_Vector* B_Vec,ARM_GP_Vector* C_Vec,ARM_GP_Vector* R_Vec);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/