/*************************************************************************
 VectorMatrix.h
 Author: M.Huq, Credit DR
 *************************************************************************/
#ifndef _SC_VECTORMATRIX_H_
#define _SC_VECTORMATRIX_H_

#include "edginc/coreConfig.hpp"
#include <vector>

CORE_BEGIN_NAMESPACE

#if 0
#include "edginc/tnt/tnt.h"
#include "edginc/tnt/vec.h"
#include "edginc/tnt/cmat.h"

// Definitions for vectors and matrices being used internally.
// Definitions for vectors and matrices based upon TNT. See tnt directory.
// TNT = Template Numerical Library from http://math.nist.gov/tnt/
#undef USE_TNT
#define USE_TNT
#undef TNT_UNROLL_LOOPS
#define TNT_UNROLL_LOOPS

typedef TNT::Vector<double> tntd_vec;
typedef TNT::Matrix<double>  tntd_mat;
typedef TNT::Vector<int>    tnti_vec;
typedef TNT::Matrix<int>  tnti_mat;
#else
#include <vector>
typedef std::vector<double> tntd_vec;
typedef std::vector<std::vector<double> >  tntd_mat;
typedef std::vector<int>    tnti_vec;
typedef std::vector<std::vector<int> >  tnti_mat;
#endif


typedef std::vector<double>  std_d_vec;

CORE_END_NAMESPACE

#endif //_SC_VECTORMATRIX_H_
