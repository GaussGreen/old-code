//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetStochGarfSolver.hpp
//
//   Description : mixing of one factor finite difference algorithm
//
//----------------------------------------------------------------------------

#ifndef EDG_FD1D_RET_STOCHGARF_SOLVER_HPP
#define EDG_FD1D_RET_STOCHGARF_SOLVER_HPP

#include "edginc/FD1DRetStochGarf.hpp"
#include "edginc/FD1DRetSolver.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------------------

/** FD1DRetStochGarfSolver algorithm class that supports FD1DRetStochGarfSolver */

//class FD1DRetStochGarfSolver : public virtual FDModel::IFDSolver{

class TREE_DLL FD1DRetStochGarfSolver : public FD1DRetSolver{

public:

    //FD1DRetStochGarfSolver(vector<FD1DRet* >engines);
    FD1DRetStochGarfSolver(FD1DRetStochGarf* engineMV);
    virtual ~FD1DRetStochGarfSolver();

    /** solve backward or forward induction through all steps. returns price*/
    virtual void roll();

    //-----------------------------------------------------------------------------------

private: 
    FD1DRetStochGarf* engineStochGarf;

    //vector<FD1DRet* > engines;
    //  ---------------------------  local ft --------------------                         

  
//----------------------------------------------------------------------------
// temporary variable for DR

//----------------------------------------------------------------------------


    int numOfEngines;
//    vector<FD1DRet* > engines; 

};

DRLIB_END_NAMESPACE
#endif
