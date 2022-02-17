#ifndef FD1DMultiSolver_HPP
#define FD1DMultiSolver_HPP

#include "edginc/FDModel.hpp"
#include "edginc/FD1DMulti.hpp"
#include "edginc/SparseMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1DMultiSolver : virtual public FDModel::IFDSolver
{
public:
    FD1DMultiSolver( const FD1DMulti & model );

private: 
    // FDModel::IFDSolver implementation
    virtual void roll();

    void rollOneStep( int step, int direction, const FDProductArray & products );

    void setBoundary(int step, 
					double pv,
					const TreeSlice & currValue,
					const TreeSlice & prevValue );

    void updateBoundary(int step, 
						bool nonNegative,
						const TreeSliceEQ & currValue );

    void interpBoundary( const TreeSlice & s, int i, int d );

    void calcCoeffPde( int step );
    void calcCoeffPdeFix( int step );
    void calcCoeffPdeVar( int step );

    void calcSource( int step );

    void euroOneStepWithSource1D(
        int step, bool nonNegative,
        const TreeSliceEQ & currValue,
        const TreeSliceEQ & prevValue );

	 /**-------------------------------------------------------------------------
		Approximate u(j, x_n + jump(l,j,x_n)) = weight (l,j,n)* u(j,index(l,j,n))
										+ (1-weight(l,j,n))*u(j,index(l,j,n)+1)
		i.e., convert jump(l,j,n) into weight(l,j,n) and index(l,j,n)         
	 --------------------------------------------------------------------------*/
	void interpJump( int step );

	// solve A X = B, where A is sparse matrix and B is vector
	void euroOneStepWithSource(double* sol_n, 
							   double* sol_nplus1);

    // temporary slices to perform DEV into
    vector< vector< TreeSliceSP > > tempSlices;

    // PDE coeffs
	DoubleArray			a[2];
	DoubleArray			c[2];
	DoubleArray			f[2];
	DoubleArray			g[2];
	DoubleArray			q[2];
	DoubleArray			jump[2];

	DoubleArray			jumpWeight[2];
	DoubleArray			jumpIdx[2];

	SparseCollection	xCoeff[2];
	DoubleArray			xSrc;
	DoubleArray			xVar;
	DoubleArray			xValue;

	int					nPdes;
	int					dim;

	IntArray			top;
	IntArray			bot;

    double				theta;
    const FD1DMulti &	model;
};

DRLIB_END_NAMESPACE
#endif
