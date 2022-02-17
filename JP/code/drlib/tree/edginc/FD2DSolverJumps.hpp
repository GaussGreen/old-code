//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSolverJumps.hpp
//
//   Description : two factor finite difference algorithm
//                  derive for FD2DSolver, treatment for jumps
//
//   Author      : Ning Shen
//                 Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#ifndef FD2DSOLVERJUMPS_HPP
#define FD2DSOLVERJUMPS_HPP

#include "edginc/FD2DSolver.hpp"

#include "edginc/FD2DSVCJ.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------------------

/** FD2DSolver algorithm class that supports FD2DSolver */
//class FD2DSolverJumps : public FD2DSolver{
class FD2DSolverJumps  {

//class TREE_DLL FD2DSolverJumps :
//            public virtual VirtualDestructorBase{

public:

    FD2DSolverJumps(FD2D* model, int xNinput, int yNinput) ;//: FD2DSolver(engine);

    virtual ~FD2DSolverJumps();

    //define the method, it's declared in the FD2DSolver as pure virtual
    //First try: 
    //we can treat Jumps part as explicit in FD2DSVCJ case.
    //with the following ft and FD2DSolver, we should be able to get a price.
    void compute_source_jumps(double* M_price, double* rhsJumps,double dx, double dy);

//-------------------------------------------------------------
//Jumps parts
//-------------------------------------------------------------

private:

    //  Index conversion
    int Pos2Index(int xn, int ym, int n, int m);

    void setup(double dx, double dy);
    
    void constructFFTGrid(double dx, double dy);

    void FFTJumpDen(double dx, double dy);

    void doInterp(double** BigSlice, int xMiddle, int yMiddle) ;

    void mJumpCompQuick(double** BigSlice, d_complex** mFFT_den,
                                    double** BigjumpPart);

    void JumpCompQuick(double** BigSlice, d_complex* FFT_den,
                                    double** BigjumpPart);

    //  The jump probability density. 
    //  moved to model SVCJ
    //  The probability that the jump hits the point (x,y)
    //  The integration of the probability density function
    double jumpPro (double x, double y, double dx, double dy);

    //  calculate the jump, ie just do the convolution
    //  bwetween the jump probability grid and
    //  the option value grid, which can be current or the next.
    void jumpComp(double* slice, double*  jumpPart, double dx, double dy);

    //  convolution of 2 double array (n by m)
    void Conv2(int n, int m, double* A, double* B, double* C);

//----------------------------------------------------------------------------

private :
    FD2D* engine;

    int xN, yN;
    int xMiddle , yMiddle ;

    //do FFT only one time for density
    int xExtend, yExtend;  // set to be xExtend =2 * xN -1,
    vector<double> xF;     //big grid for FFT, but constructed such that it includes the same points as FD grid
    vector<double> yF;

    d_complex* g_complex; //jump density
    d_complex* g_FFT;

    //price used to do the convol for jump parts
    //and also used to store the results of convolution for jumps
    d_complex* p_complex;
    d_complex* p_FFT;

    //for test 
    d_complex** mg_complex; //jump density
    d_complex** mg_FFT;

    //price used to do the convol for jump parts
    //and also used to store the results of convolution for jumps
    d_complex** mp_complex;
    d_complex** mp_FFT;


    double** BigSlice;
    double** BigjumpPart;

    double *jumpProb_slice;
    double *jumpPart;
    double *current_slice;

    int whichFFTGrid;
};

//typedef smartPtr<FD2DSolverJumps> FD2DSolverJumpsSP;
typedef refCountPtr<FD2DSolverJumps> FD2DSolverJumpsSP;


DRLIB_END_NAMESPACE
#endif
