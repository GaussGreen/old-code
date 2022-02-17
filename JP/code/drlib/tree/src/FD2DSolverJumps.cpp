
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DSolverJumps.hpp
//
//   Description : two factor finite difference algorithm
//
//   Author      : Ning Shen
//                 Xiaolan Zhang
//
//   Date        : November 29, 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD2DSolverJumps.hpp"
#include <time.h>
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

/************ for debuging tree step outputs **************/
/************ for debuging tree step outputs **************/
//-----------------------------------------------------------------------------

/** FD2DSolverJumps algorithm class that supports FD2DSolverJumps */
//FD2DSolverJumps::FD2DSolverJumps(FD2D* model, int xNinput, int yNinput): FD2DSolver(model){
FD2DSolverJumps::FD2DSolverJumps(FD2D* model, int xNinput, int yNinput):engine(model){

    int i;
    int j;
  
    xN = xNinput;
    yN = yNinput;

    xExtend = 2 * xN - 1;
    yExtend = 2 * yN - 1;

    xMiddle = (int)(xN -1) / 2;
    yMiddle = (int)(yN -1) / 2;

    xF.resize(xExtend); //big grid for FFT, but constructed such that it includes the same points as FD grid
    yF.resize(yExtend);

    g_complex = new d_complex[xExtend*yExtend];  //jump density
    g_FFT = new d_complex[xExtend*yExtend];  
    p_complex = new d_complex[xExtend*yExtend];  
    p_FFT = new d_complex[xExtend*yExtend];  

    BigSlice = new double* [xExtend];
    BigjumpPart = new double* [xExtend];
    for (i = 0; i < xExtend; i++){
        BigSlice[i] = new double [yExtend];
        BigjumpPart[i] = new double [yExtend];
    }

    jumpProb_slice = new double[xN * yN];
    jumpPart = new double[xN * yN]; 

    current_slice  = new double[xN * yN];

    for ( i = 0; i < xN; i++) {
        for ( j = 0; j < yN; j++){
            int ij = Pos2Index(i,j, xN, yN);
            current_slice[ij] = 0;
            jumpProb_slice[ij] = 0;
            jumpPart[ij] = 0;
        }
    }
}

//-----------------------------------------------------------------------------

FD2DSolverJumps::~FD2DSolverJumps(){

    delete [] g_complex; 
    delete [] g_FFT;  
    
    delete [] p_complex; 
    delete [] p_FFT;  

    int i;

    if (BigSlice !=0) {
        for (i = 0; i < xExtend; i++){
            delete [] BigSlice[i] ;
        }
        delete BigSlice;
    }
    BigSlice =0;

    if (BigjumpPart != 0) {
        for (i = 0; i < xExtend; i++){
            delete [] BigjumpPart[i] ;
        }
        delete BigjumpPart;
    }
    BigjumpPart = 0;
    
    delete[] jumpProb_slice;
    delete[] jumpPart;
    delete[] current_slice;
}


//-----------------------------------------------------------------------------
//Jumps parts
//-----------------------------------------------------------------------------    

//calculate jumps part contributing to the source term
//for now, to change in the future
void FD2DSolverJumps::compute_source_jumps(double* M_price, double* rhsJumps, 
                                           double dx, double dy){
    

    FD2DSVCJ* modelSVCJ = dynamic_cast< FD2DSVCJ*>(engine);

    int idx, idy;

    //int whichMethod = 2;

    whichFFTGrid = 1;
    int whichMethod = modelSVCJ->whichM;

    double opt ;

    switch (whichMethod) {
    case 1: {

        //only calculate one time FFT for jump density
        //only for tests
        //should be moved to FD2DSolver setup since we'll only need to call once
        //set grid for FFT, and calc FFT for density
        setup(dx, dy);

        //convert small matrix price to big matrix price        
        for (idx = 0; idx < xExtend; idx ++ ){
            for (idy = 0; idy < yExtend; idy ++ ){
                BigSlice[idx][idy] = 0.0;
                BigjumpPart[idx][idy] =0.0;
            }
        }

        if (whichFFTGrid ==1){
            for (idx = xMiddle; idx < xMiddle + xN; idx ++ ){
                for (idy = yMiddle; idy < yMiddle + yN; idy ++ ){
                    int old_idx = idx - xMiddle;
                    int old_idy = idy - yMiddle;
                    BigSlice[idx][idy] = M_price[Pos2Index(old_idx, old_idy, xN, yN)];
                }
            }

            //do linear interpolation outside?
            doInterp(BigSlice, xMiddle, yMiddle) ;
        }else{
        //interpolation to get price at new grid
            int dimX = engine->gridLevel1.size();
            int dimY = engine->gridLevel1.size();

            for (idx = 0; idx < xExtend; idx ++ ){
                for (idy = 0; idy < yExtend; idy ++ ){
                    BigSlice[idx][idy]= interpF2_linear(xF[idx], yF[idy], &*engine->gridLevel1.begin(), &*engine->gridLevel2.begin(), M_price, dimX, dimY);
                }
            }
        }

        //assumed that we have already calc g_FFT
        JumpCompQuick(BigSlice, g_FFT, BigjumpPart);

        //convert the bigjumpparts to smalll jumps parts

        FD2DSVCJ* modelSVCJ = dynamic_cast< FD2DSVCJ*>(engine);
//        opt = dt * modelSVCJ->vol->commonCrashRate;
        opt = -modelSVCJ->cCrashRatedt;

        for (idx = xMiddle; idx < xMiddle + xN; idx ++ ){
            for (idy = yMiddle; idy < yMiddle + yN; idy ++ ){
                int old_idx = idx - xMiddle;
                int old_idy = idy - yMiddle;
                int ij = Pos2Index(old_idx, old_idy, xN, yN);

                double p;
                if (whichFFTGrid !=1){
                    //need to do interpolation to get price at the original FD grid
                    double x = engine->gridLevel1[old_idx];
                    double y = engine->gridLevel2[old_idy];

                    p = interpF2_linear(x, y, xF, yF, BigjumpPart);
                }else{
                    p =BigjumpPart[idx][idy];
                }

                double temp = opt * (p  - M_price[ij]);

                //double temp = opt * (BigjumpPart[idx][idy]  - M_price[ij]);

                rhsJumps[ij] += temp ;           
            }
        }
    }
    break;
    case 2:{
        //do 2 FFT and 1 inverse at each time step
        //but, the working area is narrowed at [(xN-1)/2, (xN-1)/2 + xN] and [(yN-1)/2, (yN-1)/2 + yN], 
        //probably, not great!

        //copy  M_price to current slice
        for (idx = 0; idx < xN*yN; idx++){
            current_slice[idx] = M_price[idx];
        }

        jumpComp(current_slice, jumpPart, dx, dy);

        //should use the following lines, put in comments, due not compile
        opt = - modelSVCJ->cCrashRatedt;
       
        //Jump part is calculated now
        for (idx = 0; idx < xN*yN; idx++){
            rhsJumps[idx] += opt * (jumpPart[idx]  - current_slice[idx]); ;
        }
    }
    break;    
    }
}

//-----------------------------------------------------------------------------    
void FD2DSolverJumps::doInterp(double** BigSlice, int xMiddle, int yMiddle) {

    int idx;
    int idy;

    //constant extropolation on V direction
    for (idy = 0; idy < yMiddle ; idy ++ ){
        for (idx = xMiddle; idx < xMiddle +xN ; idx ++ ){
            BigSlice[idx][idy] =BigSlice[idx][yMiddle];
        }
    }

    //constant extropolation on V direction
    for (idy = yMiddle +yN ; idy < yExtend ; idy ++ ){
        for (idx = xMiddle; idx < xMiddle +xN ; idx ++ ){
            BigSlice[idx][idy] =BigSlice[idx][yMiddle +yN-1];
        }
    }

    //linear extropolation on Stock direction
    for (idy = 0; idy < yExtend ; idy ++ ){
        for (idx = 0; idx < xMiddle ; idx ++ ){
            double slope = (BigSlice[xMiddle+1][idy] -BigSlice[xMiddle][idy])/(xF[xMiddle+1] - xF[xMiddle]);
            BigSlice[idx][idy] = BigSlice[xMiddle][idy] + slope * (xF[xMiddle] - xF[idx]);
        }
    }

    //linear extropolation on Stock direction
    for (idy = 0; idy < yExtend ; idy ++ ){
        for (idx = xMiddle + xN; idx < xExtend; idx ++ ){
            double slope = (BigSlice[xMiddle + xN-1][idy] -BigSlice[xMiddle + xN-2][idy])/(xF[xMiddle + xN-1] - xF[xMiddle + xN-2]);
            BigSlice[idx][idy] = BigSlice[xMiddle + xN-1][idy] + slope * (xF[xMiddle + xN-1] - xF[idx]);
        }
    }
}

//-----------------------------------------------------------------------------    

// probability function, the approximation of the integration over the pdf
double FD2DSolverJumps::jumpPro(double x, double y, double dx, double dy){
    FD2DSVCJ* modelSVCJ = dynamic_cast< FD2DSVCJ*>(engine);

    return 0.25 * dx * dy * (modelSVCJ->jumpProDen(x,y + 0.5* dy) 
                + modelSVCJ->jumpProDen(x, y - 0.5 * dy) 
                + modelSVCJ->jumpProDen(x - 0.5 * dx, y) 
                + modelSVCJ->jumpProDen(x + 0.5 * dx, y));
}

//-----------------------------------------------------------------------------    

int FD2DSolverJumps::Pos2Index(int idx, int idy, int height, int width){
    //return idx * width + idy;
    return idx + idy * height;
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void FD2DSolverJumps::JumpCompQuick(double** BigSlice, d_complex* FFT_den, 
                                    double** BigjumpPart){
    // Implement the jump part using the FFT.

    //fill p_complex;
    //do mult of p_complex and g_complex
    //do inverse of the above product

    int idn, idm = 0 ;  // the index of array
    
    //for now, assuming that we have a big matrix for the slice
    for (idn = 0; idn < xExtend; idn++){
        for (idm = 0; idm < yExtend; idm++){
            int k = Pos2Index(idn, idm, xExtend, yExtend);
            p_complex[k].re = BigSlice[idn][idm];
            p_complex[k].im = 0;
        }
    }
    
    //p_FFT = imsl_z_fft_2d_complex (xExtend, yExtend, p_complex, 0);
    imsl_z_fft_2d_complex (xExtend, yExtend, p_complex, IMSL_RETURN_USER, p_FFT, 0);

    double reTemp, imTemp;

    for (idn = 0; idn < xExtend; idn++){
        for (idm = 0; idm < yExtend; idm++){
            int k = Pos2Index(idn, idm, xExtend, yExtend);
            reTemp  = p_FFT[k].re * FFT_den[k].re
                    - p_FFT[k].im * FFT_den[k].im;

            imTemp  = p_FFT[k].re * FFT_den[k].im
                    + p_FFT[k].im * FFT_den[k].re;
            p_FFT[k].re = reTemp;
            p_FFT[k].im = imTemp;
        }
    }

    //here p_complex contains the convolution of the price and the jump den.
    //it's the jump comp
    //p_complex = imsl_z_fft_2d_complex(xExtend, yExtend, p_FFT, IMSL_BACKWARD, 0);
    imsl_z_fft_2d_complex(xExtend, yExtend, p_FFT, IMSL_BACKWARD, 
                                IMSL_RETURN_USER, p_complex, 0);

    
    for (idn = 0; idn < xExtend; idn++){
        for (idm = 0; idm < yExtend; idm++){
            int k = Pos2Index(idn, idm, xExtend, yExtend);

            BigjumpPart[idn][idm] = p_complex[k].re /(xExtend * yExtend);
        }
    }

    //ask??????
    //p_FFT, and p_complex are declare global in the class, 
    //and destroy in the destructor of this class,
    // imsl ft return a point, but, I used p_FFT, p_complex,
    //is that ok at memory point of view??????????

}

//-----------------------------------------------------------------------------
//this ft need to be called before loop on time.
void FD2DSolverJumps::setup(double dx, double dy){
    //only called one time
    //at the begginnng
    constructFFTGrid(dx, dy);
    FFTJumpDen(dx, dy);
}

//-----------------------------------------------------------------------------

void FD2DSolverJumps::constructFFTGrid(double dx, double dy){

    int idx, idy;

    double xCenter = engine->gridLevel1[xMiddle];
    double yCenter = engine->gridLevel2[yMiddle];

    if (whichFFTGrid==1){
        //xF is divided by three interval
        //first [0, idxMiddle[,
        //second [idxMiddle, idxMiddle + xN -1] = GridLevel1
        //third [idxMiddle + xN, 2 * xN -1]

        //check the mid interval is correct or not,
        //in the future, to see if it worth to move to solver as a big GridLevel
        //xExtend = 2 * xN -1

        for (idx = 0; idx < xExtend; idx ++ ){
            xF[idx] = xCenter + (idx - (xN -1)) * dx; 
        }

        for (idy = 0; idy < yExtend; idy ++ ){
            yF[idy] = Maths::max(0.0, yCenter + (idy - (yN -1)) * dy); 
            //yF[idy] = idy * dy; 
        }


//    for (idx = idxMiddle; idx <= idxMiddle + xN -1; idx ++ ){
//        //xF[idx] = xCenter + (idx - (xN -1)) * dx; 
//        xF[idx] = engine->GridLevel1[idx - idxMiddle];
//    }
//
//    for (idy = idyMiddle; idy <= idyMiddle + yN -1; idy ++ ){
//        //xF[idx] = xCenter + (idx - (xN -1)) * dx; 
//        yF[idy] = engine->GridLevel1[idy - idyMiddle];
//    }

    }else{

        //variance should be positive
        for (idy = 0; idy < yExtend; idy ++ ){
            yF[idy] = Maths::max(0.0, yCenter + (idy - (yN -1)) * dy); 
            //yF[idy] = idy * dy; 
        }

        //test

        FD2DSVCJ* modelSVCJ = dynamic_cast< FD2DSVCJ*>(engine);

        double var = modelSVCJ->getJumeSizeVar();

        double maxMean = modelSVCJ->getJumpSizeMean(yF[yExtend]);
        double minMean = modelSVCJ->getJumpSizeMean(yF[0]);
        double xMax = maxMean + 3.0 * var;
        double xMin = minMean - 3.0 * var;

        for (idx = 0; idx < xExtend; idx ++ ){
    //       xF[idx] = xCenter + (idx - (xN -1)) * dx; 
            xF[idx] = xMin + idx  * dx; 
        }

    }
}

void FD2DSolverJumps::FFTJumpDen(double dx, double dy){

    double ftotal = 0.0;
    for (int idn = 0; idn < xExtend; idn++){
        for (int idm = 0; idm < yExtend; idm++){
            int k = Pos2Index(idn, idm, xExtend, yExtend);
            double pro = jumpPro(xF[idn],yF[idm],dx,dy);
            g_complex[k].re = pro;
            g_complex[k].im = 0;
            ftotal += pro;
        }
    }

    if (ftotal > 0.0) {
        ftotal  = ftotal ;
    }
    //cout <<"The probability is" <<ftotal <<"\n";
//    time_t start, end;
//    double dif;
//    time(&start);

    //g_FFT = imsl_z_fft_2d_complex (xExtend, yExtend, g_complex, 0);
    imsl_z_fft_2d_complex (xExtend, yExtend, g_complex, 
                                    IMSL_RETURN_USER, g_FFT, 0);

//    time(&end);
//    dif = difftime(end, start);
//    cout<<"the time it takes "<< dif << "seconds."<< "\n";    
}

//  computes the jump part. 
void FD2DSolverJumps::jumpComp(double* slice, double* jumpPart, double dx, double dy) {
    double ftotal=0;
//    int idxMiddle = (xN - 1) / 2;
//    int idyMiddle = (yN - 1) / 2;
    double xCenter = engine->gridLevel1[xMiddle];
    double yCenter = engine->gridLevel2[yMiddle];
    
    // Implement the jump part using the FFT.    
    int idx, idy ; // the index for the jump probability array

//    for (idx = 0; idx < xN; idx++){
//        for (idy = 0; idy < yN; idy++){
    for (idx = 0; idx < xExtend; idx++){
        for (idy = 0; idy < yExtend ; idy++){

            //assuming that dx and dy for the jumps parts are the same as FD grid
            //only work for equi dy now, 
            //need to chg if it's a variable dy

            double x = (double)(idx - xMiddle) * dx + xCenter;
            double y = (double)(idy - yMiddle) * dy + yCenter;

//            double x = (double)(idx ) * dx + xCenter;
//            double y = (double)(idy ) * dy + yCenter;


            jumpProb_slice[Pos2Index(idx,idy, xN, yN)] = jumpPro(x,y,dx,dy);
            ftotal += jumpPro(x, y, dx, dy);

            //for tests
            if (ftotal > 0.0) {
                ftotal = ftotal;
            }
        }
    }    

//    cout <<"The probability is" <<ftotal <<"\n";
//    time_t start, end;
//    double dif;
//    time(&start);

    Conv2(xN, yN, slice, jumpProb_slice, jumpPart);    

//    time(&end);
//    dif = difftime(end, start);
//    cout<<"the time it takes "<< dif << "seconds."<< "\n";    

}

//  convolution of 2 double array (n by m) 
//  return C = A con B
void FD2DSolverJumps::Conv2(int n, int m, double *A, double *B, double *C){
    // Implement the jump part using the FFT.

    int idn, idm = 0 ;  // the index of array
    int xExtend = 2 * n - 1;
    int yExtend = 2 * m - 1;

    int nMiddle, mMiddle;
    
    d_complex *A_complex, *B_complex;
    d_complex *A_FFT, *B_FFT;

    d_complex *C_complex, *C_FFT;

    A_complex = new d_complex[xExtend * yExtend];
    B_complex = new d_complex[xExtend * yExtend];

    A_FFT = new d_complex[xExtend * yExtend];
    B_FFT = new d_complex[xExtend * yExtend];
    
    C_complex = new d_complex[xExtend * yExtend];
    C_FFT = new d_complex[xExtend * yExtend];

    nMiddle = n / 2; 
    mMiddle = m / 2;    

//    nMiddle = (n -1) / 2; 
//    mMiddle = (m -1) / 2;    

//  initilization 
    for (idn = 0; idn < xExtend; idn++){
        for (idm = 0; idm < yExtend; idm++){
            int k = Pos2Index(idn, idm, xExtend, yExtend);
            A_complex[k].re = 0;
            A_complex[k].im = 0;

            B_complex[k].re = 0;
            B_complex[k].im = 0;
        }
    }

//  storing the values into the complex double array
    for (idn = 0; idn < n; idn++){
        for (idm = 0; idm < m; idm++){
            int k = Pos2Index(idn + nMiddle, idm + mMiddle, xExtend, yExtend);
//            int k = Pos2Index(idn , idm , xExtend, yExtend);

            int kk = Pos2Index(idn, idm, n, m);
            A_complex[k].re = A[kk];
            A_complex[k].im = 0;

            B_complex[k].re = B[kk];
            B_complex[k].im = 0;
        }
    }

//  convert the double real number into double complex number.
//  It is because imsl has only 2d fft for complex number, but not double.
    A_FFT = imsl_z_fft_2d_complex (xExtend, yExtend, A_complex, 0);
    B_FFT = imsl_z_fft_2d_complex (xExtend, yExtend, B_complex, 0);

//    imsl_z_fft_2d_complex (xExtend, yExtend, A_complex, IMSL_RETURN_USER, A_FFT, 0);
//    imsl_z_fft_2d_complex (xExtend, yExtend, B_complex, IMSL_RETURN_USER, B_FFT, 0);

    for (idn = 0; idn < xExtend; idn++){
        for (idm = 0; idm < yExtend; idm++){
            int k = Pos2Index(idn, idm, xExtend, yExtend);
            C_FFT[k].re 
                = A_FFT[k].re
                    * B_FFT[k].re
                    - A_FFT[k].im
                    * B_FFT[k].im;

            C_FFT[k].im 
                = A_FFT[k].re
                    * B_FFT[k].im
                    + A_FFT[k].im
                    * B_FFT[k].re;
        }
    }

    C_complex = imsl_z_fft_2d_complex(xExtend, yExtend, C_FFT, IMSL_BACKWARD, 0);
    //imsl_z_fft_2d_complex(xExtend, yExtend, C_FFT, IMSL_BACKWARD, 
    //                            IMSL_RETURN_USER, C_complex, 0);
    
    for (idn = 0; idn < n; idn++){
        for (idm = 0; idm < m; idm++) {
            int k = Pos2Index(idn, idm, n, m);
            int kk = Pos2Index(idn + nMiddle, idm + mMiddle, xExtend, yExtend);
            C[k] = C_complex[kk].re / (xExtend * yExtend);
        }
    }

    delete []A_complex;
    delete []A_FFT;
    delete []B_complex;
    delete []B_FFT;
    delete []C_complex;
    delete []C_FFT;
}

DRLIB_END_NAMESPACE
