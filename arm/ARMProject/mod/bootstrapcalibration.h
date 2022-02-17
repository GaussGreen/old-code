/*
 * $Log: bootstrapcalibration.h,v $
 * Revision 1.18  2003/06/23 12:15:21  ebenhamou
 * use ournag.h
 *
 * Revision 1.17  2003/02/11 16:46:35  mab
 * Improvements in NAg Call
 *
 * Revision 1.16  2002/11/19 11:43:31  mab
 * Formatting
 *
 * Revision 1.15  2002/09/17 16:45:40  mab
 * Improvements in Nag numerical optimization
 *
 * Revision 1.14  2002/07/02 10:40:45  mab
 * Yan LI add-on
 *
 * Revision 1.13  2002/03/01 16:32:19  mab
 * MODIF : YLI Corrections
 *
 * Revision 1.12  2001/08/14 23:29:21  smysona
 *  modif pour calage cap
 *
 * Revision 1.11  2001/04/27 09:33:32  smysona
 * Supression des sorties log
 *
 * Revision 1.10  2001/04/23 09:14:50  smysona
 * ajout #ifdef ...
 *
 * Revision 1.9  2001/02/22 19:46:51  smysona
 * Smart SetSigma pour FrmAna
 *
 * Revision 1.8  2001/01/23 09:07:31  mab
 * Rajout printf
 *
 * Revision 1.7  2001/01/19 14:08:38  smysona
 * Destruction du filteredPF remontee dans le bootstrap
 *
 * Revision 1.6  2000/11/10 14:30:16  smysona
 * Error correction
 *
 * Revision 1.5  2000/11/09 20:39:41  smysona
 * FRMANA case added in SetDataCurveToTheModel
 *
 * Revision 1.4  2000/10/25 16:28:15  ypilchen
 * Modifs ARM_HWSigVarTwoFactorModel SetSigma1 modifie par SetNewBasisSigma1
 *
 * Revision 1.3  1999/10/27 15:22:05  sgasquet
 * Modifs SetDataCurveToTheModel
 *
 * Revision 1.2  1999/05/07 14:44:58  nicolasm
 * Recherche des mktPrices et des poids dans le portfeuille et non
 * plus dans l'asset
 *
 * Revision 1.1  1998/12/23 17:19:00  ypilchen
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/
#ifndef _BOOTSTRAPCALIBRATION_H
#define _BOOTSTRAPCALIBRATION_H




#include "ipricer.h"
#include "x1fhwsigvar.h"
#include "x2fhwsigvar.h"
#include "portfolio.h"
#include "model.h"

#include "minimization.h"
#include "function.h"    
#include "frmana.h"

#ifdef WIN32

#include "ournag.h"
#include "nage04.h"
#include "nag_stdlib.h"

#endif




/*----------------------------------------------------------------------------*

  Minimization function for a model which as as an input a curve whose points
  are computed by bootstrapping
  Inputs for this function:
  - the calibration portdolio
  - the model to be calibrated
  - the curve maturities ( itsMatCurve en double )
  - the curve to be computed itsCurve
  This function computes the following value of the curve.
  If itsCurve has the same size as itsMatCurve, 
   the curve as been totally computed.

  Otherwise, the new point to be computed corresponds 
  to the maturity n if itsCurve contains n-1 plots.

*----------------------------------------------------------------------------*/

double global_powell_function(double []);



typedef void (*E04UDM)(int&, int&, int&, int&, int*, double*,
                       double*, double*,
                       int&, int*, double*);

typedef void (*GLOBAL_NAG_FUNCTION)(int& MODE, int& M,
                                    int& N, int& LDFJ, double* X,
                                    double* F, double* FJAC, int& NSTATE,
                                    int* IUSER, double* USER);

typedef void (*GLOBAL_NAG0_FUNCTION)(int& MODE, int& N,
                                     double* X, double* F, double* FJAC,
                                     int& NSTATE,
                                     int* IUSER, double* USER);

#ifdef sun

extern "C"
{
extern void e04udm_(int&, int&, int&, int&, int*, double*, double*,
                    double*,
                    int&, int*, double*);

extern void e04unf_(int&, int&, int&, int&, int&, int&, int&,
                    int&, double*, double*, double*,
                    double* Y,
                    E04UDM, GLOBAL_NAG_FUNCTION,
                    int&,
                    int*, double*, double*,
                    double*, double*, double*,
                    double&, double*, double*, int*, int&, double*,
                    int&, int*, double*, int&);

extern void x04abf_(int&, int&);
}

#endif






#ifdef WIN32
void NAG_CALL global_nag_function(long M, long N, double X[], double F[],
                                  double fjac[],
                                  long tdfjac, Nag_Comm* comm);

void NAG_CALL confun(long n, long ncnlin, long needc[], double x[],
                     double conf[], double cjac[],
                     Nag_Comm* comm);
#else
    void global_nag_function(int& MODE, int& M, int& N, int& LDFJ,
                             double* X, double* F, double* FJAC,
                             int& NSTATE,
                             int* IUSER, double* USER);

    void confun(long n, long ncnlin, long needc[], double x[],
                double conf[], double cjac[]);
#endif






/*--------------------------------------------------------------*/



class MinimizationOnMultiAssets : public RealValuedFunction
{
    private:

         int itsNbPoints;

         int itsNbFit;

         int itsNbSec;

         int* nbPoint_array;

         double* imsl_x;

         double* imsl_y;

         ARM_Vector* itsCurve;    // part of the curve already computed

         ARM_Vector* itsMatCurve;

         ARM_Model* itsModel;                    

         ARM_Portfolio* itsCalibrationPortfolio; // portfolio which has been 
                                                 // fitted to the maturity

    public:

        void Init(void)
        {
            itsNbPoints = 0;
            itsNbFit    = 0;
            itsNbSec    = 0;
            nbPoint_array = NULL;
            itsCurve      = NULL;
            itsMatCurve   = NULL;
            itsModel      = NULL;
            imsl_x        = NULL;
            imsl_y        = NULL;

            itsCalibrationPortfolio = NULL;
        }

 
        MinimizationOnMultiAssets(void)
        {
            Init();
        }

        void GetCurve(double* vals_curve_)
        {
            for (int i = 0; i < itsNbPoints; i++)
                vals_curve_[i] = (*itsCurve)[i];
        }

        MinimizationOnMultiAssets(long nb_fit_,
                                  long nb_secs_,
                                  long nb_points_,
                                  int * nb_point_array_,
                                  double* p_,
                                  double* vals_curve_, 
                                  double* mats_curve_,
                                  ARM_Model* model_,
                                  ARM_Portfolio* calibration_portfolio_)
        {
            Init();
 
            // Construction des vecteurs pour la courbe de volat.

            int i,j;
            itsNbPoints = nb_points_;
            itsNbFit = nb_fit_;
            itsNbSec = nb_secs_;
            nbPoint_array = nb_point_array_;

            if ( nb_points_ != 0 )
            {
               itsCurve    = new ARM_Vector(nb_points_);
               itsMatCurve = new ARM_Vector(nb_points_);
            }

            for ( i = 0; i < nb_points_; i++)
            {
                (*itsMatCurve)[i] = mats_curve_[i]; 
            }

            for(j=0;j<nbPoint_array[0];j++)
            {
                (*itsCurve)[j] = p_[0];
            }

            for (i = 0; i < itsNbFit-1; i++)
            {
                for (j = nbPoint_array[i]; j < nbPoint_array[i+1]; j++)
                {
                    (*itsCurve)[j] = p_[i];
                }
            }

            for (j = nbPoint_array[itsNbFit-1]; j < itsNbPoints; j++)
            {
                (*itsCurve)[j] = p_[itsNbFit-1];
            }

            itsCalibrationPortfolio = calibration_portfolio_; 

            itsModel = model_;
        }


       ~MinimizationOnMultiAssets(void)
        {
            delete itsCurve;
            delete itsMatCurve;
            delete imsl_x;
            delete imsl_y;
            //delete itsCalibrationPortfolio;
        }


    virtual double operator () ( double* p_ ) const
    {   
       int i,j;

       double* curve_in_process;

       double* mat_curve_in_process;


       if ( itsCurve == NULL )
       {
          throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,                
"MinimizationOnMultiAssets problem : no initialization or zero para. to optimized");
       }
       else
       { 
          // Rebuild the curve with the new point
          
          for (j = 0; j < itsNbFit; j++)
              p_[j] = MAX(p_[j], 0);

          curve_in_process = new double[itsCurve->GetSize()];

          mat_curve_in_process = new double[itsCurve->GetSize()];

          for (j = 0; j < nbPoint_array[0]; j++)
          {
              (*itsCurve)[j] = p_[0];
          }

          for (i = 0; i < itsNbFit-1; i++)
          {
              for (j = nbPoint_array[i]; j < nbPoint_array[i+1]; j++)
              {
                  (*itsCurve)[j] = p_[i];
              }
          }

          for (j = nbPoint_array[itsNbFit-1]; j < itsNbPoints; j++)
          {
              (*itsCurve)[j] = p_[itsNbFit-1];
          }

          for (i = 0; i < itsCurve->GetSize(); i++)
          {
              mat_curve_in_process[i]=(*itsMatCurve)[i];

              curve_in_process[i]=(*itsCurve)[i];
          }

       }

       double square_sum = 0.0;

       SetDataCurveToTheModel(itsNbPoints, mat_curve_in_process, 
                              curve_in_process);
    
       double sumMktPrice = 0.0;
       double sumtFitPrice = 0.0;

       for ( i = 0 ; i < itsCalibrationPortfolio->GetSize() ; i++)
       {
           ARM_IFPricer pricer(itsCalibrationPortfolio->GetAsset(i), itsModel);


           double theoPrice = pricer.Price();

           double marketPrice = 
                    itsCalibrationPortfolio->GetMktPrices()->Elt(i);

           double coeff = itsCalibrationPortfolio->GetWeights()->Elt(i);

           square_sum += Sqr(theoPrice-marketPrice)*coeff;

           sumMktPrice += marketPrice;
           sumtFitPrice += theoPrice;
       }    

       delete curve_in_process;
       delete mat_curve_in_process;

       return(square_sum);
    }


    void SetDataCurveToTheModel(long nbPlots, double* mat_curve_in_process,
                                double* curve_in_process) const
    {
        ARM_CLASS_NAME model_name = itsModel->GetName();

        switch (model_name)
        {
            case ARM_GYCSIGVARMODEL :

                ((ARM_GYCSigVarModel *) itsModel)->SetSigma(nbPlots, 
                                                      mat_curve_in_process,
                                                      curve_in_process); 
            break;

            case ARM_HWSIGVARTWOFACTORMODEL :
            
                ((ARM_HWSigVarTwoFactorModel *) itsModel)->SetSigma1(nbPlots, 
                                                   mat_curve_in_process,
                                                   curve_in_process);
            break;

            case ARM_FRM_ANA :

                (dynamic_cast <ARM_FrmAna *>(itsModel))->SetSigma(nbPlots,
                                                    mat_curve_in_process,
                                                    curve_in_process, true);
            break;

            default :
                   ;
        }
    }
    
#ifdef WIN32
    void nag_function(long m,
                      long n,
                      double* xx, 
                      double* ordinate,
                      double* fjac,
                      Nag_Comm* comm)
    {
        double objf;

        if ( comm->flag == 0 )
        {
           nag_function0(m, n, xx, ordinate, objf, comm);
        }
        else
        {
           nag_grad_function(m, n, xx, ordinate, fjac, objf, comm);
        }
    }

#else

   void nag_function(int MODE,
                     long m,
                     long n,
                     double* xx,
                     double* ordinate,
                     double* fjac)
    {
        double objf;

        if ( MODE == 0 )
        {
           nag_function0(m, n, xx, ordinate, objf);
        }
        else
        {
           nag_grad_function(m, n, xx, ordinate, fjac, objf);
        }
    }

#endif

    void optimizeNAG(double* x, double *f)
    {   
        long m, n, nclin, ncnlin, i, tda, tdfjac, nmax;

#ifdef WIN32

        static NagError fail;

#endif

        double *a = NULL, *y, *fjac, *xl, *xu, objf;


        m = itsNbFit;
        n = itsNbSec;
        ncnlin = 0;
        nclin = 0;

        tda = n;
        tdfjac = n;
        nmax = n + nclin + ncnlin;

        y = new double[n];
        fjac = new double[m*n];
        xl = new double[nmax];
        xu = new double[nmax];

        for (i = 0; i < m; i++)
        {
            y[i] = .0;
        }

        for (i = 0; i < n; i++)
        {
            xl[i] = 0.0;
            xu[i] = 1.5;
        }


#ifdef WIN32

       e04unc(m, n, nclin, ncnlin, NULL, tda, xl, xu,
              y, global_nag_function, confun, x, &objf, f, fjac,
              tdfjac, E04_DEFAULT, NAGCOMM_NULL, &fail);

#else

       int NCLIN, NCNLN, LDA, LDCJ, LDFJ, LDR, LIWORK, LWORK, ITER, IFAIL;
       int IUSER;
       int inM, N, *ISTATE, *IWORK;

       double A, C, CJAC, *F, *FJAC, *CLAMDA, *R, *WORK, USER;
       double* Y;

       inM = m; 
       N = n; 
       NCLIN = ncnlin;
       NCNLN = nclin;
       LDA = 1;
       LDCJ = 1;
       LDFJ = inM;
       LDR = N;
       LIWORK = 3*N+NCLIN + 2*NCNLN ;
       LWORK = 20*N+inM*(N+3);
       ITER = 50;

       ISTATE = new int[N+NCLIN+NCNLN];
       IWORK  = new int[LIWORK];

       Y = new double[inM];
       F = new double[inM];
       FJAC = new double[LDFJ*N];
       CLAMDA = new double[N+NCLIN+NCNLN];
       R = new double[LDR*N];
       WORK = new double[LWORK];


       int in, out;

       in = 1 ;
       out = -1;

       x04abf_(in, out);

       IFAIL = -1;
      
       for (int ii = 0; ii < inM; ii++)
       {
           Y[ii] = 0.0;
       }

       e04unf_(inM, N, NCLIN, NCNLN, LDA, LDCJ, LDFJ,
               LDR, &A, 
               xl, xu, 
               Y,
               e04udm_,
               global_nag_function,
               ITER, ISTATE, &C, &CJAC, F, FJAC, CLAMDA,
               *f, R, x,
               IWORK, LIWORK, WORK,
               LWORK, &IUSER, &USER, IFAIL);

       delete ISTATE;
       delete IWORK;

       delete F;
       delete FJAC;
       delete CLAMDA;
       delete R;
       delete WORK;
       delete Y;
#endif

       delete y;
       delete fjac;
       delete xl;
       delete xu;
    }

#ifdef WIN32

    void nag_function0(long m, long n, double* p_, double* ordinate, 
                       double& f, Nag_Comm* comm)

#else

    void nag_function0(long m, long n,
                       double* p_, double* ordinate,
                       double& f)

#endif

    {   
        int i;

        minfun(p_,ordinate, f);

        if (imsl_x) 
           delete imsl_x;

        imsl_x = new double[n];

        if (imsl_y) 
           delete imsl_y;

        imsl_y = new double[m];

        for (i = 0; i < n; i++)
            imsl_x[i] = p_[i];

        for (i = 0; i < m; i++)
            imsl_y[i] = ordinate[i];
    }

    void minfun(double* p_, double* ordinate, double& f)
    {
       int i,j;

       double* curve_in_process;

       double* mat_curve_in_process;


       if ( itsCurve == NULL )
       {
          throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
 "MinimizationOnMultiAssets problem : no initialization or zero para. to optimized");
       }
       else
       { 
          // Rebuild the curve with the new point
          
          curve_in_process = new double[itsCurve->GetSize()];

          mat_curve_in_process = new double[itsCurve->GetSize()];

          for (j = 0; j < nbPoint_array[0]; j++)
          {
              (*itsCurve)[j] = p_[0];
          }

          for (i = 0; i < itsNbFit-1; i++)
          {
              for (j = nbPoint_array[i]; j < nbPoint_array[i+1]; j++)
              {
                  (*itsCurve)[j] = p_[i];
              }
          }

          for (j = nbPoint_array[itsNbFit-1]; j < itsNbPoints; j++)
          {
              (*itsCurve)[j] = p_[itsNbFit-1];
          }

          for (i = 0; i < itsCurve->GetSize(); i++)
          {
              mat_curve_in_process[i]=(*itsMatCurve)[i];

              curve_in_process[i]=(*itsCurve)[i];
          }
       }

       double square_sum = 0.0;

       SetDataCurveToTheModel(itsNbPoints, mat_curve_in_process, 
                              curve_in_process);
    
       double sumMktPrice = 0.0;
       double sumtFitPrice = 0.0;
       f = 0;

       for (i = 0 ; i < itsCalibrationPortfolio->GetSize() ; i++)
       {
           ARM_IFPricer pricer(itsCalibrationPortfolio->GetAsset(i), itsModel);

           double theoPrice = pricer.Price();

           double marketPrice = 
                    itsCalibrationPortfolio->GetMktPrices()->Elt(i);

           double coeff = itsCalibrationPortfolio->GetWeights()->Elt(i);

           ordinate[i] = (theoPrice-marketPrice)*sqrt(coeff);

           f += coeff;

           sumMktPrice += marketPrice;
           sumtFitPrice += theoPrice;
       } 
       
       for (i = 0; i < itsCalibrationPortfolio->GetSize(); i++)
           ordinate[i] /= sqrt(f);

       f = 0.0;

       for ( i = 0 ; i < itsCalibrationPortfolio->GetSize(); i++)
           f += pow(ordinate[i], 2);

       f = sqrt(f);

       delete curve_in_process;
       delete mat_curve_in_process;
    }



#ifdef WIN32

    void nag_grad_function(long m, long n, 
                           double* xx, double* ordinate0,
                           double* g, double& f0, Nag_Comm* comm)

#else

    void nag_grad_function(long m, long n, 
                           double* xx, double* ordinate0,
                           double* g, double& f0)

#endif
    {
        double* ordinate = new double[m];
        double* xxx = new double[n];
        double* delta = new double[n];

        double eps = 1.e-3;

        int i = 0;
        int no_exist = 0;
        double sum;

        while (( i < n ) && ( no_exist == 0 ))
        {
            if (!imsl_x || imsl_x[i] != xx[i]) no_exist++;

            i++;
        }

        if (no_exist)
        {
#ifdef WIN32

           nag_function0(m, n, xx, ordinate0, f0, comm);

#else

           nag_function0(m, n, xx, ordinate0, f0);

#endif
        }

        for (int j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                xxx[i] = xx[i];
            }

            delta[j] = eps*MAX(fabs(xx[j]),1.0);
            xxx[j] = xx[j] + delta[j];

            minfun(xxx,ordinate,sum);

            for (i = 0; i < m; i++)
            {
#ifdef WIN32

                g[i*n+j] = (ordinate[i]-imsl_y[i])/delta[j];

#else

                g[j*m+i] = (ordinate[i]-imsl_y[i])/delta[j]; 

#endif
            }
        }

        delete xxx;
        delete ordinate;
        delete delta;
    }
};



class MinimizationOnCurveForBootStrapping : public RealValuedFunction
{
    private:

         ARM_Vector* itsCurve;    // part of the curve already computed

         ARM_Vector* itsMatCurve;

         ARM_Model* itsModel;                    

         ARM_Portfolio* itsCalibrationPortfolio; // portfolio which has been 
                                                 // fitted to the maturity

         double itsNewMatu;

    public:

        void Init(void)
        {
            itsCurve    = NULL;
            itsMatCurve = NULL;
            itsModel    = NULL;
            itsCalibrationPortfolio = NULL;

            itsNewMatu = 0.0;
        }

 
        MinimizationOnCurveForBootStrapping(void)
        {
            Init();
        }

        MinimizationOnCurveForBootStrapping(long nb_points_,
                                      double* vals_curve_, 
                                      double* mats_curve_,
                                      ARM_Model* model_,
                                      ARM_Portfolio* calibration_portfolio_)
        {
            Init();
 
            // Construction des vecteurs pour la courbe de volat.

            if ( nb_points_ != 0 )
            {
               itsCurve    = new ARM_Vector(nb_points_);

               itsMatCurve = new ARM_Vector(nb_points_);
            }

            for (long i = 0; i < nb_points_; i++)
            {
                (*itsCurve)[i] = vals_curve_[i]; 

                (*itsMatCurve)[i] = mats_curve_[i]; 
            }

            // Remember that mats_curve_ has at least nb_points_+1 plots

            itsNewMatu = mats_curve_[nb_points_];

            itsCalibrationPortfolio = calibration_portfolio_; 

            itsModel = model_;
        }


       ~MinimizationOnCurveForBootStrapping(void)
        {
            delete itsCurve;

            delete itsMatCurve;

            //delete itsCalibrationPortfolio;
        }


    virtual double operator () ( const double& new_sigma ) const
    {        
            double* curve_in_process;

            double* mat_curve_in_process;



            if ( itsCurve == NULL )
            {
               curve_in_process     = new double[1]; 
               mat_curve_in_process = new double[1]; 

               mat_curve_in_process[0] = itsNewMatu; 
 
               curve_in_process[0] = new_sigma;
            }
            else
            { 
               // Rebuild the curve with the new point

               curve_in_process = new double[itsCurve->GetSize()+1];

               mat_curve_in_process = new double[itsCurve->GetSize()+1];

               for (long i = 0; i < itsCurve->GetSize(); i++)
               {
                   mat_curve_in_process[i]=(*itsMatCurve)[i];

                   curve_in_process[i]=(*itsCurve)[i];
               }

               mat_curve_in_process[itsCurve->GetSize()] = itsNewMatu; 

               curve_in_process[itsCurve->GetSize()] = new_sigma;
            }

            double square_sum = 0.0;

            // We compute the prices of all the portfolio

            // ! function to be implemented in the models with a curve

            long nbPlots;

            if ( itsCurve == NULL )
            {
               nbPlots = 1;
            }
            else
            {
               nbPlots = itsCurve->GetSize()+1;
            } 

            SetDataCurveToTheModel(nbPlots, mat_curve_in_process, 
                                   curve_in_process);
    
            // YAN 02/2002    

            double sumMktPrice = 0.0;
            double sumtFitPrice = 0.0;

            // Fin YAN

            for (long i = 0; i < itsCalibrationPortfolio->GetSize(); i++)
            {
                ARM_IFPricer pricer(itsCalibrationPortfolio->GetAsset(i), 
                                    itsModel);


                double theoPrice = pricer.Price();

                double marketPrice = 
                         itsCalibrationPortfolio->GetMktPrices()->Elt(i);

                double coeff = itsCalibrationPortfolio->GetWeights()->Elt(i);

                square_sum += Sqr(theoPrice-marketPrice)*coeff;

                // YAN 02/2002

                sumMktPrice += marketPrice;
                sumtFitPrice += theoPrice;

                // Fin YAN
       }    

       delete curve_in_process;
       delete mat_curve_in_process;

       return(square_sum);
    }


    void SetDataCurveToTheModel(long nbPlots, double* mat_curve_in_process,
                                double* curve_in_process) const
    {
        ARM_CLASS_NAME model_name = itsModel->GetName();

        switch (model_name)
        {
            case ARM_GYCSIGVARMODEL :

               ((ARM_GYCSigVarModel *) itsModel)->SetSigma(nbPlots, 
                                                           mat_curve_in_process,
                                                           curve_in_process); 
            break;

            case ARM_HWSIGVARTWOFACTORMODEL :
            
                ((ARM_HWSigVarTwoFactorModel *) itsModel)->SetSigma1(nbPlots,
                                                       mat_curve_in_process,
                                                       curve_in_process);
            break;

            case ARM_FRM_ANA :

                (dynamic_cast <ARM_FrmAna *>(itsModel))->SetSigma(nbPlots,
                                                      mat_curve_in_process,
                                                      curve_in_process, true);
            break;

            default :
                ;
        }
    }
};



/*----------------------------------------------------------------------------*

    Creation d'une routine permettant le calage du modele de HW
    a sigma variable a partir de caps et/ou swaptions
    Cette routine renseigne la courbe de vol entree en inputs.
    On suppose que les instruments ont ete tries au prealable:
    - swaption en premier
    - cap/floor en second
    a chaque fois par ordre croissant de maturite.

*----------------------------------------------------------------------------*/

double ComputeParasCurve(ARM_Model* model,
                         ARM_Portfolio* pf,
                         long nbPoints,
                         double* matCurve,
                         double precision,
                         double min_para,
                         double max_para,
                         long max_iters,
                         double* vals_curve,
                         int type = 0);



extern double ComputeMeanRevAndParasCurve(ARM_Model* model,
                                          ARM_Portfolio* pf,
                                          long nbPoints,
                                          double* matCurve,
                                          double precision_meanRev,
                                          double precision_para,
                                          double min_meanRev,
                                          double min_para,
                                          double max_meanRev,
                                          double max_para,
                                          long max_iters,
                                          double* vals_curve);
        
extern double ComputeParasCurveDiag(ARM_Model* model,
                                    ARM_Portfolio* pf,
                                    long nbPoints,
                                    double* matCurve,
                                    double precision,
                                    double min_para,
                                    double max_para,
                                    long max_iters,
                                    double* vals_curve,
                                    int type = 0);


/*----------------------------------------------------------------------------*

   Fonction de minimisation de l'ecart des prix en fonction de a
   Permettra d'appliquer a cette fonction la routine brent    

*----------------------------------------------------------------------------*/

class MinimizationFunctionOnMeanRev: public RealValuedFunction
{
    private:

        ARM_Model* itsModel;
        double*    itsMatCurve;
        long       itsNbPoints;
 
        ARM_Portfolio* itsCalibrationPortfolio; // the whole portfolio

        double itsPrecision;
        double its_min_para;
        double its_max_para; 
        long   its_max_iters;
   

        void Init(void)
        {
            itsModel    = NULL;
            itsMatCurve = NULL;
            itsNbPoints = 0;
            itsCalibrationPortfolio = NULL;

            itsPrecision = 0.0;
            its_min_para = 0.0;
            its_max_para = 0.0;
            its_max_iters = 0;
        }

    public:

        MinimizationFunctionOnMeanRev(ARM_Model* model_,
                                      ARM_Portfolio* calibration_portfolio_,
                                      long nbPoints,
                                      double* matCurve,
                                      double precision_,
                                      double min_para,
                                      double max_para,
                                      long max_iters_)
       {
           itsModel = model_;
           itsMatCurve = new double[nbPoints];

           itsNbPoints = nbPoints;

           for (long i = 0; i < itsNbPoints; i++)
           {
               itsMatCurve[i] = matCurve[i];
           }

           itsCalibrationPortfolio = calibration_portfolio_;

           itsPrecision = precision_;
           its_min_para = min_para;
           its_max_para = max_para;
           its_max_iters = max_iters_; 
       }

      ~MinimizationFunctionOnMeanRev(void)
       {
           delete itsMatCurve;
       }

       virtual double operator () (const double& a) const
       {
           double result_a;

           double* vals_curve = new double[itsNbPoints];


           ((ARM_GYCSigVarModel *) itsModel)->SetMeanRevSpeed(a);
          
           result_a = ComputeParasCurve(itsModel,
                                        itsCalibrationPortfolio,
                                        itsNbPoints,
                                        itsMatCurve,
                                        itsPrecision,
                                        its_min_para,
                                        its_max_para,
                                        its_max_iters,
                                        vals_curve);

           delete vals_curve;

           return(result_a);
       }
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of File ----*/
