/*
 * $Log: bootstrapcalibration.cpp,v $
 * Revision 1.14  2003/06/30 09:15:42  ebenhamou
 * remove unused var
 *
 * Revision 1.13  2002/11/19 11:28:36  mab
 * Formatting
 *
 * Revision 1.12  2002/09/17 16:48:13  mab
 * Improvements in double ComputeParasCurveDiag(...) (Yan Li)
 *
 * Revision 1.11  2002/07/08 07:56:44  mab
 * 0.5 replaced by 0.1 after optimizeNAG
 *
 * Revision 1.10  2002/07/02 10:39:51  mab
 * Yan LI add-on
 *
 * Revision 1.9  2001/08/14 23:29:28  smysona
 *  modif pour calage cap
 *
 * Revision 1.8  2001/07/30 08:54:09  smysona
 * amelioration robustesse pour les modeles frm
 *
 * Revision 1.7  2001/04/27 09:33:41  smysona
 *  Supression des sorties log
 *
 * Revision 1.6  2001/01/23 09:20:17  smysona
 * Code maudit : quit_nan -> K_HUGE_DOUBLE
 *
 * Revision 1.5  2001/01/23 09:06:42  mab
 * Comment code init dans bootstrap calibration
 *
 * Revision 1.4  2001/01/19 14:08:06  smysona
 * Destruction du filteredPF remontee dans le bootstrap
 *
 * Revision 1.3  2001/01/18 14:58:18  smysona
 * Rustine localisee pour FrmFamily
 *
 * Revision 1.2  2000/12/13 22:04:15  smysona
 * Dans ComputeParasCurve, gestion des cas ou il n'y a aucune nouvelle security dans le filtered portfolio pour un timestep supplementaire
 *
 * Revision 1.1  1998/12/23 17:19:29  ypilchen
 * Initial revision
 *
 */


#include "armglob.h"
#include "portfolio.h"
#include "bootstrapcalibration.h"
#include "limits"



double ComputeParasCurve(ARM_Model* model,
                         ARM_Portfolio* pf,
                         long nbPoints, 
                         double* matCurve,              
                         double precision,
                         double min_para,
                         double max_para,
                         long max_iters,
                         double* vals_curve,
                         int FilterType)
{    
    double square_sum;
    bool frmFamily = IsFrmFamily(model->GetName());
    int size = 0;
    int j,k;

   
    if (frmFamily)
    {
       for (j = 0; j < nbPoints ; j++)
       {
           vals_curve[j] = K_HUGE_DOUBLE;
       }
    }


    for (j = 0; j < nbPoints; j++)
    {
        ARM_Portfolio* filteredPF = pf->TemporalFilter(matCurve[j], FilterType);

        // l'index j controle le nbre d'elements exactement utilises

        if ( filteredPF == (ARM_Portfolio *) NULL )
        {
           if (frmFamily)
              continue;
           else
           {
              char buf[200];

              sprintf(buf, "Filtering has given no portfolio for matu : %lf", 
                      matCurve[j]);

              throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, buf);
            }
        }

        if ( filteredPF->GetSize() > size )
        {
           MinimizationOnCurveForBootStrapping bootStrapMinFunc(j,
                                                                vals_curve,
                                                                matCurve,
                                                                model,
                                                                filteredPF);
                                                              
                                                              
           double mid_para = 0.5*(min_para+max_para);    
    
           BrentMinimization minimizator(bootStrapMinFunc, 
                                         min_para, max_para, mid_para, 
                                         precision, max_iters);

           vals_curve[j] = minimizator.Minimize();

           // back propagate to non calibrated maturities (cf SetSigma)

           for (k = j-1; ( k >= 0 ) && ( vals_curve[k] == K_HUGE_DOUBLE );
               vals_curve[k--] = vals_curve[j]);

           square_sum = minimizator.GetMinimumValueOfTheFunction();

           size = filteredPF->GetSize();
        }

        cleanit(filteredPF);
    }

    return(square_sum);
}




/*******************************************************
// Fonction qui permet de calibrer une fonction de vols
// en donnant un intervalle pour a et pour les sigmas
// Elle renvoie la valeur de a
/******************************************************/

double ComputeMeanRevAndParasCurve(ARM_Model* model,
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
                                   double* vals_curve)
{
    MinimizationFunctionOnMeanRev function_of_meanRev(model,
                                                      pf,
                                                      nbPoints,
                                                      matCurve,
                                                      precision_para,
                                                      min_para,
                                                      max_para,
                                                      max_iters);

 
    double mid_meanRev = 0.5*(min_meanRev+max_meanRev);

    BrentMinimization2 meanRev_minimizator(function_of_meanRev, 
                                           min_meanRev, max_meanRev, 
                                           mid_meanRev, 
                                           precision_meanRev, max_iters);

    double optimal_meanrev = meanRev_minimizator.Minimize();

    ((ARM_GYCSigVarModel *) model)->SetMeanRevSpeed(optimal_meanrev);

    double res = ComputeParasCurve(model,
                                   pf,
                                   nbPoints,
                                   matCurve,
                                   precision_para,
                                   min_para,
                                   max_para,
                                   max_iters,
                                   vals_curve);

    return(optimal_meanrev);
}




// STATIC DECLARATION !!!!!
static MinimizationOnMultiAssets* multiAssetsMinFunc;



double ComputeParasCurveDiag(ARM_Model* model,
                             ARM_Portfolio* pf,
                             long nbPoints, 
                             double* matCurve,              
                             double precision,
                             double min_para,
                             double max_para,
                             long max_iters,
                             double* vals_curve,
                             int FilterType)
{    
    double square_sum;

    bool frmFamily = IsFrmFamily(model->GetName());

    int size = 0;
    int i,j, *nbPoint_array, nbFit=0;

    nbPoint_array = new int[nbPoints];;
   
    if (frmFamily)
    {
       for (j = 0; j < nbPoints ; j++)
       {
           nbPoint_array[j] = 0;

           vals_curve[j] = K_HUGE_DOUBLE;
       }
    }

    ARM_Portfolio* filteredPF = NULL;;

    for (j = 0; j < nbPoints; j++)
    {
        cleanit(filteredPF);
       
        filteredPF = pf->TemporalFilter(matCurve[j], FilterType);

        // l'index j controle le nbre d'elements exactement utilises

        if ( filteredPF == (ARM_Portfolio *) NULL )
        {
           if (frmFamily)
              continue;
           else
           {
              char buf[200];

              sprintf(buf, "Filtering has given no portfolio for matu : %lf", 
                      matCurve[j]);

              throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, buf);
           }
        }
    
        if ( filteredPF->GetSize() > size )
        {
           nbPoint_array[nbFit] = j+1;

           size = filteredPF->GetSize();

           nbFit ++;
        }
       
    }

/* powell optimization

    double* p, fret, **xi;

    p =  new double  [nbFit];
    xi = new double* [nbFit];

    for (i = 0; i < nbFit;i++)
    {
        xi[i] = new double [nbFit];

        for (j = 0; j < nbFit; j++)
        {
            xi[i][j] = (i==j?1.0:0.0);
        }
    }

    double mid_para=0.5*(min_para+max_para); 

    for (j = 0; j < nbFit; j++)
        p[j] = mid_para;

    optimize minimizator;

    multiAssetsMinFunc = new MinimizationOnMultiAssets(nbFit, nbFit,
                                                       nbPoints,
                                                       nbPoint_array,
                                                       p,
                                                       vals_curve,
                                                       matCurve,
                                                       model,
                                                       filteredPF);

    minimizator.powell(p, xi, nbFit, precision, &max_iters,
                       &square_sum, global_powell_function);

    multiAssetsMinFunc->GetCurve(vals_curve);


    for (i = 0; i < nbFit; i++)
        delete [] xi[i];

    delete xi;
*/



    // start of NAG

    double* p; 
    double* y;

    p = new double[nbFit];
    y = new double[nbFit];

    double mid_para = 0.25*(min_para+max_para); 

    for (j = 0; j < nbFit; j++)
        p[j] = mid_para;

    multiAssetsMinFunc = new MinimizationOnMultiAssets(nbFit,
                                                       nbFit,
                                                       nbPoints,
                                                       nbPoint_array,
                                                       p,
                                                       vals_curve,
                                                       matCurve,
                                                       model,
                                                       filteredPF);

    multiAssetsMinFunc->optimizeNAG(p, y);

    int ii = 0, ind = 1;

    while ( ii < nbFit && ind )
    {
        if ( p[ii] <= 1.e-5 ) 
           ii++;
        else 
           ind = 0;
    }

    // Modif Yan 14/06/2002
/* OLD Code
    if (p[ii+1]-p[ii]>0.5)
        p[ii]=p[ii+1];
*/

    if ( p[ii+1]-p[ii] > 0.1 )
       p[ii] = p[ii+1];

    // Fin Modif

    int iim = MAX(5, ii);

    for (i = 0; i < iim; i++)
        p[i] = p[iim];

    ii = 0;
    while ( ii < nbFit-1 )
    {
        if ( p[ii] <= 1.e-5 ) 
        {
           if ( p[ii+1] > 1.0e-5 )
              p[ii] = p[ii+1];
           else if ( ii > 0 )
              p[ii] = p[ii-1];
           else
              p[ii] = 1.5;
        }

        ii++;
    }

    multiAssetsMinFunc->minfun(p, y, square_sum);

    cleanit(y);

// end of NAG

    cleanit(filteredPF);
    
    cleanit(p);

    cleanit(multiAssetsMinFunc);

    cleanit(nbPoint_array);

    return(square_sum);
}



double global_powell_function(double x[])
{
    double value;

    value = (*multiAssetsMinFunc)(x);

    return value;
}



#ifdef WIN32

void NAG_CALL global_nag_function(long M, long N, double X[], 
                                  double F[], double fjac[], 
                                  long tdfjac, Nag_Comm* comm)
{
    multiAssetsMinFunc->nag_function(M, N, X, F, fjac, comm);
}

#else

void global_nag_function(int& MODE, int& M, int& N, int& LDFJ,
                         double* X,
                         double* F, double* fjac,
                         int& NSTATE,
                         int* IUSER, double* USER)
{
    multiAssetsMinFunc->nag_function(MODE, M, N, X, F, fjac);
}

#endif


#ifdef WIN32

void NAG_CALL confun(long n, long ncnlin, long needc[], double x[], 
                     double conf[], double cjac[],
                     Nag_Comm* comm)
{
}

#else
void confun(long n, long ncnlin, long needc[], double x[],
            double conf[], double cjac[])
{
}

#endif





/*--------------------------------------------------------------*/
/*---- End Of File ----*/
