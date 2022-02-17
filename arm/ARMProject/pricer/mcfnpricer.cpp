/*
 * $Log: mcfnpricer.cpp,v $
 * Revision 1.8  2003/10/16 10:38:54  rguillemot
 *  Double Barrier Pricing with Mixture
 *
 * Revision 1.7  2003/06/17 14:21:26  emezzine
 *  Add   SetModelAndVariables called in try
 *
 * Revision 1.6  2003/06/05 15:52:43  ykhlif
 * ajout de SetModelVariable() dans price()
 *
 * Revision 1.5  2003/05/26 09:59:13  emezzine
 *  modif Price()
 *
 * Revision 1.4  2001/04/23 09:25:02  smysona
 * unused variable
 *
 * Revision 1.3  1999/06/09 13:51:28  nicolasm
 * Ajout CleanModel apres pricing
 *
 * Revision 1.2  1999/06/07 14:57:54  nicolasm
 * Ajout ecriture variance dans un fichier
 *
 * Revision 1.1  1998/11/20 11:08:15  nicolasm
 * Initial revision
 *
 */
#include "mcfnpricer.h"
#include "mcmodel.h"
#include "ycmodel.h"
#include "bm-gen.h"




double ARM_FNMonteCarloPricer::Price(void)
{
    itsSecurity->PrepareToPrice(itsModel->GetStartDate());
    
    ARM_MCModel *mcmod = dynamic_cast<ARM_MCModel *> (itsModel);  

    ARM_YCModel mod(mcmod->GetZeroCurve());

    itsSecurity->SetModelVariable(NULL);
    itsSecurity->SetModel(&mod);

    //itsSecurity->SetModel(mcmod);

    mcmod->BeFittedTo(itsSecurity);

    ARM_PathGenerator *gen = NULL;


    gen = mcmod->Generator();

    double moment1=0.0, moment2=0.0;
    double result;
    long n0 = SUMMATION_INIT_VALUE;
    long c = SUMMATION_GROWTH_FACTOR;
    long n=0;

    int nbIter = mcmod->NbIters();

    do
    {        
        n++;

        result=mcmod->PayoffFN(gen->Generate(), itsSecurity);        

        moment1 +=result;

        moment2 +=SQR(result);
    }                     
    while (n < nbIter);  

    double diff = fabs((moment2/n - moment1*moment1/(n*n)));
    double sdev = sqrt(diff/n);

  //   SetEstimatedSdev(sdev);
   // Ecriture de la variance dans un fichier

/*    char* homeDir;
    char  fOutName[100];
    FILE* mclog;

    homeDir = getenv("HOME");

    if (homeDir)
    {
       sprintf(fOutName, "%s/%s", homeDir, "MCLOG");
 
       mclog = fopen(fOutName, "a+");
       fprintf(mclog, "Nbre Traj : %d	Std Dev : %lf \n", nbIter, sdev);
       fclose(mclog);
    }
*/

    mcmod->CleanModel();
    itsSecurity->SetModelVariable(NULL);
    

   
 
   return moment1/n;
}

