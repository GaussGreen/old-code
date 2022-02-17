/*
 * $Log: mcrnpricer.cpp,v $
 * Revision 1.3  2000/11/20 09:44:01  mab
 * *** empty log message ***
 *
 * Revision 1.2  2000/11/20 09:35:08  mab
 * *** empty log message ***
 *
 */

#include "mcrnpricer.h"

#include "mcmodel.h"
#include "bm-gen.h"



double ARM_RNMonteCarloPricer::Price(void)
{
    itsSecurity->PrepareToPrice(itsModel->GetStartDate());

    ARM_MCModel *mcmod = dynamic_cast<ARM_MCModel *> (itsModel);    

    mcmod->BeFittedTo(itsSecurity);


    ARM_PathGenerator *gen = NULL;

    gen = mcmod->Generator();


    double moment1=0.0, moment2=0.0;
    double result;

    long n=0;


    do
    {        
        n++;


        result=itsSecurity->PayoffRN(gen->Generate(), mcmod);

        moment1 +=result;
        moment2 +=SQR(result);
    }                     
    while (n < mcmod->NbIters());


    double diff = fabs((moment2/n - moment1*moment1/(n*n)));
    double sdev = sqrt(diff/n);

    SetEstimatedSdev(sdev);

    return (moment1/n);
}

