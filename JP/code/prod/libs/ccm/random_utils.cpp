#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <supercube/InvertCumNormal.h>
#include <supercube/SCException.h>
#include <supercube/GaussianRandomSequence.h>
#include <supercube/SobolSequence.h>
extern "C"
{
#include "student.h"
#include "random_utils.h"
#include "proba_utils.h"
#include "error2.h"
}

// -------------------------------------------------------------------------
// Gaussian Random Sequence
//
extern "C" int CreateGaussianRandomSequence(
    double *Sequence,        /* (O) Ptr to object pointer        */
	long    aSeed,           /* (I) Seed for gasdev and perm of partitions   */
    int     numberDimensions,/* (I) Number of dimensions                     */
	int     numberPaths)     /* (I) Number of paths to extract               */
{
    static char routine[] = "CreateGaussianRandomSequence";
    int status = FAILURE;
    int i=0;
    int j=0;
    GaussianRandomSequence *A = 0;
    char ErrorMsg[MAXBUFF];

	try
    {    
        double *deviates = new double[numberDimensions];
        A = new GaussianRandomSequence;
        
		A->initializeObject(aSeed,
			                numberDimensions);

        for(j=0;j<numberPaths;j++)
        {
            A->populateVector();
            deviates = A->getVector();
            for(i=0;i<numberDimensions;i++)
            {
                Sequence[i+j*numberDimensions] = deviates[i];
            }
        }
	}catch(SCException &e){
        strcpy(ErrorMsg, e.what());
        goto RETURN;
    }

    delete A;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return (status);

} 

// -------------------------------------------------------------------------
// Sobol Sequence
//
extern "C" int  CreateSobolSequence(
    double  *sequence,     /* (O) Ptr to object pointer                    */
	long    aSeed,          /* (I) Seed for gasdev and perm of partitions   */
    int     numberDimensions,/* (I) Number of dimensions                    */
	int     numberPaths)    /* (I) Number of paths to extract               */
{
    static char routine[] = "CreateSobolSequence";
    int status = FAILURE;
    int i=0;
    int j=0;
    SobolSequence *A = 0;
    char errorMsg[MAXBUFF];

	try
    {
        double *deviates = new double[numberDimensions];
        A = new SobolSequence(aSeed,
			                numberDimensions);

        for(j=0;j<numberPaths;j++)
        {
            A->populateVector();
            deviates = A->getVector();
            for(i=0;i<numberDimensions;i++)
            {
                sequence[i+j*numberDimensions] = deviates[i];
            }
        }
  

	}catch(SCException &e){
        strcpy(errorMsg, e.what());
        goto RETURN;
    }

    delete A;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return (status);
}

// --------------------------------------------------------------------------
// CreateUniformSequence
//
extern "C" int  CreateUniformSequence(
    double  *sequence,     /* (O) Ptr to object pointer                    */
	long    aSeed,          /* (I) Seed for gasdev and perm of partitions   */
    int     numberDimensions,/* (I) Number of dimensions                    */
	int     numberPaths)    /* (I) Number of paths to extract               */
{
    static char routine[] = "CreateUniformSequence";
    int status = FAILURE;
    int i=0;
    int j=0;
    UniformRandomSequence *A = 0;
    char errorMsg[MAXBUFF];

	try
    {
        double *deviates = new double[numberDimensions];
        A = new UniformRandomSequence(aSeed,
			                numberDimensions);

        for(j=0;j<numberPaths;j++)
        {
            A->populateVector();
            deviates = A->getVector();
            for(i=0;i<numberDimensions;i++)
            {
                sequence[i+j*numberDimensions] = deviates[i];
            }
        }
  

	}catch(SCException &e){
        strcpy(errorMsg, e.what());
        goto RETURN;
    }

    delete A;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return (status);
}

// --------------------------------------------------------------------------
// CreateRandomGenerator
//
extern "C" void *CreateRandomGenerator(long aSeed)
{
    UniformRandomSequence *A = 0;
	try
    {
        A = new UniformRandomSequence(aSeed,
			                1);

	}catch(SCException &e){
        DR_Error(e.what());
        goto RETURN;
    }
RETURN:
    return reinterpret_cast<void *>(A);
}

// --------------------------------------------------------------------------
// RandomGeneratorFree
//
extern "C" void RandomGeneratorFree(void *A)
{
    delete reinterpret_cast<UniformRandomSequence *>(A);
}

// --------------------------------------------------------------------------
// RandomGeneratorGet
//
extern "C" double RandomGeneratorGet(void *A)
{
	try {
    UniformRandomSequence *u = reinterpret_cast<UniformRandomSequence *>(A);
    double deviate;
    u->populateVector();
    deviate = *u->getVector();
    return deviate;
    } catch (SCException &e){
        DR_Error(e.what());
        return -1.;;
    }
}
/*
extern "C" double Ran2(long aSeed)
{
    static char routine[] = "CreateUniformSequence";
    static UniformRandomSequence A(aSeed,1);
    int status = FAILURE;
    double deviate;

	try
    {  
        A.populateVector();
        deviate = *(A.getVector());
	}
    catch(SCException &e){
        DR_Error(e.what());
        goto RETURN;
    }

    status = SUCCESS;

RETURN:
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return deviate;
}
*/