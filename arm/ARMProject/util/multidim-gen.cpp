#include <math.h>   
#include <stdlib.h>       
#include "multidim-gen.h"

#define cleanit(__ptr) \
    { if(__ptr)        \
     { delete __ptr;   \
       __ptr = NULL;   \
     }                 \
    }                

////////////////////////////////////////////////////
///	Class  : ARM_ScramblingGenerator
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_ScramblingGenerator::ARM_ScramblingGenerator(long nb_path_, long time_step_,
                                                 long nb_factors_, long seed_)
                            :ARM_MultiDimGaussianGenerator(nb_path_, time_step_,
                                                           nb_factors_)
{
    current_path = 0;

    // allocation scrambledIndex
    scrambledIndex = new int** [nb_factors_];
    stratified_values = new double[nb_path_];

    for (long j=0; j< nb_factors_; j++)
        scrambledIndex[j] = InitializeIntMatrix(nb_path_, time_step_);

    itsBaseGen = new ARM_MMTGenerator(seed_);

    Stratify();
    for (long i = 0; i < nb_factors_; i++)
        Scramble(i);
}

////////////////////////////////////////////////////
///	Class  : ARM_ScramblingGenerator
///	Routine: Constructor
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////

ARM_ScramblingGenerator::~ARM_ScramblingGenerator()
{
    for (long j = 0; j< nb_factors; j++)
        FreeIntMatrix(scrambledIndex[j], nb_path);

    delete scrambledIndex;
    delete itsBaseGen;
    delete stratified_values;

}

////////////////////////////////////////////////////
///	Class  : ARM_ScramblingGenerator
///	Routine: generate
///	Returns: 
///	Action : function
////////////////////////////////////////////////////


double* ARM_ScramblingGenerator::Generate()
{
    int i, j, k = 0;

    for ( j = 0; j < nb_factors; j++)
        for (  i = 0; i < time_step; i++)
        {		
            innov[k] = GaussianVariable(current_path, i, j);
            k++;
        }

    current_path++;
    if ( current_path >= nb_path)
        current_path = 0;

    return innov;
}

////////////////////////////////////////////////////
///	Class  : ARM_ScramblingGenerator
///	Routine: Stratify
///	Returns: 
///	Action : variance reduction
////////////////////////////////////////////////////

void ARM_ScramblingGenerator::Stratify()
{
    int i;

    double x, dx;

    dx = 1.0/nb_path;
    x = 0.5*dx;
    stratified_values[0] = INV_PART_FUNC_NOR(x);

    for (i=1; i<nb_path; i++)
    {
        x+=dx;
        stratified_values[i] = INV_PART_FUNC_NOR(x);
    }
}



void ARM_ScramblingGenerator::Scramble(int nb_factor)
{
    int k, p, n0;
    int* n = new int[nb_path];

    for (k = 0; k < time_step; k++)
    {
        for (p = 0; p < nb_path; p++)
        {
            n[p] = p;
        }

        for (p = 0; p < nb_path; p++)
        {
            n0 = (int) floor( (nb_path-p)*(double)itsBaseGen->Generate());
            scrambledIndex[nb_factor][p][k] = n[n0];
            n[n0] = n[nb_path-p-1];
        }
    }

    delete n;
}



////////////////////////////////////////////////////
///	Class  : ARM_SimpleGenerator
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////



ARM_SimpleGenerator::ARM_SimpleGenerator(long nbpath, long timestep,
                                         long nbfactor,long seed,int gaussianType)
                                         :ARM_MultiDimGaussianGenerator(nbpath,timestep,nbfactor)
{
    current_path = -1;
    already_mirrored = 1;

    itsBaseGen = new ARM_MMTGenerator(seed);
    itsNormalGen = new ARM_MMTNormalGenerator(itsBaseGen,gaussianType);
    double* tmp = Generate();
}

////////////////////////////////////////////////////
///	Class  : ARM_SimpleGenerator
///	Routine: destructor
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_SimpleGenerator::~ARM_SimpleGenerator()
{
    delete itsBaseGen;
    delete itsNormalGen;
}


////////////////////////////////////////////////////
///	Class  : ARM_SimpleGenerator
///	Routine: destructor
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////


double* ARM_SimpleGenerator::Generate()
{
    current_path++;
 
    if (already_mirrored)
    {
        for (long i = 0; i < length ; i++)          
            innov[i] = itsNormalGen->Generate();         
        already_mirrored = 0;
    }
    else
    {
        for (long i = 0; i < length ; i++)          
            innov[i] *= (-1.);        
        already_mirrored = 1;
    }

    if (current_path >= nb_path)
        current_path = 0;

    return innov;
}





ARM_FaureGenerator::ARM_FaureGenerator(long nb_path_, long time_step_, 
                                       long nb_factors_)
                          :ARM_MultiDimGaussianGenerator(nb_path_, time_step_,
                                                         nb_factors_)
{

    SetName(ARM_FAUREGENERATOR);
    current_path = 0;

    // skip = 4: dim max 70 approx
    // skip = 3: dim max 1200 approx

    int dim = nb_factors_*time_step_;

    if (dim<=60)
	itsBaseGen = new ARM_GFGenerator(nb_factors_*time_step_, 4);
    else 
        itsBaseGen = new ARM_GFGenerator(nb_factors_*time_step_, 3);

    cleanit (innov);
    innov = itsBaseGen->VGenerate();
}



ARM_FaureGenerator::~ARM_FaureGenerator()
{
    delete itsBaseGen;
}



double* ARM_FaureGenerator::Generate()
{
    long k = 0;
    long current = current_path;
    current_path = 0;
  
    for (long i = 0; i < nb_factors ; i++)  
        for (long j = 0; j < time_step; j++)			  
        {      
            innov[k] = GaussianVariable(current, j, i);  
            k++;
        }

    current_path++;

    if (current_path >= nb_path)
        current_path = 0;

    return innov;
}

