/*
 * $Log: bm-gen.cpp,v $
 * Revision 1.3  2001/01/30 09:56:41  smysona
 * Supression de la macro cleanit et insertion de frmutils.h
 *
 * Revision 1.2  2000/12/08 19:39:33  smysona
 * innov n'est pas detruite si l'on est un FAURE_GENERATOR
 *
 * Revision 1.1  1998/11/19 11:09:22  nicolasm
 * Initial revision
 *
 */

#include <math.h>           
#include "bm-gen.h"
//#include "frmutils.h"
#define cleanit(__ptr) \
    { if(__ptr)        \
     { delete __ptr;   \
       __ptr = NULL;   \
     }                 \
    }


ARM_PathGenerator::ARM_PathGenerator( long path_length_ )
                     :innov(new double[path_length_]),length(path_length_)
{
}

ARM_PathGenerator::~ARM_PathGenerator()
{
   //delete [] innov;
    if (GetName() != ARM_FAUREGENERATOR)
        cleanit(innov);
}

long ARM_PathGenerator::Length() const
{
   return length;
}


ARM_RandomWalkGenerator::ARM_RandomWalkGenerator( long path_length_,
                                                 ARM_BaseGenerator*  base_gen_ )
                          :ARM_PathGenerator(path_length_),normal_gen(base_gen_)
{
}


double* ARM_RandomWalkGenerator::Generate()
{
  for (long i=0;i<length;i++)
  {
      innov[i]=normal_gen.Generate();
  }
  return innov;
}



ARM_MirrorRWGenerator::ARM_MirrorRWGenerator( long path_length_, ARM_BaseGenerator* base_gen_ )
  :ARM_RandomWalkGenerator(path_length_,base_gen_),already_mirrored(1)
{
}

double *ARM_MirrorRWGenerator::Generate()
{
  if (already_mirrored)
  {
     already_mirrored=0;
     return ARM_RandomWalkGenerator::Generate();
  }
  else
  {
     for (long i=0;i<length;i++)
        innov[i]*=-1;
     already_mirrored=1;
     return innov;
  }
}
  




ARM_BrownianPathGenerator::ARM_BrownianPathGenerator( double final_time, 
                                                     long nber_of_steps,
                                                     ARM_BaseGenerator* base_gen )
                          :sqrt_h(sqrt(final_time/nber_of_steps)),
                           ARM_PathGenerator(nber_of_steps),normal_gen(base_gen)
{
}


double* ARM_BrownianPathGenerator::Generate()
{
   for (long i=0;i<length;i++)
   {
     innov[i]=sqrt_h*normal_gen.Generate();
   }
   return innov;
}


ARM_MirrorBMPathGenerator::ARM_MirrorBMPathGenerator( double final_time, 
                                                     long nber_of_steps, 
                                                     ARM_BaseGenerator* base_gen )
  :ARM_BrownianPathGenerator(final_time,nber_of_steps,base_gen),already_mirrored(1)
{
}

double* ARM_MirrorBMPathGenerator::Generate()
{
  if (already_mirrored)
  {
     already_mirrored=0;
     return ARM_BrownianPathGenerator::Generate();
  }
  else
  {
     for (long i=0;i<length;i++)
       innov[i]*=-1;
     already_mirrored=1;
     return innov;
  }
}


