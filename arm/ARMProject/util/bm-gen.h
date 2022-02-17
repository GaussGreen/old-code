           
#ifndef RANDOM_WALK_GENERATOR_H
#define RANDOM_WALK_GENERATOR_H

#include "rand-gen.h"
#include "linalg.h"

/************************************************************************
    Generation de marches aleatoires
/***********************************************************************/

class ARM_PathGenerator : public ARM_Object
{
public:
  ARM_PathGenerator(void) {};
  ARM_PathGenerator( long path_length_ );
  virtual ~ARM_PathGenerator();

  virtual double* Generate()
  {
            throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <Generate> method");
  }

  long Length() const;

protected:
  double* innov;
  long length;
};


class ARM_RandomWalkGenerator: public ARM_PathGenerator
{
public:
  ARM_RandomWalkGenerator(void) {};
  ARM_RandomWalkGenerator( long path_length_, ARM_BaseGenerator* base_gen_ );
  virtual ~ARM_RandomWalkGenerator(){};

  virtual double* Generate();

protected:
  ARM_NormalGenerator normal_gen;
};

class ARM_MirrorRWGenerator: public ARM_RandomWalkGenerator{
public:
  ARM_MirrorRWGenerator( long path_length_, ARM_BaseGenerator* base_gen_ );

  virtual double* Generate();

protected:
  int already_mirrored;
};


class ARM_BrownianPathGenerator: public ARM_PathGenerator
{
public:
   ARM_BrownianPathGenerator(void) {};
  ARM_BrownianPathGenerator( double final_time, long nber_of_steps, 
                             ARM_BaseGenerator* base_gen );
  virtual ~ARM_BrownianPathGenerator(){};

  virtual double* Generate();

protected:
  double sqrt_h;
  ARM_NormalGenerator normal_gen;
};


class ARM_MirrorBMPathGenerator: public ARM_BrownianPathGenerator
{
public:
  ARM_MirrorBMPathGenerator(void) {};
  ARM_MirrorBMPathGenerator( double final_time, long nber_of_steps,
                             ARM_BaseGenerator* base_gen );
  virtual ~ARM_MirrorBMPathGenerator(){};

  virtual double* Generate();

protected:
  short already_mirrored;
};



#endif
