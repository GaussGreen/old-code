/*
 * $Log: rand-gen.h,v $
 * Revision 1.11  2004/05/14 09:42:47  emezzine
 * Added  function signature.
 *
 * Revision 1.10  2004/01/23 15:42:13  ebenhamou
 * change name for generic pricer
 *
 * Revision 1.9  2003/09/23 09:38:25  emezzine
 * Added MRGK5Random()
 *
 * Revision 1.8  2003/07/01 20:01:23  jpriaudel
 * suppression of include armglob
 *
 * Revision 1.7  2003/05/19 10:04:50  emezzine
 * remettre une valeur  par defaut a jour
 *
 * Revision 1.6  2003/05/16 17:37:10  emezzine
 *  Modif d'un parametre de defaut
 *
 * Revision 1.5  2003/04/28 15:17:59  emezzine
 *  Ajout d'un autre generateur gaussien
 *
 * Revision 1.4  2000/11/13 14:23:57  sgasquet
 * Retrait fonction virtuelle pure
 *
 * Revision 1.3  2000/11/13 10:50:47  sgasquet
 * Ajout generateur MMT et Faure
 *
 * Revision 1.2  1998/11/23 18:33:05  nicolasm
 * Suppression des template dans ARM_OldRandomGenerator
 *
 */
           
#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include "expt.h"

double MRGK5Random(long FirstUse);

/**************************************************
    Generation de variables aleatoires
/**************************************************/

class ARM_OldRandomGenerator
{
public:
   ~ARM_OldRandomGenerator() {};

  virtual double Generate()
  {
            throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <Generate> method");
  }

};



/**************************************************************
// 1 VARIABLES UNIFORMES
/**************************************************/

class ARM_BaseGenerator: public ARM_OldRandomGenerator
{
public:
  ARM_BaseGenerator( long seed=0,long FirstUse=1 );
  virtual double Generate();
  long ReturnGeneratorParam(){return generator_param;}
  
  void Restart();
  void SetFirstUse(long FirstUse)
  {
      itsFirstUse = FirstUse;
  
  }

protected:
  long generator_param;
  long seed;
  long itsFirstUse;
};  


/********************************************************************************/
/*               Monstruous Mersenne Twister                                    */
/*                                                                              */
/*  This 1D-RNG is quick, has a period of 2^19937 - 1, and gives a sequence     */
/*  that is 623-dimensionally equidistributed                                   */
/*                                                                              */
/*	Generate returns a double, LGenerate a long                                 */
/*                                                                              */
/********************************************************************************/
#define uint32 unsigned long
#define K_N              (624)                 // length of state vector
#define K_M              (397)                 // a period parameter
#define K_K              (0x9908B0DFU)         // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v



// Variables uniformes 1D
class ARM_MMTGenerator: public ARM_BaseGenerator
{
private :

	uint32   state[K_N+1];     // state vector + 1 extra to not violate ANSI C
	uint32   *next;          // next random value is computed from here
	int      left;      // can *next++ this many times before reloading

	void seedMT(uint32 seed);
	uint32 reloadMT(void);
	


public:
	ARM_MMTGenerator( uint32 seed = 4357U );
	
	double	Generate();
	uint32	LGenerate();
};




/********************************************************************************/
/*               Gray Faure multidimentionnal sequence                          */
/*                                                                              */
/*  Faure multidimentional low-discrepancy sequence generator                   */
/*                                                                              */
/*  The dim variable of the constructor sets the dimension.                     */
/*  Do NOT override the default skip parameter unless you know what it is,      */
/*  and what you're doing.                                                      */
/*                                                                              */
/*	Generate returns a double, LGenerate a long                                 */
/*  VGenerate return a vector of double, VLGenerate a vector of long            */
/********************************************************************************/


class ARM_GFGenerator: public ARM_BaseGenerator
{
private :

	int GFdimension;
	int GFbase;
	int GFmbit;
	int **GFadd;
	int **GFsub;
	int *GFgray;
	int *GFb_ary;
	int **GFfaure;
	int *GFplusone;
	int ***GFpascal;
	long **GFpowbase;
	long *GFoutput;
	double *DGFoutput;

	void FirstPoint(int skip);
	void FreeMemory(void);
	void Pascal(void);
	void Tables(void);
	int Base(void);

	int GrayFaureInit(int dim,int skip);
	
public:
	
	ARM_GFGenerator( int dim = 1, int skip =4);
	~ARM_GFGenerator(void);
	
	double	Generate(void);
	long	LGenerate(void);

	long *VLGenerate(void);
	double *VGenerate(void);

};



/**************************************************************
// 2 VARIABLES NORMALES
/**************************************************/
class ARM_NormalGenerator: public ARM_OldRandomGenerator
{
public:
  ARM_NormalGenerator();
  ARM_NormalGenerator( ARM_BaseGenerator* base_generator_,int gaussianType = K_BOX_MULLER);
  ~ARM_NormalGenerator();

  virtual double Generate();

private:
  ARM_BaseGenerator * base_generator;
  short own_base_generator;
  int next_available;
  double next;

  int itsGaussianType;
};


class ARM_MMTNormalGenerator: public ARM_NormalGenerator
{
public:
  ARM_MMTNormalGenerator()
	  :ARM_NormalGenerator(){}

  
  ARM_MMTNormalGenerator( ARM_MMTGenerator* base_generator_,int gaussianType)  
	  :ARM_NormalGenerator(base_generator_){}

};


/**************************************************************
// 3 VARIABLES GAUSSIENNES
/**************************************************/
class ARM_GaussianGenerator: public ARM_OldRandomGenerator
{
public:
  ARM_GaussianGenerator(double mu_, double sigma2_);
  ARM_GaussianGenerator( ARM_BaseGenerator* base_generator_, 
                         double mu_, double sigma2_);
  ~ARM_GaussianGenerator();

  virtual double Generate();

private:
  double mu,sigma; // attention on stocke l'ecart-type
  ARM_NormalGenerator ng;
};







#endif
