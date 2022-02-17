/*
 * $Log: multidim-gen.h,v $
 * Revision 1.6  2003/04/28 15:16:39  emezzine
 * Ajout d'un autre generateur gaussien
 *
 * Revision 1.5  2001/06/22 16:21:04  sgasquet
 * Correction appel Gaussian variable dans Faure
 *
 * Revision 1.4  2001/03/06 15:51:09  abizid
 * Corrections des conneries sur fonctions Generate()
 *
 * Revision 1.3  2000/12/08 19:40:51  smysona
 * Rajout d'un setname FAUREGENERATOR
 *
 * Revision 1.2  2000/10/25 14:27:54  sgasquet
 * Ajout constructeur par defaut
 *
 * Revision 1.1  2000/10/24 10:30:19  sgasquet
 * Initial revision
 *
 *
 */
           
#ifndef MULTIDIM_GENERATOR_H
#define MULTIDIM_GENERATOR_H


#include "bm-gen.h"
#include "linalg.h"
#include "gaussian.h"



class ARM_MultiDimGaussianGenerator: public ARM_PathGenerator
{
public:
    ARM_MultiDimGaussianGenerator(long nb_path_, long time_step_,
                                  long nb_factors_)
             : ARM_PathGenerator(time_step_*nb_factors_), time_step(time_step_),
               nb_factors(nb_factors_), nb_path(nb_path_)
    {
    }

    ARM_MultiDimGaussianGenerator(){}

    virtual double* Generate()
    {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <Generate> method");
    }


    virtual double GaussianVariable(int p, int t, int factor)
    {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <GaussianVariable> method");
    }

protected:
    long time_step;
    long nb_factors;
    long nb_path;
};




/************************************************************************
    Generation de marches aleatoires
    Les va sont appelees soit a l'aide de la fonction GaussianVariable
    soit avec la fonction Generate qui renvoie le vecteur de va (de longueur
    nb_timeStep*nb_factor  ) pour
    le chemin courant; on suppose que la fonction Generate est appelee
    pour des valeurs croissantes de p avec une remise a 0 des que le nombre 
    de chemins est atteint
/***********************************************************************/

class ARM_ScramblingGenerator: public ARM_MultiDimGaussianGenerator
{
public:
    ARM_ScramblingGenerator(long nb_path_, long time_step_,
                            long nb_factors_, long seed_=180);

    ARM_ScramblingGenerator(){}

    virtual ~ARM_ScramblingGenerator();

    virtual double* Generate();
  
    void Stratify();
    void Scramble(int factor);


    virtual double GaussianVariable(int p, int t, int factor)
    {
        int rank = scrambledIndex[factor][p][t];
        return (stratified_values[rank]);
    }

protected:

    ARM_MMTGenerator *itsBaseGen;
    int*** scrambledIndex;
    long current_path;
    double* stratified_values;
};




class ARM_SimpleGenerator: public ARM_MultiDimGaussianGenerator
{
public:
    ARM_SimpleGenerator(long nb_path_, long time_step_,
                        long nb_factors_, long seed_=10000, int gaussianType = K_BOX_MULLER);

    ARM_SimpleGenerator(){}

    virtual ~ARM_SimpleGenerator();

    virtual double* Generate();
  
    virtual double GaussianVariable(int p, int t, int factor)
    {
        if (current_path==p)
            return (innov[factor * time_step + t]);
        else
        {
            double* tmp = Generate();
            return (innov[ factor * time_step + t]);
        }
    }

protected:
    ARM_MMTGenerator *itsBaseGen;
    ARM_MMTNormalGenerator *itsNormalGen;
    long current_path;
    int already_mirrored; // ajout du symetrise
};





class ARM_FaureGenerator: public ARM_MultiDimGaussianGenerator
{
public:
    ARM_FaureGenerator( long nb_path_, long time_step_, long nb_factors_);

    ARM_FaureGenerator()
    {
        SetName(ARM_FAUREGENERATOR);
    }

    virtual ~ARM_FaureGenerator();

    virtual double* Generate();

    double GaussianVariable(int p, int t, int factor)
    {
        if (p == current_path)
        {
            double CurrGauss = INV_PART_FUNC_NOR(innov[(time_step)*factor+t]);	
            return (CurrGauss);
        }
        else
        {
            current_path = p;
            innov = itsBaseGen->VGenerate();
            double CurrGauss = INV_PART_FUNC_NOR(innov[(time_step)*factor+t]);

            return (CurrGauss);
        }
    }


protected:
    ARM_GFGenerator*itsBaseGen;
    long current_path;
};

#endif
