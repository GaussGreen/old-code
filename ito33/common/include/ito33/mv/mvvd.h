// ------------------------*-*-C++-*-*-----------------------------------------
// file..........: mvvd.hpp
// purpose.......: basic vector class (double precision)
// author........: ZHANG Yunzhi
// where.........: ITO 33
// credit........: MV++ Version 1.5
// created.......: 27 july 2000
// ----------------------------------------------------------------------------

#ifndef _MV_VD_H
#define _MV_VD_H

#include "mvgeneral.h"

#include "ito33/beforestd.h"

#include <iostream>
#ifdef MV_CAUTION
#  include <cstdlib>
#endif

#ifdef MV_VECTOR_BOUNDS_CHECK
#  include <cassert>
#endif

#include "ito33/afterstd.h"

#include "mvvind.h"

class CVecteurDouble
{ 
protected:
  double *m_pd;
  unsigned int m_iDim;
  int m_iRef;  // 0 or 1; does this own its own memory space?
public:                                                            
    CVecteurDouble();  
    CVecteurDouble(unsigned int _iDim);
    CVecteurDouble(unsigned int _iDim, const double& _d);
    CVecteurDouble(double* _pd, unsigned int _iDim);
    CVecteurDouble(const double* _pd, unsigned int _iDim);      

    CVecteurDouble(const CVecteurDouble & _V)
        : m_pd(new double[_V.m_iDim]), m_iDim(_V.m_iDim) , m_iRef(0)
    {
    #  if MY_TEST_DEBUG
        std::cout << "CVecteurDouble::CVecteurDouble(const CVecteurDouble & _V) " << std::endl;
    # endif
    #  if MV_CAUTION
        if (m_pd == NULL)
        {
            std::cerr << "Error:  Null pointer in CVecteurDouble(const MV_Vector&); " << std::endl;
            exit(1);
        }
    # endif
      for (unsigned int i=0; i<m_iDim; i++)
          m_pd[i] = _V.m_pd[i];
    };

    // reference d'une structure deja existe
    // ATTENTION: _iRef NE peut avoir QU'une valeur: 1
    CVecteurDouble(double* _pd, unsigned int _iDim, MV_Vector_::ref_type _iRef);
    CVecteurDouble(const CVecteurDouble &_V, MV_Vector_::ref_type _iRef);
    
    ~CVecteurDouble();        

    
//--------------------------------------------------------------------------------//
// operator()(_i) et [](_i) donne acces a l'element m_pd[_i]                      //
//--------------------------------------------------------------------------------//
//-----------//
// operators //
//-----------//

double& operator()(unsigned int _i)
{
# ifdef MV_VECTOR_BOUNDS_CHECK
    assert(_i < m_iDim);
# endif
  return m_pd[_i];
}

const double& operator()(unsigned int _i) const
{
# ifdef MV_VECTOR_BOUNDS_CHECK
    assert(_i < m_iDim);
# endif
  return m_pd[_i];
}

double& operator[](unsigned int _i)
{
# ifdef MV_VECTOR_BOUNDS_CHECK
    assert(_i < m_iDim);
# endif
  return m_pd[_i];
}


const double& operator[](unsigned int _i) const
{
# ifdef MV_VECTOR_BOUNDS_CHECK
    assert(_i < m_iDim);
# endif
  return m_pd[_i];
}


//--------------------------------------------------------------------------------//
// operator()(MV_VecIndex) retourne une sous-vecteur                              //
//--------------------------------------------------------------------------------//
    CVecteurDouble       operator()(const MV_VecIndex &_I) ;
    const CVecteurDouble operator()(const MV_VecIndex &_I) const;
    CVecteurDouble       operator()(void);
    const CVecteurDouble operator()(void) const;

//--------------------------------------------------------------------------------//
// RefSubVector() donne reference a un sous vecteur                               //
// ATTENTION A UTILISATION DE CETTE FONCTION:                                     //
//  Cette fonction ne verifie pas si m_iRef etait 0 ou pas, et l'affecte par 1    //
//  Donc si avant l'appel de cette fonction, ce vecteur a son propre memoire      //
//  (m_iRef=0), cette memoir va etre perdu.                                       //
// A utilisateur de verifier que m_iRef != 0                                      //
//--------------------------------------------------------------------------------//
    CVecteurDouble& RefSubVector(const CVecteurDouble &_V, const MV_VecIndex &_I);

//--------------------------------------------------------------------------------//
// des fonctions qui sortent les variables membres                                //
//--------------------------------------------------------------------------------//
    inline const double* get() const { return m_pd; }
    inline double* get() { return m_pd; }
    inline unsigned int  size() const { return m_iDim;}
    inline unsigned int  dim() const { return m_iDim;}
    inline int           ref() const { return  m_iRef;}
    inline int           null() const {return m_iDim == 0;}


//--------------------------------------------------------------------------------//
// newsize(_iDim) fait auto-destruction et re-initie un vecteur de taille _iDim   //
// Attention, en appelant cette fonction, vous aurrez des ennuies avec des vecteurs/
// references a le vecteur du depart                                              //
//--------------------------------------------------------------------------------//
   CVecteurDouble& newsize(unsigned int _iDim);


//--------------------------------------------------------------------------------//
// affectations                                                                    //
//--------------------------------------------------------------------------------//
    CVecteurDouble& operator=(const double& _d);   // vecteur constant
    /*
     * .si Vecteur present est une reference, il faut que m_iDim = _V.m_iDim
     *  ceci est verifie par le code qui si MV_CAUTION est active
     * .sinon, on cree un nouveau tableau pour le vecteur present. Dans ce cas
     *  vous avez le meme probleme que ce que j'ai mentionne pour newsize()
     */
    CVecteurDouble& operator=(const CVecteurDouble& _V);


//--------------------------------------------------------------------------------//
//  affectations arithmetiques
//--------------------------------------------------------------------------------//
    CVecteurDouble& operator*=(const double &_d);
    CVecteurDouble& operator+=(const CVecteurDouble &_V);
    CVecteurDouble& operator-=(const CVecteurDouble &_V);


//--------------------------------------------------------------------------------//
//  operations arithmetiques
//--------------------------------------------------------------------------------//
    CVecteurDouble operator*(const double &_d) const;
    CVecteurDouble operator+(const CVecteurDouble &_V) const;
    CVecteurDouble operator-(const CVecteurDouble &_V) const;

    friend CVecteurDouble operator*(const double &_d, const CVecteurDouble &_V);
    friend double dot(const CVecteurDouble &_x, const CVecteurDouble &_y);
    friend double norm(const CVecteurDouble &_x);


//--------------------------------------------------------------------------------//
//  affichages
//--------------------------------------------------------------------------------//
    friend std::ostream& operator<<(std::ostream &_s, const CVecteurDouble &_V);

};                                                                     
#endif  
