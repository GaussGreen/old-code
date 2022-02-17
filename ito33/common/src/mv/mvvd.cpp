// ------------------------*-*-C++-*-*-----------------------------------------
// file..........: mvvd.cpp
// purpose.......: basic vector class (double precision)
// author........: ZHANG Yunzhi
// where.........: ITO 33
// credit........: MV++ Version 1.5
// created.......: 27 july 2000
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>              
#include <memory>
#include <cmath>
#include "ito33/afterstd.h"
                
#include "ito33/mv/mvvd.h"

using namespace std;

//#define MY_TEST_DEBUG 1

//---------------------------//
// constructeurs/destructeur //
//---------------------------//
CVecteurDouble::CVecteurDouble() 
  : m_pd(0), m_iDim(0) , m_iRef(0)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble() " << endl;
# endif
}


CVecteurDouble::CVecteurDouble(unsigned int _iDim)
  : m_pd(new double[_iDim]), m_iDim(_iDim),  m_iRef(0)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(unsigned int _iDim) " << endl;
# endif
#  if MV_CAUTION
    if (m_pd == NULL)
    {
        cerr << "Error: NULL pointer in CVecteurDouble(int) constructor " << endl;
        cerr << "       Most likely out of memory... " << endl;
        exit(1);
    }
#  endif
}


CVecteurDouble::CVecteurDouble(unsigned int _iDim, const double& _d)
  : m_pd(new double[_iDim]), m_iDim(_iDim), m_iRef(0)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(unsigned int _iDim, const double& _d) " << endl;
# endif
#  if MV_CAUTION
    if (m_pd == NULL)
    {
        cerr << "Error: NULL pointer in CVecteurDouble(int) constructor " << endl;
        cerr << "       Most likely out of memory... " << endl;
        exit(1);
    }
#  endif
  for(unsigned int i = 0; i < m_iDim; m_pd[i++] = _d);
    /*
    int
      i,
      N = m_iDim,
      Nminus4 = m_iDim - 4;
    for (i=0; i < Nminus4; )
    {
        m_pd[i++] = _d;
        m_pd[i++] = _d;
        m_pd[i++] = _d;
        m_pd[i++] = _d;
    }
    for (; i < N; m_pd[i++] = _d);
    */
}


CVecteurDouble::CVecteurDouble(double* _pd, unsigned int _iDim,
  MV_Vector_::ref_type _iRef)
  : m_pd(_pd), m_iDim(_iDim), m_iRef(_iRef) 
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(double* _pd, unsigned int _iDim, MV_Vector_::ref_type _iRef) " << endl;
# endif
}


CVecteurDouble::CVecteurDouble(const CVecteurDouble &_V, MV_Vector_::ref_type _iRef)
  : m_pd(_V.m_pd), m_iDim(_V.m_iDim), m_iRef(_iRef) 
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(const CVecteurDouble &_V, MV_Vector_::ref_type _iRef) " << endl;
# endif
}

/*
CVecteurDouble::CVecteurDouble(const CVecteurDouble & _V)
  : m_pd(new double[_V.m_iDim]), m_iDim(_V.m_iDim) , m_iRef(0)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(const CVecteurDouble & _V) " << endl;
# endif
#  if MV_CAUTION
    if (m_pd == NULL)
    {
        cerr << "Error:  Null pointer in CVecteurDouble(const MV_Vector&); " << endl;
        exit(1);
    }
# endif
  for (unsigned int i=0; i<m_iDim; i++)
      m_pd[i] = _V.m_pd[i];
}
*/

CVecteurDouble::CVecteurDouble(double* _pd, unsigned int _iDim)
  : m_pd(new double[_iDim]), m_iDim(_iDim) , m_iRef(0)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(double* _pd, unsigned int _iDim) " << endl;
# endif
#  if MV_CAUTION
    if (m_pd == NULL)
    {
        cerr << "Error: Null pointer in CVecteurDouble(double*, int) " << endl;
        exit(1);
    }
#  endif
  for (unsigned int i=0; i<_iDim; i++)
    m_pd[i] = _pd[i];
}


CVecteurDouble::CVecteurDouble(const double* _pd, unsigned int _iDim)
  : m_pd(new double[_iDim]), m_iDim(_iDim) , m_iRef(0)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::CVecteurDouble(const double* _pd, unsigned int _iDim) " << endl;
# endif
#  if MV_CAUTION
    if (m_pd == NULL)
    {
        cerr << "Error: Null pointer in CVecteurDouble(double*, int) " << endl;
        exit(1);
    }
#  endif
  for (unsigned int i = 0; i < m_iDim; i++)
    m_pd[i] = _pd[i];
}


CVecteurDouble::~CVecteurDouble()
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::~CVecteurDouble() " << endl;
# endif
  if (m_pd && !m_iRef )
    delete [] m_pd;
}




/***************************************************************************
  affectation
 ***************************************************************************/
CVecteurDouble& CVecteurDouble::operator=(const double& _d) 
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator=(const double& _d) " << endl;
# endif
  for (unsigned int i = 0; i < m_iDim; m_pd[i++] = _d);
  return *this;
}

/*
 * .si Vecteur present est une reference, il faut que les deux vecteurs ayont
 *  le meme taille
 * .sinon, on cree un nouveau tableau pour le vecteur present.
 */
CVecteurDouble& CVecteurDouble::operator=(const CVecteurDouble & _V) 
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator=(const CVecteurDouble & _V) " << endl;
# endif
   unsigned int i;

  if (m_iRef )                  // is this structure just a pointer?
  {
#    if MV_CAUTION
      if (m_iDim != _V.m_iDim)     // check conformance,
      {
          cerr << "MV_VectorRef::operator=  non-conformant assignment.\n";
          exit(1);
      }
#    endif

    // handle overlapping matrix references
    if ((_V.m_pd + _V.m_iDim) >= m_pd)
    {
      // overlap case, copy backwards to avoid overwriting results
      for (i = m_iDim - 1; i < (unsigned int)-1; i--)
          m_pd[i] = _V.m_pd[i];
    }
    else
    {
      for (i = 0; i < m_iDim; i++)
          m_pd[i] = _V.m_pd[i];
    }     
  }
  else
  {
    if(m_iDim != _V.m_iDim)
    {
      if(m_pd)
        delete [] m_pd;
      m_iDim = _V.m_iDim;
      m_pd = new double [m_iDim];
    }
    for (i = 0; i< m_iDim; i++)
      m_pd[i] = _V.m_pd[i];
  }
  return *this;   
}


//--------------------------------------------------------------------------------//
// newsize(_iDim) fait auto-destruction et re-initie un vecteur de taille _iDim   //
// Attention, en appelant cette fonction, vous aurrez des ennuies avec des vecteurs/
// references a le vecteur du depart                                              //
//--------------------------------------------------------------------------------//
CVecteurDouble& CVecteurDouble::newsize(unsigned int _iDim)
{
  if(m_iRef == 0 && m_iDim != _iDim)
  {
    if(m_pd)
      delete [] m_pd;
    m_iDim = _iDim;
    m_pd = new double [m_iDim];

    memset(m_pd, 0, sizeof(double) * m_iDim);
  }
  return *this;
}


//---------------//
// sous vecteurs //
//---------------//

// remarque: sous vecteurs sont toujours donnes par reference


//------------------------------------------------------------------------------//
// operator()(void) s'agit de rentrer le vecteur entier
//------------------------------------------------------------------------------//
CVecteurDouble CVecteurDouble::operator()(void)
{
    return CVecteurDouble(m_pd, m_iDim, MV_Vector_::ref);
}


const CVecteurDouble CVecteurDouble::operator()(void) const
{
    return CVecteurDouble(m_pd, m_iDim, MV_Vector_::ref);
}


//-------------------------------------------------------------------------------//
// operator()(MV_VecIndex &I) rentre le sous vecteur defini par I
//-------------------------------------------------------------------------------//
CVecteurDouble CVecteurDouble::operator()(const MV_VecIndex &_I) 
{
  if (_I.all())
      return CVecteurDouble(m_pd, m_iDim, MV_Vector_::ref);
  else
  {
#    if MV_CAUTION
      // check that index is not out of bounds
      //
      if ( (unsigned int)_I.end() >= m_iDim)
      {
          cerr << "MV_VecIndex: (" << _I.start() << ":" << _I.end() << 
              ") too big for matrix (0:" << m_iDim - 1 << ") " << endl;
          exit(1);
      }
#    endif
    return CVecteurDouble(m_pd + _I.start(), _I.end() - _I.start() + 1,
          MV_Vector_::ref);
  }
}


const CVecteurDouble CVecteurDouble::operator()(const MV_VecIndex &_I) const
{
  if (_I.all())
      return CVecteurDouble(m_pd, m_iDim, MV_Vector_::ref);
  else
  {
#    if MV_CAUTION
      if ( (unsigned int)_I.end() >= m_iDim)
      {
        cerr << "MV_VecIndex: (" << _I.start() << ":" << _I.end() << 
                ") too big for matrix (0:" << m_iDim - 1 << ") " << endl;
        exit(1);
      }
#    endif
    return CVecteurDouble(m_pd+ _I.start(), _I.end() - _I.start() + 1,
          MV_Vector_::ref);
  }
}


//--------------------------------------------------------------------------------//
// RefSubVector() donne reference a un sous vecteur                               //
// ATTENTION A UTILISATION DE CETTE FONCTION:                                     //
//  Cette fonction ne verifie pas si m_iRef etait 0 ou pas, et l'affecte par 1    //
//  Donc si avant l'appel de cette fonction, ce vecteur a son propre memoire      //
//  (m_iRef=0), cette memoir va etre perdu.                                       //
// A utilisateur de verifier que m_iRef != 0                                      //
//--------------------------------------------------------------------------------//
CVecteurDouble& CVecteurDouble::RefSubVector(const CVecteurDouble &_V,
                                             const MV_VecIndex &_I)
{
#  if MV_CAUTION
    if(m_iRef == 0 && m_pd == NULL)
      {
      cerr << "CVecteurDouble::RefSubVector() : vecteur est occupe" << endl;
      exit(1);
      }
#  endif
  if (_I.all())
    {
    m_pd = _V.m_pd;
    m_iDim = _V.m_iDim;
    m_iRef = MV_Vector_::ref;
    }
  else
  {
#    if MV_CAUTION
      if ( (unsigned int)_I.end() >= _V.m_iDim)
      {
        cerr << "CVecteurDouble::RefSubVector() :\n"
          << "MV_VecIndex: (" << _I.start() << ":" << _I.end() << 
                ") too big for matrix (0:" << _V.m_iDim - 1 << ") " << endl;
        exit(1);
      }
#    endif
    m_pd = _V.m_pd + _I.start();
    m_iDim = _I.end() - _I.start() + 1;
    m_iRef = MV_Vector_::ref;
  }
  return *this;
}


//------------------------------------------------------------------------------
// affichages
//------------------------------------------------------------------------------
ostream& operator<<(ostream& _s, const CVecteurDouble& _V)
{
  int N = _V.size();

  for (int i=0; i< N; i++)
    _s << _V[i] << "\n";
  
  return _s;
}


//*******************************************************************************
// operateurs arithmetiques
//*******************************************************************************
CVecteurDouble& CVecteurDouble::operator*=(const double &_d)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator*=(const double &_d) " << endl;
#  endif
  for( unsigned int i = 0; i < m_iDim; m_pd[i++] *= _d);
  return *this;
}


CVecteurDouble CVecteurDouble::operator*(const double &_d) const
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator*(const double &_d) " << endl;
#  endif
  CVecteurDouble Vresult(m_iDim);
  for (unsigned int i = 0; i < m_iDim; i++)
    Vresult.m_pd[i] = m_pd[i] * _d;
  return Vresult;
}


CVecteurDouble operator*(const double &_d, const CVecteurDouble &_V)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble operator*(const double &a, const CVecteurDouble &x)" << endl;
#  endif
  CVecteurDouble Vresult(_V.m_iDim);
  for (unsigned int i=0; i < _V.m_iDim;i++)
    Vresult.m_pd[i] = _V.m_pd[i] * _d;
  return Vresult;
}


CVecteurDouble CVecteurDouble::operator+(const CVecteurDouble &_V) const
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator+(const CVecteurDouble &_V) " << endl;
#  endif
#  if MV_CAUTION
    if(m_iDim != _V.iDim)
      {
         cout << "CVecteurDouble::operator+: Incompatible vector lengths in +." << endl;
         exit(1);
      }
#  endif
  CVecteurDouble VR(m_iDim);// = new CVecteurDouble(m_iDim);
  for (unsigned int i = 0; i < m_iDim; i++)
    VR.m_pd[i] = m_pd[i] + _V.m_pd[i];
  return VR;
}
     
     
CVecteurDouble CVecteurDouble::operator-(const CVecteurDouble &_V) const
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator-() " << endl;
#  endif
#  if MV_CAUTION
    if(m_iDim != _V.iDim)
      {
         cout << "CVecteurDouble::operator-: Incompatible vector lengths in +." << endl;
         exit(1);
      }
#  endif
  CVecteurDouble VR(m_iDim);
  for (unsigned int i = 0; i < m_iDim; i++)
    VR.m_pd[i] = m_pd[i] - _V.m_pd[i];
  return VR;
}
          

CVecteurDouble& CVecteurDouble::operator+=(const CVecteurDouble &_V)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator+=() " << endl;
#  endif
#  if MV_CAUTION
    if(m_iDim != _V.iDim)
      {
         cout << "CVecteurDouble::operator+=: Incompatible vector lengths in +." << endl;
         exit(1);
      }
#  endif
  for (unsigned int i = 0; i < m_iDim; i++)
    m_pd[i] += _V.m_pd[i];
  return *this;
}
          

CVecteurDouble& CVecteurDouble::operator-=(const CVecteurDouble &_V)
{
#  if MY_TEST_DEBUG
    cout << "CVecteurDouble::operator-=() " << endl;
#  endif
#  if MV_CAUTION
    if(m_iDim != _V.iDim)
      {
         cout << "CVecteurDouble::operator-=: Incompatible vector lengths in +." << endl;
         exit(1);
      }
#  endif
  for (unsigned int i = 0; i < m_iDim; i++)
    m_pd[i] -= _V.m_pd[i];
  return *this;
}
  

double dot(const CVecteurDouble &_x, const CVecteurDouble &_y)
{
#  if MY_TEST_DEBUG
    cout << "dot(const CVecteurDouble &x, const CVecteurDouble &y)" << endl;
#  endif  
#  if MV_CAUTION
    //  Check for compatible dimensions:
    if (_x.m_iDim != _y.m_iDim)
    {
       cout << "Incompatible dimensions in dot(). " << endl;
       exit(1);
     }
#  endif

  double temp =  0;
  for (unsigned int i = 0; i < _x.m_iDim; i++)
       temp += _x.m_pd[i] * _y.m_pd[i];
  return temp;
}


double norm(const CVecteurDouble &_x)
{
  double temp =  0;
  for (unsigned int i = 0; i < _x.m_iDim; i++)
       temp += _x.m_pd[i] * _x.m_pd[i];
  return sqrt(temp);
}

