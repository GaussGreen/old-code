// preconditioner which is diagonal


#ifndef DIAGPRE_H
#define DIAGPRE_H

#include "mvvd.h"
#define VECTOR_double CVecteurDouble

#include "cmorse.hpp"

class DiagPreconditioner_double {
private:
  VECTOR_double m_Diag;

public:
  // Je ne veux pas que constructeur echoue,
  // les autres sont faits dans la fonction surchargee Init()
  DiagPreconditioner_double() {}
  DiagPreconditioner_double(unsigned int _iDim, double _d) : m_Diag(_iDim, _d) {}
  ~DiagPreconditioner_double (void) {}
  
  //============================================================================//
  // constructeur sous format init(), il retourne 0 si construction de precond. //
  // reussit, 1 si cela echoue                                                  //
  //============================================================================//

  //-----------------------------------------------------------------------------------
  //. Init(CMorseMatrix) construit la preconditionneur pour une matrice morse
  //. Demande:
  //  la matrice morse (carree bien sur) doit avoir une structure dont les indices
  //  commencent a 0 (iLineMin, iColMin = 0)
  //. Valeur retournee:
  //  =0 : construction reussite
  //  >0 : indice dont l'element diagonal est nul
  //  <0 : l'inverse est l'indice dont l'element est trop grand
  //-----------------------------------------------------------------------------------
  int Init(CMorseMatrix &_M);


  // fonctions demandees par iml++
  VECTOR_double solve (const VECTOR_double &_V) const;
  VECTOR_double trans_solve (const VECTOR_double &_V) const;
  

  // fonctions dont je ne vois pas interet
  const double& diag(int _i) const     { return m_Diag(_i); }
  double&       diag(int _i)           { return m_Diag(_i); }
};

#endif
