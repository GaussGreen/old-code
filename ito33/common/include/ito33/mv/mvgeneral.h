/*******************************-*-C++-*-**********************************
 * file...........: mvgeneral.h
 * purpose........: definition of some general data structure and paramter
 *                  for matrix/vector library
 * author.........: ZHANG Yunzhi
 * where..........: ITO 33
 * created........: 27 july 2000
 **************************************************************************/

#ifndef _MV_GENERAL_H
#define _MV_GENERAL_H

/**************************************************************************
 * quand MV_CAUTION est activee, le code va verifier les cas de defaillance
 * et sortir eventuellement des msg des erreurs
 **************************************************************************/
//#define MV_CAUTION 1

/**************************************************************************
 * quand MV_VECTOR_BOUNDS_CHECK est active, le code verifier si le bord de
 * vecteur est depasse
 **************************************************************************/
//#define MV_VECTOR_BOUNDS_CHECK 1

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

// definition de ref_type est pour eviter le msg de warning quand on compile
// le constructeur de vecteur par reference d'un pointeur.
struct MV_Vector_
{
  enum ref_type  { ref = 1};
};

# endif
