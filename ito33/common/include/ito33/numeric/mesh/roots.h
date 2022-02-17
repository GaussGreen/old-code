/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/roots.h
// Purpose:     class Root : useful for the change of variable in the meshes.
// Author:      Nabil
// Created:     2003/12/01
// RCS-ID:      $Id: roots.h,v 1.5 2006/03/22 13:10:22 yann Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/numeric/mesh/roots.h
    @brief class Root : useful for the change of variable in the meshes.
    
 */

#ifndef _ITO33_NUMERIC_MESH_ROOTS_H_
#define _ITO33_NUMERIC_MESH_ROOTS_H_

namespace ito33
{

namespace numeric
{

namespace mesh
{

enum RootType
{
  RootType_Real,
  RootType_Virtual,

  RootType_Max
};

///Declaration of the Root class
class Root
{
public:
  double
    m_dRoot;
  
  RootType
    m_iRootType;

  Root():m_iRootType(RootType_Virtual){};
  ~Root(){};

};

} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESH_ROOTS_H_
