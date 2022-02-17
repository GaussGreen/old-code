/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/cbmesh.h
// Purpose:     Some useful function for the meshes in the CB case.
// Author:      Nabil
// Created:     2003/12/01
// RCS-ID:      $Id: cbmesh.h,v 1.4 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#ifndef _ITO33_NUMERIC_MESH_CBMESH_H_
#define _ITO33_NUMERIC_MESH_CBMESH_H_

namespace ito33
{

namespace numeric
{

namespace mesh
{

/*
  given a grid _pdOld of _iNbOld nodes, the function WidenGrid widen the 
  grid to _dMinLeft and _dMaxRight.
  */
int WidenGrid(double *_pdOld, int _iNbOld,
               double _dRhoLeft, double _dDeltaMaxLeft, double _dMinLeft,
               double _dRhoRight, double _dDeltaMaxRight, double _dMaxRight,
               double **_ppdNew, int &_iNbNew, int &_iIndiceInitial);

int WidenGrid(double *_pdOld, int _iNbOld,
               double _dRhoLeft, double _dDeltaMaxLeft, double _dMinLeft,
               double _dRhoRight, double _dDeltaMaxRight, double _dMaxRight,
               double **_ppdNew, int &_iNbNew, int &_iIndiceInitial, 
               double _dXRight);


} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif
