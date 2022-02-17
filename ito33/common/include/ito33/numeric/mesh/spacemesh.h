/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/mesh/spacemesh.h
// Purpose:     Class for the space mesh
// Author:      WANG Xuewen
// RCS-ID:      $Id: spacemesh.h,v 1.6 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file ito33/numeric/mesh/spacemesh.h
   @brief Decalration for the space mesh
 */

#ifndef _ITO33_NUMERIC_MESH_SPACEMESH_H_
#define _ITO33_NUMERIC_MESH_SPACEMESH_H_

#include "ito33/array.h"

namespace ito33
{

namespace numeric
{

namespace mesh
{


/// A container for a fixed space mesh
class SpaceMesh
{

public:

  SpaceMesh(const double *pdLogS, const double *pdS, size_t nNbS)
          : m_nNbS(nNbS), m_pdS(m_nNbS), m_pdLogS(m_nNbS)
  {
    for (size_t nIdxS = 0; nIdxS < m_nNbS; nIdxS++)
    {
      m_pdLogS[nIdxS] = pdLogS[nIdxS];
      m_pdS[nIdxS] = pdS[nIdxS];
    }
  }

  const double *GetLogS() const { return m_pdLogS.Get(); }

  const double *GetS() const { return m_pdS.Get(); }

  size_t GetNbS() const
  {
    return m_nNbS;
  }

private:

  size_t m_nNbS;

  Array<double> m_pdS;

  Array<double> m_pdLogS;

}; // class SpaceMesh


} // namespace mesh

} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESH_SPACEMESH_H_

