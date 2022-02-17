/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/numeric/meshparams.h
// Purpose:     special numerical parameters class for mesh generation
// Author:      Wang
// Created:     2004/02/09
// RCS-ID:      $Id: meshparams.h,v 1.5 2004/11/10 16:04:55 afrolov Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
  @file ito33/numeric/meshparams.h
  @brief definition of the class for mesh generation Parameters 
  */
  
#ifndef _ITO33_NUMERIC_MESHPARAMS_H_
#define _ITO33_NUMERIC_MESHPARAMS_H_

#include <stddef.h>

namespace ito33
{

namespace numeric
{

/// numerical parameters class for mesh generation
class MeshParams {

public:

  /// ctor
  MeshParams() :
    m_bUniformTimeGrid(false),
    m_bUniformSpaceGrid(false),
    m_nTimeAccumulation(2),
    m_nTimeCompression(2),
    m_nSpaceAccumulation(2),
    m_nSpaceCompression(10),
    m_dTimeStretch(2.),
    m_dSpaceStretch(1.25),
    m_dGridSpan(4)
    {
    }

  /** 
    Instructs the model to use a pseudo-uniform time grid.
    
    If this method is not called a non uniform time grid will be used.
     
    @param bUniformTimeGrid if true a uniform time gride will be used.
   */
  void SetUniformTimeGrid(bool bUniformTimeGrid)
  {
    m_bUniformTimeGrid = bUniformTimeGrid;
  }
 
  /** 
    Instructs the model to use a uniform space grid.
   
    If this method is not called a non uniform space grid will be used.
    
    @param bUniformSpaceGrid if true a uniform space gride will be used.
   */
  void SetUniformSpaceGrid(bool bUniformSpaceGrid)
  {
    m_bUniformSpaceGrid = bUniformSpaceGrid;
  }

  
  /** 
    Number of fine time intervals around each time event.
   
    @param nTimeAccumulation the number of fine intervals around each
      time event.
   */  
  void SetTimeAccumulation(size_t nTimeAccumulation)
  {
    m_nTimeAccumulation = nTimeAccumulation;
  }


  /** 
    The relative size of refined and coarse intervals of the time grid.
   
    @param nTimeCompression the relative size of fine and coarse zones
      of the time grid.
   */
  void SetTimeCompression(size_t nTimeCompression)
  {
    m_nTimeCompression = nTimeCompression;
  }


  /** 
    Number of fine space intervals around each space event.
   
    @param nSpaceAccumulation the number of fine intervals around
      each space event.
   */  
  void SetSpaceAccumulation(size_t nSpaceAccumulation)
  {
    m_nSpaceAccumulation = nSpaceAccumulation;
  }


  /** 
    The relative size of refined and coarse intervals of the space grid.
   
    @param nSpaceCompression the relative size of fine and coarse
      zones of the space grid.
   */  
  void SetSpaceCompression(size_t nSpaceCompression)
  {
    m_nSpaceCompression = nSpaceCompression;
  }


  /** 
    This parameter guides the way the time grid passes from refined 
    to non refined zones.
   
    @param dTimeStretch the refinement "speed" for the time grid.
   */
  void SetTimeStretch(double dTimeStretch){ m_dTimeStretch = dTimeStretch; }
  
  /** 
    This parameter guides the way the space grid passes from refined 
    to non refined zones.
   
    @param dSpaceStretch the refinement "speed" for the space grid.
   */  
  void SetSpaceStretch(double dSpaceStretch)
  {
    m_dSpaceStretch = dSpaceStretch;
  }


  /** 
    Defines the size of the space grid.
   
    The grid will be centered on a certain point and have this number 
    of standard deviations on each side.
   
    @param dGridSpan The number of standard deviations.
   */  
  void SetGridSpan(double dGridSpan)
  {
    m_dGridSpan = dGridSpan; 
  }

  
  /**
    Get the UniformTimeGrid flag value

    @return UniformTimeGrid flag. If it is false, the model 
    makes a make a non-uniform time grid.
    */
  bool GetUniformTimeGrid() const { return m_bUniformTimeGrid; }

  /**
    Get the GetUniformSpaceGrid flag value

    @return UniformSpaceGrid flag. If it is false, the model 
    makes a make a non-uniform space grid.
    */
  bool GetUniformSpaceGrid() const { return m_bUniformSpaceGrid; }  

  /**
    Get the Number of fine time intervals around each time event.

    @return the Number of fine time intervals around each time event.
    */
  size_t GetTimeAccumulation() const { return m_nTimeAccumulation; }

  /**
    Get the relative size of refined and coarse intervals of the time grid.

    @return the relative size of refined and coarse intervals of the time grid.
    */
  size_t GetTimeCompression() const { return m_nTimeCompression; }

  /**
    Get the number of fine space intervals around each space event.
    
    @return the number of fine space intervals around each space event.
    */
  size_t GetSpaceAccumulation() const { return m_nSpaceAccumulation; }

  /** 
    Get The relative size of refined and coarse intervals of the space grid.
   
    @return the relative size of fine and coarse zones of the space grid.
   */ 
  size_t GetSpaceCompression() const { return m_nSpaceCompression; }

  /** 
    Get the refinement "speed" for the time grid
   
    @return the refinement "speed" for the time grid.
   */
  double GetTimeStretch() const { return m_dTimeStretch; }

  /**
    Get the the refinement "speed" for the space grid.

    @return the refinement "speed" for the space grid.
    */
  double GetSpaceStretch() const { return m_dSpaceStretch; }


  /**
    Get the value of grid span

    @return the value of grid span
    */
  double GetGridSpan() const { return m_dGridSpan; }

protected:
  bool 
    /// whether the time grid is uniform
    m_bUniformTimeGrid,
    /// whether the space grid is uniform
    m_bUniformSpaceGrid;
  
  size_t
    /// Maximum number of fine space intervals surrounding a time event.
    m_nTimeAccumulation,		
    /// Grid refinement ratio around a time event.
    m_nTimeCompression,		
    /// Maximum number of fine space intervals surrounding a space event.
    m_nSpaceAccumulation,
    /// Grid refinement ratio around a spot event.
    m_nSpaceCompression;
   
  double    
    /// time stretch
    m_dTimeStretch,
    /// space stretch
    m_dSpaceStretch,
    /// grid span
    m_dGridSpan;
  
}; // class MeshParams


} // namespace numeric

} // namespace ito33

#endif // #ifndef _ITO33_NUMERIC_MESHPARAMS_H_
