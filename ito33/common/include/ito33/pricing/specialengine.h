/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/specialengine.h
// Purpose:     Main PDE solver class
// Author:      David Pooley
// Created:     2003/08/13
// RCS-ID:      $Id: specialengine.h,v 1.6 2005/12/30 12:10:07 nabil Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/specialengine.h
    @brief Main PDE solver class

    Template for the engine classes used to solve PDEs (general framework).
 */

#ifndef _ITO33_PRICING_SPECIALENGINE_H_
#define _ITO33_PRICING_SPECIALENGINE_H_

#include "ito33/pricing/engine.h"

// disable MSC warning: 
// identifier was truncated to '255' characters in the debug information
#ifdef _MSC_VER
#   pragma warning(disable:4786)
#   pragma warning(push, 1)
#endif

namespace ito33{
 
namespace pricing{

/**
  \brief The compute engine
  
  It contains by reference a <em>params</em>, a <em>mesh</em>, a <em>inst_data</em>, 
  a <em>stepper</em> and a <em>num_output</em>.

  Engine solves a time discreted PDE <b>once</b>, with a <b>given</b> <em>mesh</em> made
  on a given input <em>params</em>.

  <b>REQUIRE</b>: The user of this class must garantee that the classes Params, Meshes, 
    InstData and Stepper are compatible between them.
  - Meshes class
    - It has the reference to the objet <em>params</em>. It has a mesh, but it is more
      than a mesh holder. In fact, it is a mesh manager and parameter manager.
    - Meshes class must have a constructor Meshes::Meshes(ParamsForMesh &) and Params is
      one of the derived classes of ParamsForMesh class.
    - It must offer GoAhead() function. This function go one step ahead, if possible, in its
      time mesh and update its related data. It is a kind of iterator. It is called by
      Engine::TryGoAhead().
  - InstData class
    - This class holds three kinds of data for every time step: instanous parameters to feed
      stationnal system solver Stepper (ex. the interest rate r), time independant data
      that Stepper solves (ex. the price arrays of the last and the current time period),
      data that engine needs and that will be passed to NumOutput which is the real numerical
      output container.
    - It must have a constructor InstData::InstData(ParamsForID&, MeshForID&) and Params
      is one of the dervied classes of ParamsForID class and that Meshes is one of the derived
      class of MeshForID class.
    - It must offer Init() function which does the allocation and initilization work at the
      very beginning (called by Engine::PrepareForTimestepping() function).
    - It must offer SetInitialCondtion() function which defines the start condition for
      the time stepping. Please note that it is different from the Init() function.
    - It must offer UpdateMeBeforeStep() function which is called by Engine::BeforeStep()
      at the beginning
      of every time step. This function update <em>inst_data</em> by getting from
      <em>mesh(manager)</em> and <em>params</em> the information of the current period.
    - It must offer DoEvents() function which is called by Engine::DoStep() function and after 
      Stepper::Run() in it. This function updates the data by applying the contraints.
  - Stepper class
    - This is just a (non)linear system solver.
    - It must have a constructor Stepper(ParamsForStepper&, MeshForStepper&, InstDataForStepper&)
      and Params, Meshes and InstData classes are respectively derived from ParamsForStepper
      MeshForStepper and InstDataForStepper class.
    - It offers Init() function which does the allocation and initilization work at the
      very beginning (called by Engine::PrepareForTimestepping() function).
    - It offers Run function that solve the system and update <em>inst_data</em>.
  - NumOutput class
    - This is class which store the numercial output data. That means, we don't want to give
      some of them to the end user who needs just a certain ModelOutput objet.
    - This class must offer a function UpdateMe() that takes <em>inst_data</em> as parameter.
      It is called by Engine::AfterStep().
    - Note that the object is physically allocated outside Engine. Why Engine
      doesn't create its own NumOutput and return it by Run() function? The reason is that
      for a problem composed by several parts, we may need to make an engine runs for each
      part and fill the NumOuptut one by one.
*/
template <class T, class T1, class T2, class T3, class T4>
class SpecialEngine : public Engine<T, T1, T2, T3, T4>
{
public:
  /// type of params
  typedef T Params;
  /// type of Meshmanager
  typedef T1 Meshes;
  /// type of InstData
  typedef T2 InstData;
  /// type of stepper
  typedef T3 Stepper;
  /// type of numerical output
  typedef T4 NumOutput;

  /**
    \brief constructor
  */
  SpecialEngine
      (
        Params &params,
        Meshes &mesh,
        InstData &instdata, 
        Stepper &stepper,
        NumOutput &numoutput
      ) 
    : Engine<Params, Meshes, InstData, Stepper, NumOutput>
                    (params, mesh, instdata, stepper, numoutput) 
  {
  }

  /**
    \brief destructor
  */
  ~SpecialEngine()
  {
  }

  /**
  \brief It runs the engine.

  It runs the engine and solves the problem once with the
  given params and mesh.

  REQUIRE: m_params and m_mesh must have been validated
  */
  void Run()
  {
    // allocates memory for InstData
    // Setup the initial conditions for InstData
    PrepareForTimestepping();

    // Loop until timestepping is done
    while (TryGoAhead())
    {
      if(m_meshes.IsEndOfGrid())
      {
        m_instdata.UpdateAtEndOfGrid();

        m_numoutput.UpdateMeAtEndOfGrid(m_instdata, m_meshes.GetTime());
      }
      else
        GoAhead();
    }

    m_numoutput.Finalize(m_instdata);
  }

private:
  // bring non-dependent stuff used by Run() in this scope to fix compilation
  // with compilers properly implementing 2 phase lookup (such as g++ 3.4)
  typedef Engine<Params, Meshes, InstData, Stepper, NumOutput> BaseEngine;

  using BaseEngine::PrepareForTimestepping;
  using BaseEngine::TryGoAhead;
  using BaseEngine::GoAhead;

  using BaseEngine::m_meshes;
  using BaseEngine::m_instdata;
  using BaseEngine::m_numoutput;


  SpecialEngine& operator=(const SpecialEngine&);

};


}  // namespace pricing

}  // namespace ito33


// Enable warning 4786
#ifdef _MSC_VER
#   pragma warning(pop)
#endif

#endif // #ifndef  _ITO33_PRICING_SPECIALENGINE_H_
