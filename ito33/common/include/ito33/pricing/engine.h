/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/pricing/engine.h
// Purpose:     Main PDE solver class
// Author:      David Pooley
// Created:     2003/08/13
// RCS-ID:      $Id: engine.h,v 1.26 2005/03/09 10:53:13 nabil Exp $
// Copyright:   (c) 2003-2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/pricing/engine.h
    @brief Main PDE solver class

    Template for the engine classes used to solve PDEs (general framework).
 */

#ifndef _ITO33_PRICING_ENGINE_H_
#define _ITO33_PRICING_ENGINE_H_

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
  
  It contains by reference a <em>params</em>, a <em>mesh</em>, 
  a <em>inst_data</em>, a <em>stepper</em> and a <em>num_output</em>.

  Engine solves a time discreted PDE <b>once</b>, with a <b>given</b> 
  <em>mesh</em> made on a given input <em>params</em>.

  <b>REQUIRE</b>: The user of this class must garantee that the classes Params, 
    Meshes, InstData and Stepper are compatible between them.
  - Meshes class
    - It has the reference to the objet <em>params</em>. It has a mesh, but it 
      is more than a mesh holder. In fact, it is a mesh manager and parameter 
      manager.
    - Meshes class must have a constructor Meshes::Meshes(ParamsForMesh &) and 
      Params is one of the derived classes of ParamsForMesh class.
    - It must offer GoAhead() function. This function go one step ahead, if 
      possible, in its time mesh and update its related data. It is a kind of 
      iterator. It is called by Engine::TryGoAhead().
  - InstData class
    - This class holds three kinds of data for every time step: instantaneous 
      parameters to feed stationnal system solver Stepper (ex. the interest 
      rate r), time independant data that Stepper solves (ex. the price arrays 
      of the last and the current time period), data that engine needs and that 
      will be passed to NumOutput which is the real numerical output container.
    - It must have a constructor InstData::InstData(ParamsForID&, MeshForID&) 
      and Params is one of the dervied classes of ParamsForID class and that 
      Meshes is one of the derived class of MeshForID class.
    - It must offer Init() function which does the allocation and initilization 
      work at the very beginning (called by Engine::PrepareForTimestepping() 
      function).
    - It must offer SetInitialCondtion() function which defines the start 
      condition for the time stepping. Please note that it is different from 
      the Init() function.
    - It must offer UpdateMeBeforeStep() function which is called by 
      Engine::BeforeStep() at the beginning of every time step. This function 
      update <em>inst_data</em> by getting from <em>mesh(manager)</em> and 
      <em>params</em> the information of the current period.
    - It must offer DoEvents() function which is called by Engine::DoStep() 
      function and after Stepper::Run() in it. This function updates the data 
      by applying the constraints.
  - Stepper class
    - This is just a (non)linear system solver.
    - It must have a constructor Stepper(ParamsForStepper&, MeshForStepper&, 
      InstDataForStepper&) and Params, Meshes and InstData classes are 
      respectively derived from ParamsForStepper MeshForStepper and 
      InstDataForStepper class.
    - It offers Init() function which does the allocation and initilization 
      work at the very beginning (called by Engine::PrepareForTimestepping() 
      function).
    - It offers Run function that solve the system and update 
      <em>inst_data</em>.
  - NumOutput class
    - This is the class which store the numercial output data. That means, we 
      don't want to give some of them to the end user who needs just a certain 
      ModelOutput objet.
    - This class must offer a function UpdateMe() that takes <em>inst_data</em> 
      as parameter. It is called by Engine::AfterStep().
    - Note that the object is physically allocated outside Engine. Why Engine
      doesn't create its own NumOutput and return it by Run() function? 
      The reason is that for a problem composed by several parts, we may need 
      to make an engine runs for each part and fill the NumOuptut one by one.
*/
template <class T, class T1, class T2, class T3, class T4>
class Engine
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
  Engine(Params &params, Meshes &mesh, InstData &instdata, 
         Stepper &stepper, NumOutput &numoutput) 
    : m_params(params), m_meshes(mesh),
      m_instdata(instdata),
      m_stepper(stepper),
      m_numoutput(numoutput)
  {
  }

  /**
    \brief destructor
  */
  ~Engine()
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
    // allocates memeory for InstData
    // Setup the initial conditions for InstData
    PrepareForTimestepping();

    // Loop until timestepping is done
    while (TryGoAhead())
      GoAhead();

    m_numoutput.Finalize(m_instdata);
  }

  /**
    Go one step ahead.

    The user can call directly this function for a complicated problem, like
    path dependant one. 
    
    Note that, if finally we don't have the need to call directly BeforeStep(),
    DoStep() or AfterStep() from outside, we can in the next stage remove these
    three functions.
    */
  void GoAhead()
  {
    BeforeStep();

    DoStep();

    AfterStep();
  }

  /**
  \brief prepares for the time stepping
  
  It is a part of Run() function. It sets up m_stepper the
  time independent PDE solver and sets the initial values 
  for m_instdata, that m_stepper will work on. Then it sets up 
  m_numoutput by m_instdata.

  RULE: when the user doesn't call Run(), this function must
  be called before any other public functions of Engine class
  */
  void PrepareForTimestepping()
  {
    m_instdata.Init();
    m_stepper.Init();
    m_numoutput.Init(m_instdata);

    m_meshes.SetInitialState();

    m_instdata.SetInitialValue();
    m_numoutput.UpdateMe(m_instdata, m_meshes.GetTime());
    
  }

  /**
  \brief tells whether the time stepping can go ahead

  In fact, it calls the TryGoAhead() function of m_meshes.\n
  It returns false when the time stepping reaches its end.\n
  It returns true when m_meshes succeeds in going one step ahead and update 
  its own parameters.
  */
  bool TryGoAhead()
  {
    return m_meshes.TryGoAhead();
  }

  /**
  \brief prepares the instantaneous data for each time step

  This function must be called the first at every time step.
  */
  void BeforeStep()
  {
    m_instdata.UpdateBeforeStep();
  }

  /**
     calculates the new data values for the current period

     It runs at first m_stepper the PDE solver, checks if there has events at
     current time, if yes then updates numoutput, then deals with the
     events of the current time step.

     REQUIRE: BeforeStep() must have been called before.
   */
  void DoStep()
  {
    m_stepper.Run();

    if ( m_params.HasEventsNow() )
      m_numoutput.UpdateMe(m_instdata, m_meshes.GetTime());

    m_instdata.DoEvents();
  }

  /**
  \brief save what we need to (that has just been calculated) in
  m_numoutput
  */
  void AfterStep()
  {
    m_numoutput.UpdateMe(m_instdata, m_meshes.GetTime());
  }


protected:

  /// reference to Params object
  Params &m_params;

  /// reference to MeshManager object
  Meshes &m_meshes;

  /// reference to stepper object
  Stepper &m_stepper;

  /// reference to InstData object
  InstData &m_instdata;

  /// reference to NumOutput object
  NumOutput &m_numoutput;
  

private:

  Engine& operator=(const Engine&);

};


}  // namespace pricing

}  // namespace ito33


// Enable warning 4786
#ifdef _MSC_VER
#   pragma warning(pop)
#endif

#endif // #ifndef  _ITO33_PRICING_ENGINE_H_
