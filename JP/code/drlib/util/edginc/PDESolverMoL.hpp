//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDESolverMoL.hpp
//
//   Description : Wrapper around NAG's 'method of lines' PDE solver
//
//   Date        : 20 December 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"

#ifndef EDR_PDESOLVER_MOL_HPP
#define EDR_PDESOLVER_MOL_HPP

DRLIB_BEGIN_NAMESPACE

/** Wrapper class around NAg's 'method of lines' PDE solver.
    The algorithm used by NAg is essentially the same as the 
    publicly available PDECHEB */
class UTIL_DLL PDESolverMoL: public CObject{
public:
    static CClassConstSP const TYPE;
    friend class PDESolverMoLHelper;

    /** PDE class for the evalutation of the time derivative of
        the n functions making up the system of PDEs as well as for 
        the evaluation of the boundary conditions */
    class UTIL_DLL PDE: public CObject{
    public:
        static CClassConstSP const TYPE;
        /** Given x, returns the initial values u_k */
        virtual void init(double       x,
                          DoubleArray& u) const = 0;

        /** Non-terminal exception. If thrown from one of the three methods
            declared below, an attempt will be made by the algorithm to 
            decrease the time step and re-evaluate the method */ 
        class UTIL_DLL Exception: public ModelException{
        public:
            Exception(const string& routine, const string& msg):
            ModelException(routine, msg){}
        };

        /** Given the value of the k-th function u_k and of its 1st order 
            space derivatives ux_k at the node (x, t), evaluate the quantity
            p, q and r defined by 
                Sum_j P_i,j * du_j/dt + Q_i = dR_i/dx
            where i, j = 0,.., npdes-1 */
        virtual void func(double              x,
                          double              t,
                          const CDoubleArray& u,
                          const CDoubleArray& ux,
                          CDoubleMatrix&      P,
                          CDoubleArray&       Q,
                          CDoubleArray&       R) const = 0;

        /** Evaluate the coefficients of the lower boundary condition at time t. 
            The boundary condition has to be of the form
                        beta * R = gamma         */
        virtual void lowerBound(double               t,
                                const CDoubleArray&  u,
                                const CDoubleArray&  ux,
                                CDoubleArray&        beta,
                                CDoubleArray&        gamma) const = 0;

        /** Evaluate the coefficients of the upper boundary condition at time t. 
            The boundary condition has to be of the form
                        beta * R = gamma         */
        virtual void upperBound(double               t,
                                const CDoubleArray&  u,
                                const CDoubleArray&  ux,
                                CDoubleArray&        beta,
                                CDoubleArray&        gamma) const = 0;

        int getNbPDEs() const;

        virtual ~PDE(){}

    private:
        static void load(CClassSP& clazz);

    protected:
        // For Inheritance
        PDE(const CClassConstSP& objClass);
        PDE(const CClassConstSP& objClass, int npdes);

        int npdes;     // nb of pdes
    };
    typedef smartConstPtr<PDE> PDEConstSP;
    typedef smartPtr<PDE> PDESP;
    
    PDESolverMoL(double acc,
                 int    npoly,
                 PDESP  pde);

    virtual void validatePop2Object();

    /** Given 
            - a system of PDEs with boundary conditions
              and initial values,
            - a time interval [tstart, tend], tstart < tend
            - an array of break points 'xbkpts',
        solves the systems of PDEs. 
        On output, 
            - the solutions 'ubkpts' at the break points */
    void solve(double                tstart,
               double                tend,
               const DoubleArray&    xbkpts,
               DoubleArrayArray&     ubkpts);    // must be of size npde * nbkpts

    void solve1D(double                tstart,
                 double                tend,
                 const DoubleArray&    xbkpts,
                 DoubleArray&          ubkpts);    // must be of size nbkpts

    // chebyshev method to discretize the pde 
    void solveChebyshev(double               tstart,
                       double                tend,
                       const DoubleArray&    xbkpts,
                       DoubleArrayArray&     ubkpts);    // must be of size npde * nbkpts

    void solve1DChebyshev(double                tstart,
                         double                tend,
                         const DoubleArray&    xbkpts,
                         DoubleArray&          ubkpts);    // must be of size nbkpts


    // help the compiler work out how to destruct
    // the hidden Implementation class - needs to be 
    // done from .cpp file
    virtual ~PDESolverMoL();

private:
    PDESolverMoL();

public:
    void pdedef(int        &npde, 
                double      &t, 
                double      &x, 
                double      u[], 
                double      ux[], 
                double      p[], 
                double      q[], 
                double      r[], 
                int         &ires);

    void bndary(int   &npde, 
                double  &t, 
                double  u[], 
                double  ux[], 
                int     &ibnd, 
                double  beta[], 
                double  gamma[], 
                int     &ires);

    void uinit(int             npde, 
                int             npts, 
                const double    x[], 
                double          u[]);

private:
     void resize(int npde);

    // registered vars
    double      acc;
    int         npoly;   // degree of the Chebyshev polynomial
    PDESP       pde;

    // transient variable
    DoubleArray     u_iptl;
    DoubleArray     ux_iptl;
    CDoubleMatrix   p_iptl;
    CDoubleArray    q_iptl;
    CDoubleArray    r_iptl;
    DoubleArray     beta_iptl;
    DoubleArray     gamma_iptl;
};
    
typedef smartPtr<PDESolverMoL> PDESolverMoLSP;
typedef smartConstPtr<PDESolverMoL> PDESolverMoLConstSP;
#ifndef QLIB_PDESOLVERMOL_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<PDESolverMoL>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<PDESolverMoL>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<PDESolverMoL>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<PDESolverMoL>);
#endif


DRLIB_END_NAMESPACE

#endif
