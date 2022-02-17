//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDESolverMoL.cpp
//
//   Description : Wrapper around NAG's 'method of lines' PDE solver
//
//   Date        : 20 December 2004
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_PDESOLVERMOL_CPP
#include "edginc/Class.hpp"
#include "edginc/PDESolverMoL.hpp"
#include "edginc/imslerror.hpp"
#include "edginc/Maths.hpp"
#include <nagmk19.hxx>
#include "edginc/Addin.hpp"


DRLIB_BEGIN_NAMESPACE

#define xxxDEBUG

#define IDX2(i,j) (npde * (j) + (i))
#define IDX3(i,j,k) (npde * (npde * (k) + (j)) + (i))


static PDESolverMoL* dummyPtr;

void __stdcall pdedefWRAPPER(int         &npde, 
                             double      &t, 
                             double      &x, 
                             double      u[], 
                             double      ux[], 
                             double      p[], 
                             double      q[], 
                             double      r[], 
                             int         &ires) {
    dummyPtr->pdedef(npde, t, x, u, ux, p, q, r, ires);
}

void __stdcall bndaryWRAPPER(int     &npde, 
                             double  &t, 
                             double  u[], 
                             double  ux[], 
                             int     &ibnd, 
                             double  beta[], 
                             double  gamma[], 
                             int     &ires){
    dummyPtr->bndary(npde,t,u,ux,ibnd,beta,gamma,ires);
}

void __stdcall pdedefChebyshevWRAPPER(int         &npde, 
                                     double      &t, 
                                     double      x[], 
                                     int         &nptl,
                                     double      u[], 
                                     double      ux[], 
                                     double      p[], 
                                     double      q[], 
                                     double      r[], 
                                     int         &ires) {

    DoubleArray         u_tmp(npde), ux_tmp(npde), q_tmp(npde), r_tmp(npde);
    DoubleMatrix        p_tmp(npde, npde);

    for (int i=0; i<nptl; i++){
        for (int j=0; j<npde; j++){
            u_tmp[j] = u[i+j*nptl];
            ux_tmp[j] = ux[i+j*nptl];
        }
        dummyPtr->pdedef(npde, t, x[i], &u_tmp[0], &ux_tmp[0], &p_tmp[0][0], &q_tmp[0], &r_tmp[0], ires);
        for (j=0; j<npde; j++){
            q[i+j*nptl] = q_tmp[j];
            r[i+j*nptl] = r_tmp[j];
            for (int k=0; k<npde; k++){
                p[i+j*nptl+k*npde*nptl] = p_tmp[k][j];
            }
        }
    }
}

void __stdcall uinitChebyshevWRAPPER(int     &npde, 
                                     int     &npts, 
                                     double  x[],
                                     double  u[]){
    dummyPtr->uinit(npde,npts,x,u);
}

int PDESolverMoL::PDE::getNbPDEs() const{
    return npdes;
}

PDESolverMoL::PDE::PDE(const CClassConstSP& objClass):
CObject(objClass), npdes(1) {}

PDESolverMoL::PDE::PDE(const CClassConstSP& objClass, int npdes):
CObject(objClass), npdes(npdes) {}

void PDESolverMoL::PDE::load(CClassSP& clazz){
    REGISTER(PDESolverMoL::PDE, clazz);
    SUPERCLASS(CObject);
    FIELD_INLINE(npdes, "number of PDEs");
    FIELD_MAKE_OPTIONAL(npdes);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const PDESolverMoL::PDE::TYPE = CClass::registerClassLoadMethod(
    "PDESolverMoL::PDE", typeid(PDESolverMoL::PDE), PDESolverMoL::PDE::load);


PDESolverMoL::~PDESolverMoL(){}

void PDESolverMoL::solve(double                tstart,
                         double                tend,
                         const DoubleArray&    xbkpts,
                         DoubleArrayArray&     ubkpts){
    static const string method = "PDESolverMoL::solve";

    try{
        // some validation
        if (!Maths::isPositive(tend - tstart)){
            throw ModelException(method, 
                                 "tstart must be < tend; got "
                                 + Format::toString(tstart) + " and "
                                 + Format::toString(tend) + ", resp.");
        }
        int npts = xbkpts.size();
        if (npts < 3){
            throw ModelException(method,
                                 "nb of break points ("
                                 + Format::toString(npts)
                                 + ") must be >= 3");
        }
        int npde = pde->getNbPDEs();
        int nw = (10 + 6 * npde) * npde * npts
                 + (21 + 3 * npde) * npde
                 + 7 * npts
                 + 54;
        int niw = npde * npts + 24;

        int m = 0; // cartesian coordinates
        DoubleArray w(nw);
        IntArray    iw(niw);
        int itask = 1;  // normal computation 
#ifdef DEBUG
        int itrace = 3;
#else
        int itrace = 0;    
#endif
        int ind = 0;    // start integration from scratch
        int ifail = -1;  //soft fail

        double t = tstart;
        
        resize(npde);
        
        uinit(npde,
                npts,
                &xbkpts[0],
                &ubkpts[0][0]);

        try{
            dummyPtr = const_cast<PDESolverMoL*>(this);
        
            D03PCF(npde,
                m,
                t,
                tend,
                pdedefWRAPPER, 
                bndaryWRAPPER,
                &ubkpts[0][0],
                npts,
                const_cast<double *> (&xbkpts[0]),
                acc,
                &w[0],
                nw,
                &iw[0],
                niw,
                itask,
                itrace,
                ind,
                ifail);

            dummyPtr = 0;
        
            if (0 != ifail){
                throw ModelException("NAG routine D03PCF : Failed with IFAIL = " + Format::toString(ifail));
            }      
        }
        catch(exception& e){
            dummyPtr = 0;
            throw ModelException::addTextToException(e, 
                                                     "NAG routine D03PCF : Failed at time t = "
                                                     + Format::toString(t));
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }

}

void PDESolverMoL::solve1D(double                tstart,
                           double                tend,
                           const DoubleArray&    xbkpts,
                           DoubleArray&          ubkpts){
    static const string method = "PDESolverMoL::solve1D";
    try{
        if (pde->getNbPDEs() != 1){
            throw ModelException(method,
                                "PDE must be of dimenstion 1; not "
                                + Format::toString(pde->getNbPDEs()));
        }
        DoubleArrayArray u(1);
        int nbkpts = xbkpts.size();
        u[0].resize(nbkpts);

        solve(tstart,
              tend,
              xbkpts,
              u);

        ubkpts = u[0];
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}

void PDESolverMoL::solveChebyshev(double                tstart,
                                  double                tend,
                                  const DoubleArray&    xbkpts,
                                  DoubleArrayArray&     ubkpts){
    static const string method = "PDESolverMoL::solveChebyshev";

    try{
        // some validation
        int npde = pde->getNbPDEs();

        int m = 0; // cartesian coordinates

        if (!Maths::isPositive(tend - tstart)){
            throw ModelException(method, 
                                 "tstart must be < tend; got "
                                 + Format::toString(tstart) + " and "
                                 + Format::toString(tend) + ", resp.");
        }

        int nbkpts = xbkpts.size();
        if (nbkpts < 2){
            throw ModelException(method,
                                 "nb of break points ("
                                 + Format::toString(nbkpts)
                                 + ") must be >= 2");
        }
        
        if ((npoly<1)||(npoly>49)){
            throw ModelException(method,
                                 "degree of Chebyshev polynomial ("
                                 + Format::toString(npoly)
                                 + ") must be an integer between 1 and 49");
        }
        
        int npts = (nbkpts - 1) * npoly + 1;
        DoubleArray     x(npts);
        DoubleMatrix    u(npde,npts);

        int nwkres = 3*(npoly+1)*(npoly+1) + (npoly+1)*(npde*npde+6*npde+nbkpts+1) + 13*npde + 5;
        int lenode = npde*npts*(3*npde*(npoly+1)-2);
        int nw = 11*npde*npts + 50 + nwkres + lenode;
        DoubleArray w(nw);

        int niw = npde * npts + 24;
        IntArray    iw(niw);

        int itask = 1;  // normal computation 

#ifdef DEBUG
        int itrace = 3;
#else
        int itrace = 0;    
#endif

        int ind = 0;    // start integration from scratch

        int ifail = -1;  //soft fail

        double t = tstart;
        
        resize(npde);
        
        try{
            dummyPtr = const_cast<PDESolverMoL*>(this);
        
            D03PDF(npde,
                m,
                t,
                tend,
                pdedefChebyshevWRAPPER, 
                bndaryWRAPPER,
                &u[0][0],
                nbkpts,
                const_cast<double *> (&xbkpts[0]),
                npoly,
                npts,
                &x[0],
                uinitChebyshevWRAPPER,
                acc,
                &w[0],
                nw,
                &iw[0],
                niw,
                itask,
                itrace,
                ind,
                ifail);

            dummyPtr = 0;
        
            if (0 != ifail){
                throw ModelException("NAG routine D03PDF : Failed with IFAIL = " + Format::toString(ifail));
            }    
            
            for (int i=0; i<npde; i++){
                for (int j=0; j<nbkpts; j++){
                    ubkpts[i][j] = u[i][(j*npoly)];
                }
            }
        }
        catch(exception& e){
            dummyPtr = 0;
            throw ModelException::addTextToException(e, 
                                                     "NAG routine D03PDF : Failed at time t = "
                                                     + Format::toString(t));
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void PDESolverMoL::solve1DChebyshev(double                tstart,
                                    double                tend,
                                    const DoubleArray&    xbkpts,
                                    DoubleArray&     ubkpts){
    static const string method = "PDESolverMoL::solve1D";
    try{
        if (pde->getNbPDEs() != 1){
            throw ModelException(method,
                                "PDE must be of dimenstion 1; not "
                                + Format::toString(pde->getNbPDEs()));
        }
        DoubleArrayArray u(1);
        int nbkpts = xbkpts.size();
        u[0].resize(nbkpts);

        solveChebyshev(tstart,
                      tend,
                      xbkpts,
                      u);

        ubkpts = u[0];
    }
    catch (exception& e){
        throw ModelException(e,method);
    }
}


void PDESolverMoL::uinit(int             npde, 
                         int             npts, 
                         const double    x[], 
                         double          u[]){
        DoubleArray u_ipts(npde);
        for (int ipts = 0; ipts < npts; ++ipts){            
            pde->init(x[ipts], u_ipts);
            for (int ipde = 0; ipde < npde; ++ipde){
                u[IDX2(ipde, ipts)] = u_ipts[ipde];
            }
        }
    }

void PDESolverMoL::pdedef(int        &npde, 
                        double      &t, 
                        double      &x, 
                        double      u[], 
                        double      ux[], 
                        double      p[], 
                        double      q[], 
                        double      r[], 
                        int         &ires){
        for (int ipde = 0; ipde < npde; ++ipde){
            u_iptl[ipde] = u[ipde];
            ux_iptl[ipde] = ux[ipde];
        }

        pde->func(x,
                    t, 
                    u_iptl,
                    ux_iptl,
                    p_iptl,
                    q_iptl,
                    r_iptl);   
            
        for (ipde = 0; ipde < npde; ++ipde){
            q[ipde] = q_iptl[ipde];
            r[ipde] = r_iptl[ipde];
            for (int ipde2 = 0; ipde2 < npde; ++ipde2){
                p[IDX2(ipde,ipde2)] = p_iptl[ipde][ipde2];
            }
        }
    }

void PDESolverMoL::bndary(int   &npde, 
                        double  &t, 
                        double  u[], 
                        double  ux[], 
                        int     &ibnd, 
                        double  beta[], 
                        double  gamma[], 
                        int     &ires){        
        for (int ipde = 0; ipde < npde; ++ipde){
            u_iptl[ipde] = u[ipde];
            ux_iptl[ipde] = ux[ipde];
        }

        if (ibnd == 0){
            (*pde).lowerBound(t,
                            u_iptl,
                            ux_iptl,
                            beta_iptl,
                            gamma_iptl);
        }
        else{
            (*pde).upperBound(t,
                            u_iptl,
                            ux_iptl,
                            beta_iptl,
                            gamma_iptl);
        }
 
        for (ipde = 0; ipde < npde; ++ipde){
            beta[ipde] = beta_iptl[ipde];
            gamma[ipde] = gamma_iptl[ipde];
        }
    }

void PDESolverMoL::resize(int npde){
        u_iptl.resize(npde);
        ux_iptl.resize(npde);
        p_iptl = DoubleMatrix(npde, npde);
        q_iptl.resize(npde);
        r_iptl.resize(npde);
        beta_iptl.resize(npde);
        gamma_iptl.resize(npde);
    }


PDESolverMoL::PDESolverMoL():
CObject(TYPE),
acc(0.0001),
npoly(1){}

PDESolverMoL::PDESolverMoL(double acc,
                           int    npoly,
                           PDESP  pde):
CObject(TYPE),
acc(acc),
npoly(npoly),
pde(pde) {
    validatePop2Object();
}

void PDESolverMoL::validatePop2Object(){
    Maths::checkPositive(acc, "acc");
    Maths::checkPositive(npoly, "npoly");
}

class PDESolverMoLHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PDESolverMoL, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD_INLINE(acc, "acc");
        FIELD_INLINE(npoly, "npoly");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new PDESolverMoL();
    }
};

CClassConstSP const PDESolverMoL::TYPE = CClass::registerClassLoadMethod(
    "PDESolverMoL", typeid(PDESolverMoL), PDESolverMoLHelper::load);


class PDESolverAddin:public CObject{
public:
    static CClassConstSP const TYPE;

    /**the 'addin function' */ 
    static IObjectSP PDESolver(PDESolverAddin* params){
        
        int             numCols = params->xbkpts.size();
        int             numRows = params->tend.size();
 
        PDESolverMoL    solver(params->acc,
                                params->npoly,
                                params->pde);
        
        CDoubleMatrix    ubkpts(numCols,numRows);
        CDoubleArray     u(numCols);
     
        for (int i=0; i<numRows; i++) {
            solver.solve1D(params->tstart,
                            params->tend[i],
                            params->xbkpts,
                            u);
            for (int j=0; j<numCols; j++) {
                ubkpts[j][i] = u[j];
            }
        }

        return IObjectSP(copy(&(ubkpts)));
    }

    /**the 'addin function' */ 
    static IObjectSP PDESolverChebyshev(PDESolverAddin* params){
        
        int             numCols = params->xbkpts.size();
        int             numRows = params->tend.size();
 
        PDESolverMoL    solver(params->acc,
                                params->npoly,
                                params->pde);
        
        CDoubleMatrix    ubkpts(numCols,numRows);
        CDoubleArray     u(numCols);
     
        for (int i=0; i<numRows; i++) {
            solver.solve1DChebyshev(params->tstart,
                                    params->tend[i],
                                    params->xbkpts,
                                    u);
            for (int j=0; j<numCols; j++) {
                ubkpts[j][i] = u[j];
            }
        }

        return IObjectSP(copy(&(ubkpts)));
    }

private:
    double                      tstart;
    DoubleArray                 tend;
    DoubleArray                 xbkpts;         
    PDESolverMoL::PDESP         pde;
    double                      acc;
    int                         npoly;

    /** for reflection */ 
    PDESolverAddin():CObject(TYPE), acc(0.0001), npoly(1) {}

    static IObject* defaultPDESolverAddin(){
        return new PDESolverAddin();
    }

    /**Invoke when Class is 'load' */ 
    static void load(CClassSP& clazz){
        REGISTER(PDESolverAddin,clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPDESolverAddin);
        FIELD_INLINE(tstart, "the starting point of time dimension.");
        FIELD_INLINE(tend, "the ending point of time dimension.");
        FIELD_INLINE(xbkpts, "an array of break points in space dimension");
        FIELD(pde, "PDE");
        FIELD_INLINE(acc, "accuracy");
        FIELD_MAKE_OPTIONAL(acc);
        FIELD_INLINE(npoly, "degree of the Chebyshev polynomial");
        FIELD_MAKE_OPTIONAL(npoly);

        Addin::registerClassObjectMethod("PDE_SOLVER",
                            Addin::UTILITIES,
                            "Solve a PDE",
                            TYPE,
                            false,
                            Addin::expandMulti, // returnHandle ?
                            (Addin::ObjMethod*)PDESolver);

        Addin::registerClassObjectMethod("PDE_SOLVER_CHEBYSHEV",
                            Addin::UTILITIES,
                            "Solve a PDE",
                            TYPE,
                            false,
                            Addin::expandMulti, // returnHandle ?
                            (Addin::ObjMethod*)PDESolverChebyshev);
    }
};

CClassConstSP const PDESolverAddin::TYPE = CClass::registerClassLoadMethod(
    "PDESolverAddin", typeid(PDESolverAddin), PDESolverAddin::load);

bool PDESolverAddinLinkIn() {
    return PDESolverAddin::TYPE != 0;
}

DRLIB_END_NAMESPACE
