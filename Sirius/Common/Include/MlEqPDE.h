//	MlEqPDE.h :				PDE class.
//
//	Author :				Alex Langnau
/////////////////////////////////////////////////////////////////////////////


#ifndef _MLEQPDE_H
#define _MLEQPDE_H


#include <vector>
#include "MlEqhandles.h"
#include "smart.h"
#include "localvol.h"


class tridiagonal_matrix
{
	friend void tridiagonal_matrix_multiply1(const tridiagonal_matrix& matrix, const double *Vector, double *result);
	friend void tridiagonal_matrix_scale(tridiagonal_matrix& matrix, double scalar);
	friend void tridiagonal_matrix_solve1(const tridiagonal_matrix& matrix, const CVector& rhs, double* soln, CVector& gamma );
	
	public:
		
	tridiagonal_matrix(int nrows=0);
	tridiagonal_matrix(const tridiagonal_matrix &);
	tridiagonal_matrix & operator= (const tridiagonal_matrix&);

	double & operator()(int i,int j);
	double   operator()(int i,int j) const;
		
	inline double size() const { return nrows; }
	inline double below_term() const { return *(pelements-1); }
	inline double above_term() const { return *(pelements+int(3*nrows-2)); }

	void populate(double below, double on, double above);

	void init(int nrows);
	inline void resize(int nrows_) { nrows=nrows_; elements.resize(3*nrows_); }

	CVector elements; // should be protected
	protected:

	double* pelements;

	size_t nrows;
	double zero;
	
};	


void tridiagonal_matrix_populate(tridiagonal_matrix& matrix, double below, double on, double above);	
void tridiagonal_matrix_copy(tridiagonal_matrix& dst, const tridiagonal_matrix& src);	




/************************************************/
/*                                              */
/*           PdeBoundaryConditions              */
/*                                              */
/************************************************/


class  pde_boundary_condition

{

	public:

	int    coord;					// Coordinate system to use: PDE or EXP. 
	double scale;					// Scaling to be applied to mesh coordinate system. 
	double shift;					// Shift applied to mesh coordinate system. 

	double a;						// Formula for condition a.u_xx + b.u_x + c u = d. 
	double b;
	double c;
	double d;

} ;


/************************************************/
/*                                              */
/*           PdeEnums							*/
/*                                              */
/************************************************/


class PdeEnums 
{
public :

	
    enum PdeMode
    {
	    PDE_NORMAL_MODE   = 1,
	    PDE_MODIFIED_MODE 
    };

	
	enum PdeCoordinates
	{
		PDE_COORDINATES_PDE = 1, // "Standard PDE coordinates"
		PDE_COORDINATES_EXP		//"Exponential PDE coordinates"
	};	


    enum PdeOutput
    {
        PDE_NONE = 0,
        PDE_LAST = 1,// mean last layer (CA)
        PDE_FULL = PDE_LAST << 1
    };

    enum PdeThetaDiffusion
    {
        PDE_THETA_DIFFUSION,    // discretise the diffusion generator with the theta-scheme
        PDE_STEP_DIFFUSION      // use the current time step diffusion generator, only apply theta scheme to the value function
    };
	
    enum PdeFactors
    {
        PDE_FACTORS_CONTINUOUS,    // discretise the diffusion generator with the theta-scheme
        PDE_FACTORS_DISCONTINUOUS      // use the current time step diffusion generator, only apply theta scheme to the value function
    };
	
	enum PdeDependent		//"PDE solver x independent flag"
	{
		PDE_X_INDEPENDENT = 1, //	"Independent of 'x'"
		PDE_X_DEPENDENT		   //	"Dependent on 'x'"
	};	

	
	enum PdeModel
	{
		PDE_MODEL_SAME = 1, // "PDE solver model diversity flag"
		PDE_MODEL_DIFFERENT		//"Dependent variable models are different"
	};	
	
	enum PdeDirection		// "PDE solver direction of integration flag"
	{
		PDE_DIRECTION_FORWARDS = 1, //"Integrate forwards"
		PDE_DIRECTION_BACKWARDS		//""Integrate backwards"
	};	
	
	enum PdeKeep		// "PDE solver solution grid retention flag"
	{
		PDE_KEEP_FULL = 1,  //"Keep the full solution grid"
		PDE_KEEP_LAST		//"Keep only the last entries in the solution grid"
	};	


    static inline bool isOutputModeLast(int mode);
    static inline bool isOutputModeFull(int mode);
};

inline bool PdeEnums::isOutputModeLast(int mode)
{
    return ((mode & PDE_LAST) == PDE_LAST);
}

inline bool PdeEnums::isOutputModeFull(int mode)
{
    return ((mode & PDE_FULL) == PDE_FULL);
}

void pde_interpolate(   const CVector& new_grid, const CVector& old_grid, const CVector& u,	CVector& new_u);

class pde_driver;
class MlEqPdeDriver;

class MlEqPdeHelper : public RCObject
{


public:

		double m_spot;

		virtual void pde_modify(int it, double t,  void* ptr, int nd,int nx,
								double lx, double ux, CMatrix& u, int& new_nx,
								double& new_lx, double& new_ux,CMatrix& new_u, MlEqPdeDriver* pde){};

        virtual PdeEnums::PdeMode strategy_function(int it, double t, void* ptr, 
			int nd, int nx, CMatrix& u_it, MlEqPdeDriver* pde){return PdeEnums::PDE_NORMAL_MODE;};


		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u , MlEqPdeDriver* pde){};


		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper, MlEqPdeDriver* pde){};

		virtual void greeks(double& delta,double& gamma,int id,int it,double t, void* ptr, MlEqPdeDriver* pde);


};



/************************************************/
/*                                              */
/*           PdeDriver							*/
/*                                              */
/************************************************/


class  MlEqPdeDriver : public RCObject
{
	
protected:

		friend class MlEqPdeHelper;

		MlEqPdeHelperHandle m_pHelper;

public:

		CMatrix m_saveVolGrid;//[idate][ix]
		CMatrix m_saveSpotGrid;//[idate][ix]

protected:

//		this section contains the current state for the pde

		double						m_dx;					
		double						m_dxdx;
		double						m_one_m_theta;
		bool						m_expdx_set_upper;
		bool						m_expdx_set_lower;			
		double						m_lo_expdx0;           // Mesh spacings at lower/upper boundary in transformed exp space. 
		double						m_lo_expdx1;
		double						m_hi_expdx0;
		double						m_hi_expdx1;
		
		pde_boundary_condition		m_last_lower;			// Previous boundary condition information. 
		pde_boundary_condition		m_last_upper;

		int							m_greekType;// log or spot
		
		CVector						m_convection_factors;	
		CVector						m_diffusion_factors;
		CVector						m_variable_factors;
		CVector						m_rhs_factors;
		CVector						m_residual;			// Contains the ode residual from the previous timestep. 
		
		CMatrix						m_u_prev;				
		CMatrix						m_u_interp;			// Interpolated solution provided by the user when the mesh changes. 
		CVector						m_tmp;					// Temporary space the tridiagonal solver. 
		
		tridiagonal_matrix			m_full_stencil;		// Full linear system to be solved. 
		tridiagonal_matrix			m_spatial_stencil;		// Spatial only factors, stored for next timestep. 
		
		CVector						m_known;				// Right hand side of the linear system. 
		double						m_fac_l;
		double						m_fac_u;
		
		typedef std::vector < tridiagonal_matrix > tridiagonal_vector;

		tridiagonal_vector m_old_stencils;

		CVector						m_old_fac_l;
		CVector						m_old_fac_u;
		CMatrix						m_old_rhs;

        PdeEnums::PdeThetaDiffusion m_thetaDiffusion;

		int							m_last_mode;

protected:


		int                           m_nd;						// Number of dependent variables. 

		int                           m_different;				// Whether models for dependent variables differ. 

		int                           m_max_nx;					// Maximum size nx will get to. 

		int                           m_nx;

		double                        m_lx;						// Integration range in x (lx,ux). 

		double						  m_ux;

		int                           m_nt;

		CVector                       m_t;						// Integration range in t [t[0],t[nt-1]). 

		int                           m_direction;				// Direction of integration in time. 

		int                           m_keep;						// Whether to keep all the solution or the last generated. 

		double                        m_theta;					

		void*                         m_users_data;				// Pointer to user's data. 

		int		                      m_convection_dependence;	

		int		                      m_diffusion_dependence;		

		int		                      m_factor_dependence;	

		int		                      m_constant_dependence;		

		int							  m_continuous;
				
		
		
		GVector< CMatrix> m_u_full;//	u["timestep"]["dependent variable"]["mesh point"].

		CMatrix m_u_last;//	u["dependent variable"]["mesh point"].
		
				
		int    m_curr_nx;

		double m_curr_lx;

		double m_curr_ux;	


	protected:


		void validate_pde();

		void pde_create_state( );

	    void pde_resize_state();

		void pde_spatial_factors( const int it, const double t, const int id);

	    void pde_initial_factors_prev( const int it, const double t);

		void pde_theta_timestep(const int it, const double t, const double dt, CMatrix&  u1, 
									CMatrix&  u2);


//		int pde_nt();
//		int pde_storing();


		CMatrix& pde_last_solution();



	public:

		GVector< CMatrix>& pde_full_solution();

		MlEqPdeDriver(void)  {};

        //pde_driver(const pde_driver& rhs);

		virtual ~MlEqPdeDriver(void) {};


		virtual void reinitialize_nd(int newnd,MlEqPdeHelperHandle pHelper);
		virtual double pde_integrate(CVector& result);
		virtual double pde_integrate(CVector& result,int nt);


	public:

		// diffusion dependent functions

		virtual void convection_fcn(int id, double lx, double ux, int nx, int it, double t, 
													  void* ptr, CVector& f);

		virtual void diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f);

		virtual void factor_fcn(int id, double lx, double ux, int nx, int it, double t, 
												 void* ptr, CVector& f);

		virtual void constant_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f);


	public:

		// payoff dependent functions

		virtual void pde_modify(int it, double t,  void* ptr, int nd,
												   int nx, double lx, double ux, CMatrix& u, int& new_nx, double& new_lx, double& new_ux, 
												   CMatrix& new_u);


        virtual PdeEnums::PdeMode strategy_function(int it, double t, void* ptr, 
													 int nd, int nx, CMatrix& u_it);

		virtual void initial_condition(void* ptr,int nd,int nx,CMatrix& u );


		virtual void boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper);


	public:


		virtual void pde_init(const int nd, const int different, const int nx, const int* max_nx, const double lx, 
								const double ux, const CVector& t, const int direction, const int keep,
								void* ptr, 
								double theta=0.5, 
								int continuous = PdeEnums::PDE_FACTORS_CONTINUOUS,
								int conv_dependence=PdeEnums::PDE_X_INDEPENDENT, 
								int diff_dependence=PdeEnums::PDE_X_INDEPENDENT, 
								int var_dependence=PdeEnums::PDE_X_INDEPENDENT, 
								int rhs_dependence=PdeEnums::PDE_X_INDEPENDENT,
								MlEqPdeHelperHandle pHelperOpt = NULL);


		double pde_x( double i);

        inline int getConvectionDependence(void)const {return m_convection_dependence;}

        inline int getDiffusionDependence(void) const {return m_diffusion_dependence;}

        inline int getFactorDependence(void)    const {return m_factor_dependence;}

        inline int getConstantDependence(void)  const {return m_constant_dependence;}

        inline int getDifferent(void)           const {return m_different;}

        inline void setConvectionDependence(int in){m_convection_dependence = in;}

        inline void setDiffusionDependence(int in) {m_diffusion_dependence = in;}

        inline void setFactorDependence(int in)    {m_factor_dependence = in;}

        inline void setConstantDependence(int in)  {m_constant_dependence = in;}

    	inline void setDifferent( int in)          { m_different = in;}


};	



class  pdeLocalVol  :  public MlEqPdeDriver
{


protected:

		GVector<CVector> m_Grid;//[idate][ix]
		CVector m_loanspread;//[idate]
		CVector m_rate;//[idate]
		CVector m_quantoDrift;//[idate]
		DupireLocalVolHandle m_pLv;


public:


		virtual void convection_fcn(int id, double lx, double ux, int nx, int it, double t, 
													  void* ptr, CVector& f);
		virtual void diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f);
		virtual void factor_fcn(int id, double lx, double ux, int nx, int it, double t, 
												 void* ptr, CVector& f);

		virtual void init(DupireLocalVolHandle lv,MlEqPdeHelperHandle pHelperOpt);


		virtual double priceAndGreeks( double& delta, double& gamma, double& vega );

		
};


class  pdeLocalVolEffective  :  public MlEqPdeDriver
{


protected:

		GVector<CVector> m_Grid;//[idateCalib][ix]
		CVector m_loanspread;//[idate]
		CVector m_rate;//[idate]
		CVector m_forwards;//[idate]

		GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle> m_localVolSlice;// [icalibDates]
		MlEqDateHandle m_startDate;
		GVector<long> m_calibDates;// [icalibDates]
		GVector<long> m_map;//it->icalibDates
		GVector<long> m_mapinv;//icalibDates->it

		double				getlocalVol(int idate,int ixspace);
		void				setlocalVolSlice(int idateCalib);


public:

		const GVector<long>& getMap()const{return m_map;};
		const GVector<long>& getMapInv()const{return m_mapinv;};

	  	void reinitialize(int idateCalib,MlEqAnalyticCurveWithTanhWingEdgeHandle& lv);
		void reinitialize(int idateCalib,const CVector& lvols);


public:


		virtual void convection_fcn(int id, double lx, double ux, int nx, int it, double t, 
													  void* ptr, CVector& f);
		virtual void diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f);
		virtual void factor_fcn(int id, double lx, double ux, int nx, int it, double t, 
												 void* ptr, CVector& f);




		void init(GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle>& localVolSlice,MlEqAsset& asset,
								MlEqDateHandle startDate,long maturityDate,GVector<long>& calibDates,MlEqPdeHelperHandle pHelperOpt,
								const int nd, const int nx, const int* max_nx, const double lx, int nt,
								const double ux, const int keep);
	
};


class  pdeLocalVolBump  :  public pdeLocalVol
{
protected:

		GVector<CVector> m_bumpedGrid;//[idate][ix]


		CVector m_bumpedLoanspread;//[idate]
		CVector m_bumpedRate;//[idate]


public:

		virtual void convection_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f);
		virtual void diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f);
		virtual void factor_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f);

		virtual void init(DupireLocalVolHandle lv,MlEqPdeHelperHandle pHelperOpt, double ls_bump);

		double priceAndGreeks( double& delta, double& gamma, double& vega );
};


// the following case test classes


class  pdeBlackScholes  :  public MlEqPdeDriver
{

		public:

		double m_vol;
		double m_loanspread;
		double m_discount_rate;
		double m_spot;

		virtual void convection_fcn(int id, double lx, double ux, int nx, int it, double t, 
													  void* ptr, CVector& f);
		virtual void diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f);
		virtual void factor_fcn(int id, double lx, double ux, int nx, int it, double t, 
												 void* ptr, CVector& f);

		virtual void constant_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f);

		void init(double vol,double forward,double discountRate,double spot,
				  MlEqDateHandle asOfDate,long maturityDate,
				  int numberSpacialGridPoints,int numberTimeSteps,double numberStdev,
				  MlEqPdeHelper* pHelper,double minSpotOpt,double maxSpotOpt);

};


	
#endif 
	