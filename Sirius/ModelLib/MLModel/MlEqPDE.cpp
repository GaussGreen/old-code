//	MlEqPDE.cpp :              Implementation of the PDE classes.
//
//	Author :				   Alex Langnau
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mleqpde.h"

#include "MlEqMaths.h"
#include "MlEqHandles.h"
#include "MlEqDate.h"
#include "MlEqShortModels.h"



const double epsilon = 1e-15;
#define deqz(x)				(fabs(x) < epsilon)



tridiagonal_matrix::tridiagonal_matrix(int Nrows) 

{
//	elements( -1, 3*Nrows ),
	nrows   = Nrows ;
	zero    = 0.0 ;
	elements.resize(3*Nrows);
	pelements = elements.getPtr()+1;
}


tridiagonal_matrix::tridiagonal_matrix(const tridiagonal_matrix & src)
{
	elements = src.elements;
	nrows    = src.nrows;
	zero     = src.zero;
	pelements = src.pelements;
} 


tridiagonal_matrix & 
tridiagonal_matrix::operator= (const tridiagonal_matrix & src)
{
	if ( this == &src )
		return *this;

	elements = src.elements;
	nrows    = src.nrows;
	zero     = src.zero;
	pelements = src.pelements;

	return *this;
} 


void tridiagonal_matrix::init(int Nrows)
{
	elements.resize(3*Nrows);
	nrows = Nrows;
	zero = 0.0;
	pelements = elements.getPtr()+1;

}


double & tridiagonal_matrix::operator()(int n,
							   int m)
{
	if ( !(( n == m) || (n+1 == m) || (n== m+1) ) )
	{
		throw( "Input Matrix is not tridiagonal");
	}

	return *(pelements+2*n+m);//elements[2*n+m];
}


double tridiagonal_matrix::operator()(int n,
							   int m) const
{
	if ( !(( n == m) || (n+1 == m) || (n== m+1) ) )
	{
		throw( "Input Matrix is not tridiagonal");

		return zero;
	}

	return *(pelements+2*n+m);//elements[2*n+m];
}


void tridiagonal_matrix::populate(double below, 
							 double on, 
							 double above)
{

	for (int i=0, j=-1; i<nrows; i++)
	{
		*(pelements+j) = below;
		j++;
		*(pelements+j) = on;
		j++;
		*(pelements+j) = above;
		j++;
	}


}


///////////////////////////////////////////////////////////////////////////////////////////
void tridiagonal_matrix_copy(      tridiagonal_matrix& dst, 
								  const tridiagonal_matrix& src)

{
	dst = src;
} 


void tridiagonal_matrix_populate(tridiagonal_matrix& matrix, 
									  double              below, 
									  double              on, 
									  double              above)
{
	matrix.populate(below, on, above);
}


void tridiagonal_matrix_multiply1(const tridiagonal_matrix& matrix, 
									   const double            * Vector, 
									   double                  * result)
{
	const double* matrix_ptr;
	double*       vector_ptr;
	double*       result_ptr;


	if ( Vector == NULL ){
		throw( "Second Input is NULL");
	}

	if ( result == NULL ){
		throw( "Third Vector is NULL");
	}

	matrix_ptr = matrix.pelements+2;//&(matrix.elements[0])+2;
	vector_ptr = (double*)Vector + 1;
	result_ptr = result + 1;

	result[0] =  matrix(0,0) * Vector[0];
	result[0] += matrix(0,1) * Vector[1];

	int i;
	for ( i = 1; i < matrix.nrows - 1; i++)
	{
		*result_ptr  = *(matrix_ptr++) * *(vector_ptr - 1);
		*result_ptr += *(matrix_ptr++) * *vector_ptr;
		*result_ptr += *(matrix_ptr++) * *(vector_ptr + 1);
		result_ptr++;
		vector_ptr++;
	}

	result[i]  = matrix(i,i-1) * Vector[i-1];
	result[i] += matrix(i,i  ) * Vector[i];
} 


void tridiagonal_matrix_scale(tridiagonal_matrix& matrix, 
								   double              scalar)
{
	double* ptr;
	
	ptr = matrix.pelements;//&(matrix.elements[0]);	//&(matrix(0,0));
	*ptr *= scalar;
	ptr++;
	*ptr *= scalar;
	ptr++;

	for (int i = 1; i < matrix.nrows - 1; i++)
	{
		*(ptr++) *= scalar;
		*(ptr++) *= scalar;
		*(ptr++) *= scalar;
	}

	*(ptr++) *= scalar;
	*ptr *= scalar;
} 


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//  Function solves matrix * soln = rhs
//  matrix is a n X n tridiagonal matrix ; 
//  rhs: is a vector of the right-hand side
//  soln: is the solution vector
//  tmp: temporary vector array of at least n rows or NULL
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void tridiagonal_matrix_solve1(const tridiagonal_matrix & matrix, 
									const CVector       & rhs, 
									double                   * soln, 
									CVector             & gamma)

{
	double       *gamma_ptr;
	const double *rhs_ptr;
	const double *matrix_ptr;
	double       *soln_ptr;
	
	if ( rhs.getsize() == 0 ){
		throw( "Input Vector is NULL");
	}

	if ( soln == NULL ){
		throw( "Input Vector is NULL");
	}

	int n = matrix.nrows;


	// If the user doesn't provide a workspace we need to. 

	if (gamma.getsize()==0)
	{
		// New vector for elimination factors. 

		gamma.resize(n);
	}



	// Perform the decomposition 

	double beta = matrix(0,0);

	if (deqz(beta))
	{
		throw( "Zero Pivot encountered");
	}

	soln[0] = rhs[0] / beta;


	// The following loop may appear strange but it helps the unrolling and this is a very intensive routine. 

	gamma_ptr  = gamma.getPtr()+1;;
	soln_ptr   = soln + 1;
	matrix_ptr = matrix.pelements+2;
	rhs_ptr    = rhs.getConstPtr()+1;

	for (int i = 1; i < n; i++)
	{


		*gamma_ptr = *(matrix_ptr - 1) / beta;


		beta = *(matrix_ptr + 1) - *matrix_ptr * *gamma_ptr;

		// Separate comparison for beta. 

		if (beta < -epsilon)
		{
			// soln[i] = (rhs[i] - matrix[i][i-1] * soln[i-1]) / beta; 

			*soln_ptr = (*rhs_ptr - *matrix_ptr * *(soln_ptr-1)) / beta;
		}
		else if (beta > epsilon)
		{
			// soln[i] = (rhs[i] - matrix[i][i-1] * soln[i-1]) / beta; 

			*soln_ptr = (*rhs_ptr - *matrix_ptr * *(soln_ptr-1)) / beta;
		}
		else
		{
			throw( "Zero Pivot encountered");
		}

		matrix_ptr += 3;

		gamma_ptr++;
		soln_ptr ++;
		rhs_ptr  ++;
	}


	// Perform back substitution. 

	for (int i = n-2; i >= 0; i--)
		soln[i] = soln[i] - gamma[i+1] * soln[i+1];
} 


void pde_interpolate(   const CVector& new_grid,
						const CVector& old_grid,
						const CVector& u,
						CVector& new_u)
{
	MlEqCubicSplineInterpolator interpolator;
	interpolator.initialize( old_grid, u );

	int new_nx = new_grid.getsize();
	new_u.resize(new_nx);

	for(int i=0; i<new_nx; ++i){
		new_u[i] = interpolator.getValue(new_grid[i]);
	}
}


/****************************************************************
**	Class  : pde_driver 
**	Routine: diffusion_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void MlEqPdeDriver::constant_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f)
{
	f[0] = 0.0;
} 


void MlEqPdeDriver::factor_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f)
{
	f[0] = 0.0;
} 


void MlEqPdeDriver::convection_fcn(int id, double lx, double ux, int nx, int it, double tt, void* ptr, CVector& f)
{
	f[0] = 0.0;
} 

void MlEqPdeDriver::diffusion_fcn(int id, double lx, double ux, int nx, int it, double tt, void* ptr, CVector& f)
{
	f[0] = 0.0;
} 


PdeEnums::PdeMode MlEqPdeDriver::strategy_function(int it, double t,  void* ptr, int nd, int nx, CMatrix& u_it)
{
	if ( !!m_pHelper )
	{
		return m_pHelper->strategy_function(it,t,ptr, nd,nx,u_it,this);
	}

    return PdeEnums::PDE_NORMAL_MODE;
} 


void MlEqPdeDriver::pde_modify(int it, double t,  void* ptr, int nd,

                   int nx, double lx, double ux, CMatrix& u, int& new_nx, double& new_lx, double& new_ux, 

                   CMatrix& new_u)

{

	if ( !!m_pHelper  ){
		m_pHelper->pde_modify(it, t,ptr,nd,nx, lx, ux, u, new_nx, new_lx, new_ux,new_u,this);
		return;
	}

	int id, ix;

	new_nx = nx;
	new_lx = lx;
	new_ux = ux;

	for (id = 0; id < nd; id++)

		for (ix = 0; ix < nx; ix++)

			new_u[id][ix] = u[id][ix];


} 

double testc()
{
	
	
	

//	int n = 4;
	tridiagonal_matrix x;
/*	double below = -1.0;
	double on = 0.5;
	double above = 1.0;

	x.init(n);
	
	tridiagonal_matrix_populate(x,below, on, above);

	CVector vec(n);
	CVector res(n);

	vec[0] = 1;
	vec[1] = 2;
	vec[2] = 3;
	vec[3] = 4;
*/

	int n = 2;
	double below = 0.5;
	double on = 1.0;
	double above = -0.5;

	x.init(n);
	
	tridiagonal_matrix_populate(x,below, on, above);

	CVector vec(n);
	CVector res(n);

	vec[0] = 1;
	vec[1] = 1;
//	vec[2] = 3;
//	vec[3] = 4;



//	tridiagonal_matrix_multiply1(x,vec.getPtr(), res.getPtr());
CVector sol(n);
CVector gamma;
tridiagonal_matrix_solve1(x, vec, sol.getPtr(), gamma );

	double ty;

	int ii,jj;
	for ( ii = 0 ; ii < n; ii++ )
		for ( jj = 0 ; jj < n; jj++ )
		{
	
			ty = x(ii,jj);
		}




	 CVector result(1);
//	 double xtime = pde.pde_integrate(result);

//	 return xtime;
	 return 0.0;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void MlEqPdeDriver::reinitialize_nd(int newnd,MlEqPdeHelperHandle pHelper)
{

	if ( !! pHelper )
	{
		pde_init(newnd, m_different, m_nx, &m_max_nx, m_lx, 
				 m_ux, m_t, m_direction, m_keep,NULL, m_theta, 
				 m_continuous,m_convection_dependence, m_diffusion_dependence,
				 m_factor_dependence, m_constant_dependence,pHelper );
	}
	else
	{

		pde_init(newnd, m_different, m_nx, &m_max_nx, m_lx, 
				 m_ux, m_t, m_direction, m_keep,NULL, m_theta, 
				 m_continuous,m_convection_dependence, m_diffusion_dependence,
				 m_factor_dependence, m_constant_dependence,m_pHelper );

	}

}


void MlEqPdeDriver::pde_init(const int nd, const int different, const int nx, const int* max_nx, const double lx, 
						  const double ux, const CVector& t, const int direction, const int keep,
						  void* ptr, double theta, int continuous,int conv_dependence, 
						  int diff_dependence,int var_dependence, int rhs_dependence,MlEqPdeHelperHandle pHelperOpt )

{

	m_pHelper	  = pHelperOpt;


	int            i;



	/* Check the inputs are valid. */

	if (!(nd >= 1)){
		throw("no pde entered");
	}

	if(!((different == PdeEnums::PDE_MODEL_SAME)||(different == PdeEnums::PDE_MODEL_DIFFERENT))){
		throw("pde input error");
	}

	if ( !(nx >= 3)){
		throw("enter more spacial grid points");
	}

	if ( (max_nx == NULL)|| ((max_nx != NULL)&&(*max_nx < nx)) ){
		throw(  "max_nx must be increased" );
	}


	if(!(lx < ux)){
		throw("lower bound must be lower than upper bound");
	}


	if ( t.getsize() < 2 ){
		throw("increase number of time steps");
	}


	for (i = 0 ; i < t.getsize() - 1; i++)
	{
		if ( !(t[i] < t[i+1])){
			throw("time array must be entered in increasing order");
		}
	}


	if(!((direction == PdeEnums::PDE_DIRECTION_FORWARDS)||(direction == PdeEnums::PDE_DIRECTION_BACKWARDS))){
		throw("direction must be forward or backwards");
	}


	if( !((keep == PdeEnums::PDE_KEEP_LAST)||(keep == PdeEnums::PDE_KEEP_FULL))){
		throw("either keep pde or not keep pde; that is the question!");
	}


	if(!((continuous == PdeEnums::PDE_FACTORS_CONTINUOUS) || (continuous == PdeEnums::PDE_FACTORS_DISCONTINUOUS))){
		throw("pde factors must be either continuous or non-continuous");
	}


	/* Fill in the preliminary details. */

	m_nd            = nd;

	m_different     = different;

	m_nx            = nx;

	m_lx            = lx;

	m_ux            = ux;

	m_nt            = t.getsize();

	m_t             = t;

	m_direction     = direction;

	m_keep          = keep;

	m_continuous    = continuous;


	/* These values persist after the integration to allow users to use the utility routines. */

	m_curr_nx = m_nx;

	m_curr_lx = m_lx;

	m_curr_ux = m_ux;

	if (max_nx == NULL){

		m_max_nx = nx;
	}

	else
	{
		m_max_nx = *max_nx;
	}
	
	m_theta = theta;

	m_convection_dependence = conv_dependence;

	m_diffusion_dependence = diff_dependence;

	m_factor_dependence = var_dependence;

	m_constant_dependence = rhs_dependence;




	if (keep == PdeEnums::PDE_KEEP_FULL){

		m_u_full.resize(m_nt);
		for ( int i = 0 ; i < m_nt ;i++){
			m_u_full[i].resize(m_nd, m_max_nx);
		}
	}
	else
	{
		m_u_last.resize(m_nd, m_max_nx);
	}

	validate_pde();

	m_greekType	 = PdeEnums::PDE_COORDINATES_EXP ;

} 






void MlEqPdeDriver::validate_pde()

{

	int            i;

	if ((m_nd < 1) || (m_different != PdeEnums::PDE_MODEL_SAME) && (m_different != PdeEnums::PDE_MODEL_DIFFERENT))
	{
		throw("pde input error");
	}

	if ((m_nx < 3) || (m_nt < 2) || (m_ux <= m_lx) || (m_nt != m_t.getsize()) ){

		throw("pde input error");
	}

	if (m_nx > m_max_nx){
		throw("pde input error");
	}

	for (i = 0; i < m_nt - 1; i++)
	{
		if (m_t[i] >= m_t[i+1]){
			throw("pde input error");
		}
	}

	if ((m_theta <= 0.0)||(m_theta > 1.0)){

		throw("theta must be between zero and one");
	}


} 





/* Routine to allocate space for the 'state' of the integration and setup parameters. */


void MlEqPdeDriver::pde_create_state( )

{


	int    nx     = m_curr_nx;

	double lx     = m_curr_lx;

	double ux     = m_curr_ux;

	int    nd     = m_nd;

	int          size;

	int          i;



	/* The state contains the integration object. */


	m_dx          = (ux - lx) / (double)(nx - 1);

	m_dxdx        = m_dx * m_dx;

	m_one_m_theta = (1.0 - m_theta);



	/* Allocate space for the general work arrays. */

    if (m_convection_factors.getsize() != m_max_nx)
    {

		m_convection_factors.resize(m_max_nx);

		m_diffusion_factors.resize(m_max_nx);

		m_variable_factors.resize(m_max_nx);

		m_rhs_factors.resize(m_max_nx);

		m_u_interp.resize(nd, m_max_nx);

		m_residual.resize(m_max_nx);

		m_tmp.resize(m_max_nx);
	}


	/* The tridiagonal matrices contains the full operator and spatial operator only.  The 'known' vector contains   */

	/* the right hand side of the linear system solved at each timestep - it is therefore all the known information. */


	m_full_stencil.init(nx - 2);

	m_spatial_stencil.init(nx - 2);

	m_known.resize(m_max_nx - 2);



	/* Also need to store factors from the previous timestep to multiply the current solution when solving */

	/* for the next solution.  The factors are stored rather than the residual as the current solution may */

	/* have been altered by the user in the strategy routine.                                              */

	if (m_different == PdeEnums::PDE_MODEL_SAME){
		size = 1;
	}
	else
	{
		size = nd;
	}

    if (m_old_stencils.size() != size)
    {
	    m_old_stencils.resize(size);
    	m_old_fac_u.resize(size);
	    m_old_fac_l.resize(size);
    }

	for (i = 0; i < size; i++){
		m_old_stencils[i].init(nx - 2);
	}



	m_old_fac_u.resize(size);

	m_old_fac_l.resize(size);

	m_old_rhs.resize(size, m_max_nx);


	/* Allocate space for solution at the previous timestep.  This is only needed if storing only the last time */
	/* step.                                                                                                    */

	if (m_keep == PdeEnums::PDE_KEEP_LAST){
		m_u_prev.resize(nd, m_max_nx);
	}

} 








void MlEqPdeDriver::pde_resize_state()
{

	int        nx   = m_curr_nx;

	int              size;

	int              i;

	tridiagonal_matrix new_full_stencil;

	tridiagonal_matrix new_spatial_stencil;

	tridiagonal_matrix new_stencil;



	/* Allocate space for the new sized stencils. */

	new_full_stencil.init(nx - 2);

	new_spatial_stencil.init(nx - 2);


	m_full_stencil    = new_full_stencil;

	m_spatial_stencil = new_spatial_stencil;


	/* Need to change size of the previous stencils. */

	if (m_different == PdeEnums::PDE_MODEL_SAME){
		size = 1;
	}
	else{
		size = m_nd;
	}

		
	for (i = 0; i < size; i++)
	{
		new_stencil.init(nx - 2);
		m_old_stencils[i] = new_stencil;
	}

}



/*********************************************************************************************************************/










void MlEqPdeDriver::pde_spatial_factors( const int it, const double t, const int id)
{


	int            nx                 = m_curr_nx;

	double         lx                 = m_curr_lx;

	double         ux                 = m_curr_ux;

	double         dx                 = m_dx;

	double         dxdx               = m_dxdx;


	int                  ix;

	int                  want_nx;




	/* Do we need a strip of values from the routine or just a single value. */

	if (m_convection_dependence == PdeEnums::PDE_X_INDEPENDENT){
		want_nx = 0;
	}
	else
	{
		want_nx = nx;
	}


	/* Ask the user to provide constant or x-dependent convection factors in [0] or [1,nx-1] respectively. */

	convection_fcn(id, lx, ux, want_nx, it,t,m_users_data, m_convection_factors);


	/* Do we need a strip of values from the routine or just a single value. */

	if (m_diffusion_dependence == PdeEnums::PDE_X_INDEPENDENT){

		want_nx = 0;
	}
	else
	{
		want_nx = nx;
	}



	/* Ask the user to provide constant or x-dependent diffusion factors in [0] or [1,nx-1] respectively. */

	diffusion_fcn(id,lx, ux, want_nx, it,t,m_users_data,m_diffusion_factors);
	

	/* Do we need a strip of values from the routine or just a single value. */

	if (m_factor_dependence == PdeEnums::PDE_X_INDEPENDENT)
	{
		want_nx = 0;
	}
	else
	{
		want_nx = nx;
	}


	/* Ask the user to provide constant or x-dependent variable factors in [0] or [1,nx-1] respectively. */

	factor_fcn(id, lx, ux,want_nx,it,t,m_users_data ,m_variable_factors);


	/* Setup the trigiagonal matrix. */

	if ((m_convection_dependence == PdeEnums::PDE_X_INDEPENDENT)&&(m_diffusion_dependence == PdeEnums::PDE_X_INDEPENDENT)&&

	    (m_factor_dependence == PdeEnums::PDE_X_INDEPENDENT))

	{

		double tmp1, tmp2, tmp3;


		/* Since all factors are independent of x the tridiagonal matrix has constant below, on and above diagonal entries. */

		tmp1 = m_convection_factors[0] / 2.0 * dx;

		tmp3 = m_diffusion_factors[0];

		tmp2 = tmp3 * 2.0;

		tridiagonal_matrix_populate(m_spatial_stencil, - tmp1 + tmp3, - tmp2 + m_variable_factors[0] * dxdx, tmp1 + tmp3);

	}

	else

	{

		const double dx_div_2 = dx / 2.0;

		double       tmp1, tmp2, tmp3;



		/* Factors are different for each unknown. */

		tmp1 = 0.0;

		if (m_convection_dependence == PdeEnums::PDE_X_INDEPENDENT)

			tmp1 = m_convection_factors[0] * dx_div_2;

		tmp2 = 0.0;

		tmp3 = 0.0;

		if (m_diffusion_dependence == PdeEnums::PDE_X_INDEPENDENT)

		{

			tmp3 = m_diffusion_factors[0];

			tmp2 = tmp3 * 2.0;

		}



		/* Tridiagonal matrix will have different factors on each row. */

		for (ix = 1; ix < nx - 1; ix++)

		{

			double before = 0.0;

			double on     = 0.0;

			double after  = 0.0;



			if (m_convection_dependence == PdeEnums::PDE_X_DEPENDENT)

			{

				double tmp = m_convection_factors[ix] * dx_div_2;



				before -= tmp;

				after  += tmp;

			}

			else

			{

				before -= tmp1;

				after  += tmp1;

			}



			if (m_diffusion_dependence == PdeEnums::PDE_X_DEPENDENT)

			{

				double tmp = m_diffusion_factors[ix];



				before += tmp;

				on     -= tmp * 2.0;

				after  += tmp;

			}

			else

			{

				before += tmp3;

				on     -= tmp2;

				after  += tmp3;

			}



			if (m_factor_dependence == PdeEnums::PDE_X_DEPENDENT){

				on += m_variable_factors[ix] * dxdx;
			}

			else
			{
				on += m_variable_factors[0] * dxdx;
			}


			/* Set these factors in the tridiagonal matrix. */

			if (ix - 1 > 0){
				m_spatial_stencil(ix - 1,ix - 2) = before;
			}

			m_spatial_stencil(ix - 1,ix - 1) = on;

			if (ix - 1 < nx - 3){
				m_spatial_stencil(ix - 1,ix) = after;
			}

		}

	}



	/* Compute the factors of the boundary grid points required. */

	m_fac_u = 0.0;

	m_fac_l = 0.0;


	if (m_convection_dependence == PdeEnums::PDE_X_DEPENDENT)
	{
		m_fac_l -= m_convection_factors[1] / 2.0 * dx;

		m_fac_u += m_convection_factors[nx-2] / 2.0 * dx;
	}
	else
	{
		m_fac_l -= m_convection_factors[0] / 2.0 * dx;

		m_fac_u += m_convection_factors[0] / 2.0 * dx;
	}

	if (m_diffusion_dependence == PdeEnums::PDE_X_DEPENDENT)

	{
		m_fac_l += m_diffusion_factors[1];

		m_fac_u += m_diffusion_factors[nx-2];
	}
	else
	{
		m_fac_l += m_diffusion_factors[0];

		m_fac_u += m_diffusion_factors[0];
	}



	/* Evaluate the right hand side vector.*/

	if (m_constant_dependence == PdeEnums::PDE_X_INDEPENDENT){

		want_nx = 0;
	}
	else
	{
		want_nx = nx;
	}

	/* Ask the user to provide constant or x-dependent right-hand side factors in [0] or [1,nx-1] respectively. */

	constant_fcn(id, lx, ux, nx,it,t, m_users_data,m_rhs_factors);


}



/*********************************************************************************************************************/



/* Routine to evaluate the residual of applying the spatial discretisation operators to the initial solution */

/* provided by the user.  This will then be used in the first step of the time stepping routine.             */


void MlEqPdeDriver::pde_initial_factors_prev( int it,  double t)
{


	int            nd   = m_nd;

	int            nx   = m_curr_nx;

	int            id;


	/* If all equations are the same model then we can use the same stencil. */

	if (m_different == PdeEnums::PDE_MODEL_SAME)

	{

		/* Get the stencil at time t=t0. */

		pde_spatial_factors( it, t, -1);

		/* Copy these into the previous factors. */

		tridiagonal_matrix_copy(m_old_stencils[0], m_spatial_stencil);

		m_old_fac_l[0] = m_fac_l;

		m_old_fac_u[0] = m_fac_u;

		m_old_rhs[0] =  m_rhs_factors;

	}

	else

	{

		for (id = 0; id < nd; id++)
		{

			pde_spatial_factors( it, t, id);


			/* Copy these into the previous factors. */

			tridiagonal_matrix_copy(m_old_stencils[id], m_spatial_stencil);

			m_old_fac_l[id] = m_fac_l;

			m_old_fac_u[id] = m_fac_u;

			m_old_rhs[id] =  m_rhs_factors;

		}

	}

}





void MlEqPdeDriver::pde_theta_timestep(const int it, const double t, const double dt, CMatrix&  u1, 
									CMatrix&  u2)
{


	int            nd                 = m_nd;


	double               fac_l              = 0.0;

	double               fac_u              = 0.0;

    int                  id;


	/* If dependent variables all have the same governing equation then only evaluate the stencil once. */

	if (m_different == PdeEnums::PDE_MODEL_SAME)

	{

		/* Evaluate the spatial operators stencil at 't'. */

		pde_spatial_factors(it, t, -1);

		fac_l = m_fac_l;
		fac_u = m_fac_u;

	}



	/* Take a timestep for each dependent variable but in reverse, nd-1, nd-2, .. ,0 */

	for (id = nd - 1; id >= 0; id--)

	{

		int          nx          = m_curr_nx;

		int          max_nx      = m_max_nx;

		double       dx          = m_dx;

		double       dxdx        = m_dxdx;

		double       dxdx_div_dt = dxdx / dt;


		pde_boundary_condition lower;

		pde_boundary_condition upper;

		double             uterm_1, uterm_2, rterm_1;

		double             uterm_nm2, uterm_nm3, rterm_nm2;

		int                ix;

		int                id2;



		/* Calculate spatial factors of unknown variables. Different dependent variables have different equations. */

		if (m_different == PdeEnums::PDE_MODEL_DIFFERENT) 
		{

			/* Evaluate the spatial operators stencil at 't'. */

			pde_spatial_factors( it, t, id);

			fac_l = m_fac_l;
			fac_u = m_fac_u;

		}



		/* Variables 'id2' is contains previous factor index (0 if same model, id if different. */

		if (m_different == PdeEnums::PDE_MODEL_DIFFERENT)
		{
			id2 = id;
		}
		else
		{
			id2 = 0;
		}


		/* If its continuous then re-use the old factors. */

		if (m_continuous == PdeEnums::PDE_FACTORS_CONTINUOUS)

		{

			/* Calculate the contribution from the previous timestep. */

			tridiagonal_matrix_multiply1(m_old_stencils[id2], &(u2[id][1]), &(m_residual[1]));



			/* Add in the right hand side values these do not depend on the . */

			if (m_constant_dependence == PdeEnums::PDE_X_INDEPENDENT)

			{

				double tmp = m_old_rhs[id2][0] * dxdx;

				for (ix = 1; ix < nx - 1; ix++){

					m_residual[ix] = tmp - m_residual[ix];
				}
			}

			else 
			{
				for (ix = 1; ix < nx - 1; ix++){
					m_residual[ix] = m_old_rhs[id2][ix] * dxdx - m_residual[ix];
				}
			}


			/* Add in the boundary elements. */

			m_residual[1] -= m_old_fac_l[id2] * u2[id][0];

			m_residual[nx - 2] -= m_old_fac_u[id2] * u2[id][nx - 1];

		}

		else

		{	

			/* Calculate the contribution from the previous timestep but use the latest factors. */

			tridiagonal_matrix_multiply1(m_spatial_stencil, &(u2[id][1]), &m_residual[1]);



			/* Add in the right hand side values these do not depend on the . */

			if (m_constant_dependence == PdeEnums::PDE_X_INDEPENDENT)

			{
				double tmp = m_rhs_factors[0] * dxdx;

				for (ix = 1; ix < nx - 1; ix++){
					m_residual[ix] = tmp - m_residual[ix];
				}
			}

			else 
			{

				for (ix = 1; ix < nx - 1; ix++)
				{
					m_residual[ix] = m_rhs_factors[ix] * dxdx - m_residual[ix];
				}
			}


			/* Add in the boundary elements. */

			m_residual[1] -= fac_l * u2[id][0];

			m_residual[nx - 2] -= fac_u * u2[id][nx - 1];

		}

		     

		/* Evaluate the boundary conditions for the next time step. */

		lower.coord = PdeEnums::PDE_COORDINATES_PDE;

		lower.scale = 1.0;

		lower.shift = 0.0;

		lower.a     = 0.0;

		lower.b     = 0.0;

		lower.c     = 0.0;

		lower.d     = 0.0;



		upper.coord = PdeEnums::PDE_COORDINATES_PDE;

		upper.scale = 1.0;

		upper.shift = 0.0;

		upper.a     = 0.0;

		upper.b     = 0.0;

		upper.c     = 0.0;

		upper.d     = 0.0;


		boundary_condition(id,it,t,m_users_data,lower,upper);


		/* If the transformation has changed then we'll need to recalculate the transformed mesh spacing. */

		if (lower.coord == PdeEnums::PDE_COORDINATES_EXP && ((!m_expdx_set_lower) || 

			(m_last_lower.scale != lower.scale) || (m_last_lower.shift != lower.shift)))

		{

			const double lx = m_lx;

			double       y0, y1, y2;



			/* Get the mesh spacing for the transformed exp coordinate system. */

			y0 = exp(lx * lower.scale + lower.shift);

			y1 = exp((lx + dx) * lower.scale + lower.shift);

			y2 = exp((lx + 2.0 * dx) * lower.scale + lower.shift);

			m_lo_expdx0 = y1 - y0;

			m_lo_expdx1 = y2 - y1;



			/* Keep this for next time. */

			m_last_lower = lower;

			m_expdx_set_lower = TRUE;

		}



		if (upper.coord == PdeEnums::PDE_COORDINATES_EXP && ((!m_expdx_set_upper) || 

			(m_last_upper.scale != upper.scale) || (m_last_upper.shift != upper.shift)))

		{

			const double ux = m_ux;

			double       y0, y1, y2;


			/* Get the mesh spacing for the transformed exp coordinate system. */

			y2 = exp(ux * upper.scale + upper.shift);

			y1 = exp((ux - dx) * upper.scale + upper.shift);

			y0 = exp((ux - 2.0 * dx) * upper.scale + upper.shift);

			m_hi_expdx0 = y1 - y0;

			m_hi_expdx1 = y2 - y1;


			/* Keep this for next time. */

			m_last_upper = upper;

			m_expdx_set_upper = TRUE;

		}



		/* Copy the spatial stencil (independent of boundary conditions) into the full stencil. */		    

		tridiagonal_matrix_copy(m_full_stencil, m_spatial_stencil);



		if (m_continuous == PdeEnums::PDE_FACTORS_CONTINUOUS)

		{

			/* Copy the factors obtained from spatial discretisation over to a different tridiagonal matrix   */

			/* as we will use these spatial factors on the next time step.  If using the same model then only */

			/* do it on the last timestep.                                                                    */

			if ((m_different == PdeEnums::PDE_MODEL_SAME)&&(id == 0)) 

			{

				tridiagonal_matrix_copy(m_old_stencils[0], m_spatial_stencil);

				m_old_fac_l[0] = fac_l;

				m_old_fac_u[0] = fac_u;

				m_old_rhs[0] =  m_rhs_factors;

			}

			if (m_different == PdeEnums::PDE_MODEL_DIFFERENT)

			{

				tridiagonal_matrix_copy(m_old_stencils[id], m_spatial_stencil);

				m_old_fac_l[id] = fac_l;

				m_old_fac_u[id] = fac_u;

				m_old_rhs[id] =  m_rhs_factors;

			}  

		}



		/* Multiply through by theta. */

		tridiagonal_matrix_scale(m_full_stencil, m_theta);



		/* Construct the right hand side of linear system that contains known variables and values only. */

		/* The array residual contains the result of applying the discrete operator to the solution at   */

		/* the end of the previous timestep.                                                             */

		for (ix = 1; ix < nx - 1; ix++){

			m_known[ix - 1] = m_one_m_theta * m_residual[ix];
		}


		/* Introduce the discrete time derivative operators into both unknown and known sides.      */

		for (ix = 1; ix < nx - 1; ix++)
		{

			if (m_direction == PdeEnums::PDE_DIRECTION_BACKWARDS)

			{

				/* Note that we are reversing backwards but dt is positive.  Therefore the negative    */

				/* sign is absorbed into the time derivative operator: dudt = (U(j+1,i) - U(j,i)) / dt */

				m_full_stencil(ix - 1,ix - 1) -= dxdx_div_dt;

				m_known[ix - 1] -= u2[id][ix] * dxdx_div_dt;

			}

			else

			{

				/* Forwards so same as usual. */

				m_full_stencil(ix - 1,ix - 1) += dxdx_div_dt;

				m_known[ix - 1] += u2[id][ix] * dxdx_div_dt;

			}
		}


		/* The right-hand-side factors are independent of the dependent variable and therefore appear in */

		/* the known vector only.                                                                        */

		if (m_constant_dependence == PdeEnums::PDE_X_INDEPENDENT)

		{

			double tmp = m_rhs_factors[0] * dxdx * m_theta ;

			for (ix = 1; ix < nx - 1; ix++){
				m_known[ix - 1] += tmp;
			}

		}

		else 
		{

			double dxdx_theta = dxdx * m_theta;

			for (ix = 1; ix < nx - 1; ix++){
				m_known[ix - 1] += m_rhs_factors[ix] * dxdx_theta;
			}
		}



		/* Apply the boundary conditions.  To understand this section you *have* to refer to the document */

		/* discussing this: bc1.doc.                                                                      */

		if (lower.coord == PdeEnums::PDE_COORDINATES_PDE)

		{

			const double tmp1 = lower.a / dxdx;

			const double tmp2 = lower.b / dx;

			double       lo_denom;



			lo_denom = tmp1 - tmp2 + lower.c;

			if (fabs(lo_denom) < epsilon){
				throw(  "zero boundary entered");
			}



			uterm_1 = (-2.0 * tmp1 + tmp2) / lo_denom;

			uterm_2 = tmp1 / lo_denom;

			rterm_1 = lower.d / lo_denom;

		}

		else

		{

			const double a          = lower.a;

			const double b          = lower.b;

			const double c          = lower.c;

			const double dx0        = m_lo_expdx0;

			const double dx1        = m_lo_expdx1;

			const double sum_square = dx0 * dx1 * (dx0 + dx1) / 2.0;

			double       lo_denom;



			lo_denom = a * dx1 / sum_square - b / dx0 + c;

			if (fabs(lo_denom) < epsilon){

				throw(  "zero boundary entered");
			}


			uterm_1 = (- a * (dx1 + dx0) / sum_square + b / dx0) / lo_denom;

			uterm_2 = a * dx0 / (sum_square * lo_denom);

			rterm_1 = lower.d / lo_denom;

		}

		m_known[0] -= m_theta * fac_l * rterm_1;

		m_full_stencil(0,0) -= m_theta * fac_l * uterm_1;

		m_full_stencil(0,1) -= m_theta * fac_l * uterm_2;



		if (upper.coord == PdeEnums::PDE_COORDINATES_PDE)

		{

			double tmp1 = upper.a / dxdx;

			double tmp2 = upper.b / dx;

			double       hi_denom;


			hi_denom = tmp1 + tmp2 + upper.c;

			if (fabs(hi_denom) < epsilon){

				throw(  "zero boundary entered");
			}


			uterm_nm2 = (-2.0 * tmp1 - tmp2) / hi_denom;

			uterm_nm3 = tmp1 / hi_denom;

			rterm_nm2 = upper.d / hi_denom;

		}

		else

		{

			const double a          = upper.a;

			const double b          = upper.b;

			const double c          = upper.c;

			const double dx0        = m_hi_expdx0;

			const double dx1        = m_hi_expdx1;

			const double sum_square = dx0 * dx1 * (dx0 + dx1) / 2.0;

			double       hi_denom;



			hi_denom = a * dx0 / sum_square + b / dx1 + c;

			if (fabs(hi_denom) < epsilon){

				throw(  "zero boundary entered");
			}


			uterm_nm2 = (- a * (dx1 + dx0) / sum_square - b / dx1) / hi_denom;

			uterm_nm3 = a * dx1 / (sum_square * hi_denom);

			rterm_nm2 = upper.d / hi_denom;

		}

		m_known[nx - 3] -= m_theta * fac_u * rterm_nm2;

		m_full_stencil(nx - 3,nx - 4) -= m_theta * fac_u * uterm_nm3;

		m_full_stencil(nx - 3,nx - 3) -= m_theta * fac_u * uterm_nm2;



		/* Now have a complete set of linear equations, so perform the solve. */

		tridiagonal_matrix_solve1(m_full_stencil, m_known, &(u1[id][1]), m_tmp);


		/* Fill in the boundary elements using the precomputed factors from above.                                      */

		u1[id][0] = - uterm_2 * u1[id][2] - uterm_1 * u1[id][1] + rterm_1;

		u1[id][nx - 1] = - uterm_nm3 * u1[id][nx - 3] - uterm_nm2 * u1[id][nx - 2] + rterm_nm2;

	}

} 



double MlEqPdeDriver::pde_integrate(CVector& result)
{

	int        nt     = m_nt;
	double res = pde_integrate(result,nt);
	return res;
}


double MlEqPdeDriver::pde_integrate(CVector& result,int nt)
{

	time_t start,finish;
	time(&start);


	int        max_nx = m_max_nx;


	int        nd     = m_nd;



	PdeEnums::PdeMode mode;

	CMatrix* initial_ptr;

//	double           **initial_ptr;

	int              it2;

	int              id;

	double           dt;


	/* Create all workspaces needed for the integration in a 'state' structure. */

	pde_create_state( );



	/* Variable 'initial_ptr' will point to where the initial solution is to be stored. */

	if ((m_direction == PdeEnums::PDE_DIRECTION_FORWARDS)&&(m_keep == PdeEnums::PDE_KEEP_FULL)){

		initial_ptr = &(m_u_full[0]);
	}
	else if ((m_direction == PdeEnums::PDE_DIRECTION_BACKWARDS)&&(m_keep == PdeEnums::PDE_KEEP_FULL)){

		initial_ptr = &(m_u_full[nt - 1]);
	}
	else
	{
		initial_ptr = &m_u_last;
	}


	/* Set up the initial conditions */



	initial_condition( m_users_data, m_nd, m_curr_nx, *initial_ptr);


	/* Evaluate the spatial factors at initial time. */

	if (m_direction == PdeEnums::PDE_DIRECTION_FORWARDS){
	   pde_initial_factors_prev(0, m_t[0]);
	}
	else
	{
	    pde_initial_factors_prev(nt-1, m_t[nt-1]);
	}


	/* First time through the transformed mesh spacing will have to be calculated. */

	m_expdx_set_upper = FALSE;

	m_expdx_set_lower = FALSE;



	/* Loop variable 'it2' is independent of integration direction, variable 'it' moves in the correct direction. */

	mode = PdeEnums::PDE_NORMAL_MODE;
	m_last_mode = mode;

	for (it2 = 1; it2 < nt; it2++)

	{

		int      it;

		CMatrix* curr_ptr;

		CMatrix* prev_ptr;



		/* Variable 'it' is the real loop variable. */

		if (m_direction == PdeEnums::PDE_DIRECTION_FORWARDS)

		{

			it = it2;

			dt = m_t[it] - m_t[it - 1];

		}

		else 

		{

			it = nt - it2 - 1;

			dt = m_t[it + 1] - m_t[it];

		}



		/* If we're only keeping the last timestep solution then transfer this to a temporary array before we */

		/* overwrite with the solution on this timestep.                                                      */

		if ((m_keep == PdeEnums::PDE_KEEP_LAST)&&(mode == PdeEnums::PDE_NORMAL_MODE))

		{
			for (id = 0; id < nd; id++){
				m_u_prev[id] =  m_u_last[id];
			}
		}



		/* Variables 'curr_ptr' and 'prev_ptr' will point to arrays containing solution required for this step and */

		/* solution obtained on the previous timestep.                                                             */

		if (m_keep == PdeEnums::PDE_KEEP_FULL)

		{

			curr_ptr = &(m_u_full[it]);				

			if (mode == PdeEnums::PDE_NORMAL_MODE)
			{

				if (m_direction == PdeEnums::PDE_DIRECTION_FORWARDS){

					prev_ptr = &(m_u_full[it - 1]);	
				}
				else
				{
					prev_ptr = &(m_u_full[it + 1]);	
				}

			}

			else
			{
				prev_ptr = &m_u_interp;
			}

		}

		else

		{

			curr_ptr = &m_u_last;

			if (mode == PdeEnums::PDE_NORMAL_MODE)
			{
				prev_ptr = &m_u_prev;
			}
			else
			{
				prev_ptr = &m_u_interp;
			}

		}



		/* Perform the next timestep. */

		pde_theta_timestep(it, m_t[it],dt, *curr_ptr, *prev_ptr);


		/* Run the user's strategy. */

		mode = strategy_function(it, m_t[it],m_users_data,nd,m_curr_nx,*curr_ptr); 
		m_last_mode = mode;

		if (mode == PdeEnums::PDE_NORMAL_MODE)

		{
			mode = PdeEnums::PDE_NORMAL_MODE;
		}

		else if (mode == PdeEnums::PDE_MODIFIED_MODE)
		{

			int    new_nx;

			double new_lx, new_ux;



			/* The user wishes to reset the domain of integration and [possibly] the number of mesh points. */

			 pde_modify(it, m_t[it],m_users_data, nd,m_curr_nx,
					   m_curr_lx, m_curr_ux, *curr_ptr, new_nx,
					   new_lx, new_ux,m_u_interp);



			if (new_nx > max_nx){
				throw( "° incorrect size encountered of spacial grid");
			}



			/* If the mesh size has changed then update the factors array for the current solution. */

			if (new_nx != m_curr_nx)

			{
				m_curr_nx = new_nx;
				pde_resize_state();

			}



			/* Reset the values in the integration task. */

			m_curr_lx  = new_lx;

			m_curr_ux  = new_ux;

			m_dx       = (m_curr_ux - m_curr_lx) / (double)(m_curr_nx - 1);

			m_dxdx     = m_dx * m_dx;



			/* The mesh spacing factors are now invalid independently of whether the actual boundary conditions */

			/* change.                                                                                          */

			m_expdx_set_upper = FALSE;

			m_expdx_set_lower = FALSE;



			/* The previous spatial factors we're holding are no long relevant.  Re-evaluate these. */

			pde_initial_factors_prev( it, m_t[it]);

			mode = PdeEnums::PDE_MODIFIED_MODE;

		}
		else
		{

			/* Didn't recognise the response or it was stop. */
			throw("pde error");
		} 

	}

////////////////////////////////////////////////////////////////////


//	get price:
	CVector			 x_array;

	x_array.resize(m_nx);

	for (int i = 0 ; i < m_nx; i++ )
	{
		x_array[i] = pde_x(i);
	}


//  put in some linear interpolation here


	CVector vals(m_nx);


	result.resize(m_nd);

    // if mode == PDE_MODIFIED_MODE then the last step has been modified but
    // m_u_interp has never been affected in the solution
    
    const CMatrix& pdeValues = (mode == PdeEnums::PDE_MODIFIED_MODE) ? m_u_interp :
                                (m_keep == PdeEnums::PDE_KEEP_FULL)? m_u_full[0] : m_u_last;

	for (int i = 0 ; i < m_nd; i++ ){
		result[i] =  MlEqMaths::linearInterp(x_array,pdeValues[i],0.0);
	}

    



	time(&finish);

	double xtime = difftime(finish,start);

	
	return xtime;
	return 1.0;
}





double MlEqPdeDriver::pde_x( double i)
{


	double   lx = m_curr_lx;

	double   ux = m_curr_ux;

	double   nx = m_curr_nx;

	double  z = i / (nx - 1) * (ux - lx) + lx;

	return z;

} 



/****************************************************************
**	Class  : pde_driver 
**	Routine: pde_last_solution
**	Returns: 
**	Action : 
**           
****************************************************************/


CMatrix&  MlEqPdeDriver::pde_last_solution()

{
	validate_pde();
	if ( m_keep != PdeEnums::PDE_KEEP_LAST){
		throw(  "keep variable must be PDE_KEEP_LAST in this case ");
	}
	return m_u_last;
} 



/****************************************************************
**	Class  : pde_driver 
**	Routine: pde_last_solution
**	Returns: 
**	Action : 
**           
****************************************************************/


void MlEqPdeHelper::greeks(double& delta,double& gamma,int id,int it,double t, void* ptr, MlEqPdeDriver* pde)
{

	if ( it && pde->m_keep != PdeEnums::PDE_KEEP_FULL){

		throw(" you must keep the full grid when requestion delta at date different from valuation date");
	}

    const CMatrix& pdeValues = (pde->m_last_mode == PdeEnums::PDE_MODIFIED_MODE) ? pde->m_u_interp :
                                (pde->m_keep == PdeEnums::PDE_KEEP_FULL)? pde->m_u_full[it] : pde->m_u_last;



	CVector x_array(pde->m_nx);
	for (int i = 0 ; i < pde->m_nx; i++ ){
		x_array[i] = pde->pde_x(i);
	}


	double val,valUp,valDo,x;

	if ( pde->m_greekType == PdeEnums::PDE_COORDINATES_EXP ){
		x = 0.0;
	}
	else{
		x = m_spot;
	}

	double dx = pde->m_dx;

	val		=  MlEqMaths::linearInterp(x_array,pdeValues[id],x);
	valUp	=  MlEqMaths::linearInterp(x_array,pdeValues[id],x+dx);
	valDo	=  MlEqMaths::linearInterp(x_array,pdeValues[id],x-dx);

	delta	= (valUp-valDo)/(2.0*dx);

	gamma	= (valUp+valDo-2.0*val)/(dx*dx);


	if ( pde->m_greekType == PdeEnums::PDE_COORDINATES_EXP )
	{
		gamma = ( gamma - delta ) / (m_spot*m_spot) ;
		delta /= m_spot;
	}


}


/****************************************************************
**	Class  : pde_driver 
**	Routine: pde_full_solution
**	Returns: 
**	Action : 
**           
****************************************************************/


GVector< CMatrix>& MlEqPdeDriver::pde_full_solution()

{

	validate_pde();

	if ( m_keep != PdeEnums::PDE_KEEP_FULL)
		throw(  "keep variable must be PDE_KEEP_FULL ");


	return m_u_full;
} 


/****************************************************************
**	Class  : pde_driver 
**	Routine: boundary_condition
**	Returns: 
**	Action : 
**           
****************************************************************/


void MlEqPdeDriver::boundary_condition(int id,int it,double t, void* ptr,pde_boundary_condition& lower,pde_boundary_condition& upper)
{

	if ( !!m_pHelper ){
		m_pHelper->boundary_condition(id,it,t,ptr,lower,upper,this);
		return;
	}


	lower.coord = PdeEnums::PDE_COORDINATES_EXP;
	lower.scale = 1.0;
	lower.shift = 0.0;
	lower.a = 1.0;
	lower.b = 0.0;
	lower.c = 0.0;
	lower.d = 1.0;


	upper.coord = PdeEnums::PDE_COORDINATES_EXP;
	upper.scale = 1.0;
	upper.shift = 0.0;
	upper.a = 1.0;
	upper.b = 1.0;
	upper.c = 0.0;
	upper.d = 0.0;
	
}


/****************************************************************
**	Class  : pde_driver 
**	Routine: initial_condition
**	Returns: 
**	Action : 
**           
****************************************************************/


void MlEqPdeDriver::initial_condition(void* ptr,int nd,int nx,CMatrix& u )
{

	if ( !!m_pHelper ){
		m_pHelper->initial_condition(ptr,nd,nx,u,this );
		return;
	}

}



/****************************************************************
**	Class  : pdeBlackScholes 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeBlackScholes::convection_fcn(int id, double lx, double ux, int nx, int it, double tt, void* ptr, CVector& f)

{
	f[0] = m_loanspread - 0.5*m_vol*m_vol;
} 




/****************************************************************
**	Class  : pdeBlackScholes 
**	Routine: diffusion_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/

void pdeBlackScholes::diffusion_fcn(int id, double lx, double ux, int nx, int it, double tt, 

        void* ptr, CVector& f)
{
	f[0] = 0.5 * m_vol*m_vol;
} 


/****************************************************************
**	Class  : pdeBlackScholes 
**	Routine: constant_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeBlackScholes::constant_fcn(int id, double lx, double ux, int nx, int it, double t, 

                    void* ptr, CVector& f)

{
	f[0] = 0.0;
} 



/****************************************************************
**	Class  : pdeBlackScholes 
**	Routine: factor_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeBlackScholes::factor_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f)
{
	f[0] = - m_discount_rate;
} 



/****************************************************************
**	Class  : pdeBlackScholes 
**	Routine: init
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeBlackScholes::init(double vol,double forward,double discountRate,double spot,
						   MlEqDateHandle asOfDate,long maturityDate,
						   int numberSpacialGridPoints,int numberTimeSteps,double numberStdev,
						   MlEqPdeHelper* pHelper,double minSpotOpt=-1.0,double maxSpotOpt=1e10)
{

//testc();

	m_vol				=	vol;	
	m_discount_rate		=	discountRate;
	m_spot				=	spot;

	double mat			= asOfDate->GetYearFraction(maturityDate);

	if ( mat < 1e-4 ){
		throw("incorrect maturity date entered");
	}

	double drift	= log(forward/spot)/mat;
	m_loanspread	= drift;	

	double theta = 0.5;
	void* ptr = NULL; 

	m_convection_dependence		=	PdeEnums::PDE_X_INDEPENDENT; 
	m_diffusion_dependence		=	PdeEnums::PDE_X_INDEPENDENT; 
	m_factor_dependence			=	PdeEnums::PDE_X_INDEPENDENT;
	m_constant_dependence		=	PdeEnums::PDE_X_INDEPENDENT;

	m_thetaDiffusion = PdeEnums::PDE_STEP_DIFFUSION;

	// set up timegrid

	double dt  = mat/(double)(numberTimeSteps-1);

	int i = 0;
	CVector t(numberTimeSteps);
	t[i] = 0.0;
	
	for ( i = 1 ; i < numberTimeSteps; i++ )
		t[i] = t[i-1]+dt;						


	double lx,ux;

	lx = MlEqMaths::Max(m_spot*exp(-m_vol*numberStdev*sqrt(mat)),minSpotOpt);
	ux = MlEqMaths::Min(m_spot*exp( m_vol*numberStdev*sqrt(mat)),maxSpotOpt);

	lx = log(lx/spot);
	ux = log(ux/spot);

	// reset number gridpoints to hit zero exactly

	double hitzero   = - lx*(double)(numberSpacialGridPoints - 1) / (ux - lx);
	if ( floor(hitzero + 1e-22) < hitzero )
		numberSpacialGridPoints++;


	int max_nx= numberSpacialGridPoints;

	const int nd = 1;

	int continuous = 1;

	pde_init(nd, PdeEnums::PDE_MODEL_SAME, numberSpacialGridPoints, &max_nx, lx, 
								ux, t, PdeEnums::PDE_DIRECTION_BACKWARDS, PdeEnums::PDE_KEEP_LAST,
								ptr, 
								theta, 
								continuous,
								m_convection_dependence, 
								m_diffusion_dependence, 
								m_factor_dependence, 
								m_constant_dependence,pHelper);


}



/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: init
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVol::init( DupireLocalVolHandle lv,MlEqPdeHelperHandle pHelperOpt)
{
	m_pLv			= lv;
	
	m_loanspread	= lv->getLoanSpreads();
	m_quantoDrift	= lv->getQuantoDrift();
	m_rate			= lv->getRates();
	m_Grid			= lv->getLocalVolGrid();

	double theta = 0.5;
	void* ptr = NULL; 

	m_convection_dependence		=	PdeEnums::PDE_X_DEPENDENT; 
	m_diffusion_dependence		=	PdeEnums::PDE_X_DEPENDENT; 
	m_factor_dependence			=	PdeEnums::PDE_X_INDEPENDENT;
	m_constant_dependence		=	PdeEnums::PDE_X_INDEPENDENT;

	m_thetaDiffusion = PdeEnums::PDE_THETA_DIFFUSION;

	// set up timegrid


	CVector t = lv->getTimes();

	int max_nx= m_Grid[0].getsize();

	const CVector& spots = lv->getSpotGrid();
	double spot = lv->getSpot();

	double lx = log(spots[0]/spot);
	double ux = log(spots[spots.getsize()-1]/spot);

	int numberSpacialGridPoints = max_nx;

	const int nd = 1;

	int continuous = 1;

	pde_init(nd, PdeEnums::PDE_MODEL_SAME, numberSpacialGridPoints, &max_nx, lx, 
								ux, t, PdeEnums::PDE_DIRECTION_BACKWARDS, PdeEnums::PDE_KEEP_FULL,
								ptr, 
								theta, 
								continuous,
								m_convection_dependence, 
								m_diffusion_dependence, 
								m_factor_dependence, 
								m_constant_dependence, pHelperOpt);


}


/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVol::convection_fcn(int id, double lx, double ux, int nx, int it, double tt, void* ptr, CVector& f)
{
	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	for ( int ix = 0 ; ix < nx; ix++ )
	{
		double local_vol = m_pLv->getlocalVol(it,ix);

		f[ix] = m_loanspread[itm] + local_vol * ( m_quantoDrift[itm] - 0.5 * local_vol );
	}
} 

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVol::diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f)
{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	double xvol;
	for ( int ix = 0 ; ix < nx; ix++ )
	{
		xvol = m_pLv->getlocalVol(it,ix);
		f[ix] = 0.5 * pow(xvol,2.0);
	}
}

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: factor_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVol::factor_fcn(int id, double lx, double ux, int nx, int it, double t, 

                    void* ptr, CVector& f)

{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 


	f[0] = - m_rate[itm];
} 





double pdeLocalVol::priceAndGreeks( double& delta, double& gamma, double& vega )
{
	CVector tmp(1);
	pde_integrate(tmp);
	double price = tmp[0];

	m_pHelper->greeks( delta, gamma, 0, 0, 0, 0, this );

	m_pLv->parallelVegaBump();
	init( m_pLv, m_pHelper );	
	pde_integrate(tmp);
	m_pLv->parallelVegaBump();

	vega = tmp[0] - price ;

	return price;
}






void pdeLocalVolBump::init( DupireLocalVolHandle lv,MlEqPdeHelperHandle pHelperOpt, double ls_bump)
{

	m_pLv	= lv;

	const CVector& rate		= m_pLv->getRates();
	const CVector& spread	= m_pLv->getLoanSpreads();


	m_loanspread	= spread;
	m_rate			= rate;
	m_Grid			= lv->getLocalVolGrid();

	m_pLv->parallelVegaBump();
	m_bumpedGrid	= m_pLv->getLocalVolGrid();
	m_pLv->parallelVegaBump();

	m_pLv->bumpLocalDrift( ls_bump );
	m_bumpedLoanspread = m_pLv->getLoanSpreads();
	m_pLv->bumpLocalDrift( -ls_bump );


	double theta = 0.5;
	void* ptr = NULL; 

	m_convection_dependence		=	PdeEnums::PDE_X_DEPENDENT; 
	m_diffusion_dependence		=	PdeEnums::PDE_X_DEPENDENT; 
	m_factor_dependence			=	PdeEnums::PDE_X_INDEPENDENT;
	m_constant_dependence		=	PdeEnums::PDE_X_INDEPENDENT;

	m_thetaDiffusion = PdeEnums::PDE_THETA_DIFFUSION;

	// set up timegrid


	CVector t = m_pLv->getTimes();


	const CVector& spots = m_pLv->getSpotGrid();

	int max_nx = spots.getsize();

	double spot = m_pLv->getSpot();

	double lx = log(spots[0]/spot);
	double ux = log(spots[spots.getsize()-1]/spot);

	int numberSpacialGridPoints = max_nx;

	const int nd = 3;	// you can add bumped rates...

	int continuous = 1;

	pde_init(nd, PdeEnums::PDE_MODEL_DIFFERENT, numberSpacialGridPoints, &max_nx, lx, 
								ux, t, PdeEnums::PDE_DIRECTION_BACKWARDS, PdeEnums::PDE_KEEP_FULL,
								ptr, 
								theta, 
								continuous,
								m_convection_dependence, 
								m_diffusion_dependence, 
								m_factor_dependence, 
								m_constant_dependence, pHelperOpt);


}


/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolBump::convection_fcn(int id, double lx, double ux, int nx, int it, double tt, void* ptr, CVector& f)
{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	const CVector* pvol;		
	double drift;

	switch( id )
	{
	case 0:	// regular price
		drift	= m_loanspread[itm];
		pvol	= &m_Grid[itm];
		break;
	case 1: // vega bump
		drift	= m_loanspread[itm];
		pvol	= &m_bumpedGrid[itm];
		break;
	case 2: // loan spread bump
		drift	= m_bumpedLoanspread[itm];
		pvol	= &m_Grid[itm];
		break;
	}

	for ( int ix = 0 ; ix < nx; ix++ ){
		f[ix] = drift - 0.5 * pow((*pvol)[ix],2.0);
	}
} 

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolBump::diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f)
{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	const CVector* pvol;
	double drift;

	switch( id )
	{
	case 0:	// regular price
		drift	= m_loanspread[itm];
		pvol	= &m_Grid[itm];
		break;
	case 1: // vega bump
		drift	= m_loanspread[itm];
		pvol	= &m_bumpedGrid[itm];
		break;
	case 2: // loan spread bump
		drift	= m_bumpedLoanspread[itm];
		pvol	= &m_Grid[itm];
		break;
	}

	for ( int ix = 0 ; ix < nx; ix++ ){
		f[ix] = 0.5 * pow((*pvol)[ix],2.);
	}
}

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: factor_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolBump::factor_fcn(int id, double lx, double ux, int nx, int it, double t, void* ptr, CVector& f)
{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	double rate = m_rate[itm];
	if( id == 4 )
		rate = m_bumpedRate[itm];

	f[0] = - rate;
} 

/****************************************************************
**	Class  : pdeLocalVolEffective 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


double pdeLocalVolBump::priceAndGreeks( double& delta, double& gamma, double& vega )
{
	CVector tmp(2);
	pde_integrate(tmp);
	double price = tmp[0];

	m_pHelper->greeks( delta, gamma, 0, 0, 0, 0, this );

	vega = tmp[0] - tmp[1] ;

	return price;
}

/****************************************************************
**	Class  : pdeLocalVolEffective 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/

double pdeLocalVolEffective::getlocalVol(int idate,int ixspace)
{
	return m_Grid[m_map[idate]][ixspace];
}


/****************************************************************
**	Class  : pdeLocalVolEffective 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::setlocalVolSlice(int idateCalib)
{

	for (int ix = 0; ix < m_nx; ix++){
		m_Grid[idateCalib][ix] = m_localVolSlice[idateCalib]->getValue(pde_x(ix));
	}

}

/****************************************************************
**	Class  : pdeLocalVolEffective 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::reinitialize(int idateCalib,MlEqAnalyticCurveWithTanhWingEdgeHandle& lv)
{

	m_localVolSlice[idateCalib] = lv;
	setlocalVolSlice(idateCalib);

}

/****************************************************************
**	Class  : pdeLocalVolEffective 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::reinitialize(int idateCalib,const CVector& lcoeff)
{
	CVector coeff;

	coeff = lcoeff;
	m_localVolSlice[idateCalib]->setCoeff(coeff);
	setlocalVolSlice(idateCalib);
}


/****************************************************************
**	Class  : pdeLocalVolEffective 
**	Routine: getlocalVol
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::init(GVector<MlEqAnalyticCurveWithTanhWingEdgeHandle>& localVolSlice,MlEqAsset& asset,
								MlEqDateHandle startDate,long maturityDate,GVector<long>& calibDates,MlEqPdeHelperHandle pHelperOpt,
								const int nd, const int nx, const int* max_nx, const double lx, int nt,
								const double ux, const int keep)
{


	
	m_localVolSlice	= localVolSlice;		
	m_startDate		= startDate;
	m_calibDates	= calibDates;

	m_Grid.resize(m_calibDates.getsize());

	for ( int i = 0 ; i < m_calibDates.getsize(); i++ ){
		m_Grid[i].resize(*max_nx);
	}

	int n = startDate->GetDate();;
	long dDay = (maturityDate-n )/(double) nt;

	if ( dDay < 1 ){
		dDay = 1;
	}

	int k=1;
	for ( ;; )
	{
		n += dDay;
		if ( n > maturityDate ){
			break;
		}
		else
		{
			k++;
		}
	}

	int flag=0;
	if ( n < maturityDate ){
		flag = 1;
	}

	GVector<long> mdates(k+flag);
	mdates[0] = startDate->GetDate();

	int i;
	for (i=1; i < k; i++){
		mdates[i] = mdates[i-1]+dDay;
	}

	if ( flag ){
		mdates[i] = maturityDate;
	}

	std::vector<long> in1(mdates.getsize());
	std::vector<long> in2(m_calibDates.getsize());
	std::vector<long> outPutSet;

	for ( int i = 0 ; i < mdates.getsize(); i++ ){
		in1[i] = mdates[i];
	}

	for ( int i = 0 ; i < m_calibDates.getsize(); i++ ){
		in2[i] = m_calibDates[i];
	}

	merge(outPutSet,in1,in2);
	mdates.resize(outPutSet.size());

	for ( int i = 0; i < outPutSet.size(); i++ ){
		mdates[i] = outPutSet[i];
	}

	CVector t(mdates.getsize());

	for ( int i = 0 ; i < t.getsize(); i++ ){
		t[i] = m_startDate->GetYearFraction(mdates[i]);
	}
	nt = t.getsize();

//  merge t and t_extra here

	m_map.resize(t.getsize());
	m_mapinv.resize(m_calibDates.getsize());

	int ilast=0;
	n = 0;
	for ( int k = 0 ; k < m_calibDates.getsize(); k++ )
	{
		for ( int i = ilast ; i < t.getsize(); i++ )
		{
			m_map[i] = n;

			if ( mdates[i] == m_calibDates[k] )
			{
				m_mapinv[k] = i;
				n++;
				ilast = i+1;
				break;
			}
		}
	}




////////////////////

	const MlEqZeroCurveHandle	curve = asset.GetPayZeroCurve(true);


	m_loanspread.resize(nt);
	m_rate.resize(nt);
	m_forwards.resize(nt);

	int nToday = asset.GetCurrentDate();
	
	double spot = asset.GetSpot(nToday);
	m_forwards[0] = spot;

	double dt, mu, r, disc, discprev = 1., fwd, fwdprev = spot, xt, tprev = 0.0;


	for ( int idate = 0 ; idate < nt-1; idate++ )
	{
		fwd		   = asset.GetForward(nToday,mdates[idate+1], false);
		disc	   = curve->GetDiscountFactor(mdates[idate+1]);
		xt		   = t[idate+1];

		dt  = xt - tprev;	
		mu	= fwd / fwdprev;
		r	= disc / discprev;	 	

		m_loanspread[idate] = log(mu)/dt;
		m_rate[idate]		= -log(r)/dt;

		m_forwards[idate+1] = fwd ;
	
		fwdprev = fwd;	
		discprev = disc;	
		tprev = xt;
	}

	m_loanspread[nt-1]	= m_loanspread[nt-2];
	m_rate[nt-1]		= m_rate[nt-2];

	
///////////////////////////////

	int direction = PdeEnums::PDE_DIRECTION_BACKWARDS;
	int different = PdeEnums::PDE_MODEL_SAME;
	int continuous = 1;//PdeEnums::PDE_FACTORS_CONTINUOUS;
	double theta = 0.5;

	int conv_dependence=PdeEnums::PDE_X_DEPENDENT; 
	int diff_dependence=PdeEnums::PDE_X_DEPENDENT; 
	int var_dependence=PdeEnums::PDE_X_INDEPENDENT;
	int rhs_dependence=PdeEnums::PDE_X_INDEPENDENT;

	m_thetaDiffusion = PdeEnums::PDE_THETA_DIFFUSION;


	pde_init(nd, different, nx, max_nx, lx, ux, t, direction, keep,
				NULL, theta, continuous,
				conv_dependence, diff_dependence, 
				var_dependence, rhs_dependence,pHelperOpt);


}

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::convection_fcn(int id, double lx, double ux, int nx, int it, double t, 
													  void* ptr, CVector& f)
{


	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	double xvol;
	for ( int ix = 0 ; ix < nx; ix++ )
	{
		xvol = getlocalVol(it,ix);
		f[ix] = m_loanspread[itm] - 0.5 * pow(xvol,2.0);
	}
} 

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: convection_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::diffusion_fcn(int id, double lx, double ux, int nx, int it, double t, 
													   void* ptr, CVector& f)
{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 

	double xvol;
	for ( int ix = 0 ; ix < nx; ix++ )
	{
		xvol = getlocalVol(it,ix);
		f[ix] = 0.5 * pow(xvol,2.0);
	}
}

/****************************************************************
**	Class  : pdeLocalVol 
**	Routine: factor_fcn
**	Returns: 
**	Action : 
**           
****************************************************************/


void pdeLocalVolEffective::factor_fcn(int id, double lx, double ux, int nx, int it, double t, 

                    void* ptr, CVector& f)

{

	int itm = it;	
	if ( it == m_nt-1 ){
		itm--;} 


	f[0] = - m_rate[itm];
} 




