#ifndef _FUNCTION_H
#define _FUNCTION_H

#include "linalg.h"
#include "armglob.h"


extern void Dummy_Func(void);


#define SWITCH_FOR_NUMERICAL_DERIVATIVE 1.0e-7
double Identity(double x)
{
	return x;
}
class RealValuedFunction
{
    public:
  
             RealValuedFunction(void){};

             virtual ~RealValuedFunction(){};

             virtual RealValuedFunction* Duplicate() const 
             {
                 return(NULL);
             }

             virtual double operator () ( const double& x ) const
             {
                 return(0.0);
             }		 
 
};

class RealNDValuedFunction
{
    public:
  
             RealNDValuedFunction(void){};

             virtual ~RealNDValuedFunction(){};

             virtual RealNDValuedFunction* Duplicate() const 
             {
                 return(NULL);
             }

             virtual double operator () ( const vector <double > x ) const
             {
                 return(0.0);
             }		 
 
};


class FunctionNDWithKnownDerivative: public RealNDValuedFunction
{
    public:
  
        FunctionNDWithKnownDerivative(void) {};

        // oblige de retourner un objet alloue sur le tas

        virtual RealNDValuedFunction* Derivative() const
        {
            return((RealNDValuedFunction *) 0);
        };
};

typedef double (*RealFctPtr)(double);



class SimpleFct: public RealValuedFunction
{
    public:

        SimpleFct(void) {};


        SimpleFct( RealFctPtr f_ ) :f(f_){};

        //DECLARE_DUPLICATE(RealValuedFunction,SimpleFct);
  
        virtual double operator () ( const double& x ) const{ return f(x); }

    protected:

        RealFctPtr f;
};


class FunctionWithKnownDerivative: public RealValuedFunction
{
    public:
  
        FunctionWithKnownDerivative(void) {};

        // oblige de retourner un objet alloue sur le tas

        virtual RealValuedFunction* Derivative() const
        {
            return((RealValuedFunction *) 0);
        };
};



class FunctionWithKnownPrimitive: public RealValuedFunction
{
    public:
  
        FunctionWithKnownPrimitive(void) {};

        virtual RealValuedFunction* Primitive() const 
        {
           return((RealValuedFunction *) 0);
        };
};

inline double CentredNumericalDerivative( const RealValuedFunction& f, double x, double h )
{    
	return (f(x+h)-f(x-h))/(2.0*h);
}

class DerivativeFunc : public RealValuedFunction
{
private:    
	const RealValuedFunction& itsFunction;
	double itsEpsilon;
	
public:
	virtual ~DerivativeFunc()
    {};

    DerivativeFunc(const RealValuedFunction& f,
                    double epsilon = 0.00001)
    :	
    itsFunction(f),
    itsEpsilon(epsilon)
	{}; 

    /// copy constructor
    DerivativeFunc( const DerivativeFunc& rhs )
    :
        RealValuedFunction( rhs ),
        itsFunction( rhs.itsFunction ),
        itsEpsilon( rhs.itsEpsilon )
    {}

	virtual double operator () ( const double& x ) const
    {
        return CentredNumericalDerivative(itsFunction,x,itsEpsilon);    
    }
private:
    /// assignment operator (unavailable)
    DerivativeFunc& operator=( const DerivativeFunc& rhs );

};


inline double NumericalDerivative( const RealValuedFunction& f, double x, double h )
{    
	return (fabs(x*h) < SWITCH_FOR_NUMERICAL_DERIVATIVE?(f(x+h)-f(x))/h
          :(f(x*(1+h))-f(x))/(h*x));
}

class CstStepwiseFunction : public ARM_Object
{
    private:

        ARM_Vector* X;
        ARM_Vector* Y;



        void Init(void)
        {
            X = (ARM_Vector *) 0;
            Y = (ARM_Vector *) 0;
        }

    public:

        CstStepwiseFunction(void)
        {
            Init();
        }

       ~CstStepwiseFunction(void)
        {
            if (X)
               delete X;

            if (Y)
               delete Y;
        }

        CstStepwiseFunction(long sz, double* inX, double* inY)
        {
            Init();


            X = new ARM_Vector(sz+2);

            Y = new ARM_Vector(sz+2);

            (*X)[0] = -std::numeric_limits<double>::infinity();
            (*Y)[0] = inY[0];
            long i;
            for (i = 0; i < sz; i++)
            {
                (*X)[i+1] = inX[i];
                (*Y)[i+1] = inY[i];
            } 

            (*X)[sz+1] = std::numeric_limits<double>::infinity();
            (*Y)[sz+1] = (*Y)[sz];

            // Divide volat by 100 in order to transform Exp: 5% en 0.05 !

            for (i = 0; i < (sz+2); i++)
            {
                 (*Y)[i] =  (*Y)[i]/100.0; 
            }
        };

        void BitwiseCopy(const ARM_Object* oFunc)
        {
            CstStepwiseFunction* cstFunc = (CstStepwiseFunction *) oFunc;
 
 
            if (X)
               delete X;

            X = (ARM_Vector *) cstFunc->X->Clone();

            if (Y)
               delete Y;

            Y = (ARM_Vector *) cstFunc->Y->Clone();
        }

        void Copy(const ARM_Object* f)
        {
            ARM_Object::Copy(f);

            BitwiseCopy(f);
        }
 
        ARM_Object* Clone(void)
        {
            CstStepwiseFunction* theClone = new CstStepwiseFunction();
 
            theClone->Copy(this);
 
            return(theClone);
        }

        CstStepwiseFunction& CstStepwiseFunction::operator = 
                         (const CstStepwiseFunction& oFunc)
        {
            (*this).ARM_Object::operator = (oFunc);

            this->BitwiseCopy(&oFunc);
 
            return(*this);
        }

        // la liste est stocke sous la forme (t_i,x_i) 
        // avec x_i valeur sur [t_{i-1},t_{i}[ et t_{-1}=-\infty
        // on extrapole au dela de t_N par x_N

        CstStepwiseFunction(const CstStepwiseFunction& f) : ARM_Object(f)
        {
            Init();

            BitwiseCopy(&f);
        }

        
        void ShiftShortVol(double volShift)
        {
            int sz = Y->GetSize();


            for (long i = 0; i < sz; i++)
            {
                (*Y)[i] += volShift;
            } 
        }


        void MultiplyShortVol(double volShift)
        {
            int sz = Y->GetSize();


            for (long i = 0; i < sz; i++)
            {
                (*Y)[i] *= volShift;
            } 
        }


        void Shift(double Epsilon)
        {
            int sz = Y->GetSize();
 
 
            for (long i = 0; i < sz; i++)
            {
                (*Y)[i] *= ( 1 + Epsilon );
            }
        }

        CstStepwiseFunction* Duplicate() const
        { 
            return new CstStepwiseFunction(*this); 
        }

        virtual double operator () (const double& t ) const
        {
            long last_index = X->GetSize()-1;

            if ( t >= (*X)[last_index])
               return (*Y)[last_index];

            // retourner x_{i+1} | { t_i < t <= t_{i+1} }

            long i = 0;

            while ((*X)[i] < t )
                i++;

            return((*Y)[i]);
        }
 
 
		virtual double operator [](const int& t)  const
        {
            long last_index = X->GetSize()-1;
 
            if ( t >= last_index)
               return (*X)[last_index];
 
            return((*X)[t]);
        }
 
 
        // Integrate[transform(this(x))*f(x),{x,a,b}]

        double IntegrateProductWith(const FunctionWithKnownPrimitive& f,
                                    const double& a,
                                    const double& b,
                                    const RealValuedFunction& transform ) const
        {
            double result = 0.0;

            long last_index = X->GetSize()-1;

            RealValuedFunction* primitive_f = f.Primitive();

            // i_inf = inf { i | t_i > a }

            long i_inf=0;

            while ( (*X)[i_inf] <= a )
                 i_inf++;

            // i_inf OK
    
            // Si [a,b]\in[t_i,t_{i+1}]

            if ( b <= (*X)[i_inf] )
               result = transform((*Y)[i_inf])*((*primitive_f)(b)-
                     (*primitive_f)(a));
            else
            {
               result = transform((*Y)[i_inf])
                        *((*primitive_f)((*X)[i_inf])-(*primitive_f)(a));
      
               // i_sup = sup { i | t_i < b }
               // i_sup = inf { i | t_{i+1} >= b }

               long i_sup=i_inf;
      
               while (( i_sup < last_index ) && ( (*X)[i_sup+1] < b )) 
               {
                   result += transform((*Y)[i_sup+1])
                               *((*primitive_f)((*X)[i_sup+1])-
                                (*primitive_f)((*X)[i_sup]));

                   i_sup++;
               }
      
               // i_sup OK

               long index;

               if ( b < (*X)[last_index] )
                  index = i_sup+1;
               else
                  index = last_index;
      
               result += transform((*Y)[index])*((*primitive_f)(b)-
                          (*primitive_f)((*X)[i_sup]));
            }

            delete primitive_f;
    
            return result;
        }

        ARM_Vector* GetSigmaVect(void)
        {
            return(Y);
        }
        
		 ARM_Vector* GetTimeVect(void)
        {
            return(X);
        }
        // Fonction utilitaire 

        ARM_Vector* FromSigmaToPara(double a)
        {
            long sz = Y->GetSize();

            ARM_Vector* sigVect = new ARM_Vector(sz-1);

            
            for (long i = 0; i < (sz-1); i++)
            {
                (*sigVect)[i] = (*Y)[i]*100.0;
            }
 
            (*sigVect)[0] = a;

            return(sigVect);
        }
 
        ARM_Vector* FromSigma2FToPara(double a, double b, double rho, double ratio)
        {
            long sz = Y->GetSize();
 
            ARM_Vector* sigVect = new ARM_Vector(sz+2);
 
            
            for (long i = 4; i < (sz+2); i++)
            {
                (*sigVect)[i] = (*Y)[i-3]*100.0;
            }
 
            (*sigVect)[0] = a;
            (*sigVect)[1] = b;
            (*sigVect)[2] = rho;
            (*sigVect)[3] = ratio;
 
            return(sigVect);
        }
 
 
        long NbPlots(void)
        {
            return(X->GetSize()-2);
        }

        void UpDateYValues(ARM_Vector* Yvalues)
        {
            if ( X == NULL )
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                         "X is NULL in : UpDateYValues: Can't be updated");

               return;
            }

            if (( Yvalues == NULL ) 
                || 
                ( (Yvalues->GetSize()+2) != X->GetSize() )
               )
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                "X and Yvalues has different sizes in: UpDateYValues");

               return;
            }

         /* SG
            for (long i = 0; i < Yvalues->GetSize(); i++ )
            {
                (*Y)[i+1] = (*Yvalues)[i]/100.0;
            }

            (*Y)[0] = (*Y)[1];
            (*Y)[Yvalues->GetSize()+1] = (*Y)[Yvalues->GetSize()]; 

          */

            double* XElements = X->GetElt();

            CstStepwiseFunction tmpSigma(Yvalues->GetSize(), 
                                         &XElements[1], Yvalues->GetElt());

          
            *this = tmpSigma; 
        }


		void ShiftXValues(double shiftX)
		{
			if ( X == NULL )
			{
				throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,"Unimplemented <Shift> method");
				return;			 
			}

			for (long i = 0 ; i < X->GetSize() ; i++)
			{
				(*X)[i] =(*X)[i]+shiftX; 
			}
		}

};



typedef CstStepwiseFunction RealStepFct;


class GeneralizedExp: public FunctionWithKnownPrimitive
{
    public:
 
        GeneralizedExp(void) {};

        GeneralizedExp( double a_, double x0_, double factor_=1.0 )
                        :a(a_),x0(x0_),factor(factor_)
        {};
 
        //DECLARE_DUPLICATE(RealValuedFunction,GeneralizedExp);
 
        virtual double operator () ( const double& x ) const
        {
            return factor*exp(a*(x-x0));
        }
 
        virtual RealValuedFunction* Primitive() const
        {
            return new GeneralizedExp(a,x0,factor/a);
        }
 
    protected:
 
        double a, x0, factor;
};
 
class ConstantFct: public FunctionWithKnownPrimitive
{
    public:
 
        ConstantFct(void) {};
 
        
        //DECLARE_DUPLICATE(RealValuedFunction,ConstantFct);
 
        virtual double operator () ( const double& x ) const
        {
            return 1.0;
        }
 
        virtual RealValuedFunction* Primitive() const
        {
            return new SimpleFct(Identity);
        }
 
    protected:
 
        double a, x0, factor;
};


#undef SWITCH_FOR_NUMERICAL_DERIVATIVE



#endif
/*----------------------------------------------------------------------*/
/*---- End of File ----*/
