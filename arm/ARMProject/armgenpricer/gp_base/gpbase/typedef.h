
#ifndef _INGPBASE_TYPEDEF_H
#define _INGPBASE_TYPEDEF_H

#include "removeidentifiedwarning.h"
#include "port.h"
#include "countedptr.h"
#include "functor.h"
#include "curvetypedef.h"
#include "surfacetypedef.h"
#include "gplinalgtypedef.h"
#include "enumbase.h"

#include <vector>
CC_USING_NS(std,vector)

/// forward declaration in the ARM_GP namespace
CC_BEGIN_NAMESPACE( ARM_GP )
	template<typename T> struct ReturnFunc;
	template<typename T,typename U> struct UnaryFunc;
CC_END_NAMESPACE()

/// forward declaration of no namespace object
class ARM_Object;

/// Start ARM namespace
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration in the ARM namespace
class ARM_DateStrip;
class ARM_GramFctorArgDict;
struct ARM_InterpolationType;
template <typename T,typename U> struct NumDerivativeFunc;
typedef ARM_InterpolationType::InterpolationType ARM_InterpolType;


/// typedef from the ARM_GP namespace
typedef CC_NS( ARM_GP, ReturnFunc  )<double>			ARM_VoidToDbleFunctor;
typedef CC_NS( ARM_GP, UnaryFunc  )<double,double>		ARM_DbleToDbleFunctor;
typedef CC_NS( ARM_GP, BinaryFunc )<double>				DbleBinaryFunctor;

typedef NumDerivativeFunc<double,double>				ARM_NumDerivativeDbleToDbleFunctor;

/// vector declaration
typedef vector< ARM_GP_T_Vector<double>* >              ARM_VectorVector;
typedef vector< ARM_DateStrip* >						ARM_DateStripVector;
typedef vector< ARM_GP_TriangularMatrix* >				ARM_TriangularMatrixVector;
typedef vector< ARM_GP_Matrix* >						ARM_MatrixVector;
typedef vector< string >                                ARM_StringVector;

typedef vector< ARM_StringVector >		                ARM_StringVectorVector;
typedef vector< ARM_Object* >                           ARM_ObjectVector;

/// reference counted pointor
typedef ARM_CountedPtr< std::vector<double> >			ARM_VectorPtr;
typedef ARM_CountedPtr< ARM_GP_T_Vector<double> >		ARM_GP_VectorPtr;
typedef	ARM_CountedPtr< ARM_Curve >					    ARM_GP_CurvePtr;
typedef	ARM_CountedPtr< ARM_GP_Matrix >					ARM_GP_MatrixPtr;
typedef	ARM_CountedPtr< ARM_MultiCurve >				ARM_GP_MultiCurvePtr;
typedef	ARM_CountedPtr< ARM_GP_Tensor >					ARM_GP_TensorPtr;
typedef	ARM_CountedPtr< ARM_GP_TriangularMatrix >		ARM_GP_TriangularMatrixPtr;
typedef ARM_CountedPtr< ARM_Surface >					ARM_SurfacePtr;
typedef ARM_CountedPtr< ARM_StringVector >				ARM_StringVectorPtr;
typedef ARM_CountedPtr< ARM_GP_StrVector >				ARM_GP_StrVectorPtr;


typedef CC_NS(ARM_GP,UnaryFunc)<ARM_VectorPtr,double>   ARM_VectorPtrDbleFunc;		/// here because of the necessity to be after the typedef of vectorPtr
typedef ARM_CountedPtr< ARM_VoidToDbleFunctor >			ARM_VoidToDbleFunctorPtr;
typedef ARM_CountedPtr< ARM_VectorPtrDbleFunc >			ARM_VectorPtrDbleFuncPtr;
typedef ARM_CountedPtr< ARM_DateStrip >					ARM_DateStripPtr;
typedef ARM_CountedPtr< ARM_Object >				    ARM_ObjectPtr;
typedef ARM_CountedPtr< ARM_IntVector>		        	ARM_IntVectorPtr;
typedef ARM_CountedPtr< ARM_IntMatrix>		        	ARM_IntMatrixPtr;


/// vector of counted ptr
typedef vector< ARM_GP_MatrixPtr >						ARM_MatrixPtrVector;
typedef vector< ARM_VectorPtrDbleFuncPtr >				ARM_VectorPtrDbleFuncPtrVector;
typedef vector< ARM_SurfacePtr >						ARM_SurfacePtrVector;

/// vector of counter ptr
typedef vector< ARM_GP_VectorPtr >							ARM_VectorPtrVector;
typedef ARM_CountedPtr< ARM_VectorPtrVector >			ARM_VectorPtrVectorPtr;


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
