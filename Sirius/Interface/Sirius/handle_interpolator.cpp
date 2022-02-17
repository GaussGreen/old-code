//	handle_interpolator.cpp : Implementation of CInterpolator
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_interpolator.h"

HRESULT CInterpolator::FinalConstruct(void)
{		
	return S_OK;
}

STDMETHODIMP CInterpolator::get_InterpolatorType(InterpolatorTypeEnum* pVal)
{	
	begin_function
	*pVal = (InterpolatorTypeEnum)m_h->m_type;	
	end_function
}

STDMETHODIMP CInterpolator::GetValue(double xValue, long whichData,double* pVal)
{	
	begin_function
	*pVal = m_h->getValue(xValue, whichData);
	end_function
}

STDMETHODIMP CInterpolator::get_Value(VARIANT *pVal)
{
	return CParameterMap::VectorToArray(m_vpm, pVal);
}

STDMETHODIMP CInterpolator::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IInterpolator };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CInterpolator::put_Value(VARIANT newVal)
{	
	/*parameter list is 0 - InterpolatorType
						1 - xData
						2 - yData
						3 - addTanhWings (this parameter and onwards for the CubicSpline case)
						4 - cL
						5 - cR
						6 - yPower*/
	
	HRESULT								hr;			
	std::vector<CParameterMap>			vpm;
	std::vector<CComVariant>			vv;
	
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, &vpm)) return hr;
	if (!vpm.size()) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		
	// the type of interpolator is the first parameter
	map_parameter(vpm[0], long, nInterpolatorType);
	
	switch (nInterpolatorType){

	case MonotonicSpline:
		if (vpm.size() != 7) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		{																		
			map_parameter(vpm[1], GVector<CVector>, xData);
			map_parameter(vpm[2], GVector<CVector>, yData);
			map_parameter(vpm[3], long, addTanhWings);
			map_parameter(vpm[4], double, cL);
			map_parameter(vpm[5], double, cR);
			map_parameter(vpm[6], double, yPower);			
			m_h = new MlEqMonotonicSplineInterpolator(xData, yData, addTanhWings, cL, cR, yPower);
		}
		break;
	case CubicSpline:
		if (vpm.size() != 7) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		{
			map_parameter(vpm[1], GVector<CVector>, xData);
			map_parameter(vpm[2], GVector<CVector>, yData);
			map_parameter(vpm[3], long, addTanhWings);
			map_parameter(vpm[4], double, cL);
			map_parameter(vpm[5], double, cR);
			map_parameter(vpm[6], double, yPower);			
			m_h = new MlEqCubicSplineInterpolator(xData, yData, addTanhWings, cL, cR, yPower);
		}
		break;
	case Linear:
		if (vpm.size() == 1){
			m_h = new MlEqInterpolator(CVector(), CVector());
		} else if (vpm.size() == 3){			
			map_parameter(vpm[1], GVector<CVector>, xData);
			map_parameter(vpm[2], GVector<CVector>, yData);						
			m_h = new MlEqInterpolator(xData, yData);			
		} else if ( vpm.size() == 7 ){
			map_parameter(vpm[1], GVector<CVector>, xData);
			map_parameter(vpm[2], GVector<CVector>, yData);
			map_parameter(vpm[3], long, addTanhWings);
			map_parameter(vpm[4], double, cL);
			map_parameter(vpm[5], double, cR);
			map_parameter(vpm[6], double, yPower);			
			m_h = new MlEqLinearInterpolator(xData, yData, addTanhWings, cL, cR, yPower);
		} else {			
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		}
		break;
	case Constant:
		if (vpm.size() == 7){
			map_parameter(vpm[1], GVector<CVector>, xData);
			map_parameter(vpm[2], GVector<CVector>, yData);
			map_parameter(vpm[3], long, addTanhWings);	
			map_parameter(vpm[4], double, cL);
			map_parameter(vpm[5], double, cR);
			map_parameter(vpm[6], double, yPower);			
			m_h = new MlEqConstantInterpolator(xData, yData, addTanhWings, cL, cR, yPower);
		} else if (vpm.size() == 3){
			map_parameter(vpm[1], GVector<CVector>, xData);
			map_parameter(vpm[2], GVector<CVector>, yData);			
			m_h = new MlEqConstantInterpolator(xData, yData);
		} else {
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		}
		break;


	case FitPolynomial:

		if (vpm.size() == 11){
			map_parameter(vpm[1], CVector, xData);
			map_parameter(vpm[2], CVector, yData);
			map_parameter(vpm[3], CVector, initialGuess);
			map_parameter(vpm[4], CMatrix, xValBounds);
			map_parameter(vpm[5], long,    addTanhWings);	
			map_parameter(vpm[6], double,  cL);
			map_parameter(vpm[7], double,  cR);
			map_parameter(vpm[8], double,  yPower);			
			map_optional_parameter(vpm[9], double, finalTol,1e-12);
			map_optional_parameter(vpm[10],double, stopTol,1e-12);			

			m_h = new MlEqAnalyticCurveWithTanhWingEdge(xData,yData,initialGuess,xValBounds,addTanhWings,cL,cR,yPower,finalTol,stopTol);

		} else 
		{
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		}


		break;


	case FitPolynomialExlicit:
		if (vpm.size() == 9)
		{

			map_parameter(vpm[1], CVector, coefficients);
			map_parameter(vpm[2], double, lowerEdge);
			map_parameter(vpm[3], double, upperEdge);
			map_parameter(vpm[4], long, addTanhWing);
			map_parameter(vpm[5], double, cL);
			map_parameter(vpm[6], double, cR);
			map_parameter(vpm[7], double, yPower);
			map_optional_parameter(vpm[8], long, npoints,10);

			m_h = new MlEqAnalyticCurveWithTanhWingEdge(lowerEdge,upperEdge,addTanhWing,cL,cR,yPower,npoints,coefficients);

		} else 
		{
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		}
		break;

	case TwoDimensional:
		if (vpm.size() == 3){
			map_object_parameter(vpm[1], Interpolator, InterpolatorHandleY);
			map_object_vector_parameter(vv[2], Interpolator, InterpolatorHandleX);
			m_h = new MlEq2DInterpolator(InterpolatorHandleY, InterpolatorHandleX);
		} else {
			return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);
		}
		break;
	default:
		return CParameterMap::ReturnErrorR(IDS_INTERPOLATOR_TYPE, IID_IInterpolator);
	}		
	m_vpm = vpm;
	end_function
}