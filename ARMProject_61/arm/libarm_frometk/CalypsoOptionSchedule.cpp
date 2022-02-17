/**********************************************************************
*
*	C++ %name:		CalypsoOptionSchedule.cpp %
*	Instance:		1
*	Description:	
*	%created_by:	fdegraeve %
*	%date_created:	Wed May 30 17:25:25 2007 %
*
**********************************************************************/

#include "CalypsoOptionSchedule.h"

ARM_ReferenceValue* convertCalypsoOptionSchedule(MSXML2::IXMLDOMNodePtr xmlNode) {
	ARM_ExerciseStyle* exerStyle = NULL;
	ARM_ReferenceValue* fees = NULL;
	convertCalypsoOptionSchedule(xmlNode, fees, exerStyle);
	if(exerStyle) {
		delete exerStyle;
		exerStyle = NULL;
	}
	return fees;
}

void convertCalypsoOptionSchedule(MSXML2::IXMLDOMNodePtr xmlNode,  ARM_ReferenceValue*& fees, ARM_ExerciseStyle*& exerStyle) {
	vector<double> exerciseDateVector;
	vector<double> deliveryDateVector;
	vector<double> dFees;
	ARM_Vector* exerDates = NULL;
	ARM_Vector* delivDates = NULL;

	//jclass clsOption = jni->GetObjectClass(optionSchedule); 

	//jmethodID mid= jni->GetMethodID(clsOption, "getCount", "()I");
	//jint exerciseCount = jni->CallIntMethod(optionSchedule, mid);

	//mid = jni->GetMethodID(clsOption, "getExerciseDates", "()[J");
	//jlongArray exerciseDatesArray = (jlongArray) jni->CallObjectMethod(optionSchedule, mid);
	//jlong * exerciseDates = jni->GetLongArrayElements(exerciseDatesArray, 0);

	//mid = jni->GetMethodID(clsOption, "getNoticeDates", "()[J");
	//jlongArray deliveryDatesArray = (jlongArray) jni->CallObjectMethod(optionSchedule, mid);
	//jlong * deliveryDates = jni->GetLongArrayElements(deliveryDatesArray, 0);

	//mid = jni->GetMethodID(clsOption, "getFees", "()[D");
	//jdoubleArray exerciseFeeValuesArray = (jdoubleArray) jni->CallObjectMethod(optionSchedule, mid);
	//jdouble * exerciseFeeValues = jni->GetDoubleArrayElements(exerciseFeeValuesArray, 0);

     
	//for (int i = 0; i< exerciseCount; i++) {
	//	exerciseDateVector.push_back((double) exerciseDates[i]);
	//	deliveryDateVector.push_back((double) deliveryDates[i]);
	//	dFees.push_back((double) exerciseFeeValues[i]/100);
	//}

	//jni->ReleaseLongArrayElements(exerciseDatesArray, exerciseDates, 0);
	//jni->ReleaseLongArrayElements(deliveryDatesArray, deliveryDates, 0);
	//jni->ReleaseDoubleArrayElements(exerciseFeeValuesArray, exerciseFeeValues, 0);

    long exerciseCount;
    MSXML2::IXMLDOMNodeListPtr xmlNodeList = XMLTools::selectNodes(xmlNode,"putCalldate"); 
    xmlNodeList->get_length(&exerciseCount);
    ARM_Date exerciseDate,deliveryDate;
    double exerciseFeeValue;

    for (int i = 0; i< exerciseCount; i++) {
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"expiryDate"),exerciseDate,"YYYYMMDD");
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"deliveryDate"),deliveryDate,"YYYYMMDD");
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"premium"),exerciseFeeValue);
			
        exerciseDateVector.push_back(exerciseDate.GetJulian());
		deliveryDateVector.push_back(deliveryDate.GetJulian());
		dFees.push_back(exerciseFeeValue/100);
	}

	// Create object
	exerDates = new ARM_Vector(exerciseDateVector);
	delivDates = new ARM_Vector(deliveryDateVector);
	exerStyle = new ARM_ExerciseStyle(exerDates, delivDates);
	fees = createRefValue(exerciseDateVector, dFees, K_STEPUP_RIGHT, 1, false);
	if(exerDates) {
		delete exerDates;
		exerDates = NULL;
	}
	if(delivDates) {
		delete delivDates;
		delivDates = NULL;
	}
}