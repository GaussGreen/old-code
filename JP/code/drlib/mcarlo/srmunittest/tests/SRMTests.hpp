
//include functionality needed by this suite:
#include "edginc/HyperTrigUtils.hpp"

#include <boost/test/auto_unit_test.hpp>

//to avoid poles, make sure b >= a
struct HyperTrigTests
{
    //Reference values below were computed with Scilab
    void test_kfuncexp()
    {
        double a = 1.0;
        double b = 1.0;
        double c = 1.0;
        double x = 2.0;   
        double val;    
        HyperLocalVolState localVol(a,b,c);
        localVol.evaluateKExp(x, &val);
        BOOST_CHECK_CLOSE(val, 1.078183, 0.0001); 
    }  
   
    void test_kfuncexp2()  
    {
        double a = 2.0;
        double b = 2.1;
        double c = 1.2;
        double x = 5.0;    
        HyperLocalVolState localVol(a,b,c);
        double val;  
        localVol.evaluateKExp(x, &val);
        BOOST_CHECK_CLOSE(val,  1.3200792, 0.00001); 
    }

    void test_kfuncexp3()
    {
        double a = 2.0;
        double b = 2.1;
        double c = 1.2;
        double x = 5.5;    
        double val;  
        HyperLocalVolState localVol(a,b,c);
        localVol.evaluateKExp(x, &val);
        BOOST_CHECK_CLOSE(val, 1.4182694, 0.00001); 
        //1.4182693869626448
    }

    //half way between two tabulated values: quite accurate
    void test_kfuncinverse()
    {
        double val;
        HyperLocalVolState localVol(2.0, 2.1, 1.2);     
        localVol.evaluateInverseKExp(1.4182693869626448, &val);
        BOOST_CHECK_CLOSE(val, 5.5, 0.001); 
    }  

    void test_kfuncexp4()
    {
        double a = 2.0;
        double b = 2.1;
        double c = 1.2;
        double x = 5.1;    
        double val;  
        HyperLocalVolState localVol(a,b,c);
        localVol.evaluateKExp(x, &val);
        BOOST_CHECK_CLOSE(val, 1.3397249, 0.00001);
        //1.3397248814255143
    }

    //close to a tabulated value: very accurate
    void test_kfuncinverse2()
    {
        double val;
        HyperLocalVolState localVol(2.0, 2.1, 1.2);     
        localVol.evaluateInverseKExp(1.3397248814255143, &val);
        BOOST_CHECK_CLOSE(val, 5.1, 0.00001); 
    }


    //taken from a run:
    void test_kfuncexp5()
    {
        double a = 0.008356;
        double b = 3.49164;
        double c = 3.49164;
        double x = 0.06664;    
        double val;  
        HyperLocalVolState localVol(a,b,c);
        localVol.evaluateKExp(x, &val);
        BOOST_CHECK_CLOSE(val, 0.064617885842, 0.00001); 
    }
      
    //Note with this value of c, x it's very inaccurate!!
    void test_kfuncinverse3()  
    {
        double a = 0.008356;
        double b = 3.49164;
        double c = 3.49164;
        double val;
        HyperLocalVolState localVol(a,b,c);
        localVol.evaluateInverseKExp(0.064617885842, &val);
        BOOST_CHECK_CLOSE(val, 0.06664, 1.0); 
    }      

    void test_inverter1()
    {
        double a = 0.008356;
        double b = 3.49164;
        double c = 3.49164;
        double x = 0.06664;   
        double K, result;  
        HyperLocalVolState localVol(a,b,c);
        localVol.evaluateKExp(x, &K);
        localVol.evaluateInverseKExp(K, &result);
        BOOST_CHECK_CLOSE(result, x, 1.0); 
    }
}; 

struct HyperTrigTestSuite : public test_suite
{
    HyperTrigTestSuite() : test_suite("Hyperbolic test suite")
    {
        boost::shared_ptr<HyperTrigTests> instance( new HyperTrigTests );
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncexp, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncexp2, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncexp3, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncexp4, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncexp5, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncinverse, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncinverse2, instance ));
        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_kfuncinverse3, instance ));

        add( BOOST_CLASS_TEST_CASE( &HyperTrigTests::test_inverter1, instance ));
    }
};
