#pragma once

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

namespace test
{
    class Nrlmsise00Test : public CppUnit::TestFixture
    {
        CPPUNIT_TEST_SUITE(Nrlmsise00Test);
          CPPUNIT_TEST(Test1);
          CPPUNIT_TEST(Test2);
          CPPUNIT_TEST(Test3);
          CPPUNIT_TEST(Test4);
          CPPUNIT_TEST(Test5);
          CPPUNIT_TEST(Test6);
          CPPUNIT_TEST(Test7);
          CPPUNIT_TEST(Test8);
          CPPUNIT_TEST(Test9);
          CPPUNIT_TEST(Test10);
          CPPUNIT_TEST(Test11);
          CPPUNIT_TEST(Test12);
          CPPUNIT_TEST(Test13);
          CPPUNIT_TEST(Test14);
          CPPUNIT_TEST(Test15);
          CPPUNIT_TEST(Test16);
          CPPUNIT_TEST(Test17);
        CPPUNIT_TEST_SUITE_END();

     public:
        /**
         * Method automatically called before each test by CppUnit
         */
        void setUp();
        /**
         * Method automatically called after each test CppUnit
         */
        void tearDown();
        /**
         * Check temperatures against reference
         */
        void CheckTemperatures(const std::array<double,2>& temp, const std::array<double,2>& ref_temp, const double tol=5e-7);
        /**
         * Check densities against reference
         */
        void CheckDensities(const std::array<double,9>& dens, const std::array<double,9>& ref_dens, const double tol=5e-7);
        /**
         * Test the standard Gtd7 model
         */
        void Test1();
        void Test2();
        void Test3();
        void Test4();
        void Test5();
        void Test6();
        void Test7();
        void Test8();
        void Test9();
        void Test10();
        void Test11();
        void Test12();
        void Test13();
        void Test14();
        void Test15();
        void Test16();
        void Test17();
    };
} // namespace test