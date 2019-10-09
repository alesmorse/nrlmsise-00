#include "Nrlmsise00Test.hpp"

#include "Nrlmsise00.hpp"
#include <iostream>
#include <iomanip>

namespace test
{
    CPPUNIT_TEST_SUITE_REGISTRATION(Nrlmsise00Test);

    std::array<int,24> flags;
    int doy;
    double sec, alt, g_lat, g_long, lst, f107A, f107;
    std::array<double,7> ap;

    std::array<double, 9> density, ref_density;
    std::array<double, 2> temperatures, ref_temperatures;

    void Nrlmsise00Test::setUp()
    {
        // Set flags
        flags[0]=0; // output in cm and g
        for(unsigned i=1; i<flags.size(); i++)
            flags[i] =1;

        // Set default input values
        doy = 172;
        sec = 29000.0;
        alt = 400.0;
        g_lat = 60.0;
        g_long = -70.0;
        lst = 16.0;
        f107A = 150.0;
        f107 = 150.0;

        ap[0] = 4.0;

        // Set output arrays
        density.fill(0.0);
        temperatures.fill(0.0);

        // Set reference arrays
        ref_density.fill(0.0);
        ref_temperatures.fill(0.0);
    }

    void Nrlmsise00Test::CheckTemperatures(const std::array<double,2>& temp, const std::array<double,2>& ref_temp, const double tol)
    {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, temp.at(0)/ref_temp.at(0), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, temp.at(1)/ref_temp.at(1), tol);
    }

    void Nrlmsise00Test::CheckDensities(const std::array<double,9>& dens, const std::array<double,9>& ref_dens, const double tol)
    {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(0)/ref_dens.at(0), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(1)/ref_dens.at(1), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(2)/ref_dens.at(2), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(3)/ref_dens.at(3), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(4)/ref_dens.at(4), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(5)/ref_dens.at(5), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(6)/ref_dens.at(6), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(7)/ref_dens.at(7), tol);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, dens.at(8)/ref_dens.at(8), tol);
    }

    void Nrlmsise00Test::Test1()
    {
        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.250540E+03, 1.241416E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {6.665177E+05, 1.138806E+08, 1.998211E+07, 4.022764E+05, 3.557465E+03, 4.074714E-15, 3.475312E+04, 4.095913E+06, 2.667273E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test2()
    {
        doy = 81;
        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.166754E+03, 1.161710E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {3.407293E+06, 1.586333E+08, 1.391117E+07, 3.262560E+05, 1.559618E+03, 5.001846E-15, 4.854208E+04, 4.380967E+06, 6.956682E+03};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test3()
    {
        sec = 75000.0;
        alt = 1000.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.239892E+03, 1.239891E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {1.123767E+05, 6.934130E+04, 4.247105E+01, 1.322750E-01, 2.618848E-05, 2.756772E-18, 2.016750E+04, 5.741256E+03, 2.374394E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test4()
    {
        alt = 100.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures

        ref_temperatures = {1.027318E+03, 2.068878E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {5.411554E+07, 1.918893E+11, 6.115826E+12, 1.225201E+12, 6.023212E+10, 3.584426E-10, 1.059880E+07, 2.615737E+05, 2.819879E-42};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test5()
    {
        g_lat = 0.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.212396E+03, 1.208135E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {1.851122E+06, 1.476555E+08, 1.579356E+07, 2.633795E+05, 1.588781E+03, 4.809630E-15, 5.816167E+04, 5.478984E+06, 1.264446E+03};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test6()
    {
        g_long = 0.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.220146E+03, 1.212712E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {8.673095E+05, 1.278862E+08, 1.822577E+07, 2.922214E+05, 2.402962E+03, 4.355866E-15, 3.686389E+04, 3.897276E+06, 2.667273E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test7()
    {
        lst = 4.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.116385E+03, 1.112999E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {5.776251E+05, 6.979139E+07, 1.236814E+07, 2.492868E+05, 1.405739E+03, 2.470651E-15, 5.291986E+04, 1.069814E+06, 2.667273E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test8()
    {
        f107A = 70.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.031247E+03, 1.024848E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {3.740304E+05, 4.782720E+07, 5.240380E+06, 1.759875E+05, 5.501649E+02, 1.571889E-15, 8.896776E+04, 1.979741E+06, 9.121815E+03};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test9()
    {
        f107 = 180.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.306052E+03, 1.293374E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {6.748339E+05, 1.245315E+08, 2.369010E+07, 4.911583E+05, 4.578781E+03, 4.564420E-15, 3.244595E+04, 5.370833E+06, 2.667273E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test10()
    {
        ap[0] = 40.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.361868E+03, 1.347389E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {5.528601E+05, 1.198041E+08, 3.495798E+07, 9.339618E+05, 1.096255E+04, 4.974543E-15, 2.686428E+04, 4.889974E+06, 2.805445E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test11()
    {
        alt = 0.0;
        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.027318E+03, 2.814648E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {1.375488E+14, 0.000000E+00, 2.049687E+19, 5.498695E+18, 2.451733E+17, 1.261066E-03, 0.000000E+00, 0.000000E+00, 0.000000E+00};
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(0)/ref_density.at(0), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(1), density.at(1), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(2)/ref_density.at(2), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(3)/ref_density.at(3), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(4)/ref_density.at(4), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(5)/ref_density.at(5), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(6), density.at(6), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(7), density.at(7), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(8), density.at(8), 1e-15);
    }

    void Nrlmsise00Test::Test12()
    {
        alt = 10.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.027318E+03, 2.274180E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {4.427443E+13, 0.000000E+00, 6.597567E+18, 1.769929E+18, 7.891680E+16, 4.059139E-04, 0.000000E+00, 0.000000E+00, 0.000000E+00};
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(0)/ref_density.at(0), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(1), density.at(1), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(2)/ref_density.at(2), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(3)/ref_density.at(3), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(4)/ref_density.at(4), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(5)/ref_density.at(5), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(6), density.at(6), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(7), density.at(7), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(8), density.at(8), 1e-15);
    }

    void Nrlmsise00Test::Test13()
    {
        alt = 30.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.027318E+03, 2.374389E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {2.127829E+12, 0.000000E+00, 3.170791E+17, 8.506280E+16, 3.792741E+15, 1.950822E-05, 0.000000E+00, 0.000000E+00, 0.000000E+00};
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(0)/ref_density.at(0), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(1), density.at(1), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(2)/ref_density.at(2), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(3)/ref_density.at(3), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(4)/ref_density.at(4), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(5)/ref_density.at(5), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(6), density.at(6), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(7), density.at(7), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(8), density.at(8), 1e-15);
    }

    void Nrlmsise00Test::Test14()
    {
        alt = 50.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.027318E+03, 2.795551E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {1.412184E+11, 0.000000E+00, 2.104370E+16, 5.645392E+15, 2.517142E+14, 1.294709E-06, 0.000000E+00, 0.000000E+00, 0.000000E+00};
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(0)/ref_density.at(0), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(1), density.at(1), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(2)/ref_density.at(2), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(3)/ref_density.at(3), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(4)/ref_density.at(4), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(5)/ref_density.at(5), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(6), density.at(6), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(7), density.at(7), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(8), density.at(8), 1e-15);
    }

    void Nrlmsise00Test::Test15()
    {
        alt = 70.0;

        // Initialize model
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.027318E+03, 2.190732E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {1.254884E+10, 0.000000E+00, 1.874533E+15, 4.923051E+14, 2.239685E+13, 1.147668E-07, 0.000000E+00, 0.000000E+00, 0.000000E+00};
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(0)/ref_density.at(0), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(1), density.at(1), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(2)/ref_density.at(2), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(3)/ref_density.at(3), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(4)/ref_density.at(4), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.at(5)/ref_density.at(5), 5e-7);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(6), density.at(6), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(7), density.at(7), 1e-15);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ref_density.at(8), density.at(8), 1e-15);
    }

    void Nrlmsise00Test::Test16()
    {
        for (uint i=0; i<7; i++)
            ap[i] = 100.0;

        // Initialize model
        flags.at(9) = -1; // Use array ap
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.426412E+03, 1.408608E+03};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {5.196477E+05, 1.274494E+08, 4.850450E+07, 1.720838E+06, 2.354487E+04, 5.881940E-15, 2.500078E+04, 6.279210E+06, 2.667273E+04};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::Test17()
    {
        alt = 100.0;

        for (uint i=0; i<7; i++)
            ap[i] = 100.0;

        // Initialize model
        flags.at(9) = -1; // Use array ap
        atmos::CNrlmsise00 atmo(flags);

        // Compute densities and temperatures
        atmo.gtd7(doy,sec,alt,g_lat,g_long,lst,f107A,f107,ap,density, temperatures);

        // Check temperatures
        ref_temperatures = {1.027318E+03, 1.934071E+02};
        CheckTemperatures(temperatures, ref_temperatures);

        // Check densities
        ref_density = {4.260860E+07, 1.241342E+11, 4.929562E+12, 1.048407E+12, 4.993465E+10, 2.914304E-10, 8.831229E+06, 2.252516E+05, 2.415246E-42};
        CheckDensities(density, ref_density);
    }

    void Nrlmsise00Test::tearDown()
    {

    }
} // namespace test