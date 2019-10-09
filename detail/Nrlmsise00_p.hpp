#ifndef atmos_NRLMSISE00_IMP_H
#define atmos_NRLMSISE00_IMP_H

#include <cmath>
#include <array>

namespace atmos
{
    class CNrlmsise00_p
    {
    private:
        std::array<int, 24>    a_switches;
        std::array<double, 24> a_sw;
        std::array<double, 24> a_swc;

        // PARAMB
        double d_gsurf;
        double d_re;

        // GTS3C
        double d_dd;

        // DMIX
        double d_dm04, d_dm16, d_dm28, d_dm32, d_dm40, d_dm01, d_dm14;

        // MESO7
        std::array<double,5> a_meso_tn1;
        std::array<double,4> a_meso_tn2;
        std::array<double,5> a_meso_tn3;
        std::array<double,2> a_meso_tgn1;
        std::array<double,2> a_meso_tgn2;
        std::array<double,2> a_meso_tgn3;

        double d_dfa;
        std::array<std::array<double,9>,4> a_plg;
        double d_ctloc, d_stloc;
        double d_c2tloc, d_s2tloc;
        double d_c3tloc, d_s3tloc;
        double d_apdf;
        std::array<double,4> a_apt;

        /**
         * @brief Calculate latitude variable gravity (GV) and effective radius (REFF)
         *
         * @param lat latitude value
         * @param gv variable gravity
         * @param reff effective radius
         */
        static void glatf(const double& lat, double& gv, double& reff);

        /**
         * @brief Chemistry/dissociation correction for msis models
         *
         * @param alt Altitude
         * @param r Target ratio
         * @param h1 Transition scale length
         * @param zh Altitude of 1/2 r
         * @return double
         */
        static double ccor(const double& alt, const double& r, const double h1, const double zh);

        /**
         * @brief O&O2 Chemistry/dissociation correction for msis models
         *
         * @param alt Altitude
         * @param r Target ratio
         * @param h1 Transition scale length
         * @param zh Altitude of 1/2 r
         * @param h2 Transition scale length for O2
         * @return double
         */
        static double ccor2(const double& alt, const double& r, const double h1, const double zh, const double h2);
        double scalh(const double& alt, const double xm, const double temp) const;

        /**
         * @brief Turbopause correction for msis models
         *
         * @param dd diffusive density
         * @param dm full mixed density
         * @param zhm transition scale length
         * @param xmm full mixed molecular weight
         * @param xm species molecular weight
         * @return double combined density
         */
        double dnet(double& dd, const double& dm, const double& zhm, const double& xmm, const double xm) const;

        double zeta(const double& zz, double z) const;

        /**
         * @brief Calculate Temperature and Density Profiles for lower atmos
         *
         * @param alt altitude above surface (km)
         * @param d0
         * @param xm
         * @param tz
         * @param mn3
         * @param zn3
         * @param tn3
         * @param tgn3
         * @param mn2
         * @param zn2
         * @param tn2
         * @param tgn2
         * @return double
         */
        double densm (const double& alt, const double& d0, const double xm, double *tz,
                 const int mn3, const double *zn3, const double *tn3, const double *tgn3,
                 const int mn2, const double *zn2, const double *tn2, const double *tgn2) const;

        /**
         * @brief Calculate Temperature and Density Profiles for MSIS models
         *
         * @param alt altitude above surface (km)
         * @param dlb
         * @param tinf
         * @param tlb
         * @param xm
         * @param alpha
         * @param tz
         * @param zlb
         * @param s2
         * @param mn1
         * @param zn1
         * @param tn1
         * @param tgn1
         * @return double
         */
        double densu (const double& alt, const double& dlb, const double& tinf, const double& tlb, const double& xm,
                 const double alpha, double *tz, const double zlb, const double& s2,
                 const int mn1, const double *zn1, double *tn1, double *tgn1) const;

        /**
         * @brief Calculate G(L) function
         *
         * @param p coefficients
         * @param doy day of year
         * @param sec seconds in day
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @return double G(L) value
         */
        double globe7(const std::array<double,150>& p, const int doy, const double sec,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap);

        /**
         * @brief Calculate G(L) function for lower atmosphere
         *
         * @param p coefficients
         * @param doy day of year
         * @param g_long geodetic longitude
         * @return double G(L) value
         */
        double glob7s(const std::array<double,100>& p, const int doy,const double& g_long);

        /**
         * @brief Thermospheric portion of NRLMSISE-00
         *
         * @param doy day of year
         * @param sec seconds in day
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @param d density array
         * @param t temperature array
         */
        void gts7(const int doy, const double sec, const double& alt,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t);

    public:

        /**
         * @brief Construct a new CNrlmsise00_p object
         *
         * @param switches array with nrlmsise-00 switches
         *
         * The 24 switches have the following meaning:
         * -  0: output in meters and kilograms instead of centimeters and grams
         * -  1: F10.7 effect on mean
         * -  2: time independent
         * -  3: symmetrical annual
         * -  4: symmetrical semiannual
         * -  5: asymmetrical annual
         * -  6: asymmetrical semiannual
         * -  7: diurnal
         * -  8: semidiurnal
         * -  9: daily ap [when this is set to -1 (!) the pointer ap_a in struct nrlmsise_input must point to a struct ap_array]
         * - 10: all UT/long effects
         * - 11: longitudinal
         * - 12: UT and mixed UT/long
         * - 13: mixed AP/UT/LONG
         * - 14: terdiurnal
         * - 15: departures from diffusive equilibrium
         * - 16: all TINF var
         * - 17: all TLB var
         * - 18: all TN1 var
         * - 19: all S var
         * - 20: all TN2 var
         * - 21: all NLB var
         * - 22: all TN3 var
         * - 23: turbo scale height var
         */
        CNrlmsise00_p(const std::array<int, 24>& switches);

        /**
         * @brief Neutral Atmosphere Empircial Model from the surface to lower exosphere.
         *
         * @param doy day of year
         * @param sec seconds in day
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @param d density array
         * @param t temperature array
         */
        void gtd7(const int doy, const double sec, const double& alt,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t);

        /**
         * @brief Include the anomalous oxygen contribution.
         *
         * This subroutine provides Effective Total Mass Density for output
         *   d[5] which includes contributions from "anomalous oxygen" which can
         *   affect satellite drag above 500 km.
         *
         * @param doy day of year
         * @param sec seconds in day
         * @param alt altitude above surface (km)
         * @param g_lat geodetic latitude
         * @param g_long geodetic longitude
         * @param lst local apparent solar time (hours)
         * @param f107A 81 day average of F10.7 flux (centered on doy)
         * @param f107 daily F10.7 flux for previous day
         * @param ap magnetic index array
         * @param d density array
         * @param t temperature array
         */
        void gtd7d(const int doy, const double sec, const double& alt,
                  const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                  std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t);
    };

} // namespace atmos
#endif // atmos_NRLMSISE00_IMP_H
