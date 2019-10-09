#include "Nrlmsise00.hpp"

namespace atmos
{
    CNrlmsise00::CNrlmsise00(const std::array<int, 24>& flags)
    {
        p_detail.reset(new CNrlmsise00_p(flags));
    }

    void CNrlmsise00::gtd7(const int doy, const double sec, const double& alt,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t)
    {
        p_detail->gtd7(doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t);
    }

    void CNrlmsise00::gtd7d(const int doy, const double sec, const double& alt,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t)
    {
        p_detail->gtd7d(doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, d, t);
    }

    double CNrlmsise00::density(const int doy, const double sec,
                  const double& alt, const double& g_lat, const double& g_long,
                  const double f107A, const double f107, std::array<double,7>& ap)
    {
        std::array<double,9> dens;
        std::array<double,2> temp;
        double lst = sec/3600.0 + g_long/15.0;
        p_detail->gtd7d(doy, sec, alt, g_lat, g_long, lst, f107A, f107, ap, dens, temp);

        return dens.at(5);
    }
}