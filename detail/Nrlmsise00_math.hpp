#ifndef math_NRLMSISE00_MATH_H
#define math_NRLMSISE00_MATH_H

#include <cmath>

namespace atmos
{
    /**
     * @brief Namespace that contains mathematical utilities
     *
     */
    namespace math
    {
        /**
         * @brief Integrate cubic spline function from xa(1) to x
         *
         * @param xa Array of interpolation nodes
         * @param ya Arrays of tabulated function in ascending order by xa with ya = f(xa)
         * @param y2a Array of second derivatives
         * @param n Array size
         * @param x Ascissa endpoint of integration
         * @return double integrated value
         */
        double splini(const double *xa, const double *ya, const double *y2a, const int n, const double& x);
        /**
         * @brief Calculate cubic spline interp value
         *
         * @param xa Array of interpolation nodes
         * @param ya Arrays of tabulated function values in ascending xa order
         * @param y2a Arrays of second derivatives
         * @param n Array size
         * @param x Abscissa of interpolation
         * @return double cubic spline interp value
         */
        double splint(const double *xa, const double *ya, const double *y2a, const int n, const double& x);
        /**
         * @brief Calculate 2nd derivatives of cubic spline interp function
         *
         * @param x Array of interpolation nodes
         * @param y Arrays of tabulated function in ascending order by x with y = f(x)
         * @param n Array size
         * @param yp1 Specified derivatives at x(1)
         * @param ypn Specified derivatives at x(n)
         * @param y2 Output array of second derivatives
         */
        void spline(const double *x, const double *y, const int n, const double& yp1, const double& ypn, double *y2);

    } // namespace math
} // namespace atmos

#endif // math_NRLMSISE00_MATH_H