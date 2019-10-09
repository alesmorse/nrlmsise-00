#include "Nrlmsise00_math.hpp"
#include <stdexcept>

namespace atmos
{
    namespace math
    {
        double splini(const double *xa, const double *ya, const double *y2a, const int n, const double& x)
        {
            double yi=0;
            int klo(0), khi(1);
            double xx, h, a, b, a2, b2;
            while ((x>xa[klo]) and (khi<n))
            {
                xx=x;
                if (khi<(n-1))
                {
                    if (x<xa[khi])
                    {
                        xx=x;
                    }
                    else
                    {
                        xx=xa[khi];
                    }
                }
                h = xa[khi] - xa[klo];
                a = (xa[khi] - xx)/h;
                b = (xx - xa[klo])/h;
                a2 = a*a;
                b2 = b*b;
                yi += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 + ((-(1.0+a2*a2)/4.0 + a2/2.0) * y2a[klo] + (b2*b2/4.0 - b2/2.0) * y2a[khi]) * h * h / 6.0) * h;
                klo++;
                khi++;
            }
            return yi;
        }
        double splint(const double *xa, const double *ya, const double *y2a, const int n, const double& x)
        {
            int klo(0), khi=(n-1);
            int k;
            double h;
            double a, b;
            while ((khi-klo)>1)
            {
                k=(khi+klo)/2;
                if (xa[k]>x)
                {
                    khi=k;
                }
                else
                {
                    klo=k;
                }
            }
            h = xa[khi] - xa[klo];
            if (h==0.0)
            {
                throw std::runtime_error("CNrlmsise00_p::splint: bad XA input to splint");
            }
            a = (xa[khi] - x)/h;
            b = (x - xa[klo])/h;
            return (a * ya[klo] + b * ya[khi] + ((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h * h/6.0);
        }
        void spline(const double *x, const double *y, const int n, const double& yp1, const double& ypn, double *y2)
        {
            double u[n];
            double sig, p, qn, un;
            int i, k;
            // Values > 1E30 signals second derivative zero
            if (yp1>0.99E30)
            {
                y2[0]=0;
                u[0]=0;
            }
            else
            {
                y2[0]=-0.5;
                u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
            }
            for (i=1;i<(n-1);i++)
            {
                sig = (x[i]-x[i-1])/(x[i+1] - x[i-1]);
                p = sig * y2[i-1] + 2.0;
                y2[i] = (sig - 1.0) / p;
                u[i] = (6.0 * ((y[i+1] - y[i])/(x[i+1] - x[i]) -(y[i] - y[i-1]) / (x[i] - x[i-1]))/(x[i+1] - x[i-1]) - sig * u[i-1])/p;
            }
            // Values > 1E30 signals second derivative zero
            if (ypn>0.99E30)
            {
                qn = 0;
                un = 0;
            }
            else
            {
                qn = 0.5;
                un = (3.0 / (x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
            }
            y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
            for (k=n-2;k>=0;k--)
            {
                y2[k] = y2[k] * y2[k+1] + u[k];
            }
        }
    } // namespace math
} // namespace atmos
