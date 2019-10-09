#include "Nrlmsise00_p.hpp"
#include "Nrlmsise00_math.hpp"
#include "Nrlmsise00_data.cpp"
#include <iostream>

// Implementation
namespace atmos
{
    CNrlmsise00_p::CNrlmsise00_p(const std::array<int, 24>& switches)
    {
        for (size_t i = 0; i < 24; i++)
        {
            a_switches.at(i) = switches.at(i);
            if(i!=9)
            {
                a_sw.at(i) = (switches.at(i)==1) ? 1.0 : 0.0;
                a_swc.at(i) = (switches.at(i)>0) ? 1.0 : 0.0;
            }
            else
            {
                a_sw.at(i) = switches.at(i);
                a_swc.at(i) = switches.at(i);
            }
        }
    }

    void CNrlmsise00_p::glatf(const double& lat, double& gv, double& reff)
    {
        double c2 = cos(2.0*data::d_DGTR*lat);
        gv = 980.616*(1.0-0.0026373*c2);
        reff = 2.0*gv/(3.085462e-6+2.27e-9*c2)*1e-5;
    }

    double CNrlmsise00_p::ccor(const double& alt, const double& r, const double h1, const double zh)
    {
        double e = (alt-zh)/h1;
        if(e>70.0)
        {
            return exp(0);
        }
        if(e<-70.0)
        {
            return exp(r);
        }
        return exp(r/(1.0+exp(e)));
    }


    double CNrlmsise00_p::ccor2(const double& alt, const double& r, const double h1, const double zh, const double h2)
    {
        double e1 = (alt - zh) / h1;
        double e2 = (alt - zh) / h2;
        if ((e1 > 70.0) or (e2 > 70.0))
        {
            return exp(0);
        }
        if ((e1 < -70.0) and (e2 < -70.0))
        {
            return exp(r);
        }
        return exp(r/(1.0 + 0.5*(exp(e1) + exp(e2))));
    }


    double CNrlmsise00_p::scalh(const double& alt, const double xm, const double temp) const
    {
        double g = d_gsurf/(pow((1.0+alt/d_re),2.0));
        return data::d_RGAS*temp/(g*xm);
    }


    double CNrlmsise00_p::dnet(double& dd, const double& dm, const double& zhm, const double& xmm, const double xm) const
    {
        double a;
        double ylog;
        a  = zhm / (xmm-xm);
        if (not ((dm>0) and (dd>0)))
        {
            std::cout << "dnet log: error dm: " << dm << ", dd: " << dd << ", xm: " << xm << std::endl;
            if ((dd==0) and (dm==0))
            {
                dd=1;
            }
            if (dm==0)
            {
                return dd;
            }
            if (dd==0)
            {
                return dm;
            }
        }
        ylog = a * log(dm/dd);
        if (ylog<-10)
        {
            return dd;
        }
        if (ylog>10)
        {
            return dm;
        }
        a = dd*pow((1.0 + exp(ylog)),(1.0/a));
        return a;
    }



    inline double CNrlmsise00_p::zeta(const double& zz, double zl) const
    {
        return ((zz-zl)*(d_re+zl)/(d_re+zz));
    }


    double CNrlmsise00_p::densm(const double& alt, const double& d0, const double xm, double *tz,
                                 const int mn3, const double *zn3, const double *tn3, const double *tgn3,
                                 const int mn2, const double *zn2, const double *tn2, const double *tgn2) const
    {
        std::array<double,10> xs, ys, y2out;
        double z;
        double z1, z2;
        double t1, t2, zg, zgdif;
        double yd1, yd2;
        double x, y, yi;
        double expl, gamm, glb;
        double densm_tmp;
        int mn;
        int k;
        densm_tmp=d0;
        if (alt>zn2[0])
        {
            if (xm==0.0)
            {
                return *tz;
            }
            else
            {
                return d0;
            }
        }

        /* STRATOSPHERE/MESOSPHERE TEMPERATURE */
        if (alt>zn2[mn2-1])
        {
            z=alt;
        }
        else
        {
            z=zn2[mn2-1];
        }
        mn=mn2;
        z1=zn2[0];
        z2=zn2[mn-1];
        t1=tn2[0];
        t2=tn2[mn-1];
        zg = zeta(z, z1);
        zgdif = zeta(z2, z1);

        /* set up spline nodes */
        for (k=0;k<mn;k++)
        {
            xs[k]=zeta(zn2[k],z1)/zgdif;
            ys[k]=1.0 / tn2[k];
        }
        yd1=-tgn2[0] / (t1*t1) * zgdif;
        yd2=-tgn2[1] / (t2*t2) * zgdif * (pow(((d_re+z2)/(d_re+z1)),2.0));

        /* calculate spline coefficients */
        math::spline(xs.data(), ys.data(), mn, yd1, yd2, y2out.data());
        x = zg/zgdif;
        y = math::splint(xs.data(), ys.data(), y2out.data(), mn, x);

        /* temperature at altitude */
        *tz = 1.0 / y;
        if (xm!=0.0)
        {
            /* calaculate stratosphere / mesospehere density */
            glb = d_gsurf / (pow((1.0 + z1/d_re),2.0));
            gamm = xm * glb * zgdif / data::d_RGAS;

            /* Integrate temperature profile */
            yi = math::splini(xs.data(), ys.data(), y2out.data(), mn, x);
            expl=gamm*yi;
            if (expl>50.0)
            {
                expl=50.0;
            }

            /* Density at altitude */
            densm_tmp = densm_tmp * (t1 / *tz) * exp(-expl);
        }

        if (alt>zn3[0])
        {
            if (xm==0.0)
            {
                return *tz;
            }
            else
            {
                return densm_tmp;
            }
        }

        /* troposhere / stratosphere temperature */
        z = alt;
        mn = mn3;
        z1=zn3[0];
        z2=zn3[mn-1];
        t1=tn3[0];
        t2=tn3[mn-1];
        zg=zeta(z,z1);
        zgdif=zeta(z2,z1);

        /* set up spline nodes */
        for (k=0;k<mn;k++)
        {
            xs[k] = zeta(zn3[k],z1) / zgdif;
            ys[k] = 1.0 / tn3[k];
        }
        yd1=-tgn3[0] / (t1*t1) * zgdif;
        yd2=-tgn3[1] / (t2*t2) * zgdif * (pow(((d_re+z2)/(d_re+z1)),2.0));

        /* calculate spline coefficients */
        math::spline(xs.data(), ys.data(), mn, yd1, yd2, y2out.data());
        x = zg/zgdif;
        y = math::splint(xs.data(), ys.data(), y2out.data(), mn, x);

        /* temperature at altitude */
        *tz = 1.0 / y;
        if (xm!=0.0)
        {
            /* calaculate tropospheric / stratosphere density */
            glb = d_gsurf / (pow((1.0 + z1/d_re),2.0));
            gamm = xm * glb * zgdif / data::d_RGAS;

            /* Integrate temperature profile */
            yi = math::splini(xs.data(), ys.data(), y2out.data(), mn, x);
            expl=gamm*yi;
            if (expl>50.0)
            {
                expl=50.0;
            }

            /* Density at altitude */
            densm_tmp = densm_tmp * (t1 / *tz) * exp(-expl);
        }
        if (xm==0.0)
        {
            return *tz;
        }
        else
        {
            return densm_tmp;
        }
    }


    double CNrlmsise00_p::densu(const double& alt, const double& dlb, const double& tinf, const double& tlb, const double& xm,
                                 const double alpha, double *tz, const double zlb, const double& s2,
                                 const int mn1, const double *zn1, double *tn1, double *tgn1) const
    {
        double yd2, yd1, x=0, y;
        double densu_temp=1.0;
        double za, z1=0.0, z2;
        double z, zg2, tt, ta;
        double dta, t1=0.0, t2, zg, zgdif=0.0;
        int mn=0;
        int k;
        double glb;
        double expl;
        double yi;
        double densa;
        double gamma, gamm;
        std::array<double,5> xs, ys, y2out;
        /* joining altitudes of Bates and spline */
        za=zn1[0];
        if (alt>za)
        {
            z=alt;
        }
        else
        {
            z=za;
        }

        /* geopotential altitude difference from ZLB */
        zg2 = zeta(z, zlb);

        /* Bates temperature */
        tt = tinf - (tinf - tlb) * exp(-s2*zg2);
        ta = tt;
        *tz = tt;
        densu_temp = *tz;

        if (alt<za)
        {
            /* calculate temperature below ZA
             * temperature gradient at ZA from Bates profile */
            dta = (tinf - ta) * s2 * pow(((d_re+zlb)/(d_re+za)),2.0);
            tgn1[0]=dta;
            tn1[0]=ta;
            if (alt>zn1[mn1-1])
            {
                z=alt;
            }
            else
            {
                z=zn1[mn1-1];
            }
            mn=mn1;
            z1=zn1[0];
            z2=zn1[mn-1];
            t1=tn1[0];
            t2=tn1[mn-1];
            /* geopotental difference from z1 */
            zg = zeta(z, z1);
            zgdif = zeta(z2, z1);
            /* set up spline nodes */
            for (k=0;k<mn;k++)
            {
                xs[k] = zeta(zn1[k], z1) / zgdif;
                ys[k] = 1.0 / tn1[k];
            }
            /* end node derivatives */
            yd1 = -tgn1[0] / (t1*t1) * zgdif;
            yd2 = -tgn1[1] / (t2*t2) * zgdif * pow(((d_re+z2)/(d_re+z1)),2.0);
            /* calculate spline coefficients */
            math::spline(xs.data(), ys.data(), mn, yd1, yd2, y2out.data());
            x = zg / zgdif;
            y = math::splint(xs.data(), ys.data(), y2out.data(), mn, x);
            /* temperature at altitude */
            *tz = 1.0 / y;
            densu_temp = *tz;
        }
        if (xm==0)
        {
            return densu_temp;
        }

        /* calculate density above za */
        glb = d_gsurf / pow((1.0 + zlb/d_re),2.0);
        gamma = xm * glb / (s2 * data::d_RGAS * tinf);
        expl = exp(-s2 * gamma * zg2);
        if (expl>50.0)
        {
              expl=50.0;
        }
        if (tt<=0)
        {
            expl=50.0;
        }

        /* density at altitude */
        densa = dlb * pow((tlb/tt),((1.0+alpha+gamma))) * expl;
        densu_temp=densa;
        if (alt>=za)
        {
            return densu_temp;
        }

        /* calculate density below za */
        glb = d_gsurf / pow((1.0 + z1/d_re),2.0);
        gamm = xm * glb * zgdif / data::d_RGAS;

        /* integrate spline temperatures */
        yi = math::splini(xs.data(), ys.data(), y2out.data(), mn, x);
        expl = gamm * yi;
        if (expl>50.0)
        {
            expl=50.0;
        }
        if (*tz<=0)
        {
            expl=50.0;
        }

        /* density at altitude */
        densu_temp = densu_temp * pow((t1 / *tz),(1.0 + alpha)) * exp(-expl);
        return densu_temp;
    }


    // 3hr Magnetic activity functions
    // Eq. A24c
    static inline double sumex(const double& ex)
    {
        return (1.0 + (1.0 - pow(ex,19.0)) / (1.0 - ex) * pow(ex,0.5));
    }

    // 3hr Magnetic activity functions
    // Eq. A24d
    static inline double g0(const double a, const std::array<double,150>& p)
    {
        return (a - 4.0 + (p[25] - 1.0) * (a - 4.0 + (exp(-sqrt(p[24]*p[24]) * (a - 4.0)) - 1.0) / sqrt(p[24]*p[24])));
    }

    // 3hr Magnetic activity functions
    // Eq. A24a
    static inline double sg0(const double ex, const std::array<double,150>& p, std::array<double,7>& ap)
    {
        return (g0(ap[1],p) + (g0(ap[2],p)*ex + g0(ap[3],p)*ex*ex + \
                g0(ap[4],p)*pow(ex,3.0)    + (g0(ap[5],p)*pow(ex,4.0) + \
                g0(ap[6],p)*pow(ex,12.0))*(1.0-pow(ex,8.0))/(1.0-ex)))/sumex(ex);
    }


    double CNrlmsise00_p::globe7(const std::array<double,150>& p, const int doy, const double sec,
                 const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                 std::array<double,7>& ap)
    {
        double t[15];
        int i,j;
        double apd;
        double tloc;
        double c, s, c2, c4, s2;
        double cd32, cd18, cd14, cd39;
        double df;
        double f1, f2;
        double tinf;

        tloc=lst;
        for (j=0;j<14;j++)
        {
            t[j]=0;
        }

        /* calculate legendre polynomials */
        c = sin(g_lat * data::d_DGTR);
        s = cos(g_lat * data::d_DGTR);
        c2 = c*c;
        c4 = c2*c2;
        s2 = s*s;

        a_plg[0][1] = c;
        a_plg[0][2] = 0.5*(3.0*c2 -1.0);
        a_plg[0][3] = 0.5*(5.0*c*c2-3.0*c);
        a_plg[0][4] = (35.0*c4 - 30.0*c2 + 3.0)/8.0;
        a_plg[0][5] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c)/8.0;
        a_plg[0][6] = (11.0*c*a_plg[0][5] - 5.0*a_plg[0][4])/6.0;
        //a_plg[0][7] = (13.0*c*a_plg[0][6] - 6.0*a_plg[0][5])/7.0;
        a_plg[1][1] = s;
        a_plg[1][2] = 3.0*c*s;
        a_plg[1][3] = 1.5*(5.0*c2-1.0)*s;
        a_plg[1][4] = 2.5*(7.0*c2*c-3.0*c)*s;
        a_plg[1][5] = 1.875*(21.0*c4 - 14.0*c2 +1.0)*s;
        a_plg[1][6] = (11.0*c*a_plg[1][5]-6.0*a_plg[1][4])/5.0;
        //a_plg[1][7] = (13.0*c*a_plg[1][6]-7.0*a_plg[1][5])/6.0;
        //a_plg[1][8] = (15.0*c*a_plg[1][7]-8.0*a_plg[1][6])/7.0;
        a_plg[2][2] = 3.0*s2;
        a_plg[2][3] = 15.0*s2*c;
        a_plg[2][4] = 7.5*(7.0*c2 -1.0)*s2;
        a_plg[2][5] = 3.0*c*a_plg[2][4]-2.0*a_plg[2][3];
        a_plg[2][6] =(11.0*c*a_plg[2][5]-7.0*a_plg[2][4])/4.0;
        a_plg[2][7] =(13.0*c*a_plg[2][6]-8.0*a_plg[2][5])/5.0;
        a_plg[3][3] = 15.0*s2*s;
        a_plg[3][4] = 105.0*s2*s*c;
        a_plg[3][5] =(9.0*c*a_plg[3][4]-7.*a_plg[3][3])/2.0;
        a_plg[3][6] =(11.0*c*a_plg[3][5]-8.*a_plg[3][4])/3.0;

        if (!(((a_sw[7]==0)&&(a_sw[8]==0))&&(a_sw[14]==0)))
        {
            d_stloc = sin(data::d_HR*tloc);
            d_ctloc = cos(data::d_HR*tloc);
            d_s2tloc = sin(2.0*data::d_HR*tloc);
            d_c2tloc = cos(2.0*data::d_HR*tloc);
            d_s3tloc = sin(3.0*data::d_HR*tloc);
            d_c3tloc = cos(3.0*data::d_HR*tloc);
        }

        cd32 = cos(data::d_DR*(doy-p[31]));
        cd18 = cos(2.0*data::d_DR*(doy-p[17]));
        cd14 = cos(data::d_DR*(doy-p[13]));
        cd39 = cos(2.0*data::d_DR*(doy-p[38]));

        /* F10.7 EFFECT */
        df = f107 - f107A;
        d_dfa = f107A - 150.0;
        t[0] =  p[19]*df*(1.0+p[59]*d_dfa) + p[20]*df*df + p[21]*d_dfa + p[29]*pow(d_dfa,2.0);
        f1 = 1.0 + (p[47]*d_dfa +p[19]*df+p[20]*df*df)*a_swc[1];
        f2 = 1.0 + (p[49]*d_dfa+p[19]*df+p[20]*df*df)*a_swc[1];

        /*  TIME INDEPENDENT */
        t[1] = (p[1]*a_plg[0][2]+ p[2]*a_plg[0][4]+p[22]*a_plg[0][6]) + \
              (p[14]*a_plg[0][2])*d_dfa*a_swc[1] +p[26]*a_plg[0][1];

        /*  SYMMETRICAL ANNUAL */
        t[2] = p[18]*cd32;

        /*  SYMMETRICAL SEMIANNUAL */
        t[3] = (p[15]+p[16]*a_plg[0][2])*cd18;

        /*  ASYMMETRICAL ANNUAL */
        t[4] =  f1*(p[9]*a_plg[0][1]+p[10]*a_plg[0][3])*cd14;

        /*  ASYMMETRICAL SEMIANNUAL */
        t[5] =    p[37]*a_plg[0][1]*cd39;

        /* DIURNAL */
        if (a_sw[7])
        {
            double t71, t72;
            t71 = (p[11]*a_plg[1][2])*cd14*a_swc[5];
            t72 = (p[12]*a_plg[1][2])*cd14*a_swc[5];
            t[6] = f2*((p[3]*a_plg[1][1] + p[4]*a_plg[1][3] + p[27]*a_plg[1][5] + t71) * \
                   d_ctloc + (p[6]*a_plg[1][1] + p[7]*a_plg[1][3] + p[28]*a_plg[1][5] \
                        + t72)*d_stloc);
        }

        /* SEMIDIURNAL */
        if (a_sw[8])
        {
            double t81, t82;
            t81 = (p[23]*a_plg[2][3]+p[35]*a_plg[2][5])*cd14*a_swc[5];
            t82 = (p[33]*a_plg[2][3]+p[36]*a_plg[2][5])*cd14*a_swc[5];
            t[7] = f2*((p[5]*a_plg[2][2]+ p[41]*a_plg[2][4] + t81)*d_c2tloc +(p[8]*a_plg[2][2] + p[42]*a_plg[2][4] + t82)*d_s2tloc);
        }

        /* TERDIURNAL */
        if (a_sw[14])
        {
            t[13] = f2 * ((p[39]*a_plg[3][3]+(p[93]*a_plg[3][4]+p[46]*a_plg[3][6])*cd14*a_swc[5])* d_s3tloc +(p[40]*a_plg[3][3]+(p[94]*a_plg[3][4]+p[48]*a_plg[3][6])*cd14*a_swc[5])* d_c3tloc);
        }

        /* magnetic activity based on daily ap */
        if (a_sw[9]==-1)
        {
            if (p[51]!=0)
            {
                double exp1;
                exp1 = exp(-10800.0*sqrt(p[51]*p[51])/(1.0+p[138]*(45.0-sqrt(g_lat*g_lat))));
                if (exp1>0.99999)
                {
                    exp1=0.99999;
                }
//                if (p[24]<1.0E-4)
//                    p[24]=1.0E-4;
                a_apt[0]=sg0(exp1,p,ap);
                /* apt[1]=sg2(exp1,p,ap);
                   apt[2]=sg0(exp2,p,ap);
                   apt[3]=sg2(exp2,p,ap);
                */
                if (a_sw[9])
                {
                    t[8] = a_apt[0]*(p[50]+p[96]*a_plg[0][2]+p[54]*a_plg[0][4]+ \
                            (p[125]*a_plg[0][1]+p[126]*a_plg[0][3]+p[127]*a_plg[0][5])*cd14*a_swc[5]+ \
                            (p[128]*a_plg[1][1]+p[129]*a_plg[1][3]+p[130]*a_plg[1][5])*a_swc[7]* \
                               cos(data::d_HR*(tloc-p[131])));
                }
            }
        }
        else
        {
            double p44, p45;
            apd=ap[0]-4.0;
            p44=p[43];
            p45=p[44];
            if (p44<0)
            {
                p44 = 1.0E-5;
            }
            d_apdf = apd + (p45-1.0)*(apd + (exp(-p44 * apd) - 1.0)/p44);
            if (a_sw[9])
            {
                t[8]=d_apdf*(p[32]+p[45]*a_plg[0][2]+p[34]*a_plg[0][4]+ \
                        (p[100]*a_plg[0][1]+p[101]*a_plg[0][3]+p[102]*a_plg[0][5])*cd14*a_swc[5]+
                        (p[121]*a_plg[1][1]+p[122]*a_plg[1][3]+p[123]*a_plg[1][5])*a_swc[7]*
                        cos(data::d_HR*(tloc-p[124])));
            }
        }

        if ((a_sw[10])&&(g_long>-1000.0))
        {
            /* longitudinal */
            if (a_sw[11])
            {
                t[10] = (1.0 + p[80]*d_dfa*a_swc[1])* \
                        ((p[64]*a_plg[1][2]+p[65]*a_plg[1][4]+p[66]*a_plg[1][6]\
                         +p[103]*a_plg[1][1]+p[104]*a_plg[1][3]+p[105]*a_plg[1][5]\
                         +a_swc[5]*(p[109]*a_plg[1][1]+p[110]*a_plg[1][3]+p[111]*a_plg[1][5])*cd14)* \
                             cos(data::d_DGTR*g_long) \
                         +(p[90]*a_plg[1][2]+p[91]*a_plg[1][4]+p[92]*a_plg[1][6]\
                         +p[106]*a_plg[1][1]+p[107]*a_plg[1][3]+p[108]*a_plg[1][5]\
                         +a_swc[5]*(p[112]*a_plg[1][1]+p[113]*a_plg[1][3]+p[114]*a_plg[1][5])*cd14)* \
                         sin(data::d_DGTR*g_long));
            }

            /* ut and mixed ut, longitude */
            if (a_sw[12])
            {
                t[11]=(1.0+p[95]*a_plg[0][1])*(1.0+p[81]*d_dfa*a_swc[1])*\
                    (1.0+p[119]*a_plg[0][1]*a_swc[5]*cd14)*\
                    ((p[68]*a_plg[0][1]+p[69]*a_plg[0][3]+p[70]*a_plg[0][5])*\
                    cos(data::d_SR*(sec-p[71])));
                t[11]+=a_swc[11]*\
                    (p[76]*a_plg[2][3]+p[77]*a_plg[2][5]+p[78]*a_plg[2][7])*\
                    cos(data::d_SR*(sec-p[79])+2.0*data::d_DGTR*g_long)*(1.0+p[137]*d_dfa*a_swc[1]);
            }

            /* ut, longitude magnetic activity */
            if (a_sw[13])
            {
                if (a_sw[9]==-1)
                {
                    if (p[51])
                    {
                        t[12]=a_apt[0]*a_swc[11]*(1.+p[132]*a_plg[0][1])*\
                            ((p[52]*a_plg[1][2]+p[98]*a_plg[1][4]+p[67]*a_plg[1][6])*\
                             cos(data::d_DGTR*(g_long-p[97])))\
                            +a_apt[0]*a_swc[11]*a_swc[5]*\
                            (p[133]*a_plg[1][1]+p[134]*a_plg[1][3]+p[135]*a_plg[1][5])*\
                            cd14*cos(data::d_DGTR*(g_long-p[136])) \
                            +a_apt[0]*a_swc[12]* \
                            (p[55]*a_plg[0][1]+p[56]*a_plg[0][3]+p[57]*a_plg[0][5])*\
                            cos(data::d_SR*(sec-p[58]));
                    }
                }
                else
                {
                    t[12] = d_apdf*a_swc[11]*(1.0+p[120]*a_plg[0][1])*\
                        ((p[60]*a_plg[1][2]+p[61]*a_plg[1][4]+p[62]*a_plg[1][6])*\
                        cos(data::d_DGTR*(g_long-p[63])))\
                        +d_apdf*a_swc[11]*a_swc[5]* \
                        (p[115]*a_plg[1][1]+p[116]*a_plg[1][3]+p[117]*a_plg[1][5])* \
                        cd14*cos(data::d_DGTR*(g_long-p[118])) \
                        + d_apdf*a_swc[12]* \
                        (p[83]*a_plg[0][1]+p[84]*a_plg[0][3]+p[85]*a_plg[0][5])* \
                        cos(data::d_SR*(sec-p[75]));
                }
            }
        }

        /* parms not used: 82, 89, 99, 139-149 */
        tinf = p[30];
        for (i=0;i<14;i++)
        {
            tinf = tinf + fabs(a_sw[i+1])*t[i];
        }
        return tinf;
    }


    double CNrlmsise00_p::glob7s(const std::array<double,100>& p, const int doy, const double& g_long)
    {
        double t[14];
        double tt;
        double cd32, cd18, cd14, cd39;
        int i,j;

        for (j=0;j<14;j++)
        {
            t[j]=0.0;
        }

        cd32 = cos(data::d_DR*(doy-p[31]));
        cd18 = cos(2.0*data::d_DR*(doy-p[17]));
        cd14 = cos(data::d_DR*(doy-p[13]));
        cd39 = cos(2.0*data::d_DR*(doy-p[38]));

        /* F10.7 */
        t[0] = p[21]*d_dfa;

        /* time independent */
        t[1]=p[1]*a_plg[0][2] + p[2]*a_plg[0][4] + p[22]*a_plg[0][6] + p[26]*a_plg[0][1] + p[14]*a_plg[0][3] + p[59]*a_plg[0][5];

        /* SYMMETRICAL ANNUAL */
        t[2]=(p[18]+p[47]*a_plg[0][2]+p[29]*a_plg[0][4])*cd32;

        /* SYMMETRICAL SEMIANNUAL */
        t[3]=(p[15]+p[16]*a_plg[0][2]+p[30]*a_plg[0][4])*cd18;

        /* ASYMMETRICAL ANNUAL */
        t[4]=(p[9]*a_plg[0][1]+p[10]*a_plg[0][3]+p[20]*a_plg[0][5])*cd14;

        /* ASYMMETRICAL SEMIANNUAL */
        t[5]=(p[37]*a_plg[0][1])*cd39;

        /* DIURNAL */
        if (a_sw[7])
        {
            double t71, t72;
            t71 = p[11]*a_plg[1][2]*cd14*a_swc[5];
            t72 = p[12]*a_plg[1][2]*cd14*a_swc[5];
            t[6] = ((p[3]*a_plg[1][1] + p[4]*a_plg[1][3] + t71) * d_ctloc + (p[6]*a_plg[1][1] + p[7]*a_plg[1][3] + t72) * d_stloc) ;
        }

        /* SEMIDIURNAL */
        if (a_sw[8])
        {
            double t81, t82;
            t81 = (p[23]*a_plg[2][3]+p[35]*a_plg[2][5])*cd14*a_swc[5];
            t82 = (p[33]*a_plg[2][3]+p[36]*a_plg[2][5])*cd14*a_swc[5];
            t[7] = ((p[5]*a_plg[2][2] + p[41]*a_plg[2][4] + t81) * d_c2tloc + (p[8]*a_plg[2][2] + p[42]*a_plg[2][4] + t82) * d_s2tloc);
        }

        /* TERDIURNAL */
        if (a_sw[14])
        {
            t[13] = p[39] * a_plg[3][3] * d_s3tloc + p[40] * a_plg[3][3] * d_c3tloc;
        }

        /* MAGNETIC ACTIVITY */
        if (a_sw[9])
        {
            if (a_sw[9]==1)
            {
                t[8] = d_apdf * (p[32] + p[45] * a_plg[0][2] * a_swc[2]);
            }
            if (a_sw[9]==-1)
            {
                t[8]=(p[50]*a_apt[0] + p[96]*a_plg[0][2] * a_apt[0]*a_swc[2]);
            }
        }

        /* LONGITUDINAL */
        if (!((a_sw[10]==0) || (a_sw[11]==0) || (g_long<=-1000.0)))
        {
            t[10] = (1.0 + a_plg[0][1]*(p[80]*a_swc[5]*cos(data::d_DR*(doy-p[81]))\
                    +p[85]*a_swc[6]*cos(2.0*data::d_DR*(doy-p[86])))\
                +p[83]*a_swc[3]*cos(data::d_DR*(doy-p[84]))\
                +p[87]*a_swc[4]*cos(2.0*data::d_DR*(doy-p[88])))\
                *((p[64]*a_plg[1][2]+p[65]*a_plg[1][4]+p[66]*a_plg[1][6]\
                +p[74]*a_plg[1][1]+p[75]*a_plg[1][3]+p[76]*a_plg[1][5]\
                )*cos(data::d_DGTR*g_long)\
                +(p[90]*a_plg[1][2]+p[91]*a_plg[1][4]+p[92]*a_plg[1][6]\
                +p[77]*a_plg[1][1]+p[78]*a_plg[1][3]+p[79]*a_plg[1][5]\
                )*sin(data::d_DGTR*g_long));
        }
        tt=0;
        for (i=0;i<14;i++)
        {
            tt+=fabs(a_sw[i+1])*t[i];
        }
        return tt;
    }


    void CNrlmsise00_p::gts7(const int doy, const double sec, const double& alt,
                const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t)
    {
        double za;
        int i, j;
        double z;
        std::array<double,5> zn1 = {120.0, 110.0, 100.0, 90.0, 72.5};
        double tinf;
        int mn1 = 5;
        double g0;
        double tlb;
        double s;
        double db01, db04, db14, db16, db28, db32, db40;
        double zh28, zh04, zh16, zh32, zh40, zh01, zh14;
        double zhm28, zhm04, zhm16, zhm32, zhm40, zhm01, zhm14;
        double xmd;
        double b28, b04, b16, b32, b40, b01, b14;
        double tz;
        double g28, g4, g16, g32, g40, g1, g14;
        double zhf, xmm;
        double zc04, zc16, zc32, zc40, zc01, zc14;
        double hc04, hc16, hc32, hc40, hc01, hc14;
        double hcc16, hcc32, hcc01, hcc14;
        double zcc16, zcc32, zcc01, zcc14;
        double rc16, rc32, rc01, rc14;
        double rl;
        double g16h, db16h, zsho;
        double tho, zsht, zmho;
        std::array<double,9> alpha = {-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0};
        std::array<double,8> altl = {200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0};
        double dd;
        double hc216, hcc232;
        za = data::pdl[1][15];
        zn1[0] = za;
        for (j=0;j<9;j++)
        {
            d[j]=0;
        }

        /* TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1) */
        if (alt>zn1[0])
        {
            tinf = data::ptm[0]*data::pt[0]*(1.0+a_sw[16]*globe7(data::pt,doy,sec,g_lat,g_long,lst,f107A,f107,ap));
        }
        else
        {
            tinf = data::ptm[0]*data::pt[0];
        }
        t[0]=tinf;

        /*  GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5) */
        if (alt>zn1[4])
        {
            g0 = data::ptm[3]*data::ps[0]*(1.0+a_sw[19]*globe7(data::ps,doy,sec,g_lat,g_long,lst,f107A,f107,ap));
        }
        else
        {
            g0 = data::ptm[3]*data::ps[0];
        }
        tlb = data::ptm[1] * (1.0 + a_sw[17]*globe7(data::pd[3],doy,sec,g_lat,g_long,lst,f107A,f107,ap))*data::pd[3][0];
        s = g0 / (tinf - tlb);

        /* Lower thermosphere temp variations not significant for
         *  density above 300 km */
        if (alt<300.0)
        {
            a_meso_tn1[1]=data::ptm[6]*data::ptl[0][0]/(1.0-a_sw[18]*glob7s(data::ptl[0],doy,g_long));
            a_meso_tn1[2]=data::ptm[2]*data::ptl[1][0]/(1.0-a_sw[18]*glob7s(data::ptl[1],doy,g_long));
            a_meso_tn1[3]=data::ptm[7]*data::ptl[2][0]/(1.0-a_sw[18]*glob7s(data::ptl[2],doy,g_long));
            a_meso_tn1[4]=data::ptm[4]*data::ptl[3][0]/(1.0-a_sw[18]*a_sw[20]*glob7s(data::ptl[3],doy,g_long));
            a_meso_tgn1[1]=data::ptm[8]*data::pma[8][0]*(1.0+a_sw[18]*a_sw[20]*glob7s(data::pma[8],doy,g_long))*a_meso_tn1[4]*a_meso_tn1[4]/(pow((data::ptm[4]*data::ptl[3][0]),2.0));
        }
        else
        {
            a_meso_tn1[1]=data::ptm[6]*data::ptl[0][0];
            a_meso_tn1[2]=data::ptm[2]*data::ptl[1][0];
            a_meso_tn1[3]=data::ptm[7]*data::ptl[2][0];
            a_meso_tn1[4]=data::ptm[4]*data::ptl[3][0];
            a_meso_tgn1[1]=data::ptm[8]*data::pma[8][0]*a_meso_tn1[4]*a_meso_tn1[4]/(pow((data::ptm[4]*data::ptl[3][0]),2.0));
        }

        /* N2 variation factor at Zlb */
        g28=a_sw[21]*globe7(data::pd[2],doy,sec,g_lat,g_long,lst,f107A,f107,ap);

        /* VARIATION OF TURBOPAUSE HEIGHT */
        zhf=data::pdl[1][24]*(1.0+a_sw[5]*data::pdl[0][24]*sin(data::d_DGTR*g_lat)*cos(data::d_DR*(doy-data::pt[13])));
        t[0]=tinf;
        xmm = data::pdm[2][4];
        z = alt;


        /**** N2 DENSITY ****/

        /* Diffusive density at Zlb */
        db28 = data::pdm[2][0]*exp(g28)*data::pd[2][0];
        /* Diffusive density at Alt */
        d[2]=densu(z,db28,tinf,tlb,28.0,alpha[2],&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[2];
        /* Turbopause */
        zh28=data::pdm[2][2]*zhf;
        zhm28=data::pdm[2][3]*data::pdl[1][5];
        xmd=28.0-xmm;
        /* Mixed density at Zlb */
        b28=densu(zh28,db28,tinf,tlb,xmd,(alpha[2]-1.0),&tz,data::ptm[5],s,mn1, zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        if ((a_sw[15])&&(z<=altl[2]))
        {
            /*  Mixed density at Alt */
            d_dm28=densu(z,b28,tinf,tlb,xmm,alpha[2],&tz,data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            /*  Net density at Alt */
            d[2]=dnet(d[2],d_dm28,zhm28,xmm,28.0);
        }


        /**** HE DENSITY ****/

        /*   Density variation factor at Zlb */
        g4 = a_sw[21]*globe7(data::pd[0],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
        /*  Diffusive density at Zlb */
        db04 = data::pdm[0][0]*exp(g4)*data::pd[0][0];
            /*  Diffusive density at Alt */
        d[0]=densu(z,db04,tinf,tlb, 4.,alpha[0],&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[0];
        if ((a_sw[15]) && (z<altl[0]))
        {
            /*  Turbopause */
            zh04=data::pdm[0][2];
            /*  Mixed density at Zlb */
            b04=densu(zh04,db04,tinf,tlb,4.0-xmm,alpha[0]-1.,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            /*  Mixed density at Alt */
            d_dm04=densu(z,b04,tinf,tlb,xmm,0.0,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            zhm04=zhm28;
            /*  Net density at Alt */
            d[0]=dnet(d[0],d_dm04,zhm04,xmm,4.0);
            /*  Correction to specified mixing ratio at ground */
            rl=log(b28*data::pdm[0][1]/b04);
            zc04=data::pdm[0][4]*data::pdl[1][0];
            hc04=data::pdm[0][5]*data::pdl[1][1];
            /*  Net density corrected at Alt */
            d[0]=d[0]*ccor(z,rl,hc04,zc04);
        }


        /**** O DENSITY ****/

        /*  Density variation factor at Zlb */
        g16= a_sw[21]*globe7(data::pd[1],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
        /*  Diffusive density at Zlb */
        db16 =  data::pdm[1][0]*exp(g16)*data::pd[1][0];
            /*   Diffusive density at Alt */
        d[1]=densu(z,db16,tinf,tlb, 16.,alpha[1],&t[1],data::ptm[5],s,mn1, zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[1];
        if ((a_sw[15]) && (z<=altl[1]))
        {
            /*   Turbopause */
            zh16=data::pdm[1][2];
            /*  Mixed density at Zlb */
            b16=densu(zh16,db16,tinf,tlb,16.0-xmm,(alpha[1]-1.0), &t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            /*  Mixed density at Alt */
            d_dm16=densu(z,b16,tinf,tlb,xmm,0.0,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            zhm16=zhm28;
            /*  Net density at Alt */
            d[1]=dnet(d[1],d_dm16,zhm16,xmm,16.0);
            rl=data::pdm[1][1]*data::pdl[1][16]*(1.0+a_sw[1]*data::pdl[0][23]*(f107A-150.0));
            hc16=data::pdm[1][5]*data::pdl[1][3];
            zc16=data::pdm[1][4]*data::pdl[1][2];
            hc216=data::pdm[1][5]*data::pdl[1][4];
            d[1]=d[1]*ccor2(z,rl,hc16,zc16,hc216);
            /*   Chemistry correction */
            hcc16=data::pdm[1][7]*data::pdl[1][13];
            zcc16=data::pdm[1][6]*data::pdl[1][12];
            rc16 =data::pdm[1][3]*data::pdl[1][14];
            /*  Net density corrected at Alt */
            d[1]=d[1]*ccor(z,rc16,hcc16,zcc16);
        }


        /**** O2 DENSITY ****/

        /*   Density variation factor at Zlb */
        g32= a_sw[21]*globe7(data::pd[4],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
        /*  Diffusive density at Zlb */
        db32 = data::pdm[3][0]*exp(g32)*data::pd[4][0];
        /*   Diffusive density at Alt */
        d[3]=densu(z,db32,tinf,tlb, 32.,alpha[3],&t[1],data::ptm[5],s,mn1, zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[3];
        if (a_sw[15])
        {
            if (z<=altl[3])
            {
                /* Turbopause */
                zh32=data::pdm[3][2];
                /* Mixed density at Zlb */
                b32=densu(zh32,db32,tinf,tlb,32.-xmm,alpha[3]-1., &t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
                /* Mixed density at Alt */
                d_dm32=densu(z,b32,tinf,tlb,xmm,0.,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
                zhm32=zhm28;
                /* Net density at Alt */
                d[3]=dnet(d[3],d_dm32,zhm32,xmm,32.0);
                /* Correction to specified mixing ratio at ground */
                rl=log(b28*data::pdm[3][1]/b32);
                hc32=data::pdm[3][5]*data::pdl[1][7];
                zc32=data::pdm[3][4]*data::pdl[1][6];
                d[3]=d[3]*ccor(z,rl,hc32,zc32);
            }
            /*  Correction for general departure from diffusive equilibrium above Zlb */
            hcc32 =data::pdm[3][7]*data::pdl[1][22];
            hcc232=data::pdm[3][7]*data::pdl[0][22];
            zcc32 =data::pdm[3][6]*data::pdl[1][21];
            rc32  =data::pdm[3][3]*data::pdl[1][23]*(1.+a_sw[1]*data::pdl[0][23]*(f107A-150.0));
            /*  Net density corrected at Alt */
            d[3]=d[3]*ccor2(z,rc32,hcc32,zcc32,hcc232);
        }


        /**** AR DENSITY ****/

        /*   Density variation factor at Zlb */
        g40= a_sw[21]*globe7(data::pd[5],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
        /*  Diffusive density at Zlb */
        db40 = data::pdm[4][0]*exp(g40)*data::pd[5][0];
        /*   Diffusive density at Alt */
        d[4]=densu(z,db40,tinf,tlb, 40.,alpha[4],&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[4];
        if ((a_sw[15]) && (z<=altl[4]))
        {
            /*   Turbopause */
            zh40=data::pdm[4][2];
            /*  Mixed density at Zlb */
            b40=densu(zh40,db40,tinf,tlb,40.0-xmm,alpha[4]-1.0,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            /*  Mixed density at Alt */
            d_dm40=densu(z,b40,tinf,tlb,xmm,0.0,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            zhm40=zhm28;
            /*  Net density at Alt */
            d[4]=dnet(d[4],d_dm40,zhm40,xmm,40.0);
            /*   Correction to specified mixing ratio at ground */
            rl=log(b28*data::pdm[4][1]/b40);
            hc40=data::pdm[4][5]*data::pdl[1][9];
            zc40=data::pdm[4][4]*data::pdl[1][8];
            /*  Net density corrected at Alt */
            d[4]=d[4]*ccor(z,rl,hc40,zc40);
        }


        /**** HYDROGEN DENSITY ****/

        /*   Density variation factor at Zlb */
        g1 = a_sw[21]*globe7(data::pd[6],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
        /*  Diffusive density at Zlb */
        db01 = data::pdm[5][0]*exp(g1)*data::pd[6][0];
        /*   Diffusive density at Alt */
        d[6]=densu(z,db01,tinf,tlb,1.,alpha[6],&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[6];
        if ((a_sw[15]) && (z<=altl[6]))
        {
            /*   Turbopause */
            zh01=data::pdm[5][2];
            /*  Mixed density at Zlb */
            b01=densu(zh01,db01,tinf,tlb,1.0-xmm,alpha[6]-1., &t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            /*  Mixed density at Alt */
            d_dm01=densu(z,b01,tinf,tlb,xmm,0.,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            zhm01=zhm28;
            /*  Net density at Alt */
            d[6]=dnet(d[6],d_dm01,zhm01,xmm,1.0);
            /*   Correction to specified mixing ratio at ground */
            rl=log(b28*data::pdm[5][1]*sqrt(data::pdl[1][17]*data::pdl[1][17])/b01);
            hc01=data::pdm[5][5]*data::pdl[1][11];
            zc01=data::pdm[5][4]*data::pdl[1][10];
            d[6]=d[6]*ccor(z,rl,hc01,zc01);
            /*   Chemistry correction */
            hcc01=data::pdm[5][7]*data::pdl[1][19];
            zcc01=data::pdm[5][6]*data::pdl[1][18];
            rc01 =data::pdm[5][3]*data::pdl[1][20];
            /*  Net density corrected at Alt */
            d[6]=d[6]*ccor(z,rc01,hcc01,zcc01);
        }


        /**** ATOMIC NITROGEN DENSITY ****/

        /*   Density variation factor at Zlb */
        g14 = a_sw[21]*globe7(data::pd[7],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
            /*  Diffusive density at Zlb */
        db14 = data::pdm[6][0]*exp(g14)*data::pd[7][0];
            /*   Diffusive density at Alt */
        d[7]=densu(z,db14,tinf,tlb,14.0,alpha[7],&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        dd=d[7];
        if ((a_sw[15]) && (z<=altl[7]))
        {
            /*   Turbopause */
            zh14=data::pdm[6][2];
            /*  Mixed density at Zlb */
            b14=densu(zh14,db14,tinf,tlb,14.0-xmm,alpha[7]-1., &t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            /*  Mixed density at Alt */
            d_dm14=densu(z,b14,tinf,tlb,xmm,0.0,&t[1],data::ptm[5],s,mn1,zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
            zhm14=zhm28;
            /*  Net density at Alt */
            d[7]=dnet(d[7],d_dm14,zhm14,xmm,14.0);
            /*   Correction to specified mixing ratio at ground */
            rl=log(b28*data::pdm[6][1]*sqrt(data::pdl[0][2]*data::pdl[0][2])/b14);
            hc14=data::pdm[6][5]*data::pdl[0][1];
            zc14=data::pdm[6][4]*data::pdl[0][0];
            d[7]=d[7]*ccor(z,rl,hc14,zc14);
            /*   Chemistry correction */
            hcc14=data::pdm[6][7]*data::pdl[0][4];
            zcc14=data::pdm[6][6]*data::pdl[0][3];
            rc14 =data::pdm[6][3]*data::pdl[0][5];
            /*  Net density corrected at Alt */
            d[7]=d[7]*ccor(z,rc14,hcc14,zcc14);
        }


        /**** Anomalous OXYGEN DENSITY ****/

        g16h = a_sw[21]*globe7(data::pd[8],doy,sec,g_lat,g_long,lst,f107A,f107,ap);
        db16h = data::pdm[7][0]*exp(g16h)*data::pd[8][0];
        tho = data::pdm[7][9]*data::pdl[0][6];
        dd=densu(z,db16h,tho,tho,16.,alpha[8],&t[1],data::ptm[5],s,mn1, zn1.data(),a_meso_tn1.data(),a_meso_tgn1.data());
        zsht=data::pdm[7][5];
        zmho=data::pdm[7][4];
        zsho=scalh(zmho,16.0,tho);
        d[8]=dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1.));


        /* total mass density */
        d[5] = 1.66E-24*(4.0*d[0]+16.0*d[1]+28.0*d[2]+32.0*d[3]+40.0*d[4]+ d[6]+14.0*d[7]);


        /* temperature */
        z = sqrt(alt*alt);
        densu(z,1.0, tinf, tlb, 0.0, 0.0, &t[1], data::ptm[5], s, mn1, zn1.data(), a_meso_tn1.data(), a_meso_tgn1.data());
        if (a_sw[0])
        {
            for(i=0;i<9;i++)
            {
                d[i]=d[i]*1.0E6;
            }
            d[5]=d[5]/1000;
        }
    }


    void CNrlmsise00_p::gtd7(const int doy, const double sec, const double& alt,
                const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t)
    {
        double xlat;
        double xmm;
        int mn3 = 5;
        std::array<double,5> zn3 = {32.5,20.0,15.0,10.0,0.0};
        int mn2 = 4;
        std::array<double,4> zn2 = {72.5,55.0,45.0,32.5};
        double altt;
        double zmix=62.5;
        double dm28m;
        double tz;
        double dmc;
        double dmr;
        double dz28;
        std::array<double,9> sdens;
        std::array<double,2> stemp;
        int i;

        /* Latitude variation of gravity (none for sw[2]=0) */
        xlat=g_lat;
        if (a_sw[2]==0)
            xlat=45.0;
        glatf(xlat, d_gsurf, d_re);

        xmm = data::pdm[2][4];

        /* THERMOSPHERE / MESOSPHERE (above zn2[0]) */
        if (alt>zn2[0])
        {
            altt=alt;
        }
        else
        {
            altt=zn2[0];
        }

        gts7(doy, sec, altt, g_lat, g_long, lst, f107A, f107, ap, sdens, stemp);
        altt=alt;
        if (a_sw[0])   /* metric adjustment */
        {
            dm28m=d_dm28*1.0E6;
        }
        else
        {
            dm28m=d_dm28;
        }
        t[0]=stemp[0];
        t[1]=stemp[1];
        if (alt>=zn2[0])
        {
            for (i=0;i<9;i++)
            {
                d[i]=sdens[i];
            }
            return;
        }

        /* LOWER MESOSPHERE/UPPER STRATOSPHERE (between zn3[0] and zn2[0])
         *   Temperature at nodes and gradients at end nodes
         *   Inverse temperature a linear function of spherical harmonics
         */
        a_meso_tgn2[0]=a_meso_tgn1[1];
        a_meso_tn2[0]=a_meso_tn1[4];
        a_meso_tn2[1]=data::pma[0][0]*data::pavgm[0]/(1.0-a_sw[20]*glob7s(data::pma[0],doy,g_long));
        a_meso_tn2[2]=data::pma[1][0]*data::pavgm[1]/(1.0-a_sw[20]*glob7s(data::pma[1],doy,g_long));
        a_meso_tn2[3]=data::pma[2][0]*data::pavgm[2]/(1.0-a_sw[20]*a_sw[22]*glob7s(data::pma[2],doy,g_long));
        a_meso_tgn2[1]=data::pavgm[8]*data::pma[9][0]*(1.0+a_sw[20]*a_sw[22]*glob7s(data::pma[9],doy,g_long))*a_meso_tn2[3]*a_meso_tn2[3]/(pow((data::pma[2][0]*data::pavgm[2]),2.0));
        a_meso_tn3[0]=a_meso_tn2[3];

        if (alt<=zn3[0])
        {
            /* LOWER STRATOSPHERE AND TROPOSPHERE (below zn3[0])
             *   Temperature at nodes and gradients at end nodes
             *   Inverse temperature a linear function of spherical harmonics
             */
            a_meso_tgn3[0]=a_meso_tgn2[1];
            a_meso_tn3[1]=data::pma[3][0]*data::pavgm[3]/(1.0-a_sw[22]*glob7s(data::pma[3],doy,g_long));
            a_meso_tn3[2]=data::pma[4][0]*data::pavgm[4]/(1.0-a_sw[22]*glob7s(data::pma[4],doy,g_long));
            a_meso_tn3[3]=data::pma[5][0]*data::pavgm[5]/(1.0-a_sw[22]*glob7s(data::pma[5],doy,g_long));
            a_meso_tn3[4]=data::pma[6][0]*data::pavgm[6]/(1.0-a_sw[22]*glob7s(data::pma[6],doy,g_long));
            a_meso_tgn3[1]=data::pma[7][0]*data::pavgm[7]*(1.0+a_sw[22]*glob7s(data::pma[7],doy,g_long))*a_meso_tn3[4]*a_meso_tn3[4]/(pow((data::pma[6][0]*data::pavgm[6]),2.0));
        }

        /* LINEAR TRANSITION TO FULL MIXING BELOW zn2[0] */

        dmc=0;
        if (alt>zmix)
        {
            dmc = 1.0 - (zn2[0]-alt)/(zn2[0] - zmix);
        }
        dz28=sdens[2];

        /**** N2 density ****/
        dmr=sdens[2] / dm28m - 1.0;
        d[2]=densm(alt,dm28m,xmm, &tz, mn3, zn3.data(), a_meso_tn3.data(), a_meso_tgn3.data(), mn2, zn2.data(), a_meso_tn2.data(), a_meso_tgn2.data());
        d[2]=d[2] * (1.0 + dmr*dmc);

        /**** HE density ****/
        dmr = sdens[0] / (dz28 * data::pdm[0][1]) - 1.0;
        d[0] = d[2] * data::pdm[0][1] * (1.0 + dmr*dmc);

        /**** O density ****/
        d[1] = 0;
        d[8] = 0;

        /**** O2 density ****/
        dmr = sdens[3] / (dz28 * data::pdm[3][1]) - 1.0;
        d[3] = d[2] * data::pdm[3][1] * (1.0 + dmr*dmc);

        /**** AR density ***/
        dmr = sdens[4] / (dz28 * data::pdm[4][1]) - 1.0;
        d[4] = d[2] * data::pdm[4][1] * (1.0 + dmr*dmc);

        /**** Hydrogen density ****/
        d[6] = 0;

        /**** Atomic nitrogen density ****/
        d[7] = 0;

        /**** Total mass density */
        d[5] = 1.66E-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7]);

        if (a_sw[0])
        {
            d[5]=d[5]/1000;
        }

        /**** temperature at altitude ****/
        d_dd = densm(alt, 1.0, 0, &tz, mn3, zn3.data(), a_meso_tn3.data(), a_meso_tgn3.data(), mn2, zn2.data(), a_meso_tn2.data(), a_meso_tgn2.data());
        t[1]=tz;
    }


    void CNrlmsise00_p::gtd7d(const int doy, const double sec, const double& alt,
                const double& g_lat, const double& g_long, const double& lst, const double f107A, const double f107,
                std::array<double,7>& ap, std::array<double,9>& d, std::array<double,2>& t)
    {
        gtd7(doy, sec, alt,g_lat, g_long, lst, f107A, f107, ap, d, t);
        d[5] = 1.66E-24 * (4.0 * d[0] + 16.0 * d[1] + 28.0 * d[2] + 32.0 * d[3] + 40.0 * d[4] + d[6] + 14.0 * d[7] + 16.0 * d[8]);
        if (a_sw[0])
        {
            d[5]=d[5]/1000.0;
        }
    }
}