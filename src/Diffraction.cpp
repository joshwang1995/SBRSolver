#include "Diffraction.hpp"

void wedge_diff_coeff
(
	double r,
	double ph, 
	double php, 
	double bo, 
	double fn, 
	std::complex<double>& ds, 
	std::complex<double>& dh, 
	std::complex<double>& dps, 
	std::complex<double>& dph
)
{
	double betn = ph - php;
	
	std::complex<double> din = di(r,betn,bo,fn);
    std::complex<double> dpn = dpi(r,betn,bo,fn);
	
	if((abs(php) > 2.5e-4) && (abs(php) < fn*180.0 - 2.5e-4))
	{
		double betp = ph + php;
		std::complex<double> dip = di(r,betp,bo,fn);
        std::complex<double> dpp = dpi(r,betp,bo,fn);
        ds = din - dip;
        dh = din + dip;
        dps = dpn + dpp;
        dph = dpn - dpp;
	}
	else
	{
        ds = 0.0;
        dh = din;
        dps = dpn;
        dph = 0.0;
	}
}

void dielec_wedge_diff_coeff
(
	double r,
	double ph, 
	double php, 
	double bo, 
	double fn,
	std::complex<double> gammah,
	std::complex<double> gammas,
	std::complex<double>& ds, 
	std::complex<double>& dh, 
	std::complex<double>& dps, 
	std::complex<double>& dph
)
{
	double betn = ph - php;
	
	std::complex<double> din = di(r,betn,bo,fn);
    std::complex<double> dpn = dpi(r,betn,bo,fn);
	
	if((abs(php) > 2.5e-4) && (abs(php) < fn*180.0 - 2.5e-4))
	{
		double betp = ph + php;
		std::complex<double> dip = di(r,betp,bo,fn);
        std::complex<double> dpp = dpi(r,betp,bo,fn);
        ds = din + gammas*dip;
        dh = din + gammah*dip;
        dps = dpn - gammas*dpp;
        dph = dpn - gammah*dpp;
	}
	else
	{
        ds = 0.0;
        dh = din;
        dps = dpn;
        dph = 0.0;
	}
}


std::complex<double> di (double r, double bet, double bo, double fn)
{
	double ang = bet * PI/180.0;
	std::complex<double> com = -exp(-j*PI/4.0)/(2.0 * TWOPI * fn * sin(bo*PI/180.0));
	double sqr = sqrt(TWOPI*r);
	
	double dns=(PI+ang)/(2.0*fn*PI);
	double sgn = (dns >= 0) ? 1 : -1;
	double n=floor(abs(dns)+0.5);
	double dn = sgn * n;
	double a = abs(1.0+cos(ang-2.0*fn*PI*dn));
    std::pair<double,double> frnels = frnels_int(2.0*sqrt(abs(r*a)));
    double c=sqrt(PI/2.0)*(0.5-frnels.first);
    double s=sqrt(PI/2.0)*(frnels.second-0.5);
    std::complex<double> fa = 2.0*j*sqr*exp(j*TWOPI*r*a)*(c+j*s);
	double rag = (PI+ang)/(2.0*fn);
    double tsin = sin(rag);
    double cota;
    if(abs(tsin) > 1.e-5)
    {
        cota=sqrt(a)*cos(rag)/tsin;
    }
    else
    {
        cota=-sqrt(2.0)*fn*sin(ang/2.0-fn*PI*dn);
    }
    std::complex<double> uppi = com*cota*fa;
    
	dns=(-PI+ang)/(2.0*fn*PI);
	sgn = (dns >= 0) ? 1 : -1;
	n=floor(abs(dns)+0.5);
	dn = sgn * n;
	a = abs(1.0+cos(ang-2.0*fn*PI*dn));
    frnels = frnels_int(2.0*sqrt(abs(r*a)));
    c=sqrt(PI/2.0)*(0.5-frnels.first);
    s=sqrt(PI/2.0)*(frnels.second-0.5);
    fa = 2.0*j*sqr*exp(j*TWOPI*r*a)*(c+j*s);
	rag = (PI-ang)/(2.0*fn);
    tsin = sin(rag);
    if(abs(tsin) > 1.e-5)
    {
        cota=sqrt(a)*cos(rag)/tsin;
    }
    else
    {
        cota=sqrt(2.0)*fn*sin(ang/2.0-fn*PI*dn);
        if (cos(ang/2.0-fn*PI*dn) < 0.0) 
        {
            cota=-cota;
        }
    }
    std::complex<double> unpi = com*cota*fa;
    std::complex<double> dir = uppi + unpi;
	return dir;
}

std::complex<double> dpi (double r, double bet, double bo, double fn)
{
	double ang = bet * PI/180.0;
	double sbo = sin(bo* PI/180.0);
	std::complex<double> com = exp(-j*PI/4.0)/(4.0*TWOPI*fn*fn*sbo*sbo);
	
	double dns=(PI+ang)/(2.0*fn*PI);
	double sgn = (dns >= 0) ? 1 : -1;
	double n=floor(abs(dns)+0.5);
	double dn = sgn * n;
	double a = abs(1.0+cos(ang-2.0*fn*PI*dn));
    std::pair<double,double> frnels = frnels_int(2.0*sqrt(abs(r*a)));
    double c=sqrt(PI/2.0)*(0.5-frnels.first);
    double s=sqrt(PI/2.0)*(frnels.second-0.5);
    std::complex<double> fpa = TWOPI*r*(2.0*j+4.0*sqrt(abs(TWOPI*r*a))*exp(j*TWOPI*r*a)*(c+j*s));
	double rag = (PI+ang)/(2.0*fn);
    double ts = sin(rag)*sin(rag);
    double csca;
    if(ts > 1.e-5)
    {
        csca=a/ts;
    }
    else
    {
        csca = -2.0*fn*fn*cos(ang-TWOPI*fn*dn)/cos((PI*ang)/fn);
    }
    std::complex<double> uppi = com*csca*fpa;
    
	dns=(-PI+ang)/(2.0*fn*PI);
	sgn = (dns >= 0) ? 1 : -1;
	n=floor(abs(dns)+0.5);
	dn = sgn * n;
	a = abs(1.0+cos(ang-2.0*fn*PI*dn));
    frnels = frnels_int(2.0*sqrt(abs(r*a)));
    c=sqrt(PI/2.0)*(0.5-frnels.first);
    s=sqrt(PI/2.0)*(frnels.second-0.5);
    fpa = TWOPI*r*(2.0*j+4.0*sqrt(abs(TWOPI*r*a))*exp(j*TWOPI*r*a)*(c+j*s));
	rag = (PI-ang)/(2.0*fn);
    ts = sin(rag)*sin(rag);
    if(ts > 1.e-5)
    {
        csca = a/ts;
    }
    else
    {
        csca = -2.0*fn*fn*cos(ang-TWOPI*fn*dn)/cos((PI*ang)/fn);
    }
    std::complex<double> unpi = com*csca*fpa;
    std::complex<double> dpir = uppi - unpi;
	return dpir;
}



std::pair<double,double> frnels_int(double xs)
{
    std::vector<double> a = {1.595769140 , -0.000001702 , -6.808568854 , -0.000576361 , 6.920691902 , -0.016898657 , -3.050485660 , -0.075752419 , 0.850663781 , -0.02563901 , -0.150230960 , 0.034404779};
    std::vector<double> b = {-0.000000033 , 4.255387524 , -0.000092810 , -7.780020400 , -0.009520895 , 5.075161298 , -0.138341947 , -1.363729124 , -0.403349276 , 0.70222206 , -0.216195929 , 0.019547031};
    std::vector<double> cc= {0.0 , -0.024933975 , 0.000003936 , 0.005770956 , 0.000689892 , -0.009497136 , 0.011948809 , -0.006748873 , 0.000246420 , 0.002102967 , -0.00121730 , 0.000233939};
    std::vector<double> d = {0.199471140 , 0.000000023 , -0.009351341 , 0.000023006 , 0.004851466 , 0.001903218 , -0.017122914 , 0.029064067 , -0.027928955 , 0.016497308 , -.005598515 , 0.000838386};
    
	double c,s;
    if(xs < 0)
    {
        c = -0.0;
        s = -0.0;
    }
    else
    {
        double x = xs;
        x=PI*x*x/2.0;
        
        double fr=0.0;
        double fi=0.0;
        double y;
        
        if((x-4.0) < 0)
        {
            y = x/4.0;
            for(int k = 11; k>=1; --k)
            {
                fr=(fr+a[k])*y;
                fi=(fi+b[k])*y;
            }
            fr=fr+a[0];
            fi=fi+b[0];
            c = (fr*cos(x)+fi*sin(x))*sqrt(y);
            s = (fr*sin(x)-fi*cos(x))*sqrt(y);
        }
        else
        {
            y = 4.0/x;
            for(int k = 11; k>=1; --k)
            {
                fr=(fr+cc[k])*y;
                fi=(fi+d[k])*y;
            }
            fr=fr+cc[0];
            fi=fi+d[0];
            c =0.5+(fr*cos(x)+fi*sin(x))*sqrt(y);
            s =0.5+(fr*sin(x)-fi*cos(x))*sqrt(y);
        }
    }
	return std::make_pair(c,s);
}