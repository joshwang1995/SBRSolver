#include "Interpolation.hpp"

void InterpPattern
(
	double theta, 
	double phi, 
	const std::vector<std::vector<double>>& Gxy, 
	const std::vector<std::vector<double>>& Gxz, 
	const std::vector<std::vector<double>>& Gyz, 
	double & GBS,
	double & GBP
)
{
	double theta1, theta2;
	if(std::islessequal(theta,PI/2.0))
	{
		theta1 = theta;
		theta2 = PI/2.0 - theta;
	}
	else
	{
		theta1 = PI-theta;
		theta2 = theta-(PI/2.0);
	}
	
	double phi1,phi2;
	if(phi == TWOPI)
	{
		phi = 0;
	}
	
	if(std::islessequal(phi, PI))
	{
		if(std::islessequal(phi, PI/2.0))
		{
			phi1 = phi;
			phi2 = (PI/2.0) - phi;
		}
		else
		{
			phi1 = phi-(PI/2.0);
			phi2 = PI - phi;
		}	
	}
	else
	{
		if(std::islessequal(phi,3.0*PI/2.0))
		{
			phi1 = phi - PI;
			phi2 = (3.0*PI/2.0) - phi;
		}
		else
		{
			phi1= phi-(3.0*PI/2.0);
			phi2= TWOPI - phi;
		}
	}
	
	std::complex<double> gt1;
	std::complex<double> gt2;
	std::complex<double> gp1;
	std::complex<double> gp2;
	InterpLine(theta,phi,Gxy,Gxz,Gyz,gt1,gt2,gp1,gp2);
	
	std::complex<double> GBC;
	double facttheta, factphi, denom;
	if(theta==0 || theta == PI)
	{
		GBC = gt1;
	}
	else
	{
		facttheta=(theta1*theta2)/((theta1+theta2)*(theta1+theta2));
		factphi=(phi1*phi2)/((phi1+phi2)*(phi1+phi2));
		denom=((phi1+phi2)*facttheta)+((theta1+theta2)*factphi);
		GBC=((((phi1*gp2) + (phi2*gp1))*facttheta)+((theta1*gt2+theta2*gt1)*factphi))/denom;
	}
	GBS=10*log10(abs(GBC));
	GBP=180*arg(GBC)/PI;
}

void InterpLine
(
	double theta, 
	double phi,
	const std::vector<std::vector<double>>& Gxy, 
	const std::vector<std::vector<double>>& Gxz, 
	const std::vector<std::vector<double>>& Gyz,
	std::complex<double>& gt1,
	std::complex<double>& gt2,
	std::complex<double>& gp1,
	std::complex<double>& gp2
)
{
	std::vector<double> p;
	std::vector<std::complex<double>> Gt2;
	
	for (int k = 0; k < Gxy.size(); k++) 
	{
		if(Gxy[k][0] >= 0)
		{
			p[k] = PI*Gxy[k][0]/180.0;
		}
		else
		{
			p[k] = PI*(Gxy[k][0]+360.0)/180.0;
		}
		Gt2[k] = pow(10.0,(Gxy[k][1]/10.0))*exp(j*PI*Gxy[k][2]/180.0);
	}
	reorder_vector(p,sort_indices(p));
	reorder_vector(Gt2,sort_indices(p));
	
	std::vector<double> t1;
	std::vector<double> t3;
	std::vector<std::complex<double>> Gp1;
	std::vector<std::complex<double>> Gp3;
	for (int k = 0; k < Gxz.size(); k++) 
	{
		if(Gxz[k][0] >= 0)
		{
			t1.push_back(PI*Gxz[k][0]/180.0);
			Gp1.push_back(pow(10.0,(Gxz[k][1]/10.0))*exp(j*PI*Gxz[k][2]/180.0));
		}
		else
		{
			t3.push_back(PI*abs(Gxz[k][0])/180.0);
			Gp3.push_back(pow(10.0,(Gxz[k][1]/10.0))*exp(j*PI*Gxz[k][2]/180.0));
		}
	}
	reorder_vector(t1,sort_indices(t1));
	reorder_vector(Gp1,sort_indices(t1));
	reorder_vector(t3,sort_indices(t3));
	reorder_vector(Gp3,sort_indices(t3));
	
	int n = t1.size()-1;
	t3.push_back(t1[n]);
	t3.insert(t3.begin(),t1[0]);
	Gp3.push_back(Gp1[n]);
	Gp3.insert(Gp3.begin(),Gp1[0]);
	
	std::vector<double> t2;
	std::vector<double> t4;
	std::vector<std::complex<double>> Gp2;
	std::vector<std::complex<double>> Gp4;
	for (int k = 0; k < Gyz.size(); k++) 
	{
		if(Gxz[k][0] >= 0)
		{
			t2.push_back(PI*Gyz[k][0]/180.0);
			Gp2.push_back(pow(10.0,(Gyz[k][1]/10.0))*exp(j*PI*Gyz[k][2]/180.0));
		}
		else
		{
			t4.push_back(PI*abs(Gyz[k][0])/180.0);
			Gp4.push_back(pow(10.0,(Gyz[k][1]/10.0))*exp(j*PI*Gyz[k][2]/180.0));
		}
	}
	reorder_vector(t2,sort_indices(t2));
	reorder_vector(Gp2,sort_indices(t2));
	reorder_vector(t4,sort_indices(t4));
	reorder_vector(Gp4,sort_indices(t4));
	
	n = t2.size()-1;
	t4.push_back(t2[n]);
	t4.insert(t4.begin(),t2[0]);
	Gp4.push_back(Gp2[n]);
	Gp4.insert(Gp4.begin(),Gp2[0]);
	
	if(std::islessequal(theta, PI/2.0))
	{
		gt1 = Gp1[0];
	}
	else
	{
		gt1 = Gp1[n];
	}
	
	gt2 = interp1(p,Gt2,phi);
	if(std::islessequal(phi, PI))
	{
		if(std::islessequal(phi, PI/2.0))
		{
			gp1=interp1(t1,Gp1,theta);
			gp2=interp1(t2,Gp2,theta);
		}
		else
		{
			gp1=interp1(t1,Gp2,theta);
			gp2=interp1(t2,Gp3,theta);
		}
	}
	else
	{
		if(std::islessequal(phi, 3.0*PI/2.0))
		{
			gp1=interp1(t1,Gp3,theta);
			gp2=interp1(t2,Gp4,theta);
		}
		else
		{
			gp1=interp1(t1,Gp4,theta);
			gp2=interp1(t2,Gp1,theta);
		}
	}
}

std::vector<size_t> sort_indices(const std::vector<double> &v) 
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

void reorder_vector(std::vector<std::complex<double>>& v, const std::vector<size_t>& order)
{
	for(int s = 1, d; s < order.size(); ++ s ) 
	{
		for ( d = order[s]; d < s; d = order[d] ) ;
		if ( d == s ) while ( d = order[d], d != s ) std::swap( v[s], v[d] );
	}
}

void reorder_vector(std::vector<double>& v, const std::vector<size_t>& order)
{
	for(int s = 1, d; s < order.size(); ++ s ) 
	{
		for ( d = order[s]; d < s; d = order[d] ) ;
		if ( d == s ) while ( d = order[d], d != s ) std::swap( v[s], v[d] );
	}
}

std::complex<double> interp1(std::vector<double> &x, std::vector<std::complex<double>> &y, double &x_new )
{
	std::complex<double> y_new;
	double dx;
	std::complex<double> dy;
	std::complex<double> slope;
	std::complex<double> intercept;
	
	int idx = findNearestNeighbourIndex(x_new,x);
	if(idx >= x.size()-1)
	{
		idx = x.size() - 2;
	}
	
	dx = x[idx+1] - x[idx];
	dy = y[idx+1] - y[idx];
	slope = dy/dx;
	intercept = y[idx] - (x[idx]*slope);
	y_new = (slope * x_new) + intercept;
	return y_new;
}

int findNearestNeighbourIndex(double value, std::vector< double > &x )
{
    double dist = INF;
    int idx = -1;
    for( int i = 0; i < x.size(); ++i ) 
	{
        double newDist = value - x[i];
        if ( newDist > 0 && newDist < dist ) 
		{
            dist = newDist;
            idx = i;
        }
    }
    return idx;
}

