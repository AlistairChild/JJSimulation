//property of Alistair child 
//alistair@mtoto.org

#define _USE_MATH_DEFINES

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <typeinfo>
#include <iostream>
#include <cmath>
#include <fstream>
//#include <stdio.h>   
#include <algorithm>   
#include <complex>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <typeinfo>
#include <string>
//#include <iterator>
using namespace std;


void generic_args(po::options_description &desc,po::positional_options_description &p) {
    desc.add_options()
        ("help,h", "produce help message")
        ("lowfield", po::value<double>()->default_value(0), "lowest field used (default 0)")
        ("highfield", po::value<double>()->default_value(0.15), "highest field used (default 0.15)")
        ("magres", po::value<int>()->default_value(60), "n points per angle (default 100)")
        ("intres", po::value<int>()->default_value(60), "n trapezoid integral! (default 100)")
        ("length", po::value<double>()->default_value(500E-9), "length default (500E-9)")
        ("width", po::value<double>()->default_value(500E-9), "width default (500E-9)")
        ("height", po::value<double>()->default_value(1E-9), "height default (1E-9)")
        ("angleincrement", po::value<double>()->default_value(0.5), "how much to inrement angle by 0-90 (default 1).")
        ("step_height", po::value<double>()->default_value(10), "fraction of max critica current. (default 0.1)")
        ("write_to",   po::value<std::string>()->default_value("tsv"),      "csv or tsv ")
        ("distribution",   po::value<std::string>()->default_value("2d"),      "current density 1d or 2d")
        ("view",   po::value<std::string>()->default_value("critical_current"),      "profile or critical_current ")
        ("step_thickness", po::value<double>()->default_value(0.1), "fraction of length.(default 0.1) ");
    p.add("input", -1);
}


struct Results
{
    vector<vector<double> > critical_current;
    vector<vector<double> > magnetic_flux;
    vector<vector<double> > current_density_profile;
};

struct Options
{
    double lowfield;
    double highfield;  
    int N;
    int n;
    double L;
    double w;
    double d;
    double angleincrement;
    double stepheight;
    double stepthickness;
    string write_to;
    string view;
    string distribution;
};

void resolvedfields(const Options &opts, double k, double theta, double field, double *Bx, double *By, double *kx, double *ky)
{
            
	*Bx = sin((theta * M_PI)/180)*(field);
	*By = cos((theta * M_PI)/180)*(field);
	*kx = k * *By;// + ((4*M_PI)/opts.L)* ((a / (1- (b*exp(-K * (*By - C) ) ) ) ) - (a*0.5));
	*ky = k * *Bx; //+ ((4*M_PI)/opts.w)* ((a / (1- (b*exp(-K * (*Bx - C) ) ) ) ) - (a*0.5));
}

double geometry(const Options &opts, double x, double y, double thickness, double L, double w, double height)
{
    double J0;
    if ((x <= (- L/2 + thickness) || x >= (L/2 - thickness)) && (opts.distribution == "2d"))
    {
        //J0 = height * cos((2 * M_PI * x) / L - M_PI) + height;
        J0 = height;
    }
    else if ((y <= (- w/2 + thickness) || y >= (w/2 - thickness)))
    {
        J0 = height;
        //J0 = height * cos((2 * M_PI * y) / L - M_PI) + height;
    }
    else
    {
        J0 = 0;
        //J0 = height * cos((4 * M_PI * x ) / L) + height;
    }
    return J0;
}

complex<double> current_density(const Options &opts,double Bx,double By,double x, double y, complex<double> Kx, complex<double> Ky, double thickness,double L, double w,double height)
{

    double a =100;
    double b = 0.33341;
    double C = 0.1;
    double K = 0.12089;
    
    double J0 = geometry(opts, x, y, thickness, L, w, height);
	return 1E7 * J0 * exp( Kx * x - Ky * y);// + ((((4*M_PI)/opts.L)* ((a / (1- (b * exp(-K * (2 - C) ) ) ) ))) - (a*0.5))+ (((4*M_PI)/opts.L)* ((a / (1- (b * exp(-K * (0 - C) ) ) ) ) - (a*0.5))));//complex exponential (Phase relations)
}

double sumdoubleintegral(const Options &opts, double Bx, double By, complex<double> Kx, complex<double> Ky,double lowbound1, double lowbound2 ,int n, double dy, double dx, double thickness,double L, double w,double height)
{
	complex<double> cumbigsum (0,0);//define and initialise a complex cumsum for outer integral
	
	for(int i=0; i<n ;i++)//set up for loop for outer
	{	
		double xi = lowbound1 + i*dx;//define the trapezoid width
		complex<double>  cumsmallsum (0,0);//inner cumsum initialised
		
		for (int j=0; j<n; j++)//inner loop
		{
			double yi = lowbound2 + dy * j;//define trapeziod width 
			complex<double> funvalue = current_density(opts, Bx ,By,xi, yi, Kx, Ky, thickness, L, w, height);//call function to evaluate return complex double
			complex<double> rectanglearea = funvalue * dy;//multiply width by height
			cumsmallsum += rectanglearea;//add to inner cumsum
		}
		
		complex<double> secondrectanglearea = cumsmallsum*dx;//use total inner cumsum as cross section
		cumbigsum += secondrectanglearea;//add up sliced area to get total
	}
		
	return abs(cumbigsum);//return the absolute value of complex number
}

void make_profile(const Options &opts, Results &result)
{
    double lowbound1 = -opts.L/2, lowbound2 = -opts.w/2;
	double upbound1 = opts.L/2, upbound2 = opts.w/2;
    double dx = (double)(upbound2-lowbound2)/opts.n;//dx
	double dy = (double)(upbound1-lowbound1)/opts.n;//dy
    
    for(int xaxis = 0; xaxis<opts.n ; ++xaxis){
        double xi = lowbound1 + xaxis * dx;
        result.current_density_profile.push_back(vector<double>());
        for(int yaxis = 0 ; yaxis<opts.n ; ++yaxis){
            double yi = lowbound2 + dy * yaxis;
            double level = geometry(opts, xi, yi, opts.stepthickness*opts.L, opts.L, opts.w, opts.stepheight);
            
            result.current_density_profile[xaxis].push_back(level);
        }
    }
}

void write_csv(const vector<string> &labels, const vector< vector<double> > &data)
{
    // output labels
    for(size_t i = 0; i != labels.size(); ++i)    
    {
        if (i != 0) cout << ",";
        cout << labels[i];
    }
    cout << endl;
    
    // output data
    for(const auto& datum: data) {
        for(size_t i = 0; i != datum.size(); ++i)    
        {
            if (i != 0) cout << ",";
            cout << datum[i];
        }
        cout << endl;    
    }
}

void write_tsv(const vector< vector<double> > &data)
{
    // output data
    for(const auto& datum: data) {
        for(size_t i = 0; i != datum.size(); ++i)    
        {
            if (i != 0) cout << "\t";
            cout << datum[i];
        }
        cout << endl;    
    }
}

void file_output(const Options &opts, Results &result)
{
    if (opts.view == "profile")
    {
        for(int xaxis = 0; xaxis < opts.n; ++xaxis)
	    {
	        for (int yaxis = 0; yaxis < opts.n; ++yaxis)
	        {
                if (opts.write_to == "csv"){
        	        cout << xaxis << "," << yaxis << "," << result.current_density_profile[xaxis][yaxis]  << endl;
	            }
                else if(opts.write_to == "tsv"){
                    cout << xaxis << "\t" << yaxis << "\t" << result.current_density_profile[xaxis][yaxis]  << endl;
                }
            }
	    } 
    }    
    else if (opts.view == "critical_current")
    {
        for(int angle = 0; angle <= 90/opts.angleincrement; ++angle)
	    {
            double max = *max_element(result.critical_current[angle].begin(), result.critical_current[angle].end());
            for (int flux=0; flux < opts.N; ++flux)
	        {
                if(opts.write_to == "csv"){
        	        cout << angle*opts.angleincrement << "," << result.magnetic_flux[angle][flux] << "," << result.critical_current[angle][flux]/max << endl;
	            }
                else if(opts.write_to == "tsv"){
                    cout << angle*opts.angleincrement << "\t" << result.magnetic_flux[angle][flux]<< "\t" << result.critical_current[angle][flux]/max << endl;
                }
            }   
	    } 
    } 
    else
    {
        cerr << "Error: i do not understand the option view. Please choose from the view options --help for info" << endl;
    }
    
}

void getargs(Options &opts, int ac, const char * av[])
{   
    //Options opts;
    try {
        po::options_description desc("Allowed options");
        po::positional_options_description p;
        generic_args(desc, p);//calling function
        po::variables_map param;
        po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(), param);
        po::notify(param);
        
        if (param.count("help")) {
            cout << desc << "\n";
            exit(0);
        }
        
        opts.lowfield = param["lowfield"].as<double>();
        opts.highfield = param["highfield"].as<double>();
        opts.N = param["magres"].as<int>();
        opts.n = param["intres"].as<int>();
        opts.L = param["length"].as<double>();
        opts.w = param["width"].as<double>();
        opts.d = param["height"].as<double>();
        opts.stepheight = param["step_height"].as<double>();
        opts.stepthickness = param["step_thickness"].as<double>();
        opts.angleincrement = param["angleincrement"].as<double>();
        opts.write_to = param["write_to"].as<string>();
        opts.view = param["view"].as<string>();
        opts.distribution = param["distribution"].as<string>();
    }
    catch(exception &e) {
        cerr << "error: " << e.what() << "\n";
        throw e;
    }
}

void make_results(const Options &opts, Results &result)
{
    //define Physics constants to double precision!
    double fluxquantum = 2.06783383E-15, lambdax = 90E-9;
	complex <double> im(0,1);
    double k = 2 * M_PI * (2 * lambdax + opts.d) / fluxquantum;
    
    //define the upper/lower bounds of double integral
    double lowbound1 = -opts.L/2, lowbound2 = -opts.w/2;
	double upbound1 = opts.L/2, upbound2 = opts.w/2;
    double dx = (double)(upbound2-lowbound2)/opts.n;//dx
	double dy = (double)(upbound1-lowbound1)/opts.n;//dy
	double step = (opts.highfield - opts.lowfield) / (opts.N-1);
	int angleblock = 90 / opts.angleincrement;
    
    //create the data.

    for (int angle = 0; angle <= angleblock; ++angle)
	{	
		
        //local values set to global values.
		float loc_lowfield = opts.lowfield;
		float loc_highfield = opts.highfield;
		double loc_step = step;
        
        result.critical_current.push_back(vector<double>());
        result.magnetic_flux.push_back(vector<double>());
        
		for (int flux = 0; flux < opts.N; ++flux)
		{			
			float field = loc_lowfield;
			loc_lowfield += loc_step; 
			//int theta = angle * opts.angleincrement;
            double theta = angle * opts.angleincrement;

			double Bx, By ,kx, ky;
			resolvedfields(opts, k, theta, field, &Bx, &By , &kx , &ky);

            
            double appendff = sumdoubleintegral(opts, Bx, By,kx * im, ky * im, lowbound1, lowbound2, opts.n, dy, dx,opts.stepthickness*opts.L,opts.L, opts.w, opts.stepheight);
            double mflux =( ( ( (Bx * opts.w) + (By * opts.L) ) * (2 * lambdax + opts.d))  / fluxquantum );//+ ( ((opts.L)* ((a / (1- (b*exp(-K * (2- C) ) ) ) ) - (a*0.5))) / fluxquantum )+ ( ((1)* ((a / (1- (b*exp(-K * (2 - C) ) ) ) ) - (a*0.5))) / fluxquantum );
            result.critical_current[angle].push_back(appendff);
            result.magnetic_flux[angle].push_back(mflux);
		}	
	}  
    
    
}

int main(int ac, const char* av[])
{
    //Options opts;
    try {
        //Define the user options
        Options opts;  
        
        //parse arguments  
        getargs(opts, ac, av);
        
        //Results(opts) output defined;
        Results result;
        
        //create results call Math
        make_results(opts, result);
    
        //make the current density profile
        make_profile(opts, result);
        
        //write results to a file csv, tsv, 1d 2d.
        file_output(opts, result);
    }
    catch(exception &e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
