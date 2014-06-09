/* -------------------------------------------------
 *				TEST HaarTools class
 * -------------------------------------------------
 *
 * Authors : 
 *		Paul Cazeaux
 * 
 * History :
 * 	 2014-06-07 - 1st version
 */

  // Compiler
 #include <iostream>
 #include <cassert>
 #include <cstdlib>
 #include <HaarTools.h>

  // Norm 
double NormDiff(std::vector<double>& u, std::vector<double>& v)
{
	double norm2=0;
	assert(u.size()==v.size());
	for (int i=0; i<u.size(); i++)
		norm2 += (u[i]-v[i])*(u[i]-v[i]);
	return std::sqrt(norm2);
}


  // Main
 int main(int argc, char *argv[])
 {
 	std::cout << "Start of t-HaarTools" << std::endl;

 	// Instanciation
 	int size = 16;
 	int maxdepth = 4;
 	int mindepth = 1;
 	std::vector<double> Ht(size), H = {0.,0.,0.,0.,0.,0.0326,0.1577,0.8482,2.0444,2.1641,1.0005,0.1740,0.0272,0.,0.,0.};
  	std::vector<double> St(size), S = {0.,0.,1.0385,5.4102,0.,1.0385,5.3830,0.0272,0.,0.,0.0326,1.0059,4.2085,1.1745,0.0272,0.};
  	std::vector<double> Dt(size), D = {0.,0.,-1.0385,5.3558,0.,-0.9733,3.0340,0.0272,0.,0.,-0.0326,-0.6905,-0.1197,0.8265,0.0272,0.};


 	HaarTools::ForwardTransform(H, St, Dt, maxdepth, mindepth);

 	assert(NormDiff(St,S)<1e-6);
 	assert(NormDiff(Dt,D)<1e-6);
 	std::cout << "\t> Test Forward transform ...... ok" << std::endl;

 	HaarTools::InverseTransform(Ht, S, D, maxdepth, mindepth);
 	
 	assert(NormDiff(Ht,H)<1e-6);
 	std::cout << "\t> Test Inverse transform ...... ok" << std::endl;
 	
  	std::vector<double> Ds = {0.,0,-1.0385,5.3558,0.,-0.9733,3.0340,0.0272,0.,0.,0.,0.,0.,0.,0.,0.};
  	std::vector<bool> Pt(size), P = {false,false,true,true,false,true,true,true,false,false,false,false,false,false,false,false};

  	HaarTools::PUREShrink(S, D, Pt, maxdepth, mindepth);

 	assert(Pt == P);
 	assert(NormDiff(Ds,D)<1e-6);
 	std::cout << "\t> Test PUREShrink ...... ok" << std::endl;

 	std::cout << "End ot t-HaarTools" << std::endl;
 }