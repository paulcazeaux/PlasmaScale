#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>
#include <eigen3/Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;


void test_Poisson(double L, double epsilon, int n)
{
	for (int grid_end = 10; grid_end <= std::pow(10, n); grid_end *= 10)
	{
		
		double dx = L/grid_end;

		double scaling = epsilon/(dx*dx);

		std::vector<Eigen::Triplet<double>> M;
		SpMat _M;
		Eigen::SimplicialCholesky<SpMat> _solver;


		for (int i=0; i<grid_end-1; i++)
		{
			if (i>0) M.push_back(Eigen::Triplet<double>(i, i-1, -scaling));
			M.push_back(Eigen::Triplet<double>(i, i,  2*scaling));
			M.push_back(Eigen::Triplet<double>(i, i+1, -scaling));
		}
		/* Neumann condition at x = _length */
		M.push_back(Eigen::Triplet<double>(grid_end-1, grid_end-1, scaling));
		M.push_back(Eigen::Triplet<double>(grid_end-1, grid_end-2, -scaling));

		_M = SpMat(grid_end, grid_end);
		_M.setFromTriplets(M.begin(), M.end());
		_solver.compute(_M);

		std::vector<double> charge = std::vector<double>(grid_end+1);
		std::vector<double> potential = std::vector<double>(grid_end+1);
		for (int i=0; i<=grid_end; i++)
		{
			double y = 1. - std::pow(i*dx-1., 2.);
			charge.at(i) = 2.*(1.-2.*y)*std::exp(y);
		}

		Eigen::Map<Eigen::VectorXd> Charge(charge.data()+1, grid_end);
		Eigen::Map<Eigen::VectorXd> Potential(potential.data()+1, grid_end);

		potential.at(0) = 0; charge.at(grid_end) *= .5;
		Potential = _solver.solve(Charge);
		charge.at(grid_end) *= 2.;

		double err = 0., norm = 0.;
		for (int i=0; i<=grid_end; i++)
		{
			double y = 1. - std::pow(i*dx-1., 2.);
			err += std::pow(potential.at(i) - (1.-std::exp(y)), 2.);
			norm += std::pow(1.-std::exp(y), 2.);
		}
		std::cout << grid_end << "\t" << sqrt(err/norm) << std::endl;
	}
}

void test_Nonlin_Poisson(double L, double epsilon, double occupancy, double n)
{
	std::vector<double> ion_density, exp_potential;
	std::cout << "Residuals: ";

	int N = 100;
	for (int grid_end = N; grid_end <= N*std::pow(10, n-1); grid_end *= 10)
	{
		double debye_scaling = 1.;
		int bin = occupancy*grid_end;
		double s = occupancy*grid_end - bin;
		double dx = L/grid_end;

		ion_density.resize(grid_end+1);
		exp_potential.resize(grid_end+1);

		Eigen::Map<Eigen::VectorXd> Ion_density(ion_density.data(), grid_end+1);

		for (int i=0; i<bin; i++) Ion_density(i) = 1;
		Ion_density(bin) = s;
		Ion_density *= grid_end/Ion_density.sum();
		

		double tol2 = 1e-30, iter_max = 100;
		double scaling = debye_scaling/std::pow(dx, 2.);

		Eigen::Map<Eigen::VectorXd> Exp_potential(exp_potential.data(), grid_end+1);
		Eigen::VectorXd Residual(grid_end+1), Potential(grid_end+1), Update(grid_end+1);

		std::vector<Eigen::Triplet<double> > J_T;
		Eigen::SparseMatrix<double> J;
		Eigen::SimplicialCholesky<SpMat> solver;

		/* Neumann condition at x = 0 */
		J_T.push_back(Eigen::Triplet<double>(0,0, scaling));
		J_T.push_back(Eigen::Triplet<double>(0,1, -scaling));
		for (int i=1; i<grid_end; i++)
		{
			J_T.push_back(Eigen::Triplet<double>(i, i-1, -scaling));
			J_T.push_back(Eigen::Triplet<double>(i, i,  2*scaling));
			J_T.push_back(Eigen::Triplet<double>(i, i+1, -scaling));
		}
		/* Neumann condition at x = _length */
		J_T.push_back(Eigen::Triplet<double>(grid_end, grid_end, scaling));
		J_T.push_back(Eigen::Triplet<double>(grid_end, grid_end-1, -scaling));

		J = SpMat(grid_end+1, grid_end+1);
		J.setFromTriplets(J_T.begin(), J_T.end());

		Potential = (Ion_density.array()+1e-15).log();
		Exp_potential = Potential.array().exp();
		Residual = Exp_potential - Ion_density + J*Potential;

		/* Newton's method loop */
		for (int count=0; count<iter_max; count++)
		{
			SpMat Je = J;
			Je.coeffRef(0,0) += .5* Exp_potential(0); 
			for (int i=1; i<grid_end; i++) Je.coeffRef(i,i) += Exp_potential(i); 
			Je.coeffRef(grid_end, grid_end) += .5* Exp_potential(grid_end);
			
			solver.compute(Je);

			Residual(0) *= .5; 
			Residual(grid_end) *= .5;
			Update = solver.solve(Residual);

			Potential -= Update;
			Exp_potential = Potential.array().exp();
			Residual = Exp_potential - Ion_density + J*Potential;
			std::cout << Residual.squaredNorm() << "\t";
		
			if (Residual.squaredNorm() < tol2)
				break;
		}
		std::cout << std::endl;
		std::cout << "Ion density: ";
		for (double & n : ion_density) std::cout << n << "\t";
		std::cout << std::endl;
		std::cout << "Sum: " << Ion_density.sum() << std::endl;
		std::cout << "Potential: ";
		for (int i=0; i<=grid_end; i++) std::cout << Potential(i) << "\t";
		std::cout << std::endl;
		std::cout << "Electron density: ";
		for (double & n : exp_potential) std::cout << n << "\t";
		std::cout << std::endl;
		std::cout << "Sum: " << Exp_potential.sum() << std::endl;
	}
}

int main(int argc, char **argv)
{
	double L = 1.;
	double epsilon = 1;

	test_Poisson(L, epsilon, 5);

	L = 10000;
	double occupancy = .02;

	test_Nonlin_Poisson(L, epsilon, occupancy, 2);
}
