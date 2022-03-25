#include <iostream>
#include <random>

#include "Utils.h"
#include "Simulation.h"
#include "Number.h"



template <class T>
T f(T x)
{
    
    Number exp1 = 2*x * x + log(x);

    exp1.propogateAdjoints();

    std::cout << x.adjoint() << std::endl;



    return exp1;
}


int main()
{
 ///*   std::vector<std::vector<double>> A;
 //   int n;
 //   double H;
 //   double T;

 //   std::random_device rd;
 //   std::mt19937 mt(rd());
 //   std::normal_distribution<double> norm(0, 1);

 //   H = 0.12;
 //   n = 1000;
 //   T = 10;

 //   A = covMatrix(H, T, n);

 //   A = cholesky(A);
 //   
 //   std::vector<double> gaussian(1000, 0);

 //   for (size_t i = 0; i < gaussian.size(); i++)
 //   {
 //       gaussian[i] = norm(mt);
 //   }


 //   gaussian = vectorMult(A, gaussian);

 //   printVector(gaussian);*/


    Number a(5);
    Number c(10);

    Number b = f(c);


    std::cout << b.getNode()->result() << std::endl;

    b = f(a);

    std::cout << b.getNode()->result();

}
