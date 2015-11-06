#include <iostream>
#include <fstream>
#include "Matrix.h"

using namespace std;
#include "r8mat_expm1.h"

int main() {

    Matrix testMatrix(3);
    //std::cin >> testMatrix; [ 1 2 3 ;
    //                          2 6 1 ;
    //                          2 5 4 ] \n;

    testMatrix[0][0] = 1;
    testMatrix[0][1] = 2;
    testMatrix[0][2] = 3;
    testMatrix[1][0] = 2;
    testMatrix[1][1] = 6;
    testMatrix[1][2] = 1;
    testMatrix[2][0] = 2;
    testMatrix[2][1] = 5;
    testMatrix[2][2] = 4;

    Matrix temp = expMatrix(testMatrix);
    std::cout << "Approximation \n";
    std::cout.precision(15);
    std::cout << temp;

    double a[] = {1,2,3,2,6,1,2,5,4};
    double* reference = r8mat_expm1(3,a);
    string fileName = "example.txt";
    temp.printMatrix(fileName);

    Matrix refMatrix(3);
    for(std::size_t i = 0; i < 3;i++){
        for(std::size_t j = 0; j<3 ; j++){
            refMatrix[i][j] = reference[3*i+j];
        }
    }
    std::cout << "\n\n Reference matrix \n";
    std::cout << refMatrix << "\n";

    std::cout << "norm of error\n";
    Matrix difference = temp - refMatrix;
    std::cout << difference.norm() << "\n";

    return 0;
}