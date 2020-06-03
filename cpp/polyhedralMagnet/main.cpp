#include <iostream>
#include "magnet.h"

using namespace std;

int main()
{
    cout << "Hello world!" << endl;

    magnet myMag("magnet1.csv");

    float *mag = myMag.getMagn();
    float **vertices = myMag.getVertices();
    int numVertices = myMag.getNumVertices();

    std::cout << "The magnetisation is [" << mag[0] << ", " << mag[1] << ", " << mag[2] << "].\n";

    std::cout << "The magnet has " << numVertices << " vertices. These are:\n";
    for (int ii = 0; ii < numVertices; ii++) {
        for (int jj = 0; jj < 2; jj++) {
            std::cout << vertices[ii][jj] << ",\t";
        }
        std::cout << vertices[ii][2] << std::endl;
    }

    return 0;
}
