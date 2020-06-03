#include "magnet.h"
#include <iostream>
#include <fstream>
#include <string>

magnet::magnet(std::string filename)
{
    std::string line;
    bool readMag = false;
    bool readVerts = false;
    int ii;
    int strStart;
    int strLen;
    int ctr = 0;
    int vertctr = 0;
    int row = 0;

    // Open the magnet file
    std::ifstream magFile (filename);

    // Calculate the number of vertices
    while (getline(magFile,line)) {
        if (readVerts) ctr++;
        if (!line.compare("vertices {")) readVerts = true;
        if (readVerts && !line.compare("}")) readVerts = false;
    }
    numVertices = ctr-1;

    // Initialise 2D vertices array
    vertices = new float*[numVertices];
    for (ii = 0; ii < numVertices; ii++) {
        vertices[ii] = new float[3];
    }

    // Start again at top of file
    magFile.clear();
    magFile.seekg(0);

    // Read magnet data
    while (getline(magFile,line)) {
        // Set flags
        if (readMag && !line.compare("}")) readMag = false;
        if (readVerts && !line.compare("}")) readVerts = false;

        // Read parameters if necessary
        if (readMag) {
            strStart = 0;
            ctr = 0;
            strLen = line.length();

            for (ii = 0; ii < strLen+1; ii++) {
                if (line[ii] == ' ' || line[ii] == ',' || line[ii] == '\0') {
                    char *myString = new char[ii-strStart];
                    line.copy(myString,ii-strStart,strStart);
                    strStart = ii+1;
                    magnetisation[ctr] = std::stof(myString);
                    delete myString;
                    ctr++;
                }
            }
        }

        if (readVerts) {
            strLen = line.length();
            int col = 0;
            strStart = 0;

            for (ii = 0; ii < strLen+1; ii++) {
                if (line[ii] == ' ' || line[ii] == ',' || line[ii] == '\0') {
                    char *myString = new char[ii-strStart+1];
                    line.copy(myString,ii-strStart,strStart);
                    myString[ii-strStart] = '\0';
                    strStart = ii+1;
                    vertices[row][col] = std::stof(myString);
                    delete myString;
                    col++;
                }
            }
            row++;
        }

        // Set flags
        if (!line.compare("magnetisation {")) readMag = true;
        if (!line.compare("vertices {")) readVerts = true;
    }
}

float *magnet::getMagn() {
    return magnetisation;
}

float **magnet::getVertices() {
    return vertices;
}

int magnet::getNumVertices() {
    return numVertices;
}

magnet::~magnet()
{
    //dtor
}
