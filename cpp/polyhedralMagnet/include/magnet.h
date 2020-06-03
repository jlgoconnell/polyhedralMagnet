#ifndef MAGNET_H
#define MAGNET_H

#include <string>

class magnet
{
    public:
        magnet(std::string filename);
        ~magnet();

        float *getMagn();
        float **getVertices();
        int getNumVertices();

    private:
        float magnetisation[3];
        float **vertices;
        int numVertices;
};

#endif // MAGNET_H
