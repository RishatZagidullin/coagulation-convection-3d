#pragma once
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <cmath>

template<typename T>
std::vector<T> split(const std::string& line)
{
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is),
                          std::istream_iterator<T>());
}


enum parse {header, verts, inds};

class domain_obj {
public:
    int v_size, i_size;
    float * vertices;
    unsigned int * indices;
    domain_obj(std::string filename)
    {
        parse p = header;
        std::ifstream in;
        in.open(filename);
        std::string line;
        getline(in, line);
        if (line != "OFF")
        {
           throw std::invalid_argument("File not in OFF format");
        }
        int c = 0;
        while(getline(in, line))
        {
            if (line == "") continue;
            if (p == header)
            {
                std::vector<int> vec = split<int>(line);
                v_size = vec[0];
                i_size = vec[1];
                p = verts;
                vertices = new float [v_size*3];
                indices = new unsigned int [i_size*3];
            }
            else if (p == verts)
            {
                std::vector<float> vec = split<float>(line);
                vertices[c*3] = vec[2];
                vertices[c*3+1] = vec[1];
                vertices[c*3+2] = vec[0];
                c++;
                if (c == v_size)
                {
                    p = inds;
                    c = 0;
                }
            }
            else if (p == inds)
            {
                std::vector<unsigned int> vec = split<unsigned int>(line);
                indices[c*3] = vec[1];
                indices[c*3+1] = vec[2];
                indices[c*3+2] = vec[3];
                c++;
                if (c == i_size)
                {
                    c = 0;
                    break;
                }
            }
            
        }
        in.close();
    }
    ~domain_obj();
};

domain_obj::~domain_obj()
{
    delete [] vertices;
    delete [] indices;
}
