#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <vector>
#include <iterator>
#include <cmath>
#include "vector3d.h"

template<typename T>
std::vector<T> split(const std::string& line)
{
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is),
                          std::istream_iterator<T>());
}

std::string LED_format(double value)
{
    std::string str = std::to_string(value);
    str.erase ( str.find_last_not_of('0') + 1, std::string::npos );
    if (str[str.size()-1] == '.')
        str.append("0");
    return str;
}

template<int n_velos>
class parse_foam{
public:
    Vector3d<double> * velos [n_velos];
    Vector3d<double> * cell_centers;
    int n_points;

    parse_foam(std::string filename, size_t n_p) : n_points(n_p)
    {
        cell_centers = new Vector3d<double> [n_p];
        for (int i = 0; i < n_velos; i++)
            velos[i] = new Vector3d<double> [n_p];

        std::ifstream in;
        std::string f = filename+"/cc.txt";
        in.open(f);
        //std::cout << "file: " << f << "\n";

        std::string line;
        int c = 0;
        while(getline(in, line))
        {
            std::vector<float> vec = split<float>(line);
            cell_centers[c].x = vec[0];
            cell_centers[c].y = vec[1];
            cell_centers[c].z = vec[2];
            c++;
        }
        in.close();

        for (int i = 0; i < n_velos; i++)
        {
            c = 0;
            std::string f = filename + "/new/";
            std::stringstream stream;
            stream << std::fixed << std::setprecision(2) << LED_format((i+1)*0.01);
            std::string s = stream.str();
            f = f + s + "/u.txt";

            //std::cout << "file: " << f << "\n";
            in.open(f);
            while(getline(in, line))
            {
                std::vector<float> vec = split<float>(line);
                velos[i][c].x = vec[0];
                velos[i][c].y = vec[1];
                velos[i][c].z = vec[2];
                c++;
            }
            in.close();
        }
    }
    ~parse_foam();
};

template< int n_velos>
parse_foam<n_velos>::~parse_foam()
{
    delete [] cell_centers;
    for (int i = 0; i < n_velos; i++)
        delete [] velos[i];
}
