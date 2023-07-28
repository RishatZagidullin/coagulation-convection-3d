#pragma once
#include <time.h>
#include <sys/time.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>

typedef CGAL::Surface_mesh_default_triangulation_3 Tr2d;
typedef Tr2d::Geom_traits GT;
typedef GT::Point_3 Point_3;


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

float sphere(Point_3 p, Point_3 p0, float radius = 1.)
{
    float x2 = pow(p.x() - p0.x(), 2);
    float y2 = pow(p.y() - p0.y(), 2);
    float z2 = pow(p.z() - p0.z(), 2);
    return x2 + y2 + z2 - radius*radius;
}

float curve_off(Point_3 p)
{
    std::ofstream out;
    out.open("offs/thin_foil.off");
    int max_y = 20;
    int max_x = 40;
    out << "OFF\n"<< max_x*max_y*2 << " " << 
           4*(max_x-1)*(max_y-1)+4*(max_x-1)+4*(max_y-1) << " 0\n\n";

    for (int i = 0; i < max_y; i++)
    {
        for (int j = 0; j < max_x; j++)
        {
            out << j/(max_x-1.)*5-1  << " " 
                << 0.5*sin(j/(max_x-1.)*5-1)-1. << " " 
                << i/(max_y-1.)*2. - 1 << std::endl;
        }
    }
    for (int i = 0; i < max_y; i++)
    {
        for (int j = 0; j < max_x; j++)
        {
            out << j/(max_x-1.)*5-1  << " " 
                << 0.5*sin(j/(max_x-1.)*5-1)-0.8 << " " 
                << i/(max_y-1.)*2. - 1 << std::endl;
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        int offset = i*max_x;
        for (int j = 0; j < max_x-1; j++)
        {
            out << "3  " << j+offset << " " << j+1+offset << " " 
                << j+max_x+offset << "\n";
            out << "3  " << j+max_x+offset << " " << j+1+offset 
                << " " << j+1+max_x+offset << "\n";
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        int offset = i*max_x+max_x*max_y;
        for (int j = 0; j < max_x-1; j++)
        {
            out << "3  " << j+1+offset << " " << j+offset << " " 
                << j+max_x+offset << "\n";
            out << "3  " << j+1+offset << " " << j+max_x+offset 
                << " " << j+1+max_x+offset << "\n";
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        out << "3  " << i*max_x << " " << i*max_x+max_x 
            << " " << i*max_x+max_x*max_y << "\n";
        out << "3  " << i*max_x+max_x*max_y << " " << i*max_x+max_x
            << " " << i*max_x+max_x*max_y+max_x << "\n";
        out << "3  " << i*max_x+max_x+max_x-1 << " " << i*max_x+max_x-1
            << " " << i*max_x+max_x-1+max_x*max_y << "\n";
        out << "3  " << i*max_x+max_x+max_x-1 << " "
            << i*max_x+max_x-1+max_x*max_y  << " " 
            << i*max_x+max_x-1+max_x*max_y+max_x << "\n";
    }
    for (int j = 0; j < max_x-1; j++)
    {
        out << "3  " << j+1 << " " << j << " " << j+max_x*max_y <<"\n";
        out << "3  " << j+1 << " " << j+max_x*max_y << " "
            << j+max_x*max_y+1 << "\n";
        out << "3  " << j+max_x*(max_y-1) << " "<< j+1+max_x*(max_y-1)
            << " " << j+max_x*(max_y-1)+max_x*max_y << "\n";
        out << "3  " << j+max_x*(max_y-1)+max_x*max_y << " " 
            << j+1+max_x*(max_y-1) << " " 
            << j+max_x*(max_y-1)+max_x*max_y+1 << "\n";
    }
    out.close();
    return 0.;
}

float sphere_tet1(Point_3 p)
{
    Point_3 p0(0., 1.5, 0.);
    return sphere(p, p0);
}

float sphere_tet2(Point_3 p)
{
    Point_3 p0(3., 1.5, 0.);
    return sphere(p, p0);
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
