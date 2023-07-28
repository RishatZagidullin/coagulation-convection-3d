#pragma once
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <cmath>

//d->y, w->z, h->x (for grid)
///1->x, 2->z, 3->y (for camera)
struct grid_params{
    int d, w, h;
    double left_d, left_w, left_h, dd;
    grid_params(int d, int w, int h,
                double left_d, double left_w,
                double left_h, double dd) : d(d), w(w), h(h),
                    left_d(left_d), left_w(left_w), left_h(left_h),
                    dd(dd){}
};

class projection_reg {
public:
    //grid sections parameters
    grid_params grid;
    int max_v, max_i;
    //vertices and indices of the grid sections
    float * vertices;
    unsigned int * indices;
    //vertices and indices of the vector field
    //vectors are located at the center of each grid cell
    float * vect_verts;
    unsigned int * vect_inds;
    //vertices and indices of additional objects in the computational domain
    projection_reg(grid_params & grid) : grid(grid)
    {
        init_data();
    }
    ~projection_reg();
private:
    void init_data();
};


void projection_reg::init_data()
{
    max_v = (grid.h*grid.w+grid.h*grid.d+grid.w*grid.d);
    max_i = ((grid.h-1)*(grid.w-1)+(grid.h-1)*(grid.d-1)
                                  +(grid.w-1)*(grid.d-1));

    vect_inds = new unsigned int [max_i*2];
    for (int i = 0; i < max_i; i++)
    {
        vect_inds[2*i] = i;
        vect_inds[2*i+1] = i+max_i;
    }

    vertices = new float [max_v*6];
    indices = new unsigned int [max_i*6];

    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.w; j++)
        {
            int ind = 6*(i*grid.w+j);
            vertices[ind] = grid.left_h+i*grid.dd; 
            vertices[ind+1] = grid.left_w+j*grid.dd;
            vertices[ind+2] = 0.0;

            vertices[ind+3] = 0.0;
            vertices[ind+4] = 0.0;
            vertices[ind+5] = 0.0;
        }

    for (int i = 0; i < grid.h-1; i++)
    {
        int offset = i*grid.w;
        for (int j = 0; j < grid.w-1; j++)
        {
            indices[6*(i*(grid.w-1)+j)] = j+offset;
            indices[6*(i*(grid.w-1)+j)+1] = j+1+offset;
            indices[6*(i*(grid.w-1)+j)+2] = j+grid.w+offset;
            indices[6*(i*(grid.w-1)+j)+4] = j+1+offset;
            indices[6*(i*(grid.w-1)+j)+3] = j+grid.w+offset;
            indices[6*(i*(grid.w-1)+j)+5] = j+grid.w+1+offset;
        }
    }
    int global_offset = grid.h*grid.w;
    for (int i = 0; i < grid.w; i++)
        for (int j = 0; j < grid.d; j++)
        {
            vertices[6*(i*grid.d+j+global_offset)] = 0.0; 
            vertices[6*(i*grid.d+j+global_offset)+1] = 
                                        grid.left_w+i*grid.dd;
            vertices[6*(i*grid.d+j+global_offset)+2] =
                                        grid.left_d+j*grid.dd;
            vertices[6*(i*grid.d+j+global_offset)+3] = 0.0;
            vertices[6*(i*grid.d+j+global_offset)+4] = 0.0;
            vertices[6*(i*grid.d+j+global_offset)+5] = 0.0;
        }
    global_offset = (grid.h-1)*(grid.w-1);
    for (int i = 0; i < grid.w-1; i++)
    {
        int offset = i*grid.d + grid.h*grid.w;
        for (int j = 0; j < grid.d-1; j++)
        {
            indices[6*(i*(grid.d-1)+j+global_offset)] = j+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+1] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+2] = 
                                                     j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+4] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+3] = 
                                                     j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+5] = 
                                                   j+grid.d+1+offset;
        }
    }
    global_offset = grid.h*grid.w + grid.w*grid.d;
    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.d; j++)
        {
            vertices[6*(i*grid.d+j+global_offset)] = 
                                        grid.left_h+i*grid.dd; 
            vertices[6*(i*grid.d+j+global_offset)+1] = -3.0;
            vertices[6*(i*grid.d+j+global_offset)+2] =
                                        grid.left_d+j*grid.dd;
            vertices[6*(i*grid.d+j+global_offset)+3] = 0.0;
            vertices[6*(i*grid.d+j+global_offset)+4] = 0.0;
            vertices[6*(i*grid.d+j+global_offset)+5] = 0.0;
        }
    global_offset = (grid.h-1)*(grid.w-1) + (grid.w-1)*(grid.d-1);
    for (int i = 0; i < grid.h-1; i++)
    {
        int offset = i*grid.d + grid.h*grid.w + grid.w*grid.d;
        for (int j = 0; j < grid.d-1; j++)
        {
            indices[6*(i*(grid.d-1)+j+global_offset)] = j+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+1] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+2] = 
                                                    j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+4] = j+1+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+3] = 
                                                    j+grid.d+offset;
            indices[6*(i*(grid.d-1)+j+global_offset)+5] = 
                                                  j+grid.d+1+offset;
        }
    }

    return;
}

projection_reg::~projection_reg()
{
    delete [] vertices;
    delete [] indices;
    delete []  vect_inds;
}
