#pragma once
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <iterator>
#include <cmath>

//d->x, w->y, h->z
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
        gen_rot_mats(45./180.*3.14, 45/180.*3.14, 45/180.*3.14);
        init_data();
        
    }
    void out_to_file(std::string vfile, std::string ffile);
    ~projection_reg();
private:
    float rot_mat_1 [9];
    float rot_mat_2 [9];
    float rot_mat_3 [9];
    float trans_vect [3] = {0., 0., 1.};
    void gen_rot_mats(float th_1, float th_2, float th_3);
    void init_data();
    void matvec(float * mat, float * vec);
};

void projection_reg::matvec(float * mat, float * vec)
{
    float res [3] = {0, 0, 0};
    res[0] = mat[0]*vec[0]+mat[1]*vec[1]+mat[2]*vec[2];
    res[1] = mat[3]*vec[0]+mat[4]*vec[1]+mat[5]*vec[2];
    res[2] = mat[6]*vec[0]+mat[7]*vec[1]+mat[8]*vec[2];
    vec[0] = res[0];
    vec[1] = res[1];
    vec[2] = res[2];
    return;
}

void projection_reg::gen_rot_mats(float th_1, float th_2, float th_3)
{
    rot_mat_1[0] = 1.;
    rot_mat_1[1] = 0.;
    rot_mat_1[2] = 0.;
    rot_mat_1[3] = 0.;
    rot_mat_1[4] = cos(th_1);
    rot_mat_1[5] = -sin(th_1);
    rot_mat_1[6] = 0.;
    rot_mat_1[7] = sin(th_1);
    rot_mat_1[8] = cos(th_1);

    rot_mat_2[0] = cos(th_2);
    rot_mat_2[1] = 0.;
    rot_mat_2[2] = sin(th_2);
    rot_mat_2[3] = 0.;
    rot_mat_2[4] = 1.;
    rot_mat_2[5] = 0.;
    rot_mat_2[6] = -sin(th_2);
    rot_mat_2[7] = 0;
    rot_mat_2[8] = cos(th_2);

    rot_mat_3[0] = cos(th_3);
    rot_mat_3[1] = -sin(th_3);
    rot_mat_3[2] = 0.;
    rot_mat_3[3] = sin(th_3);
    rot_mat_3[4] = cos(th_3);
    rot_mat_3[5] = 0.;
    rot_mat_3[6] = 0.;
    rot_mat_3[7] = 0.;
    rot_mat_3[8] = 1.;

    return;
}

void projection_reg::out_to_file(std::string vfile, std::string ffile)
{
    std::ofstream v_reader;
    v_reader.open(vfile);
    v_reader << grid.d << " " << grid.w << " " << grid.h << " " << grid.dd << "\n";
    std::ofstream f_reader;
    f_reader.open(ffile);
    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.w; j++)
            v_reader << grid.left_h+i*grid.dd << " " 
                     << grid.left_w+j*grid.dd << " 0.0\n";

    for (int i = 0; i < grid.h-1; i++){
        int offset = i*grid.w;
        for (int j = 0; j < grid.w-1; j++){
            f_reader << j+offset << " " << j+1+offset << " "
                     << j+grid.w+offset << "\n";
            f_reader << j+1+offset << " " << j+grid.w+offset
                     << " " << j+grid.w+1+offset << "\n";
        }
    }
    for (int i = 0; i < grid.w; i++)
        for (int j = 0; j < grid.d; j++)
            v_reader << "0.0 " << grid.left_w+i*grid.dd << " "
                     << grid.left_d+j*grid.dd << std::endl;

    for (int i = 0; i < grid.w-1; i++){
        int offset = i*grid.d + grid.h*grid.w;
        for (int j = 0; j < grid.d-1; j++){
            f_reader << j+offset << " " << j+1+offset << " "
                     << j+grid.d+offset << "\n";
            f_reader << j+1+offset << " " << j+grid.d+offset
                     << " " << j+grid.d+1+offset << "\n";
        }
    }

    for (int i = 0; i < grid.h; i++)
        for (int j = 0; j < grid.d; j++)
            v_reader << grid.left_h+i*grid.dd << " 0.0 "
                     << grid.left_d+j*grid.dd << std::endl;

    for (int i = 0; i < grid.h-1; i++){
        int offset = i*grid.d + grid.h*grid.w + grid.w*grid.d;
        for (int j = 0; j < grid.d-1; j++){
            f_reader << j+offset << " " << j+1+offset << " "
                     << j+grid.d+offset << "\n";
            f_reader << j+1+offset << " " << j+grid.d+offset
                     << " " << j+grid.d+1+offset << "\n";
        }
    }

    v_reader.close();
    f_reader.close();
}

void projection_reg::init_data()
{
    max_v = (grid.h*grid.w+grid.h*grid.d+grid.w*grid.d);
    max_i = ((grid.h-1)*(grid.w-1)+(grid.h-1)*(grid.d-1)
                                  +(grid.w-1)*(grid.d-1));

    vect_verts = new float [max_i*6];
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

            //matvec(rot_mat_1, &vertices[ind]);
            //matvec(rot_mat_2, &vertices[ind]);
            //matvec(rot_mat_3, &vertices[ind]);
            //vertices[ind] += trans_vect[0];
            //vertices[ind+1] += trans_vect[1];
            //vertices[ind+2] += trans_vect[2];

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

            vect_verts[3*(i*(grid.w-1)+j)] = grid.left_h+i*grid.dd+0.5*grid.dd; 
            vect_verts[3*(i*(grid.w-1)+j)+1] = grid.left_w+j*grid.dd+0.5*grid.dd;
            vect_verts[3*(i*(grid.w-1)+j)+2] = 0.0;

            vect_verts[3*(i*(grid.w-1)+j)+max_i*3] = grid.left_h+i*grid.dd+1.2*grid.dd; 
            vect_verts[3*(i*(grid.w-1)+j)+1+max_i*3] = grid.left_w+j*grid.dd+0.5*grid.dd;
            vect_verts[3*(i*(grid.w-1)+j)+2+max_i*3] = 0.0;
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

            vect_verts[3*(i*(grid.d-1)+j+global_offset)] = 0.0;
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+1] = grid.left_w+i*grid.dd+0.5*grid.dd;
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+2] = grid.left_d+j*grid.dd+0.5*grid.dd;

            vect_verts[3*(i*(grid.d-1)+j+global_offset)+max_i*3] = 0.0+0.7*grid.dd; 
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+1+max_i*3] = grid.left_w+i*grid.dd+0.5*grid.dd;
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+2+max_i*3] = grid.left_d+j*grid.dd+0.5*grid.dd;
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

            vect_verts[3*(i*(grid.d-1)+j+global_offset)] = grid.left_h+i*grid.dd+0.5*grid.dd;
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+1] = 0.0;
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+2] = grid.left_d+j*grid.dd+0.5*grid.dd;

            vect_verts[3*(i*(grid.d-1)+j+global_offset)+max_i*3] = grid.left_h+i*grid.dd+1.2*grid.dd; 
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+1+max_i*3] = 0.0;
            vect_verts[3*(i*(grid.d-1)+j+global_offset)+2+max_i*3] = grid.left_d+j*grid.dd+0.5*grid.dd;
        }
    }

    return;
}

projection_reg::~projection_reg()
{
    delete [] vertices;
    delete [] indices;
    delete []  vect_verts;
    delete []  vect_inds;
}
