#pragma once
#include "../geometry/vector3d.h"
#include "../geometry/mesh.h"
#include "../geometry/projection_reg.h"
#include "../solvers/advection_tet.h"

class AnimationData_tet {
public:
    float * colors_mesh;
    float * colors_proj;
    std::vector<int> * inds;
    projection_reg &proj;
    SpaceHeter3d_tet<double, int> &eqn;
    AnimationData_tet(projection_reg & proj, 
        SpaceHeter3d_tet<double, int> & eqn) : proj(proj), eqn(eqn)
    {
       init_color_data();
    }
    void update_proj_colors();
    void update_mesh_colors();
    ~AnimationData_tet();
private:
    void init_color_data();
};

void AnimationData_tet::update_proj_colors()
{
    #pragma omp parallel for
    for (int i = 0; i < proj.max_v; i++)
    {
        colors_proj[i] = 0.0;
        for (int j = 0; j < inds[i].size(); j++)
        {
            colors_proj[i] += eqn.f[inds[i][j]];
        }
        if (inds[i].size() > 0)
            colors_proj[i] /= inds[i].size();
    }
}

void AnimationData_tet::update_mesh_colors()
{
    int c = 0;
    for (int c = 0; c < eqn.data.v_size; c++)
    {
        for (int i = 0; i < eqn.data.incident_cells[c].size(); i++)
            colors_mesh[c] += eqn.f[eqn.data.incident_cells[c][i]];
        colors_mesh[c] /= ((float) eqn.data.incident_cells[c].size());
    }
}

void AnimationData_tet::init_color_data()
{
    inds = new std::vector<int> [proj.max_v];
    colors_mesh = new float[eqn.data.v_size];
    colors_proj = new float[proj.max_v];

    for (int ind=0; ind<eqn.data.c_size; ind++)
    {
        float val = 0.1;
        float x = eqn.data.cell_centroids[ind].x;
        float y = eqn.data.cell_centroids[ind].y;
        float z = eqn.data.cell_centroids[ind].z;
        if (x*x < val)
        {
            int ind_y = (int) ((y-proj.grid.left_w)/proj.grid.dd);
            int ind_z = (int) ((z-proj.grid.left_h)/proj.grid.dd);
            int proj_ind = ind_z*proj.grid.w+ind_y;
            inds[proj_ind].push_back(ind);
        }
        if (z*z < val)
        {
            int ind_x = (int) ((x-proj.grid.left_d)/proj.grid.dd);
            int ind_y = (int) ((y-proj.grid.left_w)/proj.grid.dd);
            int offset = proj.grid.w*proj.grid.h;
            int proj_ind = ind_y*proj.grid.d+ind_x+offset;
            inds[proj_ind].push_back(ind);   
        }
        if (y*y < val)
        {
            int ind_x = (int) ((x-proj.grid.left_d)/proj.grid.dd);
            int ind_z = (int) ((z-proj.grid.left_h)/proj.grid.dd);
            int offset = proj.grid.w*proj.grid.h + 
                         proj.grid.w*proj.grid.d;
            int proj_ind = ind_z*proj.grid.d+ind_x+offset;
            inds[proj_ind].push_back(ind);
        }
    }
    return;
}

AnimationData_tet::~AnimationData_tet()
{
    std::cout << "animation tet destructor\n";
    delete [] colors_mesh;
    delete [] colors_proj;
    delete [] inds;
}
