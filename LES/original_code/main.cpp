#include "glad/glad.h"
#include "utils/utils.cuh"
//#include "utils/utils.h"
#include "solvers/advection_tet.h"
//#include "solvers/advection_reg.h"
//#include "solvers/diffusion_reg.h"
#include "geometry/mesh.h"
#include "geometry/projection_reg.h"

#include "gl/viewer.h"
#include "gl/grid_viewer.h"
	
#include "gl/cuda_gl.cuh"
#include "gl/animation_tet.h"
//#include "gl/animation_reg.h"
#include <time.h>
#include <sys/time.h>

#include "geometry/parse_foam.h"

#include "wrappers.h"

//solve first Smoluchowski integral
double L1(const int &N, const int &i, const double *n)
{
    double l1 = 0;
    for (int i1 = 0; i1 < i; i1++)
    {
        l1 += n[i1] * n[i - i1 - 1] * wrappers::K((i - i1 - 1), i1, 1);
    }
    return l1;
}

//solve second Smoluchowski integral
double L2(const int &N, const int &i, const double *n)
{
    double l2 = 0;
    for (int i1 = 0; i1 < N; i1++)
    {
        l2 += n[i1] * wrappers::K(i, i1, 1);
    }
    return l2;
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void init_data(double * data, double dd, int N, int M, int K, int S)
{
    #pragma omp parallel for collapse(2)
    for (int m = 0; m < S; m++)
    {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < 1; j++)
            for (int k = 0; k < K; k++)
            {
                if (m == 0)
                    data[m*N*M*K+i*M*K+j*K+k] = 0.3*exp(-30*pow((m)*0.0625,2)-100*pow((j)*dd,2)-100*pow((i-K/2)*dd,2)-100*pow((k-K/2)*dd,2));
                if (m > 0) data[m*N*M*K+i*M*K+j*K+k] = 0.0;
            }
    }
    return;
}

//hardcodes:
//here, init_data K/2 twice
//animation_reg.h update_proj_colors N/4 and 0/2
//projection_reg.h init_data -3.0
//advection_reg.cpp find_foam +2.0 +2.0 +3.0
//think that's it but not sure

int main(int argc, char ** argv)
{
    double start = get_wall_time();
    //int n_vort = 288000;
    //parse_foam<100> foam("./grid", n_vort);

    //domain_obj obj("offs/sphere.off");

    short SCR_W = 1600;
    short SCR_H = 1200;
    
    geometry3d geom("./offs/sphere.off", "./offs/boundary.off");
    geometry_data<double, int> grid(geom);
    geometry_vis vis(geom);
    geometry_solver_utils<double, int> solver_utils(grid);

    //grid_params proj_grid = grid_params(32,32,80,-1.0,-1.0,-3.0,0.0625);

    grid_params proj_grid = grid_params(40, 40, 80, -2.0, -3.0, -2.0, 0.1);
    projection_reg proj(proj_grid);

    int N = proj_grid.h;
    int M = proj_grid.w;
    int K = proj_grid.d;

    int S = 8;
    double *n_k;
    n_k = new double [S];

    double dd = proj_grid.dd;

    double dt = 0.001;

    SpaceHeter3d_tet<double, int> eqn(grid, solver_utils, dt);
    //solvers::Advection3d_reg<double> adv(N, M, K, dd, dt);
    //solvers::Diffusion3d_reg<double> dif(0.01, N, M, K, dd, dt);
    //adv.find_foam(foam.cell_centers, n_vort);
    //return 0;
    //double * data = new double [N*M*K*S];
    //double * data_new = new double [N*M*K*S];
    //init_data(data, dd, N, M, K, S);

    AnimationData_tet colors(proj, eqn);
    //AnimationData_reg colors(N, M, K);

    //true is to hide window
    //constexpr int num_objs = 3;
    //int obj_lens [num_objs] = {proj.max_i*6, obj.i_size*3, proj.max_i*2};
    //gl::viewer<num_objs, gl::grid_viewer> v(SCR_W, SCR_H, false);

    constexpr int num_objs = 1;
    int obj_lens [num_objs] = {proj.max_i*6};
    gl::viewer<num_objs, gl::grid_viewer> v(SCR_W, SCR_H, true);


    //viewer.buffer_vbo(vis.v_size*4, vis.vertices);
    //viewer.buffer_ebo(vis.i_size*3, vis.indices);

    //v.buffer<GL_ARRAY_BUFFER>(obj.v_size*3, obj.vertices, 2);
    //v.buffer<GL_ELEMENT_ARRAY_BUFFER>(obj.i_size*3, obj.indices, 3);
    //v.buffer<GL_ARRAY_BUFFER>(proj.max_i*6, proj.vect_verts, 4);
    //v.buffer<GL_ELEMENT_ARRAY_BUFFER>(proj.max_i*2,proj.vect_inds,5);
    v.buffer<GL_ARRAY_BUFFER>(proj.max_v*6, proj.vertices, 0);
    v.buffer<GL_ELEMENT_ARRAY_BUFFER>(proj.max_i*6, proj.indices, 1);

    //custom_gl::cuda_gl shape_interop(&viewer.VBO2,
    //                           grid.vert_data_size);
    gl::interop<float> scalar_interop(&v.buffer_objects[0],
                                      proj.max_v);

    int time = 0;
    int TIME_MAX = 1000;//200*8;//2000;
    std::cout << "Preprocessing time: " << get_wall_time()-start << "\n";
    start = get_wall_time();
    int ccount = 0;
    while (!glfwWindowShouldClose(v.window))
    {
        //do not pass the time argument 
        //if you do not want any images to be saved
        if (time % 5 == 1)
        {
            //colors.update_mesh_colors();
            //shape_interop.update_gpu_data(colors.colors_mesh);

            //colors.update_proj_colors();
            colors.update_proj_colors(data+N*M*K*1);
            scalar_interop.update_gpu_data(colors.colors_proj, 0);
            colors.update_proj_colors(data+N*M*K*4);
            scalar_interop.update_gpu_data(colors.colors_proj, 1);
            colors.update_proj_colors(data+N*M*K*7);
            scalar_interop.update_gpu_data(colors.colors_proj, 2);
            double maxi = 0.0;
            int iii = -1;
            for (int i = 0; i < N*M*K; i++)
            {
                double val = data[i+N*M*K*1];//(log10(data[i+N*M*K*1]+1e-15)+15)/15;
                if (maxi < val)
                {
                    maxi = val;
                    iii = i;
                }
            }
            std::cout << "1: " << maxi << " " << iii << "\n";

            maxi = 0.0;
            iii = -1;
            for (int i = 0; i < N*M*K; i++)
            {
                double val = data[i+N*M*K*4];//(log10(data[i+N*M*K*4]+1e-15)+15)/15;
                if (maxi < val)
                {
                    maxi = val;
                    iii = i;
                }
            }
            std::cout << "3: " << maxi << " " << iii << "\n";

            maxi = 0.0;
            iii = -1;
            for (int i = 0; i < N*M*K; i++)
            {
                double val = data[i+N*M*K*7];//(log10(data[i+N*M*K*7]+1e-15)+15)/15;
                if (maxi < val)
                {
                    maxi = val;
                    iii = i;
                }
            }
            std::cout << "7: " << maxi << " " << iii << "\n\n";

            //viewer.view(vis.i_size*3, 0, 0);
            v.view(obj_lens, time);
        }

        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                for (int k = 0; k < K; k++)
                {
                    int x = i*M*K+j*K+k;
                    for (int m = 0; m < S; m++)
                    {
                        int ind = x+M*N*K*m;
                        n_k[m] = data[ind];
                    }
                    /*if (j < M/2)
                    {
                        for (int m = 0; m < S; m++)
                        {
                            int ind = x+M*N*K*m;
                            data[ind] = n_k[m];
                        }
                        continue;
                    }
                    else*/
                    {
                        //#pragma omp parallel for
                        for (int m = 0; m < S; m++)
                        {
            	            int ind = x+M*N*K*m;
                            data_new[ind] = ( L1(S,m,n_k)*0.5-n_k[m]*L2(S,m,n_k) )*(dt)+n_k[m];
                            if (data_new[ind] < 0.0) data_new[ind] = 0.0;
                        }
                    }
                }

        for (int gg = 0; gg < S; gg++)
        {
            dif.iteration(data_new+N*M*K*gg, data+N*M*K*gg, 1, (gg==0) );

            if (time % 10 == 0)
                adv.set_foam_vels(foam.velos[ccount]);
            adv.find_src();
            adv.iteration(data+N*M*K*gg, data_new+N*M*K*gg, 1, false);

            //adv.upd_downward_velo(gg);
            adv.find_src_reg();
            for (int kk = 0; kk < 3; kk++)
            {
                adv.iteration(data_new+N*M*K*gg, data+N*M*K*gg, 1, false);
                adv.iteration(data+N*M*K*gg, data_new+N*M*K*gg, 1, false);
            }

        }
        if (time % 10 == 0)
            ccount++;

        std::swap(data, data_new);
        //init_data(data, dd, N, M, K, S);

        if (time==TIME_MAX)
        {
            time=0;
            break;
        }
        time++;
        //printProgress((double)time/TIME_MAX);
    }
    std::cout <<"\nComputation time: "<< get_wall_time()-start<<"\n";
    return 0;
}
