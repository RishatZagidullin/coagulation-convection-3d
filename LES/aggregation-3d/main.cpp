#include "glad/glad.h"
#include "utils/utils.cuh"
#include "solvers/advection_tet.h"
#include "geometry/mesh.h"
#include "geometry/projection_reg.h"

#include "gl/viewer.h"
#include "gl/grid_viewer.h"
	
#include "gl/cuda_gl.cuh"
#include "gl/animation_tet.h"
#include <time.h>
#include <sys/time.h>

#include "geometry/parse_foam.h"

#include "wrappers.h"

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
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
    int n_vort = 288000;
    parse_foam<100> foam("./grid", n_vort);

    //domain_obj obj("offs/sphere.off");

    short SCR_W = 1600;
    short SCR_H = 1200;
    
    geometry3d geom("./offs/sphere.off", "./offs/boundary.off");
    geometry_data<double, int> grid(geom);
    geometry_vis vis(geom);
    geometry_solver_utils<double, int> solver_utils(grid);


   //grid_params proj_grid=grid_params(32,32,80,-1.0,-1.0,-3.0,0.0625);
    grid_params proj_grid = grid_params(
                                20, 30, 60, -2.0, -3.0, -2.0, 0.2);

    projection_reg proj(proj_grid);

    double dd = proj_grid.dd;
    double dt = 0.001;

    SpaceHeter3d_tet<double, int> eqn(grid, solver_utils, dt);
    eqn.find_inds(foam.cell_centers, n_vort);
    AnimationData_tet colors(proj, eqn);

    //true is to hide window
    constexpr int num_objs = 2;
    int obj_lens [num_objs] = {proj.max_i*6, vis.i_size*3};//{proj.max_i*6, obj.i_size*3, proj.max_i*2};
    gl::viewer<num_objs, gl::grid_viewer> v(SCR_W, SCR_H, true);
 
    v.buffer<GL_ARRAY_BUFFER>(proj.max_v*6, proj.vertices, 0);
    v.buffer<GL_ELEMENT_ARRAY_BUFFER>(proj.max_i*6, proj.indices, 1);

    v.buffer<GL_ARRAY_BUFFER>(vis.v_size*6, vis.vertices, 2);
    v.buffer<GL_ELEMENT_ARRAY_BUFFER>(vis.i_size*3, vis.indices, 3);




    //v.buffer<GL_ARRAY_BUFFER>(obj.v_size*3, obj.vertices, 2);
    //v.buffer<GL_ELEMENT_ARRAY_BUFFER>(obj.i_size*3, obj.indices, 3);
    //v.buffer<GL_ARRAY_BUFFER>(proj.max_i*6, proj.vect_verts, 4);
    //v.buffer<GL_ELEMENT_ARRAY_BUFFER>(proj.max_i*2,proj.vect_inds,5);

    gl::interop<float> scalar_interop(&v.buffer_objects[0],
                                      proj.max_v);

    gl::interop<float> shape_interop(&v.buffer_objects[2],
                               vis.v_size);

    int time = 0;
    int TIME_MAX = 5000;//200*8;//2000;
    std::cout << "Preprocessing time: " << get_wall_time()-start << "\n";
    start = get_wall_time();
    int ccount = 0;
    while (!glfwWindowShouldClose(v.window))
    {
        //do not pass the time argument 
        //if you do not want any images to be saved
        if (time % 50 == 0)
        {
            colors.update_mesh_colors();
            shape_interop.update_gpu_data(colors.colors_mesh, 0);

            colors.update_proj_colors();
            scalar_interop.update_gpu_data(colors.colors_proj, 0);

            vis.upd_indices(geom, eqn.f, eqn.data.incident_cells);
            obj_lens[1] = vis.i_size*3;
            //std::cout << "obj len: " << obj_lens[0] << "\n";
            v.buffer<GL_ELEMENT_ARRAY_BUFFER>(vis.i_size*3, 
                                              vis.indices, 3);

            v.view(obj_lens, time);
        }

        if (time % 50 == 0)
        {
            eqn.find_vels(foam.velos[ccount]);
            ccount++;
            //if (time % 50 == 0)
                eqn.initialize_again(eqn.f, Vector3d<double>{0.0,-3.0,0.0});
        }
        //eqn.iteration_diffusion(eqn.f);
        eqn.iteration_advection(eqn.u, eqn.f);

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
