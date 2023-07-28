#include "glad/glad.h"
#include "utils/utils.h"
#include "geometry/mesh.h"
#include "geometry/polygons.h"
#include "gl/viewer.h"
#include "gl/grid_viewer.h"

int main(int argc, char ** argv)
{
    double start = get_wall_time();

    short SCR_W = 1600;
    short SCR_H = 1200;
    
    
    polyhedron_factory pf;
    auto polyhedron = pf.gen_polygon({&sphere_tet1, &sphere_tet2, &curve_off});
    geometry3d geom(polyhedron);
    geometry_vis vis(geom);

    constexpr int num_objs = 1;
    int obj_lens [num_objs] = {vis.i_size*3};
    gl::viewer<num_objs, gl::grid_viewer> v(SCR_W, SCR_H, false);
 
    v.buffer<GL_ARRAY_BUFFER>(vis.v_size*6, vis.vertices, 0);
    v.buffer<GL_ELEMENT_ARRAY_BUFFER>(vis.i_size*3, vis.indices, 1);

    int time = 0;
    int TIME_MAX = 5000;
    std::cout << "Preprocessing time: " << get_wall_time()-start << "\n";
    start = get_wall_time();

    while (!glfwWindowShouldClose(v.window))
    {
        if (time % 50 == 0)
        {
            if (time == 0)
                v.view(obj_lens, time);
            else
                v.view(obj_lens);
        }

        //if (time==TIME_MAX)
        //{
        //    time=0;
        //    break;
        //}
        time++;
    }
    std::cout <<"Computation time: "<< get_wall_time()-start<<"\n";
    return 0;
}
