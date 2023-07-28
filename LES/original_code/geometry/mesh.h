#pragma once
#include "vector3d.h"
#include <cmath>
#include <vector>

#include <CGAL/centroid.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/centroid.h>
#include <CGAL/Mesh_vertex_base_3.h>
#include <CGAL/Mesh_cell_base_3.h>
#include <CGAL/make_mesh_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef CGAL::Point_3<K> Point_3;
typedef C3t3::Cell_handle Cell_handle;
typedef C3t3::Vertex_handle Vertex_handle;
typedef Tr::Finite_cells_iterator Cell_iterator;

using namespace CGAL::parameters;

class geometry3d {
public:
    Tr tet;
    int c_size;
    int v_size;
    geometry3d(std::string smin, std::string smout)
    {
        Polyhedron inside;
        std::ifstream input(smin);
        input >> inside;
        input.close();
        Polyhedron outside;
        input.open(smout);
        input >> outside;
        input.close();
        Mesh_domain domain(inside, outside);

        domain.detect_features();
        Mesh_criteria crits(edge_size = 0.14, facet_size = 0.14,
                            cell_size = 0.15);

        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, crits,
                                            no_lloyd(), odt(),
                                            no_perturb(), no_exude());
        this->tet = c3t3.triangulation();
        std::cout <<"Number of vertices: ";
        std::cout << tet.number_of_vertices() << "\n";
        std::cout <<"Cells in complex/finite cells: ";
        std::cout <<  c3t3.number_of_cells_in_complex()<< "/";
        std::cout << tet.number_of_finite_cells() << "\n";
        std::cout << "Facets in complex/finite facets: ";
        std::cout << c3t3.number_of_facets_in_complex()<< "/";
        std::cout << tet.number_of_finite_facets() << "\n";

        int integer = 0;
        for (auto it(tet.finite_vertices_begin());
                  it!=tet.finite_vertices_end(); it++)
        {
            double x = it->point().x();
            double y = it->point().y();
            double z = it->point().z()+1.7;
            if (x*x+y*y+z*z < 0.089)
                it->set_index(-1);
            else
                it->set_index(integer++);
        }
        v_size = integer;
        integer = 0;
        for (Cell_iterator it(tet.finite_cells_begin());
                       it!=tet.finite_cells_end(); it++)
            for (int i = 0; i < 4; i ++)
                it->neighbor(i)->set_subdomain_index(-1);

        for (auto it(tet.finite_cells_begin());
                  it!=tet.finite_cells_end(); it++)
        {
            bool bad_cell = false;
            for (int i = 0; i < 4; i++)
            {
                if (it->vertex(i)->index() == -1)
                {
                    bad_cell = true;
                    break;
                }
            }
            if (bad_cell) continue;
            it->set_subdomain_index(integer++);
        }
        c_size = integer;
    }
};

class geometry_vis {
public:
    float * vertices;
    unsigned int * indices;
    int v_size;
    int i_size;
  
    geometry_vis(geometry3d const & g)
    {
        v_size = g.v_size;
        i_size = g.tet.number_of_finite_facets();

        vertices = new float [v_size*4];
        int c = 0;
        for (auto it(g.tet.finite_vertices_begin());
                  it!=g.tet.finite_vertices_end(); it++)
        {
            if (it->index() == -1)
            {
                continue;
            }
            vertices[c*4] = it->point().z();
            vertices[c*4+1] = it->point().y();
            vertices[c*4+2] = it->point().x();
            vertices[c*4+3] = 0.0;
            c++;
        }
        indices = new unsigned int [i_size*3];
        c = 0;
        for (auto it(g.tet.finite_facets_begin());
                  it != g.tet.finite_facets_end(); it++)
        {
            bool bad_face = false;
            for (int i = 0; i < 4; i++) if (i!=it->second)
            {
                if (it->first->vertex(i)->index() == -1)
                {
                    bad_face = true;
                    break;
                }
            }
            if (bad_face) continue;
            for (int i = 0; i < 4; i++) if (i!=it->second)
            {
                indices[c] = it->first->vertex(i)->index();
                c++;
            }
        }
        i_size = (int) (c/3);
        indices = (unsigned int*) realloc(indices, i_size*3*sizeof(unsigned int));
        return;
    }
    ~geometry_vis()
    {
        delete [] vertices;
        delete [] indices;
    }
};

template<typename F, typename I>
class geometry_data {
public:
    I c_size;
    I v_size;
    //inner cell has 4 facet neighbors
    I * cell_neighbor_ids;
    //cell's stencil neighbors share a facet, line or point
    std::vector<I> * cell_stencil_ids;
    //cell has 4 facets, here are it's areas
    F * tr_areas;
    //cell has a volume
    F * tet_volume;
    //normals point out from 4 facets
    Vector3d<F> * normals;
    //facet centroids
    Vector3d<F> * facet_centroids;
    //cell centroids
    Vector3d<F> * cell_centroids;

    Vector3d<F> * vertices;
    std::vector<I> * adjacent_verts_ids;
    std::vector<I> * incident_cells;

    geometry_data(geometry3d const & g);
    ~geometry_data();
};

template<typename F, typename I>
class geometry_solver_utils {
public:
    I size;

    F * beta1_p;
    F * beta2_p;
    F * beta3_p;
    F * beta1_m;
    F * beta2_m;
    F * beta3_m;

    I * a_p;
    I * b_p;
    I * c_p;
    I * a_m;
    I * b_m;
    I * c_m;

    F * dist_center_face;
    F * dist_inter_p;
    F * dist_inter_m;

    std::vector<I> * sorted_ids;

    geometry_solver_utils(geometry_data<F, I> const &g);
    ~geometry_solver_utils();
};
