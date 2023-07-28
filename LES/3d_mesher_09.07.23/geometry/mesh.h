#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
//#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
//typedef CGAL::Mesh_domain_with_polyline_features_3<CGAL::Polyhedral_mesh_domain_with_features_3<K> > Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

//typedef std::vector<K::Point_3> Polyline_3;
//typedef std::list<Polyline_3> Polylines;

using namespace CGAL::parameters;

class geometry3d {
public:
    Tr tet;
    int c_size;
    int v_size;
    geometry3d(Polyhedron & obj)
    {
        Polyhedron outside;
        std::ifstream input("offs/boundary.off");
        input >> outside;

        //Mesh_domain domain(outside);
        Mesh_domain domain(obj, outside);
        domain.detect_features();
        Mesh_criteria crits(edge_size = 0.2,
                         facet_size = 0.2,
                         cell_size = 0.2);

        /*Polylines polylines1 (1);
        Polyline_3& polyline1 = polylines1.front();
        int max_x = 40;
        int max_y = 20;
        for(int j = 0; j < max_x; ++j)
        {
            K::Point_3 p (j/(max_x-1.)*5-1, 0.5*sin(j/(max_x-1.)*5-1)-1, -1);
            polyline1.push_back(p);
        }
        for(int j = max_x-1; j >= 0; --j)
        {
            K::Point_3 p (j/(max_x-1.)*5-1, 0.5*sin(j/(max_x-1.)*5-1)-0.8, -1);
            polyline1.push_back(p);
        }
        polyline1.push_back(polyline1.front());
        domain.add_features(polylines1.begin(), polylines1.end());

        Polylines polylines2 (1);
        Polyline_3& polyline2 = polylines2.front();
        for(int j = 0; j < max_x; ++j)
        {
            K::Point_3 p (j/(max_x-1.)*5-1, 0.5*sin(j/(max_x-1.)*5-1)-1, 1);
            polyline2.push_back(p);
        }
        for(int j = max_x-1; j >= 0; --j)
        {
            K::Point_3 p (j/(max_x-1.)*5-1, 0.5*sin(j/(max_x-1.)*5-1)-0.8, 1);
            polyline2.push_back(p);
        }
        polyline2.push_back(polyline2.front());
        domain.add_features(polylines2.begin(), polylines2.end());

        Polylines polylines3 (1);
        Polyline_3& polyline3 = polylines3.front();
        for(int i = 0; i < max_y; ++i)
        {
            K::Point_3 p (0/(max_x-1.)*5-1, 0.5*sin(0/(max_x-1.)*5-1)-1, i/(max_y-1.)*2. - 1);
            polyline3.push_back(p);
        }
        for(int i = max_y-1; i >= 0; --i)
        {
            K::Point_3 p (0/(max_x-1.)*5-1, 0.5*sin(0/(max_x-1.)*5-1)-0.8, i/(max_y-1.)*2. - 1);
            polyline3.push_back(p);
        }
        polyline3.push_back(polyline3.front());
        domain.add_features(polylines3.begin(), polylines3.end());

        Polylines polylines4 (1);
        Polyline_3& polyline4 = polylines4.front();
        for(int i = 0; i < max_y; ++i)
        {
            K::Point_3 p ((max_x-1.)/(max_x-1.)*5-1, 0.5*sin((max_x-1.)/(max_x-1.)*5-1)-1, i/(max_y-1.)*2. - 1);
            polyline4.push_back(p);
        }
        for(int i = max_y-1; i >= 0; --i)
        {
            K::Point_3 p ((max_x-1.)/(max_x-1.)*5-1, 0.5*sin((max_x-1.)/(max_x-1.)*5-1)-0.8, i/(max_y-1.)*2. - 1);
            polyline4.push_back(p);
        }
        polyline4.push_back(polyline4.front());
        domain.add_features(polylines4.begin(), polylines4.end());*/

        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, crits,
                                            no_lloyd(), no_odt(),
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
            double z = it->point().z();
            if (x*x+(y-1.5)*(y-1.5)+z*z <= 1.04)
                it->set_index(integer++);
            else if ((x-3.)*(x-3.)+(y-1.5)*(y-1.5)+z*z <= 1.04)
                it->set_index(integer++);
            else if ((abs(y - (0.5*sin(x)-1.)) <= 0.04) && (z <= 1.04) && (z >= -1.04) && (x >= -1.04) && (x <= 4.04))
                it->set_index(integer++);
            else if ((abs(y - (0.5*sin(x)-0.8)) <= 0.04) && (z <= 1.04) && (z >= -1.04) && (x >= -1.04) && (x <= 4.04))
                it->set_index(integer++);
            else
                it->set_index(-1);
        }
        v_size = integer;
        integer = 0;
        for (auto it(tet.finite_cells_begin());
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
    /*int max_y = 20;
    int max_x = 40;
    v_size = max_x*max_y*2;
    i_size = 4*(max_x-1)*(max_y-1)+4*(max_x-1)+4*(max_y-1);
    vertices = new float [v_size*6];
    indices = new unsigned int [i_size*3];
    int c = 0;
    for (int i = 0; i < max_y; i++)
    {
        for (int j = 0; j < max_x; j++)
        {
            vertices[c*6] = j/(max_x-1.)*7-2;
            vertices[c*6+1] = 0.5*sin(j/(max_x-1.)*7-2)-1.;
            vertices[c*6+2] = i/(max_y-1.)*4. - 2;
            vertices[c*6+3] = 0.0;
            vertices[c*6+4] = 0.0;
            vertices[c*6+5] = 0.0;
            c++;
        }
    }
    for (int i = 0; i < max_y; i++)
    {
        for (int j = 0; j < max_x; j++)
        {
            vertices[c*6] = j/(max_x-1.)*7-2;
            vertices[c*6+1] = 0.5*sin(j/(max_x-1.)*7-2)-0.8;
            vertices[c*6+2] = i/(max_y-1.)*4. - 2;
            vertices[c*6+3] = 0.0;
            vertices[c*6+4] = 0.0;
            vertices[c*6+5] = 0.0;
            c++;
        }
    }
    c = 0;
    for (int i = 0; i < max_y-1; i++)
    {
        int offset = i*max_x;
        for (int j = 0; j < max_x-1; j++)
        {
            indices[6*c] = j+offset;
            indices[6*c+1] = j+1+offset;
            indices[6*c+2] = j+max_x+offset;
            indices[6*c+3] = j+max_x+offset;
            indices[6*c+4] = j+1+offset;
            indices[6*c+5] = j+1+max_x+offset;
            c++;
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        int offset = i*max_x+max_x*max_y;
        for (int j = 0; j < max_x-1; j++)
        {
            indices[6*c] = j+1+offset;
            indices[6*c+1] = j+offset;
            indices[6*c+2] = j+max_x+offset;
            indices[6*c+3] = j+1+offset;
            indices[6*c+4] = j+max_x+offset;
            indices[6*c+5] = j+1+max_x+offset;
            c++;
        }
    }
    for (int i = 0; i < max_y-1; i++)
    {
        indices[6*c] = i*max_x;
        indices[6*c+1] = i*max_x+max_x;
        indices[6*c+2] = i*max_x+max_x*max_y;
        indices[6*c+3] = i*max_x+max_x*max_y;
        indices[6*c+4] = i*max_x+max_x;
        indices[6*c+5] = i*max_x+max_x*max_y+max_x;
        c++;
        indices[6*c] = i*max_x+max_x+max_x-1;
        indices[6*c+1] = i*max_x+max_x-1;
        indices[6*c+2] = i*max_x+max_x-1+max_x*max_y;
        indices[6*c+3] = i*max_x+max_x+max_x-1;
        indices[6*c+4] = i*max_x+max_x-1+max_x*max_y;
        indices[6*c+5] = i*max_x+max_x-1+max_x*max_y+max_x;
        c++;
    }
    for (int j = 0; j < max_x-1; j++)
    {
        indices[6*c] = j+1;
        indices[6*c+1] = j;
        indices[6*c+2] = j+max_x*max_y;
        indices[6*c+3] = j+1;
        indices[6*c+4] = j+max_x*max_y;
        indices[6*c+5] = j+max_x*max_y+1;
        c++;
        indices[6*c] = j+max_x*(max_y-1);
        indices[6*c+1] = j+1+max_x*(max_y-1);
        indices[6*c+2] = j+max_x*(max_y-1)+max_x*max_y;
        indices[6*c+3] = j+max_x*(max_y-1)+max_x*max_y;
        indices[6*c+4] = j+1+max_x*(max_y-1);
        indices[6*c+5] = j+max_x*(max_y-1)+max_x*max_y+1;
        c++;
    }*/


        v_size = g.v_size;
        i_size = g.tet.number_of_finite_facets();

        vertices = new float [v_size*6];
        int c = 0;
        for (auto it(g.tet.finite_vertices_begin());
                  it!=g.tet.finite_vertices_end(); it++)
        {
            if (it->index() == -1)
            {
                continue;
            }
            vertices[c*6] = it->point().x();
            vertices[c*6+1] = it->point().y();
            vertices[c*6+2] = it->point().z();
            vertices[c*6+3] = 0.0;
            vertices[c*6+4] = 0.0;
            vertices[c*6+5] = 0.0;
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
