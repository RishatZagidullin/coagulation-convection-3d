#pragma once
#include <initializer_list>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Surface_mesh_default_triangulation_3 Tr2d;
typedef CGAL::Complex_2_in_triangulation_3<Tr2d> C2t3;
typedef Tr2d::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

class polyhedron_factory {
public:
    Polyhedron gen_polygon(std::initializer_list<float (*)(Point_3)> funcs);
};

Polyhedron polyhedron_factory::gen_polygon(std::initializer_list<float (*)(Point_3)> funcs)
{
    std::vector<float (*)(Point_3)> funcs_v(funcs);
    Tr2d tr1;
    C2t3 c2t3_1(tr1);
    Surface_3 surface_1(funcs_v[0], Sphere_3(Point_3(0., 1.5, 0.), 2.));
    CGAL::Surface_mesh_default_criteria_3<Tr2d> criteria(30., 0.1, 0.1);
    CGAL::make_surface_mesh(c2t3_1, surface_1, criteria, CGAL::Non_manifold_tag());
    std::ofstream out("offs/temp1.off");
    CGAL::output_surface_facets_to_off(out, c2t3_1);
    out.close();
    
    Tr2d tr2;
    C2t3 c2t3_2(tr2);
    Surface_3 surface_2(funcs_v[1], Sphere_3(Point_3(3., 1.5, 0.), 2.));
    CGAL::make_surface_mesh(c2t3_2, surface_2, criteria, CGAL::Non_manifold_tag());
    out.open("offs/temp2.off");
    CGAL::output_surface_facets_to_off(out, c2t3_2);
    out.close();

    funcs_v[2](Point_3(0.,0.,0.));
    /*Tr2d tr3;
    C2t3 c2t3_3(tr2);
    Surface_3 surface_3(funcs_v[2], Sphere_3(CGAL::ORIGIN, 7.));
    CGAL::make_surface_mesh(c2t3_3, surface_3, criteria, CGAL::Non_manifold_tag());
    out.open("offs/temp3.off");
    CGAL::output_surface_facets_to_off(out, c2t3_3);
    out.close();*/

    Polyhedron polyhedron;
    std::ifstream in("offs/temp1.off");
    in >> polyhedron;
    in.close();
    in.open("offs/temp2.off");
    in >> polyhedron;
    in.close();
    in.open("offs/thin_foil.off");
    in >> polyhedron;
    in.close();
    return polyhedron;
}

