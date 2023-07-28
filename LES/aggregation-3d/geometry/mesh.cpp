#include "mesh.h"

template<typename F>
void find_plane_equation(Vector3d<F> const &normal,
                         Vector3d<F> const &point,
                         std::vector<F> &output)
{
    output.push_back(normal.x);
    output.push_back(normal.y);
    output.push_back(normal.z);
    output.push_back(
         -(normal.x*point.x+normal.y*point.y+normal.z*point.z)
    );
}

template <typename T, typename Compare>
std::vector<int> sort_permutation(const std::vector<T>& vec,
                                  Compare&& compare)
{
    std::vector<int> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
              [&](int &i, int &j){
                  return compare(vec[i], vec[j]);
              }
    );
    return p;
}

template <typename T>
void apply_permutation_in_place(std::vector<T>& vec,
                                const std::vector<int>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done[i]) continue;
        done[i] = true;
        int prev_j = i;
        int j = p[i];
        while (i != j)
        {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = p[j];
        }
    }
}

template<typename F, typename I>
void sort_neighbors(std::vector<Vector3d<F>> const &cell_centers, 
                    Vector3d<F> const &B, Vector3d<F> const &M,
                    std::vector<I> &output, I &crit)
{
    std::vector<F> values(cell_centers.size());
    Vector3d<F> a{M.x-B.x, M.y-B.y, M.z-B.z};
    F dist = pow(a.x*a.x+a.y*a.y+a.z*a.z, 0.5);
    a.x/=dist;
    a.y/=dist;
    a.z/=dist;
    #pragma omp parallel for
    for (int i = 0; i < values.size(); i++)
    {	
        Vector3d<F> b{B.x-cell_centers[i].x,
                      B.y-cell_centers[i].y,
                      B.z-cell_centers[i].z
        };
        F dist = pow(b.x*b.x+b.y*b.y+b.z*b.z, 0.5);
        b.x/=dist;
        b.y/=dist;
        b.z/=dist;
        values[i] = a.x*b.x + a.y*b.y + a.z*b.z;
    }
    //figure out the allocations timings
    output = sort_permutation(values,
                  [&](F const &a, F const &b){return a > b;});
    //std::cout << values[output[0]] << " ";
    //std::cout << values[output[1]] << " "
    //std::cout << values[output[2]] << "\n";
    for (int i = 0; i < values.size(); i++)
        if (values[output[i]] < 0.0)
        {
             crit = i;
            break;
        }
    //std::cout << "H: " << values[output[crit-1]] << " ";
    //std::cout << values[output[crit]] << " ";
    //std::cout << values[output[crit+1]]<< "\n\n";
}

//I will rely on one of the coordinates to be zero
template<typename F>
F sign (Vector3d<F> p1, Vector3d<F> p2, Vector3d<F> p3)
{
    //Vector3d first {(p1.x==0.0 ? p1.z:p1.x),
    //                (p1.y==0.0 ? p1.z:p1.y), 0.0}; 
    //Vector3d secon {(p2.x==0.0 ? p2.z:p1.x),
    //                (p2.y==0.0 ? p2.z:p2.y), 0.0};
    //Vector3d thirs {(p3.x==0.0 ? p3.z:p1.x),
    //                (p3.y==0.0 ? p3.z:p3.y), 0.0};
    return (p1.x - p3.x) * (p2.y - p3.y) -
           (p2.x - p3.x) * (p1.y - p3.y);
}

template<typename F>
bool if_in_tr_2d (Vector3d<F> pt, Vector3d<F> v1,
                  Vector3d<F> v2, Vector3d<F> v3)
{
    F d1, d2, d3;
    bool has_neg, has_pos;

    d1 = sign(pt, v1, v2);
    d2 = sign(pt, v2, v3);
    d3 = sign(pt, v3, v1);

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}

template<typename F>
bool if_inside_triangle(std::vector<Vector3d<F>> const &tr,
                        Vector3d<F> const &p)
{
    bool a = if_in_tr_2d(Vector3d<F>{p.x, p.y, 0.0}, 
                         Vector3d<F>{tr[0].x, tr[0].y, 0.0},
                         Vector3d<F>{tr[1].x, tr[1].y, 0.0},
                         Vector3d<F>{tr[2].x, tr[2].y, 0.0}); 
    bool b = if_in_tr_2d(Vector3d<F>{p.x, p.z, 0.0},
                         Vector3d<F>{tr[0].x, tr[0].z, 0.0},
                         Vector3d<F>{tr[1].x, tr[1].z, 0.0},
                         Vector3d<F>{tr[2].x, tr[2].z, 0.0}); 
    bool c = if_in_tr_2d(Vector3d<F>{p.z, p.y, 0.0},
                         Vector3d<F>{tr[0].z, tr[0].y, 0.0},
                         Vector3d<F>{tr[1].z, tr[1].y, 0.0},
                         Vector3d<F>{tr[2].z, tr[2].y, 0.0}); 
    return a && b && c;
}

template<typename F>
void find_basis_vectors(Vector3d<F> const &O, Vector3d<F> const &B, 
                        std::vector<F> const &plane,
                        std::vector<Vector3d<F>> &output)
{
    Vector3d<F> a{B.x-O.x, B.y-O.y, B.z-O.z};
    F dist = pow(a.x*a.x+a.y*a.y+a.z*a.z, 0.5);
    a.x/=dist;
    a.y/=dist;
    a.z/=dist;
    output.push_back(a);
    Vector3d<F> b{0.0, 0.0, 0.0};
    int counter = 0;
    if (plane[0] != 0.0) {
        b.y = 1.0;
        b.z = (a.x*plane[1]/plane[0] - a.y)/
              (a.z - a.x*plane[2]/plane[0]);
        b.x = -(plane[2]*b.z + plane[1])/plane[0];
    }
    else if (plane[1] != 0.0) {
        b.z = 1.0;
        b.x = (a.y*plane[2]/plane[1] - a.z)/
              (a.x - a.y*plane[0]/plane[1]);
        b.y = -(plane[0]*b.x + plane[2])/plane[1];
    }
    else if (plane[2] != 0.0) {
		b.x = 1.0;
		b.y = (a.z*plane[0]/plane[2] - a.x)/
		      (a.y - a.z*plane[1]/plane[2]);
		b.z = -(plane[1]*b.y + plane[0])/plane[2];
    }
    else std::cout << "\n\nTHIS SHOULD NEVER HAPPEN\n\n";
    dist = pow(b.x*b.x+b.y*b.y+b.z*b.z, 0.5);
    b.x/=dist;
    b.y/=dist;
    b.z/=dist;
    output.push_back(b);
}

template<typename F>
void find_projection(Vector3d<F> const &point,
                     std::vector<F> const &plane,
                     Vector3d<F> &output)
{
    F alpha = -(plane[3] + plane[0]*point.x + plane[1]*point.y + 
                plane[2]*point.z)/(plane[0]*plane[0] + 
                plane[1]*plane[1] + plane[2]*plane[2]);
    output.x = alpha*plane[0]+point.x;
    output.y = alpha*plane[1]+point.y;
    output.z = alpha*plane[2]+point.z;
}

template<typename F>
void find_intersection(Vector3d<F> const &point,
                       Vector3d<F> const &vec,
                       std::vector<F> const &plane,
                       Vector3d<F> &output)
{
    F alpha = -(plane[3] + plane[0]*point.x + plane[1]*point.y +
                plane[2]*point.z)/
               (plane[0]*vec.x + plane[1]*vec.y + plane[2]*vec.z);
    output.x = alpha*vec.x+point.x;
    output.y = alpha*vec.y+point.y;
    output.z = alpha*vec.z+point.z;
}

template<typename F, typename I>
void find_u_face(std::vector<I> const &sorted_ids, 
                std::vector<Vector3d<F>> const &stencil_cell_centers, 
                Vector3d<F> const &grid,
                Vector3d<F> const &face_centroid,
                I const &int1, I const &int2, I const &int3, 
                Vector3d<F> &intersection,
                std::vector<Vector3d<F>> &tri)
{
    Vector3d<F> Normal = face_centroid - grid;
    Vector3d<F> FirstCentroid = 
                       stencil_cell_centers[sorted_ids[int1]];
    Vector3d<F> SecondCentroidProj = 
                       stencil_cell_centers[sorted_ids[int2]];
    Vector3d<F> ThirdCentroidProj = 
	               stencil_cell_centers[sorted_ids[int3]];

    std::vector<F> plane_tr;
    //plane equation by three points
    F coef1_x = (SecondCentroidProj.y - FirstCentroid.y) *
                (ThirdCentroidProj.z - FirstCentroid.z);
    F coef2_x = (ThirdCentroidProj.y - FirstCentroid.y) *
                (SecondCentroidProj.z - FirstCentroid.z);
    F coef1_y = (SecondCentroidProj.z - FirstCentroid.z) *
                (ThirdCentroidProj.x - FirstCentroid.x);
    F coef2_y = (ThirdCentroidProj.z - FirstCentroid.z) * 
                (SecondCentroidProj.x - FirstCentroid.x);
    F coef1_z = (SecondCentroidProj.x - FirstCentroid.x) * 
                (ThirdCentroidProj.y - FirstCentroid.y);
    F coef2_z = (ThirdCentroidProj.x - FirstCentroid.x) * 
                (SecondCentroidProj.y - FirstCentroid.y);
    plane_tr.push_back(coef1_x-coef2_x);
    plane_tr.push_back(coef1_y-coef2_y);
    plane_tr.push_back(coef1_z-coef2_z);
    plane_tr.push_back(
              (coef2_x-coef1_x)*FirstCentroid.x +
              (coef2_y-coef1_y)*FirstCentroid.y +
              (coef2_z-coef1_z)*FirstCentroid.z
    );
    find_intersection(grid, Normal, plane_tr, intersection);
    tri.push_back(FirstCentroid);
    tri.push_back(SecondCentroidProj);
    tri.push_back(ThirdCentroidProj);
}

template<typename F>
F area(Vector3d<F> const &a, Vector3d<F> const &b,
       Vector3d<F> const &c)
{
    F val1 = pow((a.x-b.x) * (a.x-b.x) +
                 (a.y-b.y) * (a.y-b.y) +
                 (a.z-b.z) * (a.z-b.z), 0.5);
    F val2 = pow((b.x-c.x) * (b.x-c.x) +
                 (b.y-c.y) * (b.y-c.y) +
                 (b.z-c.z)*(b.z-c.z), 0.5);
    F res = 0.0;
    if (val1!=0.0 && val2!=0.0)
             res = 0.5 * pow(1.0 - pow(((a.x-b.x)*(b.x-c.x) + 
                                        (a.y-b.y)*(b.y-c.y) +
                                        (a.z-b.z)*(b.z-c.z))/
                                        (val1*val2), 2),
                             0.5) * val1 * val2;
    if (std::isnan(res))
        std::cout << "NAN INSIDE AREA FUNCTION\n";
    return res;
}

template<typename F, typename I>
void get_optimal_cells(std::vector<I> const &sorted_ids, 
        std::vector<Vector3d<F>> const &stencil_cell_centers,
        Vector3d<F> const &grid, Vector3d<F> const &face_centroid,
        std::vector<I> &output, Vector3d<F> &intersection,
        std::vector<Vector3d<F>> &tri, bool left_to_right, I crit)
{
    int i = left_to_right ? 0 : sorted_ids.size()-1;
    bool condition = false;
    int sign = left_to_right ? 1 : -1;
    std::vector<int> orig;
    while(output.size()!=3)
    {
        if (!(stencil_cell_centers[sorted_ids[i]] ==
              Vector3d<F>{100., 100., 100.}))
        {
            output.push_back(i);
            orig.push_back(i);
        }
        left_to_right ? i++ : i--;
    }
    bool first = true;
    do {
        std::vector<Vector3d<F>>().swap(tri);
        if (!first)
        {
            if (output[2]-output[1]==sign)
            {
                output[2]+=sign;
                output[1]=output[0]+sign;
            }
            else output[1]+=sign;
        }
        else first = false;
        Vector3d<F> dummy{100.,100.,100.};
        condition = 
            (stencil_cell_centers[sorted_ids[output[1]]] == dummy)
            || 
            (stencil_cell_centers[sorted_ids[output[2]]] == dummy);
        find_u_face(sorted_ids, stencil_cell_centers, grid, 
                    face_centroid, output[0], output[1],
                    output[2], intersection, tri);
        if (output[2]*sign==sign*crit)
        {
            output[0] = orig[0];
            output[1] = orig[1];
            output[2] = orig[2];
            std::vector<Vector3d<F>>().swap(tri);
            find_u_face(sorted_ids, stencil_cell_centers,
                        grid, face_centroid, output[0],
                        output[1], output[2], intersection, tri);
            break;
        }
    }
    while(condition || !if_inside_triangle(tri, intersection)
          /*&& (output[2]*sign < sign*crit )*/);
}

void recursive_push_adjacent_verts(std::vector<int> & adjacent_verts, Tr const &tet, Vertex_handle const cur_vert, int const &depth)
{
    if (depth == 0) return;

    std::vector<Vertex_handle> finite_adjacent_verts;
    tet.finite_adjacent_vertices(cur_vert, std::back_inserter(finite_adjacent_verts));
    for (int j = 0; j < finite_adjacent_verts.size(); j++)
    {
        auto adj_vert = finite_adjacent_verts[j];
        int ind = adj_vert->index();
        if (ind==-1) continue;
        bool is_copy = false;
        for (int k = 0; k < adjacent_verts.size(); k++)
        {
            if (adjacent_verts[k] == ind)
            {
                is_copy = true;
                break;
            }
        }
        if (!is_copy) adjacent_verts.push_back(ind);
        recursive_push_adjacent_verts(adjacent_verts, tet, adj_vert, depth-1);
    }
}

template<typename F, typename I>
geometry_data<F,I>::geometry_data(geometry3d const & g)
{
    c_size = g.c_size;;
    v_size = g.v_size;
    cell_neighbor_ids = new I [c_size*4];
    cell_stencil_ids = new std::vector<I> [c_size];
    tr_areas = new F [c_size*4];
    tet_volume = new F [c_size];
    normals = new Vector3d<F> [c_size*4];
    facet_centroids = new Vector3d<F> [c_size*4];
    cell_centroids = new Vector3d<F> [c_size];
    vertices = new Vector3d<F> [v_size];
    adjacent_verts_ids = new std::vector<I> [v_size];
    incident_cells = new std::vector<I> [v_size];

    int c = 0;
    for (auto it(g.tet.finite_vertices_begin());
              it!=g.tet.finite_vertices_end(); it++)
    {
        if (it->index() == -1) continue;
        recursive_push_adjacent_verts(adjacent_verts_ids[c], g.tet, it, 3);
        vertices[c] = Vector3d<F>{(F) it->point().x(), (F) it->point().y(), (F) it->point().z()};

        std::vector<Cell_handle> cells;
        g.tet.incident_cells(it, std::back_inserter(cells));
        for (int i = 0; i < cells.size(); i++)
        {
            int ind = cells[i]->subdomain_index();
            if (ind != -1) incident_cells[c].push_back(ind);
        }

        c++;
    }
    c = 0;
    for (Cell_iterator it(g.tet.finite_cells_begin());
         it!=g.tet.finite_cells_end(); it++)
    {
        if (it->subdomain_index() == -1) continue;
        auto point = CGAL::centroid(g.tet.tetrahedron(it));
        cell_centroids[c] = Vector3d<F>{(F)point.x(), 
                                 (F)point.y(), (F)point.z()};
		
        c++;
    }
    c = 0;
    for (Cell_iterator it(g.tet.finite_cells_begin());
         it!=g.tet.finite_cells_end(); it++)
    {
        if (it->subdomain_index() == -1) continue;
        std::vector<I> ids;
        for (int i = 0; i < 4; i ++)
        {
            I val = it->neighbor(i)->subdomain_index();
            ids.push_back(val);
            cell_neighbor_ids[c*4+i] = val;
        }
        for (int i = 0; i < 4; i++)
        {
            auto ver = it->vertex(i);
            std::vector<Cell_handle> finite_incident_cells;
            g.tet.finite_incident_cells(ver, 
                        std::back_inserter(finite_incident_cells));
            int id_size = ids.size();
            int size = finite_incident_cells.size();
            for (int j = 0; j < size; j++)
            {
                I val = finite_incident_cells[j]->subdomain_index();
                for (int k = 0; k < id_size; k++)
                {
                    if (val == ids[k] || val == c) break;
                    if (k == id_size-1) ids.push_back(val);
                }
            }
        }
        cell_stencil_ids[c] = ids;
        for (int i = 0; i < 4; i++)
        {
            F x0 = g.tet.triangle(it, i)[0].x();
            F x1 = g.tet.triangle(it, i)[1].x();
            F x2 = g.tet.triangle(it, i)[2].x();
            F y0 = g.tet.triangle(it, i)[0].y();
            F y1 = g.tet.triangle(it, i)[1].y();
            F y2 = g.tet.triangle(it, i)[2].y();
            F z0 = g.tet.triangle(it, i)[0].z();
            F z1 = g.tet.triangle(it, i)[1].z();
            F z2 = g.tet.triangle(it, i)[2].z();
            F A = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
            F B = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);
            F C = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);
            F coef = 1.0/pow(A*A+B*B+C*C,0.5);			

            F dist1 = pow(CGAL::centroid(g.tet.tetrahedron(it)).x()
                 -(-coef*A+CGAL::centroid(g.tet.triangle(it,i)).x()),2)
                 +pow(CGAL::centroid(g.tet.tetrahedron(it)).y()
                 -(-coef*B+CGAL::centroid(g.tet.triangle(it,i)).y()),2)
                 +pow(CGAL::centroid(g.tet.tetrahedron(it)).z()
                 -(-coef*C+CGAL::centroid(g.tet.triangle(it,i)).z()),2);
            F dist2 = pow(CGAL::centroid(g.tet.tetrahedron(it)).x()
                 -(coef*A+CGAL::centroid(g.tet.triangle(it,i)).x()),2)
                 +pow(CGAL::centroid(g.tet.tetrahedron(it)).y()
                 -(coef*B+CGAL::centroid(g.tet.triangle(it,i)).y()),2)
                 +pow(CGAL::centroid(g.tet.tetrahedron(it)).z()
                 -(coef*C+CGAL::centroid(g.tet.triangle(it,i)).z()),2);
            if (dist1>dist2) coef *= -1.0;
            normals[c*4+i] = Vector3d<F>{coef*A, coef*B, coef*C};
        }

        for (int i = 0; i < 4; i++)
        {
            auto point = CGAL::centroid(g.tet.triangle(it, i));
            tr_areas[c*4+i] = pow(g.tet.triangle(it, i).squared_area(), 
                          0.5);
            facet_centroids[c*4+i] = Vector3d<F>{(F)point.x(), 
                                     (F)point.y(), (F)point.z()};
        }
        tet_volume[c] = g.tet.tetrahedron(it).volume();
        c++;
    }
}

template<typename F, typename I>
geometry_data<F,I>::~geometry_data() 
{
    std::cout << "geometry data destructor\n";
    delete [] cell_neighbor_ids;
    delete [] cell_stencil_ids;
    delete [] tr_areas;
    delete [] tet_volume;
    delete [] normals;
    delete [] facet_centroids;
    delete [] cell_centroids;
    delete [] vertices;
    delete [] adjacent_verts_ids;
    delete [] incident_cells;
}

template class geometry_data<double,int>;

template<typename F, typename I>
geometry_solver_utils<F,I>::geometry_solver_utils(geometry_data<F,I> const &g)
{
    size = g.c_size;
    beta1_p = new F [size*4];
    beta2_p = new F [size*4];
    beta3_p = new F [size*4];
    beta1_m = new F [size*4];
    beta2_m = new F [size*4];
    beta3_m = new F [size*4];
    a_p = new I [size*4];
    b_p = new I [size*4];
    c_p = new I [size*4];
    a_m = new I [size*4];
    b_m = new I [size*4];
    c_m = new I [size*4];
    dist_center_face = new F [size*4];
    dist_inter_p = new F [size*4];
    dist_inter_m = new F [size*4];

    sorted_ids = new std::vector<I> [size*4];

    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        std::vector<Vector3d<F>> stencil_cell_centers(
                                  g.cell_stencil_ids[i].size());
        for (int k = 0; k < g.cell_stencil_ids[i].size(); k++)
        {
            if (g.cell_stencil_ids[i][k] !=-1)
                stencil_cell_centers[k] = 
                         g.cell_centroids[g.cell_stencil_ids[i][k]];
             else stencil_cell_centers[k] =
                Vector3d<F>{100., 100., 100.};
        }
        for (int j = 0; j < 4; j++)
        {
            if (g.cell_neighbor_ids[i*4+j]!=-1)
            {
                std::vector<I> s_ids;
                int crit;
                sort_neighbors(stencil_cell_centers,
                               g.cell_centroids[i],
                               g.facet_centroids[i*4+j],
                               s_ids, crit);
                sorted_ids[i*4+j] = s_ids;
                Vector3d<F> intersection_plus;
                std::vector<Vector3d<F>> tri_plus;
                std::vector<I> optimal_cells_p;
                get_optimal_cells(s_ids,
                         stencil_cell_centers, g.cell_centroids[i], 
                         g.facet_centroids[i*4+j], optimal_cells_p, 
                         intersection_plus, tri_plus, false, crit);
                F denom = area(tri_plus[0], tri_plus[1], tri_plus[2]);

		beta1_p[i*4+j] = area(intersection_plus,
		                     tri_plus[1], tri_plus[2])/denom;
                beta2_p[i*4+j] = area(intersection_plus,
                                     tri_plus[0], tri_plus[2])/denom;
                beta3_p[i*4+j] = area(intersection_plus,
                                     tri_plus[0], tri_plus[1])/denom;

                a_p[i*4+j] = optimal_cells_p[0];
                b_p[i*4+j] = optimal_cells_p[1];
                c_p[i*4+j] = optimal_cells_p[2];

                Vector3d<F> intersection_minus;
                std::vector<Vector3d<F>> tri_minus;
                std::vector<I> optimal_cells_m;
                get_optimal_cells(s_ids, stencil_cell_centers, 
                        g.cell_centroids[i], g.facet_centroids[i*4+j], 
                        optimal_cells_m, intersection_minus,
                        tri_minus, true, crit);

                denom = area(tri_minus[0], tri_minus[1], tri_minus[2]);
                a_m[i*4+j] = optimal_cells_m[0];
                b_m[i*4+j] = optimal_cells_m[1];
                c_m[i*4+j] = optimal_cells_m[2];

                beta1_m[i*4+j] = area(intersection_minus,
                                 tri_minus[1], tri_minus[2])/denom;
                beta2_m[i*4+j] = area(intersection_minus,
                                 tri_minus[0], tri_minus[2])/denom;
                beta3_m[i*4+j] = area(intersection_minus,
                                 tri_minus[0], tri_minus[1])/denom;

                dist_center_face[i*4+j] =
                          pow(pow(g.cell_centroids[i].x -
                                  g.facet_centroids[i*4+j].x, 2)
                            + pow(g.cell_centroids[i].y -
                                  g.facet_centroids[i*4+j].y, 2) 
                            + pow(g.cell_centroids[i].z -
                                  g.facet_centroids[i*4+j].z, 2),
                           0.5);
                dist_inter_p[i*4+j] = 
                          pow(pow(g.cell_centroids[i].x -
                                      intersection_plus.x, 2) 
                            + pow(g.cell_centroids[i].y -
                                      intersection_plus.y, 2) 
                            + pow(g.cell_centroids[i].z -
                                      intersection_plus.z, 2),
                            0.5);
                dist_inter_m[i*4+j] =
                          pow(pow(g.cell_centroids[i].x -
                                      intersection_minus.x, 2) 
                            + pow(g.cell_centroids[i].y -
                                      intersection_minus.y, 2) 
                            + pow(g.cell_centroids[i].z -
                                      intersection_minus.z, 2),
                             0.5);
            }
            else {}
                /*essentially you don't need to do anything,
                  but if main algo changes,
                  this might lead to problems
                */
        }
    }
}

template<typename F, typename I>
geometry_solver_utils<F,I>::~geometry_solver_utils() 
{
    std::cout << "solver utils destructor\n";
    delete [] beta1_p;
    delete [] beta2_p;
    delete [] beta3_p;
    delete [] beta1_m;
    delete [] beta2_m;
    delete [] beta3_m;
    delete [] a_p;
    delete [] b_p;
    delete [] c_p;
    delete [] a_m;
    delete [] b_m;
    delete [] c_m;
    delete [] dist_center_face;
    delete [] dist_inter_p;
    delete [] dist_inter_m;
    delete [] sorted_ids;
}

template class geometry_solver_utils<double,int>;
