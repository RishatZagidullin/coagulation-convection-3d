#include "advection_tet.h"

template<typename F>
F sigma(F const &r, F const &dr, int const &PML, int const &size)
{
    if (r>(dr*(size-PML))) return pow((r-(size-PML)*dr)/(PML*dr),2)*3.0*log(10.0)*13.0/(PML*dr);
    //&& dr!=0.05 must be temporary, remove later
    else if (r < dr*(PML) && dr!=0.05) return pow((PML*dr-r)/(PML*dr),2)*3.0*log(10.0)*13.0/(PML*dr);
    else return 0;
}

template<typename F>
F limiter(F const &r_factor, F const & n_plus, F const &n_minus)
{
    F phi_superbee = std::max(std::min(n_minus*r_factor, (F)1.0),
                              std::min(r_factor, n_plus));
    phi_superbee = std::max(phi_superbee, (F)0.0);
    return phi_superbee;
}

template<typename F>
void find_gradients_eigen(std::vector<std::vector<F>> &res,
                          std::vector<F> const &u,
                          int const * const &num_dums,
                          Vector3d<F> const * const cell_centroids,
                          std::vector<int> const * const neighbor_ids)
{
    #pragma omp parallel for
    for (int i = 0; i < u.size(); i++)
    {
        int s = neighbor_ids[i].size();
        Eigen::VectorXd T(s-num_dums[i]);
        Eigen::MatrixXd d(s-num_dums[i], 3);
        int val = 0;
        for (int j = 0; j < s; j++)
        {
            if ((neighbor_ids[i][j]!=-1)&&(neighbor_ids[i][j]!=i))
            {
                d(val,0) = (cell_centroids[neighbor_ids[i][j]].x
                          - cell_centroids[i].x);
                d(val,1) = (cell_centroids[neighbor_ids[i][j]].y
                          - cell_centroids[i].y);
                d(val,2) = (cell_centroids[neighbor_ids[i][j]].z
                          - cell_centroids[i].z);
                F w = 1.0/pow(d(val,0)*d(val,0)+d(val,1)*d(val,1)
                             +d(val,2)*d(val,2), 0.5);
                d(val,0) *= w;
                d(val,1) *= w;
                d(val,2) *= w;
                T[val] = (u[neighbor_ids[i][j]] - u[i])*w;
                //if (abs(T[val]) < 1e-4) T[val] = 0.0;
                val++;
            }

            //Eigen::Vector3d output = mats[i].bdcSvd(Eigen::ComputeThinU 
            //                     | Eigen::ComputeThinV).solve(T);
            //Eigen::Vector3d output = mats[i].householderQr().solve(T);
            Eigen::Vector3d output = (d.transpose() * d)
                     .llt().solve(d.transpose() * T);
            res[0][i] = output(0);
            res[1][i] = output(1);
            res[2][i] = output(2);
        }
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::initialize(std::vector<F> &val, Vector3d<F> const  &source)
{
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        F x_sq = (data.cell_centroids[i].x -source.x) * 
                 (data.cell_centroids[i].x-source.x);
        F y_sq = (data.cell_centroids[i].y - source.y) * 
                 (data.cell_centroids[i].y-source.y);
        F z_sq = (data.cell_centroids[i].z - source.z) * 
                 (data.cell_centroids[i].z-source.z);
        val[i] = 1.0*exp(-4.0*(z_sq+y_sq+x_sq));
        //F x = data.cell_centroids[i].x;
        //F y = data.cell_centroids[i].y;
        //F z = data.cell_centroids[i].z;
        //val[i] = 2.*exp(-(z*z+y*y+x*x))*(2.*(z*z+y*y+x*x)-3);
        //val[i] = 20.0*(20*x*x+20*y*y+20*z*z-3)*exp(-10.0*(z*z+y*y+x*x));
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::interpolate_to_faces(
                          std::vector<F> const &centered, 
                          std::vector<std::vector<F>> &faced)
{
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (data.cell_neighbor_ids[i*4+j]!=-1)
            {
                F result_p_x = utils.beta1_p[i*4+j] * 
                               centered[data.cell_stencil_ids[i]
                               [utils.sorted_ids[i*4+j][utils.a_p[i*4+j]]]]
                             + utils.beta2_p[i*4+j] * 
                               centered[data.cell_stencil_ids[i]
                               [utils.sorted_ids[i*4+j][utils.b_p[i*4+j]]]]
                            + utils.beta3_p[i*4+j] * 
                               centered[data.cell_stencil_ids[i]
                               [utils.sorted_ids[i*4+j][utils.c_p[i*4+j]]]];

                if (fabs(utils.beta1_p[i*4+j] + utils.beta2_p[i*4+j] + utils.beta3_p[i*4+j] - 1.0) > 1e-4)
                    result_p_x = 0.0;

                F result_m_x = utils.beta1_m[i*4+j] * 
                               centered[data.cell_stencil_ids[i]
                               [utils.sorted_ids[i*4+j][utils.a_m[i*4+j]]]]
                             + utils.beta2_m[i*4+j] * 
                               centered[data.cell_stencil_ids[i]
                               [utils.sorted_ids[i*4+j][utils.b_m[i*4+j]]]]
                             + utils.beta3_m[i*4+j] * 
                               centered[data.cell_stencil_ids[i]
                               [utils.sorted_ids[i*4+j][utils.c_m[i*4+j]]]];

                if (fabs(utils.beta1_m[i*4+j] + utils.beta2_m[i*4+j] + utils.beta3_m[i*4+j] - 1.0) > 1e-4)
                    result_m_x = 0.0;

                F slope_plus_x = (result_p_x - centered[i])
                             /utils.dist_inter_p[i*4+j];
                F slope_minus_x = (centered[i] - result_m_x)
                             /utils.dist_inter_m[i*4+j];
                F n_plus = utils.dist_inter_p[i*4+j]
                             /utils.dist_center_face[i*4+j];
                F n_minus = utils.dist_inter_m[i*4+j]
                             /utils.dist_center_face[i*4+j];
                F small_cor = 1e-8;
                faced[i][j] = -(centered[i] + slope_plus_x * 
                                limiter((slope_minus_x+small_cor)/(slope_plus_x+small_cor), n_plus, n_minus) * 
                                utils.dist_center_face[i*4+j]);
                if (result_m_x == 0.0 || result_p_x == 0.0)
                    faced[i][j] = 0.0;
                if (std::isnan(faced[i][j]))
                {
                    std::cout << i << " is nan" << std::endl;
                    std::cout << utils.beta1_m[i*4+j] << " ";
                    std::cout << utils.beta2_m[i*4+j] << " ";
                    std::cout << utils.beta3_m[i*4+j] << std::endl;
                    std::cout << utils.beta1_p[i*4+j] << " ";
                    std::cout << utils.beta2_p[i*4+j] << " ";
                    std::cout << utils.beta3_p[i*4+j] << std::endl;
                }
            }
            else
                faced[i][j] = 0.0;
        }
    }
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::iteration_advection(std::vector<std::vector<F>> & vel, std::vector<F> & val)
{                       
    for (int j = 0; j < 3; j++)
        interpolate_to_faces(vel[j], temp_face_vec[j]);
    interpolate_to_faces(val, temp_face);  
    iteration_muscl(val, temp_face, temp_face_vec);
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::iteration_diffusion(std::vector<F> & val)
{
    find_gradients_eigen(temp_x, val, num_dums, data.cell_centroids, data.cell_stencil_ids);
    iteration_advection(temp_x, val);
}

bool is_vert_boundary(double x, double y, double z){
    if (fabs(x+1) < 1e-2 || fabs(y+1) < 1e-2 || fabs(z+3) < 1e-2)
        return true;
    if (fabs(x-1) < 1e-2 || fabs(y-1) < 1e-2 || fabs(z-2) < 1e-2)
        return true;
    return false;
}

bool is_on_circle(double x, double y, double z){
    float x2=x*x, y2=y*y, z2=(z+1.7)*(z+1.7);
    if (fabs(x2+y2+z2 - 0.09) < 1e-2)
    {
        return true;
    }
    else return false;
}


template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::iteration_muscl(std::vector<F> &f, 
                std::vector<std::vector<F>> const &f_face,
                std::vector<std::vector<std::vector<F>>> const &vel)
{
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        F temp = 0.0;
        bool boundary = false;

        for (int j = 0; j < 4; j++)
        {
            if (data.cell_neighbor_ids[i*4+j]==-1)
            {
                boundary = true;
                continue;
            }
            F sc_pr = (vel[0][i][j]*data.normals[i*4+j].x 
                    + vel[1][i][j]*data.normals[i*4+j].y 
                    + vel[2][i][j]*data.normals[i*4+j].z);
            int neigh_id = -1;
            for (int m = 0; m < 4; m ++)
            {
	        if (data.cell_neighbor_ids[data.cell_neighbor_ids[i*4+j]*4+m] == i)
                {
                    neigh_id = m;
                    break;
                }
            }
            if(neigh_id != -1) 
                 temp -= data.tr_areas[i*4+j]*(std::max(sc_pr, 
                    (F) 0.0) * f_face[i][j] 
                 + std::min(sc_pr,
                    (F) 0.0) * 
                    f_face[data.cell_neighbor_ids[i*4+j]][neigh_id]);
            else
            {
                std::cout << "THIS HAPPENED\n";
                temp -= data.tr_areas[i*4+j] * f_face[i][j];
            }
        }
        temp/=data.tet_volume[i];
        if (std::isnan(f[i]))
        {
            std::cout << i << " is nan" << std::endl;
            std::cout << data.cell_centroids[i].x << " " 
                      << data.cell_centroids[i].y << " " 
                      << data.cell_centroids[i].z << std::endl;
        }
        if (fabs(f[i]) > 150)
        {
            std::cout << i << " " << f[i] << " f is bigger than 150\n";
            std::cout << data.cell_centroids[i].x << " "
                      << data.cell_centroids[i].y << " " 
                      << data.cell_centroids[i].z << std::endl;
        }
        f[i] -= dt*(temp + 1./3.*f[i] * sigma(data.cell_centroids[i].z+3.0, (F) 0.05, 10, 100) +
                           1./3.*f[i] * sigma(data.cell_centroids[i].y+1.0, (F) 0.02, 10, 100) +
                           1./3.*f[i] * sigma(data.cell_centroids[i].x+1.0, (F) 0.02, 10, 100)
                   )/
                   (1.+1./3.*dt*(sigma(data.cell_centroids[i].z+3.0, (F) 0.05, 10, 100) +
                                 sigma(data.cell_centroids[i].y+1.0, (F) 0.02, 10, 100) +
                                 sigma(data.cell_centroids[i].x+1.0, (F) 0.02, 10, 100)
                                 )
                   );
    }
}

double W(double v)
{
    return v < 0.2 ? 1. - pow(v/0.2, 0.4) : 0.0;
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::init_mats()
{
    int siz = v_size;
    mat = SpMat(siz, siz);
    std::vector<Triplet> tripletList;
    tripletList.reserve(100*siz);
    
    std::vector<Triplet> small_triplets[4];

    #pragma omp parallel for
    for (int worker = 0; worker<4; worker++)
    {
        small_triplets[worker].reserve(25*siz);

    for (int i = (int) (siz/4.*worker); i < (int) (siz*((worker+1)/4.)); i++)
    {
        bool is_boundary = is_vert_boundary(data.vertices[i].x, data.vertices[i].y, data.vertices[i].z);

        if (!is_boundary && !is_on_circle(data.vertices[i].x, data.vertices[i].y, data.vertices[i].z))
        {
            Eigen::MatrixXd B_inv(3, 3);
            B_inv << 0.,0.,0.,0.,0.,0.,0.,0.,0.;
            Eigen::VectorXd o(3);
            o << 0.,0.,0.;
            Eigen::Vector3d xi(data.vertices[i].x, data.vertices[i].y, data.vertices[i].z);
            int s = data.adjacent_verts_ids[i].size();
            double total_psis = 0.0;
            for (int j = 0; j < s; j++)
            {
                I new_i = data.adjacent_verts_ids[i][j];
                if (new_i!=i && new_i!=-1)
                {
                    Eigen::Vector3d xj(data.vertices[new_i].x, 
                                       data.vertices[new_i].y, 
                                       data.vertices[new_i].z);
                    F psi = W((xj-xi).norm());
                    F omega = 0.0;
                    for (int k = 0; k < s; k++)
                    {
                        I fin_i = data.adjacent_verts_ids[i][k];
                        if (fin_i!=i && fin_i!=-1)
                        {
                            Eigen::Vector3d xk(data.vertices[fin_i].x, data.vertices[fin_i].y, data.vertices[fin_i].z);
                            omega+= W((xk-xi).norm());
                        }
                    }
                    psi /= omega;
                    total_psis += psi;
                    o += psi * (xj-xi);
                    B_inv +=  psi * (xj-xi)*(xj-xi).transpose();
                }
            }
            assert(fabs(total_psis-1.) < 1e-4);
            total_psis = 0.0;
            F denominator = 0;
            for (int j = 0; j < s; j++)
            {
                I new_i = data.adjacent_verts_ids[i][j];
                if (new_i!=i && new_i!=-1)
                {
                    Eigen::Vector3d xj(data.vertices[new_i].x, data.vertices[new_i].y, data.vertices[new_i].z);
                    F term = 1. - (xj-xi).transpose()
                                  * (B_inv.inverse()*o);
                    
                    F psi = W((xj-xi).norm());
                    F omega = 0.0; 
                    for (int k = 0; k < s; k++)
                    {
                        I fin_i = data.adjacent_verts_ids[i][k];
                        if (fin_i!=i && fin_i!=-1)
                        {
                            Eigen::Vector3d xk(data.vertices[fin_i].x, data.vertices[fin_i].y, data.vertices[fin_i].z);
                            omega+= W((xk-xi).norm());
                        }
                    }
                    psi /= omega;
                    total_psis += psi;
                    denominator += psi * (xj-xi).norm() *
                                   (xj-xi).norm() * term;
                }
            }
            assert(fabs(total_psis-1.) < 1e-4);
            total_psis = 0.0;
            F ith_val = 0;
            for (int j = 0; j < s; j++)
            {
                I new_i = data.adjacent_verts_ids[i][j];
                if (new_i!=i && new_i!=-1)
                {
                    Eigen::Vector3d xj(data.vertices[new_i].x, data.vertices[new_i].y, data.vertices[new_i].z);
                    F term = 1 - (xj-xi).transpose()
                                 *(B_inv.inverse()*o);
                    F psi = W((xj-xi).norm());
                    F omega = 0.0; 
                    for (int k = 0; k < s; k++)
                    {
                        I fin_i = data.adjacent_verts_ids[i][k];
                        if (fin_i!=i && fin_i!=-1)
                        {
                            Eigen::Vector3d xk(data.vertices[fin_i].x, data.vertices[fin_i].y, data.vertices[fin_i].z);
                            omega+= W((xk-xi).norm());
                        }
                    }
                    psi /= omega;
                    total_psis+=psi;
                    F value = 6*psi*term/denominator;
                    ith_val -= value;
                    small_triplets[worker].push_back(Triplet(i, new_i, value));
                }
            }
            assert(fabs(total_psis-1.) < 1e-4);
            small_triplets[worker].push_back(Triplet(i, i, ith_val));
        }
        else {
            small_triplets[worker].push_back(Triplet(i, i, 1.0));
        }
    }
    }
    for (int worker = 0; worker<4; worker++)
    {
        tripletList.insert(tripletList.end(), std::make_move_iterator(small_triplets[worker].begin()), 
                    std::make_move_iterator(small_triplets[worker].end()));
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    num_dums = new I [size];
    for (int i = 0; i < size; i++)
    {
        int s = data.cell_stencil_ids[i].size();
        int num_dummies = 0;
        for (int j = 0; j < s; j++) 
            if ((data.cell_stencil_ids[i][j]==-1)||
                (data.cell_stencil_ids[i][j]==i))
                num_dummies++;
        num_dums[i] = num_dummies;
    }
    return;
}

template<typename F, typename I>
void SpaceHeter3d_tet<F,I>::init_help_vectors()
{
    size = data.c_size;
    v_size = data.v_size;
    f.resize(size);
    
    u.resize(3);
    temp_x.resize(3);
    temp_y.resize(3);
    temp_z.resize(3);
    for (int i = 0; i < 3; i++) 
    {
        u[i] = std::vector<F>(size);
        temp_x[i] = std::vector<F>(size);
        temp_y[i] = std::vector<F>(size);
        temp_z[i] = std::vector<F>(size);
    }
    initialize(f, Vector3d<F>{0.0,0.0,-3.0});
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        F x = data.cell_centroids[i].x;
        F y = data.cell_centroids[i].y;
        //1.7 shift hardcoded to sphere position
        F z = data.cell_centroids[i].z+1.7;
        F r = pow(x*x+y*y+z*z, 0.5);
        F t = acos(z/r);
        F p = atan2(y, x);
        //0.3 is hardcoded to sphere radius
        F coef = 0.3;//0.4;
        F u_r = -coef*(1.-3./2.*0.3/r+1./2.*pow(0.3/r,3))*cos(t);
        F u_t = coef*(1.-3./4.*0.3/r-1./4.*pow(0.3/r,3))*sin(t);

        u[2][i] = u_r*cos(t) - r*sin(t)*u_t;
        if (u[2][i] > 0.0) std::cout << "vel is positive\n";
        u[0][i] = u_r*sin(t)*cos(p)+cos(t)*cos(p)*u_t;
        p = atan2(x, y);
        u[1][i] = u_r*sin(t)*cos(p)+cos(t)*cos(p)*u_t;
    }
    
    temp_face.resize(size);
    #pragma omp parallel for
    for (int i = 0; i < size; i++)
        temp_face[i] = std::vector<F>(4);

    temp_face_vec.resize(3);
    for (int i = 0; i < 3; i++)
    {
        temp_face_vec[i] = std::vector<std::vector<F>>(size);
        #pragma omp parallel for
        for (int j = 0; j < size; j++)
            temp_face_vec[i][j] = std::vector<F>(4);
    }
}

template<typename F, typename I>
SpaceHeter3d_tet<F,I>::~SpaceHeter3d_tet()
{
    std::cout << "eqn destructor\n";
    delete [] num_dums;
}

template class SpaceHeter3d_tet<double, int>;

