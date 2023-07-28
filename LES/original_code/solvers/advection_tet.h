#pragma once
#include <iostream>
#include <algorithm>
#include <math.h>
#include <numeric>
#include <omp.h>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include "../geometry/vector3d.h"
#include "../geometry/mesh.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;

template<typename F, typename I>
class SpaceHeter3d_tet
{
protected:
public:
    geometry_solver_utils<F, I> & utils;
    geometry_data<F, I> & data;
    I size;
    I v_size;
    F dt;
    std::vector<F> f;
    std::vector<std::vector<F>> u;
    
    std::vector<std::vector<F>> temp_x;
    std::vector<std::vector<F>> temp_y;
    std::vector<std::vector<F>> temp_z;

    std::vector<std::vector<std::vector<F>>> temp_face_vec;
    std::vector<std::vector<F>> temp_face;

    SpMat mat;
    I * num_dums;
    SpaceHeter3d_tet(geometry_data<F, I> & data, geometry_solver_utils<F, I> & utils, F dt) : data(data), utils(utils), dt(dt)
    {
        init_help_vectors();
        init_mats();
    }
    ~SpaceHeter3d_tet();

    void initialize(std::vector<F> &input, Vector3d<F> const & source);
    void init_help_vectors();
    void init_mats();

    void interpolate_to_faces(std::vector<F> const &centered, 
                              std::vector<std::vector<F>> &faced);

    void iteration_muscl(std::vector<F> &f,
                 std::vector<std::vector<F>> const &f_face,
                 std::vector<std::vector<std::vector<F>>> const &vel);
    void iteration_advection(std::vector<std::vector<F>> & vel, std::vector<F> & val);
    void iteration_diffusion(std::vector<F> & val);
};

