#pragma once
#include <iostream>
#include <vector>
#include "../geometry/vector3d.h"

namespace solvers
{
    template <typename T>
    class Advection3d_reg
    {
    private:
        //N->x, M->z, K->y
        int N, M, K;
        T dd, dt;
        std::vector<T> * weights;
        std::vector<int> * inds;

        Vector3d<T> * vels_foam;
        std::vector<T> * weights_foam;
        std::vector<int> * inds_foam;

        std::pair<T, int> find_inds(T src_pos, int len);
        void init_data();
        
        void set_inds(std::pair<T, int> pair_h,
                      std::pair<T, int> pair_w,
                      std::pair<T, int> pair_d, int idx);
    public:
        void upd_downward_velo(int p);
        void find_src();
        void set_foam_vels(Vector3d<T> * vel);
        void find_foam(Vector3d<T> * cc, int n, Vector3d<T> st_point);
        void iteration(T * data, T * data_new);
        Advection3d_reg(int N, int M, int K, T dd, T dt):
                    N(N), M(M), K(K), dd(dd), dt(dt)
                    {
                        init_data();
                    };
        ~Advection3d_reg();
    };
}
