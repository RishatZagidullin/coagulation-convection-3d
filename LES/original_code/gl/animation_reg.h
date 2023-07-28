#pragma once
#include "../geometry/vector3d.h"
#include "../geometry/projection_reg.h"
#include "../solvers/advection_reg.h"

class AnimationData_reg {
public:
    float * colors_proj;
    int N, M, K;
    AnimationData_reg(int N, int M, int K) : N(N), M(M), K(K)
    {
       init_color_data();
    }
    void update_proj_colors(double const * data);
    ~AnimationData_reg();
private:
    void init_color_data();
};

void AnimationData_reg::update_proj_colors(double const * data)
{
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            colors_proj[i*M+j] = data[i*M*K+j*K+K/2];//(log10(data[i*M*K+j*K+K/2]+1e-15)+15)/15;
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < M; j++)
        for (int k = 0; k < K; k++)
            colors_proj[N*M+j*K+k] = data[N/4*M*K+j*K+k];//(log10(data[N/4*M*K+j*K+k]+1e-15)+15)/15;
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++)
        for (int k = 0; k < K; k++)
            colors_proj[N*M+M*K+i*K+k] = data[i*M*K+0/2*K+k];//(log10(data[i*M*K+0/2*K+k]+1e-15)+15)/15;
}

void AnimationData_reg::init_color_data()
{
    colors_proj = new float[N*M+M*K+N*K];
    return;
}

AnimationData_reg::~AnimationData_reg()
{
    delete [] colors_proj;
}
