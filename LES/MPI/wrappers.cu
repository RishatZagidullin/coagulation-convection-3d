#include "wrappers.h"

namespace wrappers
{
    __global__ void calc_smoluch(const double *L1, const double *L2, 
                                 double *n_k, int N, double dt)
    {
        const int numThreads = blockDim.x * gridDim.x;
        const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        for (int i = threadID; i < N; i += numThreads)
        {
            //mod is cause you don't call bubble from smol_conv_discrete
            double temp = n_k[i];
            n_k[i] = ( L1[(N+i-1)%N] * 0.5 - temp * L2[i] ) * dt + temp;
        }
    }
}

