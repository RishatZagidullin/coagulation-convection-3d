#pragma once
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "cuda_gl.cuh"	

namespace gl
{
    __global__ void section(unsigned int * i_out, unsigned int * inds, int i_size,
                            float * v_out, float * verts, int v_size, float val)
    {
        int numThreads = blockDim.x * gridDim.x;
        int global_id = threadIdx.x + blockIdx.x * blockDim.x;

        for (int id = global_id; id < i_size; id+=numThreads)
        {
            int ind = inds[3*id];
            int ind1 = inds[3*id+1];
            int ind2 = inds[3*id+2];
            float p = verts[4*ind+2];
            float p1 = verts[4*ind1+2];
            float p2 = verts[4*ind2+2];
            if ((p>val && p1>val && p2>val) || (p<val && p1<val && p2<val))
            {
                i_out[3*id+1] = ind;
                i_out[3*id+2] = ind;
            }
        }
        return;
    }

    __global__ void find_inds(unsigned int * output, unsigned int * inds, float * verts, 
                              int i_size, int v_size, int lim0, int lim1, int lim2)
    {
        int numThreads = blockDim.x * gridDim.x;
        int global_id = threadIdx.x + blockIdx.x * blockDim.x;

        for (int id = global_id; id < i_size; id+=numThreads)
        {
            int ind = inds[3*id];
            int ind1 = inds[3*id+1];
            int ind2 = inds[3*id+2];

            /*if ((verts[4*ind] > lim0+0.1 || verts[4*ind] < lim0-0.1) &&
                (verts[4*ind1]  > lim0+0.1 || verts[4*ind1]  < lim0-0.1) &&
                (verts[4*ind2]  > lim0+0.1 || verts[4*ind2]  < lim0-0.1))
            {
                output[3*id+1] = ind;
                output[3*id+2] = ind;
            }
            else if ((verts[4*ind+1] > lim1+0.1 || verts[4*ind+1] < lim1-0.1) &&
                (verts[4*ind1+1]  > lim1+0.1 || verts[4*ind1+1]  < lim1-0.1) &&
                (verts[4*ind2+1]  > lim1+0.1 || verts[4*ind2+1]  < lim1-0.1))
            {
                output[3*id+1] = ind;
                output[3*id+2] = ind;
            }*/
            /*else */if ((verts[4*ind+2] > lim2+0.1 || verts[4*ind+2] < lim2-0.1) &&
                (verts[4*ind1+2]  > lim2+0.1 || verts[4*ind1+2]  < lim2-0.1) &&
                (verts[4*ind2+2]  > lim2+0.1 || verts[4*ind2+2]  < lim2-0.1))
            {
                output[3*id+1] = ind;
                output[3*id+2] = ind;
            }
        }
        return;
    }

    class slicer_gl
    {
        private:
        cuda_gl<float> vert_gl;
        cuda_gl<unsigned int> ind_gl;

        public:
        slicer_gl(unsigned int * VBO, int vsize, float * v_data,
                  unsigned int * EBO, int isize, unsigned int * i_data)
                  : vert_gl{VBO, vsize}, ind_gl{EBO, isize}
        {
            //CHANGE TO CUDAMEMCPY
            #pragma omp parallel for
            for (int i = 0; i < vsize; i++)
                vert_gl.data_device[i] = v_data[i];
            #pragma omp parallel for
            for (int i = 0; i < isize; i++)
                ind_gl.data_device[i] = i_data[i];
        }

        void update_gpu_data(float lim0, float lim1, float lim2)
        {
            section<<<ind_gl.grid, ind_gl.block>>>(ind_gl.gl_data, 
                                                 ind_gl.data_device,
                                                 ind_gl.size/3,
                                                 vert_gl.gl_data,
                                                 vert_gl.data_device, 
                                                 vert_gl.size/4,
                                                 lim2);
            return;
        }
    };
}
