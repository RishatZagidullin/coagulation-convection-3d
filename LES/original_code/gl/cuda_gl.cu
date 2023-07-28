#include "cuda_gl.cuh"
#include <iostream>
#include <algorithm>

namespace gl
{
    //all that needs to change is the function bellow and ctor args
    //to fully variate the interop behaviour I think.
    //The rest can stay the same.
    template<typename T>
    __global__ void map_to_gpu(T * output, T * input, int size, int ind)
    {
        int numThreads = blockDim.x * gridDim.x;
        int global_id = threadIdx.x + blockIdx.x * blockDim.x;

        for (int id = global_id; id < size; id+=numThreads)
        {
            output[6*(id)+3+ind] = input[id];
        }
        return;
    }

    template<typename T>
    void interop<T>::update_gpu_data(T * data, int ind)
    {
        cudaMemcpy(data_device, data, sizeof(T)*size,
                   cudaMemcpyHostToDevice);

        map_to_gpu<<<grid, block>>>(gl_data, data_device, size, ind);
        cudaDeviceSynchronize();
    }

    template<typename T>
    interop<T>::interop(unsigned int * BO, int size)
    {
        this->size = size;
        block = 256;
        grid = (size + block.x-1)/block.x;
        cudaMalloc(&data_device, sizeof(T)*size);
        cudaGraphicsGLRegisterBuffer(&gl_buffer, *BO,
                                 cudaGraphicsMapFlagsNone);
        cudaGraphicsMapResources(1, &gl_buffer, 0);
        cudaGraphicsResourceGetMappedPointer((void**)&gl_data,
                   nullptr, gl_buffer);
    }

    template<typename T>
    interop<T>::~interop()
    {
        cudaGraphicsUnmapResources(1, &gl_buffer, 0);
        cudaFree(data_device);
    }

    template class interop<float>;
}
