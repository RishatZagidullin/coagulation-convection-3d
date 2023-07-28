#pragma once
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

namespace gl
{
    template<typename T>
    class interop
    {
        private:
        cudaGraphicsResource * gl_buffer;    
        dim3 block;
        dim3 grid;

        public:
        int size;
        T * gl_data;
        T * data_device;
        void update_gpu_data(T * data, int ind);
        interop(unsigned int * BO, int size);
        ~interop();
    };
}
