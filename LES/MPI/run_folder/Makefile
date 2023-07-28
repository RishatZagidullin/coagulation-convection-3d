GCC = mpicxx
NVCC = nvcc

LIB  = -lgomp -lm -lgmp -lmpfr -lm -ldl -lglfw -lfftw3 -lmpi -lmpi_cxx -lblas -llapacke

GCC_COPT = -fopenmp

VPATH = ../:../glad:../utils:../solvers:../gl:../geometry:../coagulation

SOURCES_SOL = advection_reg.cpp diffusion_reg.cpp vector3d.cpp wrappers.cpp
SOURCES_GL =  glad.c
SOURCES_MAIN = blas_double.cpp blas_int.cpp cross.cpp matrix.cpp maxvol.cpp parallel_cross_omp.cpp tensor.cpp tensor_train.cpp convolution.cpp main.cpp
SOURCES_CU = cuda_gl.cu parallel_cross_omp.cu convolution.cu wrappers.cu# main.cu

OBJECTS_SOL = $(SOURCES_SOL:%.cpp=$(OBJDIR)/%.o)
OBJECTS_GL = $(SOURCES_GL:%.c=$(OBJDIR)/%.o)
OBJECTS_MAIN = $(SOURCES_MAIN:%.cpp=$(OBJDIR)/%.o)
OBJECTS_CU = $(SOURCES_CU:%.cu=$(OBJDIR)/%_cu.o)

FFTW_INSTALL_PATH = /usr
FFTWCFLAGS=-I${FFTW_INSTALL_PATH}/include
FFTWLDFLAGS=-L${FFTW_INSTALL_PATH}/lib
FFTW_LIB = -lfftw3

CUDA_INSTALL_PATH = /usr/local/cuda-11.6
CUDACFLAGS=-I${CUDA_INSTALL_PATH}/include
CUDALDFLAGS=-L${CUDA_INSTALL_PATH}/lib64
CUDA_LIB = -lcufft -lcublas

MPI_INSTALL_PATH = /usr/lib/x86_64-linux-gnu/openmpi
MPICFLAGS=-I${MPI_INSTALL_PATH}/include
MPILDFLAGS=-L${MPI_INSTALL_PATH}/lib

COPT = -DCUDA_FFT -DFFTW

OBJDIR = build

all: tet.exe

tet.exe: $(OBJECTS_MAIN) $(OBJECTS_GL) $(OBJECTS_SOL) $(OBJECTS_CU)
	$(NVCC) $^ -o $@ $(LIB) $(FFTWLDFLAGS) $(FFTW_LIB) $(CUDA_LIB)

$(OBJDIR)/%.o: %.cpp
	$(GCC) -c $(CUDACFLAGS) $(COPT) $(FFTWCFLAGS) $(GCC_COPT) $< -o $@

$(OBJDIR)/%.o: %.c
	$(GCC) -c $< -o $@ $(CUDACFLAGS) $(COPT) $(FFTWCFLAGS) $(GCC_COPT)

$(OBJDIR)/%_cu.o: %.cu
	$(NVCC) -c $< -o $@ $(COPT) $(FFTWCFLAGS) $(MPICFLAGS)

clean:
	rm -rf $(OBJDIR)/*.o *.exe
