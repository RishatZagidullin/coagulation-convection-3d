#include "glad/glad.h"
#include "utils/utils.cuh"
#include "solvers/mpi_solver.h"
#include "geometry/projection_reg.h"
#include "gl/viewer.h"
#include "gl/grid_viewer.h"
#include "gl/cuda_gl.cuh"
#include "gl/animation_reg.h"
#include <time.h>
#include <sys/time.h>
#include "geometry/parse_foam.h"
#include "wrappers.h"
#include <mpi.h>

void smol_iteration_fft_c(solvers::mpi_solver<double> & solver, int S,
          double *& n_k, double *& L1, double *& L2,
          fftw_complex *& ub, fftw_complex *& vb, double dt,
          TCross_Parallel_v1 & crossed_kernel,
          fftw_plan *& plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inv)
{
    for (int i = 0; i < solver.N_; i++)
        for (int j = 0; j < solver.M_; j++)
            for (int k = 0; k < solver.K_; k++)
    {
        //double st = MPI_Wtime();
        int x = i*solver.M_*solver.K_+j*solver.K_+k;
        for (int m = 0; m < S; m++)
        {
            int ind = x+solver.M_*solver.N_*solver.K_*m;
            n_k[m] = solver.data[ind];
        }
        L2 = crossed_kernel.matvec(n_k);
        L1 = crossed_kernel.smol_conv_discrete(n_k, ub, vb, 
                                             plan_v, plan_u, plan_inv);

        for (int m = 0; m < S; m++)
        {
            int ind = x+solver.M_*solver.N_*solver.K_*m;
            solver.data_new[ind] = ( L1[m] * 0.5 -
                                 n_k[m] * L2[m]) * dt + n_k[m];
            if (solver.data_new[ind] < 0.0) solver.data_new[ind] = 0.0;
        }
        delete [] L2;
        delete [] L1;
    }
}

void smol_iteration_fft_g(solvers::mpi_solver<double> & solver, int S,
          double *& n_k, double *& L1, double *& L2,
          cublasHandle_t & handle, cufftHandle & plan, double dt,
          TCross_Parallel_v1 & crossed_kernel)
{
    dim3 bl(128);
    dim3 gr((S+bl.x-1)/bl.x);
    for (int i = 0; i < solver.N_; i++)
        for (int j = 0; j < solver.M_; j++)
            for (int k = 0; k < solver.K_; k++)
    {
        
        int x = i*solver.M_*solver.K_+j*solver.K_+k;
        for (int m = 0; m < S; m++)
        {
            int ind = x+solver.M_*solver.N_*solver.K_*m;
            n_k[m] = solver.data[ind];
        }

        crossed_kernel.matvec(n_k, handle, L2);
        crossed_kernel.smol_conv_discrete(n_k, plan, handle, L1);

        checkCudaErrors(cudaDeviceSynchronize());
        wrappers::calc_smoluch<<<gr, bl>>>(L1, L2, n_k, S, dt);
        checkCudaErrors(cudaDeviceSynchronize());

        for (int m = 0; m < S; m++)
        {
            int ind = x+solver.M_*solver.N_*solver.K_*m;
            solver.data_new[ind] = n_k[m];
            if (solver.data_new[ind] < 0.0) solver.data_new[ind] = 0.0;
        }
    }
}

//from no on: x (N) - progragation, z (M) - up/down

//hardcodes:
//animation_reg.h update_proj_colors N/6 and 0/2
//projection_reg.h init_data -3.0
//advection_reg.cpp find_foam +2.0 +2.0 +3.0
//advection_reg.cpp find_foam foam_dx,foam_dy...
//advection_reg.cpp upd_downward_velo N/6 and K/2
//diffusion_reg.cpp st_point.z
//diffusion_reg.cpp source location (0.,0.)
//think that's it but not sure

int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    double start = MPI_Wtime();
    int n_vort = 288000;
    parse_foam<100> foam("./grid", n_vort);

    grid_params proj_grid = grid_params(10, 15, 30,
                                 -2.0, -3.0, -2.0, 0.4);
    short SCR_W = 1600;
    short SCR_H = 1200;
    int N = proj_grid.h;
    int M = proj_grid.w;
    int K = proj_grid.d;
    double dd = proj_grid.dd;
    int S = 8;
    double dt = 0.001;
    //assert that size == X*Y*Z
    int X = 1;
    int Y = 1;
    int Z = 1;
    if (rank == 0)
    {
        std::cout << "SIZE: " << size << "\n";
        std::cout << "X: " << X << "\n";
        std::cout << "Y: " << Y << "\n";
        std::cout << "Z: " << Z << "\n";
    }

    //***************GPU STUFF START*******************
    double tolerance = 1e-4;    			
    TCross_Parallel_v1 crossed_kernel = wrappers::default_crossed_kernel(tolerance, S, 1);		

    double R_value = crossed_kernel.get_rank();
    double V_value = crossed_kernel.get_columns_number();
    std::cout << "R value: " << R_value << "\n";
    std::cout << "V value: " << V_value << "\n";

    cublasHandle_t handle;
    cublasCreate(&handle);
    int Rnum = crossed_kernel.get_rows_number();
    int n[] = { Rnum }; 
    int inembed[] = { 0 };
    int onembed[] = { 0 }; 
    cufftHandle plan;
    cufftPlanMany(&plan, 1, n, inembed, 1, Rnum,
                  onembed, 1, Rnum, CUFFT_Z2Z, R_value);

    double *n_k;
    gpuErrchk(cudaMallocManaged(&n_k, S*sizeof(double)));
    //n_k = new double [S];
    double *L1_res_g, *L2_res_g;
    gpuErrchk(cudaMallocManaged(&L1_res_g, S*sizeof(double)));
    gpuErrchk(cudaMallocManaged(&L2_res_g, S*sizeof(double)));
    //***************GPU STUFF END*********************

        //****************LOW RANK STUFF START*************
    fftw_complex *ub = (fftw_complex *) fftw_malloc(R_value * V_value * sizeof(fftw_complex));
    fftw_complex *vb = (fftw_complex *) fftw_malloc(R_value * V_value * sizeof(fftw_complex));
    fftw_plan * plan_v = (fftw_plan *) fftw_malloc(R_value * sizeof(fftw_plan));
    fftw_plan * plan_u = (fftw_plan *) fftw_malloc(R_value * sizeof(fftw_plan));
    fftw_plan * plan_inv = (fftw_plan *) fftw_malloc(R_value * sizeof(fftw_plan));
    for (int i = 0; i < R_value; i++)
    {
        plan_v[i] = fftw_plan_dft_1d(S, vb+i*S, vb+i*S, FFTW_FORWARD, FFTW_ESTIMATE);
        plan_u[i] = fftw_plan_dft_1d(S, ub+i*S, ub+i*S, FFTW_FORWARD, FFTW_ESTIMATE);
        plan_inv[i] = fftw_plan_dft_1d(S, ub+i*S, ub+i*S, FFTW_BACKWARD, FFTW_ESTIMATE);	
    }
    double *L1_res_c, *L2_res_c;
    //****************LOW RANK STUFF END***************


    MPI_Barrier(MPI_COMM_WORLD);
    double D = 0.01;
    solvers::mpi_solver<double> solver(foam, rank, dd, dt, D,
                                      N/X, M/Z, K/Y, S, X, Y, Z);

    double * data, * data_recv;
    if (rank == 0)
    {
        data = new double [N*M*K*S];
        data_recv = new double [N/X*M/Z*K/Y*S];
    }
    else
    {
        data = new double [N/X*M/Z*K/Y*S];
    }

    AnimationData_reg * colors;
    constexpr int num_objs = 1;
    gl::viewer<num_objs, gl::grid_viewer> * v;
    gl::interop<float> * scalar_interop;
    projection_reg * proj;
    if (rank == 0)
    {
        proj = new projection_reg(proj_grid);
        //true is to hide window
        v = new gl::viewer<num_objs, gl::grid_viewer>(SCR_W, SCR_H, true);
        colors = new AnimationData_reg(N,M,K);
        v->buffer<GL_ARRAY_BUFFER>(proj->max_v*6,proj->vertices,0);
        v->buffer<GL_ELEMENT_ARRAY_BUFFER>(proj->max_i*6, 
                                           proj->indices, 1);
        scalar_interop = new gl::interop<float>(&v->buffer_objects[0],
                                                proj->max_v);
    }

    //int time = 0;
    int TIME_MAX = 1000;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "Preprocessing time: " << MPI_Wtime()-start;
        std::cout << "\n";
    }
    start = MPI_Wtime();
    //int ccount = 0;
    for (int time = 0; time < TIME_MAX; time++)
    {
        if (time % 5 == 1)
        {

            if (rank == 0)
            {
                for (int s = 0; s < S; s++)
                    for (int i = 0; i < N/X; i++)
                        for (int j = 0; j < M/Z; j++)
                            for (int k = 0; k < K/Y; k++)
                {
                    int ind = k+1 + (j+1)*solver.K_ +
                              (i+1)*solver.K_*solver.M_ +
                              s*solver.K_*solver.M_*solver.N_;

                    data[k+j*K+i*K*M+s*K*M*N] = solver.data[ind];
                }
                int n = 0;
                for (int x = 0; x < X; x++)
                    for (int z = 0; z < Z; z++)
                        for (int y = 0; y < Y; y++)
                {
                    if (n > 0)
                    {
                        MPI_Recv(data_recv, N/X*M/Z*K/Y*S, MPI_DOUBLE,
                             n, n, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for (int s = 0; s < S; s++)
                            for (int i = 0; i < N/X; i++)
                                for (int j = 0; j < M/Z; j++)
                                    for (int k = 0; k < K/Y; k++)
                        {
                            int ind = k + y*K/Y + (j + z*M/Z) * K + (i + x*N/X) * M*K + s * N*M*K;
                            data[ind] = data_recv[k+j*K/Y+i*K/Y*M/Z+s*K/Y*M/Z*N/X];
                        }
                    }
                    n++;
                }
                colors->update_proj_colors(data+N*M*K*1);
                scalar_interop->update_gpu_data(colors->colors_proj, 0);
                colors->update_proj_colors(data+N*M*K*3);
                scalar_interop->update_gpu_data(colors->colors_proj, 1);
                colors->update_proj_colors(data+N*M*K*7);
                scalar_interop->update_gpu_data(colors->colors_proj, 2);

                int obj_lens [num_objs] = {proj->max_i*6};
                //do not pass the time argument 
                //if you do not want any images to be saved
                v->view(obj_lens, time);
            }
            else
            {
                for (int s = 0; s < S; s++)
                    for (int i = 0; i < N/X; i++)
                        for (int j = 0; j < M/Z; j++)
                            for (int k = 0; k < K/Y; k++)
                {
                    int ind = k+1 + (j+1)*solver.K_ +
                                    (i+1)*solver.K_*solver.M_ +
                                    s*solver.K_*solver.M_*solver.N_;

                    data[k+j*K/Y+i*K/Y*M/Z+s*K/Y*M/Z*N/X] = 
                                    solver.data[ind];
                }
                MPI_Send(data, N/X*M/Z*K/Y*S, MPI_DOUBLE, 0, rank,
                         MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        smol_iteration_fft_g(solver, S, n_k, L1_res_g, 
                           L2_res_g, handle, plan, dt, 
                           crossed_kernel);

        //smol_iteration_fft_c(solver, S, n_k, L1_res_c, 
        //                   L2_res_c, ub, vb, dt, crossed_kernel,
        //                   plan_v, plan_u, plan_inv);

        cudaDeviceSynchronize();
        solver.iteration(time);
        if (size>1)
            solver.exchange();

        //printProgress((double)time/TIME_MAX);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "Processing time: " << MPI_Wtime()-start << endl;
    }
    cudaFree(n_k);
    cudaFree(L1_res_g);
    cudaFree(L2_res_g);
    cublasDestroy(handle);

    for (int i = 0; i < R_value; i++)
    {
        fftw_destroy_plan(plan_v[i]);
        fftw_destroy_plan(plan_u[i]);
        fftw_destroy_plan(plan_inv[i]);	
    }
    fftw_free(vb);
    fftw_free(ub);
    fftw_free(plan_u);
    fftw_free(plan_v);
    fftw_free(plan_inv);

    MPI_Finalize();
    return 0;
}
