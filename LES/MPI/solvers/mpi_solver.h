#pragma once

#include "advection_reg.h"
#include "diffusion_reg.h"
#include "../geometry/parse_foam.h"
#include <mpi.h>

namespace solvers
{
    template <typename T>
    class mpi_solver
    {
    public:
        // size is equal to X*Y*Z
        int size;
        //ranking first along K, then M, then N
        int rank;
        int ccount = 0;
        //X <-> N
        //Z <-> M
        //Y <-> K
        int X, Y, Z;
        //sizes without ghost cells
        int N, M, K, S;
        int N_, M_, K_;

        T dd, dt, D;
        Vector3d<T> st_point = {-2.,-2.,-3.};

        bool * neighs;
        int * neigh_ids;
        int * sizes;
        //you will need to allocate this with ghost cells
        T * data;
        T * data_new;

        //length should be 6
        //1 - front(-) x, 2 - back(+) x
        //3 - left(-) y, 4 - right(+) y
        //5 - down(-) z, 6 - up(+) z
        std::vector<T> * buffer_send;
        std::vector<T> * buffer_recv;

        parse_foam<100> & foam;
        Diffusion3d_reg<T> dif;
        Advection3d_reg<T> adv;

        void transfer_in_buffer();
        void transfer_out_of_buffer();
        void calc_mpi();
        void iteration(int time);
        void exchange();

        mpi_solver(parse_foam<100> & foam, int rank, T dd, T dt, T D,
                   int N, int M, int K, int S, int X, int Y, int Z):
            foam(foam), rank(rank), dd(dd), dt(dt), D(D),
            N(N), M(M), K(K), S(S), X(X), Y(Y), Z(Z),
            N_(N+2), M_(M+2), K_(K+2), size(X*Y*Z),
            dif(D, N_, M_, K_, dd, dt),
            adv(N_, M_, K_, dd, dt)
        {
            calc_mpi();
        }
        ~mpi_solver();
    };


}

template<typename T>
solvers::mpi_solver<T>::~mpi_solver()
{
    delete [] data;
    delete [] data_new;
    delete [] neighs;
    delete [] neigh_ids;
    delete [] sizes;
    delete [] buffer_send;
    delete [] buffer_recv;
}

template <typename T>
void solvers::mpi_solver<T>::calc_mpi()
{
    st_point.x += (rank / (Y*Z))*N*dd;
    st_point.y += (rank%Y)*K*dd;
    st_point.z += ((rank%(Y*Z)) /Y)*M*dd;
    adv.find_foam(foam.cell_centers, foam.n_points, st_point);

    data = new T [K_*M_*N_*S];
    data_new = new T [K_*M_*N_*S];

    neighs = new bool [6];
    neigh_ids = new int [6];
    sizes = new int [6];

    //front
    if (rank>=Y*Z) neighs[0] = true;
    else neighs[0] = false;
    neigh_ids[0] = rank - Y*Z;
    sizes[0] = M*K*S;
    //printf("rank: %d, has front: %d, n_id: %d\n", rank, neighs[0], neigh_ids[0]);
    //back
    if (rank<Y*Z*(X-1)) neighs[1] = true;
    else neighs[1] = false;
    neigh_ids[1] = rank + Y*Z;
    sizes[1] = M*K*S;
    //printf("rank: %d, has back: %d, n_id: %d\n", rank, neighs[1], neigh_ids[1]);

    //left
    if (rank%Y != 0) neighs[2] = true;
    else neighs[2] = false;
    neigh_ids[2] = rank - 1;
    sizes[2] = N*M*S;
    //printf("rank: %d, has left: %d, n_id: %d\n", rank, neighs[2], neigh_ids[2]);

    //right
    if (rank%Y != Y-1) neighs[3] = true;
    else neighs[3] = false;
    neigh_ids[3] = rank + 1;
    sizes[3] = N*M*S;
    //printf("rank: %d, has right: %d, n_id: %d\n", rank, neighs[3], neigh_ids[3]);

    //down
    if ( (rank%(Y*Z)) / Y != 0) neighs[4] = true;
    else neighs[4] = false;
    neigh_ids[4] = rank - Y;
    sizes[4] = N*K*S;
    //printf("rank: %d, has down: %d, n_id: %d\n", rank, neighs[4], neigh_ids[4]);

    //up
    if ( (rank%(Y*Z)) / Y != Z-1) neighs[5] = true;
    else neighs[5] = false;
    neigh_ids[5] = rank + Y;
    sizes[5] = N*K*S;
    //printf("rank: %d, has up: %d, n_id: %d\n", rank, neighs[5], neigh_ids[5]);

    buffer_send = new std::vector<T> [6];
    buffer_recv = new std::vector<T> [6];

    for (int i = 0; i < 6; i++)
    {
        if (neighs[i])
        {
            buffer_send[i].resize(sizes[i]);
            buffer_recv[i].resize(sizes[i]);
        }
    }
}

template <typename T>
void solvers::mpi_solver<T>::transfer_out_of_buffer()
{
    std::vector<T> * p = buffer_recv;

    if (neighs[0])
    {
        for (int s = 0; s < S; s++)
            for (int j = 0; j < M; j++)
                for (int k = 0; k < K; k++)
        {
            //front
            data[(k+1)+(j+1)*K_+0*M_*K_+s*M_*K_*N_] = p[0][k+j*K+s*M*K];
            //std::cout << rank << "\n";
            //std::cout << data[(k+1)+(j+1)*K_+0*M_*K_+s*M_*K_*N_] << " " << buffer_recv[0][k+j*K+s*M*K] << "\n";
        }
    }
    if (neighs[1])
    {
        for (int s = 0; s < S; s++)
            for (int j = 0; j < M; j++)
                for (int k = 0; k < K; k++)
        {
            //back
            data[(k+1)+(j+1)*K_+(N_-1)*M_*K_+s*M_*K_*N_] = p[1][k+j*K+s*M*K];
        }
    }
    if (neighs[2])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
        {
            //left
            data[0+(j+1)*K_+(i+1)*M_*K_+s*M_*K_*N_] = p[2][j+i*M+s*N*M];
        }
    }
    if (neighs[3])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
        {
            //right
            data[(K_-1)+(j+1)*K_+(i+1)*M_*K_+s*M_*K_*N_] = p[3][j+i*M+s*N*M];
        }
    }
    if (neighs[4])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int k = 0; k < K; k++)
        {
            //down
            data[(k+1)+0*K_+(i+1)*M_*K_+s*M_*K_*N_] = p[4][k+i*K+s*N*K];
        }
    }
    if (neighs[5])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int k = 0; k < K; k++)
        {
            //up
            data[(k+1)+(M_-1)*K_+(i+1)*M_*K_+s*M_*K_*N_] = p[5][k+i*K+s*N*K];
        }
    }

    return;
}

template <typename T>
void solvers::mpi_solver<T>::transfer_in_buffer()
{
    std::vector<T> * p = buffer_send;

    if (neighs[0])
    {
        for (int s = 0; s < S; s++)
            for (int j = 0; j < M; j++)
                for (int k = 0; k < K; k++)
        {
            //front
            p[0][k+j*K+s*M*K] = data[(k+1)+(j+1)*K_+1*M_*K_+s*M_*K_*N_];
        }
    }
    if (neighs[1])
    {
        for (int s = 0; s < S; s++)
            for (int j = 0; j < M; j++)
                for (int k = 0; k < K; k++)
        {
            //back
            p[1][k+j*K+s*M*K] = data[(k+1)+(j+1)*K_+(N)*M_*K_+s*M_*K_*N_];
        }
    }
    if (neighs[2])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
        {
            //left
            p[2][j+i*M+s*N*M] = data[1+(j+1)*K_+(i+1)*M_*K_+s*M_*K_*N_];
        }
    }
    if (neighs[3])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < M; j++)
        {
            //right
            p[3][j+i*M+s*N*M] = data[(K)+(j+1)*K_+(i+1)*M_*K_+s*M_*K_*N_];
        }
    }
    if (neighs[4])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int k = 0; k < K; k++)
        {
            //down
            p[4][k+i*K+s*N*K] = data[(k+1)+1*K_+(i+1)*M_*K_+s*M_*K_*N_];
        }
    }
    if (neighs[5])
    {
        for (int s = 0; s < S; s++)
            for (int i = 0; i < N; i++)
                for (int k = 0; k < K; k++)
        {
            //up
            p[5][k+i*K+s*N*K] = data[(k+1)+(M)*K_+(i+1)*M_*K_+s*M_*K_*N_];
        }
    }

    return;
}

//MPI_DOUBLE HARDCODED
template<typename T>
void solvers::mpi_solver<T>::exchange()
{
    this->transfer_in_buffer();
    MPI_Request recv_req[6];
    MPI_Request send_req[6];
    for (int i = 0; i < 6; i++)
    {
        if (neighs[i])
        {
            MPI_Isend(buffer_send[i].data(), sizes[i], MPI_DOUBLE, 
                      neigh_ids[i], rank,
                      MPI_COMM_WORLD, &send_req[i]);
            MPI_Irecv(buffer_recv[i].data(), sizes[i], MPI_DOUBLE, 
                      neigh_ids[i], neigh_ids[i],
                      MPI_COMM_WORLD, &recv_req[i]);
        }
    }
    for (int i = 0; i < 6; i++)
    {
        if (neighs[i])
        {
            MPI_Waitall(1, &send_req[i], MPI_STATUS_IGNORE);
            MPI_Waitall(1, &recv_req[i], MPI_STATUS_IGNORE);
        }
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    this->transfer_out_of_buffer();
}

template<typename T>
void solvers::mpi_solver<T>::iteration(int time)
{
    for (int gg = 0; gg < S; gg++)
    {
        dif.iteration(data_new+N_*M_*K_*gg, data+N_*M_*K_*gg, (gg==0), st_point);
        if (time % 10 == 0)
        { 
            adv.set_foam_vels(foam.velos[ccount]);
            adv.upd_downward_velo(gg);
            adv.find_src();
        }
        for (int kk = 0; kk < 1; kk++)
        {
            adv.iteration(data+N_*M_*K_*gg, data_new+N_*M_*K_*gg);
            adv.iteration(data_new+N_*M_*K_*gg, data+N_*M_*K_*gg);
        }
    }
    if (time % 10 == 0)
        ccount++;
}
