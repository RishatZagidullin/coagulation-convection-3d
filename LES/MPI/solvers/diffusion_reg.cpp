#include "diffusion_reg.h"

namespace solvers
{

    //keep in mind that a,b,c are altered after the execution
    template<typename T>
    void Diffusion3d_reg<T>::solve_matrix (int n, T *a, T *c, 
                                       T *b, T *f, T *x)
    {
        T m;
        for (int i = 1; i < n; i++)
        {
            m = a[i] / c[i-1];
            c[i] = c[i] - m * b[i-1];
            f[i] = f[i] - m * f[i-1];
        }
        x[n-1] = f[n-1]/c[n-1];

        for (int i = n - 2; i >= 0; i--)
        {
            x[i] = ( f[i] - b[i] * x[i+1] ) / c[i];
        }
    }

    template<typename T>
    void Diffusion3d_reg<T>::iteration(T* data, T* data_new, bool monomer, Vector3d<T> st_point)
    {
        //=======================N=================================
        T * a = new T [N];
        T * b = new T [N];
        T * c = new T [N];
        T * rhs = new T [N];

        T * output = new T [N];

        for (int j = 1; j < M-1; j++)
        {
            for (int k = 1; k < K-1; k++)
            {
                for (int i = 1; i < N-1; i++)
                {
                    a[i] = -D*0.5/dd/dd*dt;
                    b[i] = -D*0.5/dd/dd*dt;
                    c[i] = 1.+D/dd/dd*dt;
                    rhs[i] = 1.*data[k+j*K+i*K*M]
                             +0.5*dt*D/dd/dd*
                             (data[k+j*K+(i+1)*K*M]
                             -2.*data[k+j*K+(i)*K*M]
                             +data[k+j*K+(i-1)*K*M])
                             +dt*D/dd/dd*
                             (data[k+(j+1)*K+i*K*M]
                             -2.*data[k+(j)*K+i*K*M]
                             +data[k+(j-1)*K+i*K*M])
                             +dt*D/dd/dd*
                             (data[(k+1)+j*K+i*K*M]
                             -2.*data[(k)+j*K+i*K*M]
                             +data[(k-1)+j*K+i*K*M]);
                }
                a[0] = 0.0;
                b[0] = 0.0;
                c[0] = 1.0;
                a[N-1] = 0.0;
                b[N-1] = 0.0;
                c[N-1] = 1.0;

                rhs[0] = data[k+j*K+0*K*M];
                rhs[N-1] = data[k+j*K+(N-1)*K*M];
                solve_matrix(N, a, c, b, rhs, output);
                for (int i = 0; i < N; i++)
                {
                    data_new[k+j*K+i*K*M] = output[i];
                }
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;
        //=======================M=================================
        a = new T [M];
        b = new T [M];
        c = new T [M];
        rhs = new T [M];

        output = new T [M];

        for (int i = 1; i < N-1; i++)
        {
            for (int k = 1; k < K-1; k++)
            {
                for (int j = 1; j < M-1; j++)
                {
                    a[j] = -D*0.5/dd/dd*dt;
                    b[j] = -D*0.5/dd/dd*dt;
                    c[j] = 1.+D/dd/dd*dt;
                    rhs[j] = 1.*data_new[k+j*K+i*K*M]
                             +0.5*dt*D/dd/dd*
                             (-data[k+(j+1)*K+i*K*M]
                             +2.*data[k+(j)*K+i*K*M]
                             -data[k+(j-1)*K+i*K*M]);
                }
                a[0] = 0.0;
                b[0] = 0.0;
                c[0] = 1.0;
                a[M-1] = 0.0;
                b[M-1] = 0.0;
                c[M-1] = 1.0;
                double r_sq = pow(st_point.x + i*dd, 2) + pow(st_point.y + k*dd, 2);
                if (monomer && st_point.z == -3.0 && r_sq<0.09) 
                    rhs[0] = 1.0;
                else
                    rhs[0] = data[k+0*K+i*K*M];
                rhs[M-1] = data[k+(M-1)*K+i*K*M];
                solve_matrix(M, a, c, b, rhs, output);
                for (int j = 0; j < M; j++)
                {
                    data_new[k+j*K+i*K*M] = output[j];
                }
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;
        //=======================K=================================
        a = new T [K];
        b = new T [K];
        c = new T [K];
        rhs = new T [K];

        output = new T [K];

        for (int i = 1; i < N-1; i++)
        {
            for (int j = 1; j < M-1; j++)
            {
                for (int k = 1; k < K-1; k++)
                {
                    a[k] = -D*0.5/dd/dd*dt;
                    b[k] = -D*0.5/dd/dd*dt;
                    c[k] = 1.+D/dd/dd*dt;
                    rhs[k] = 1.*data_new[k+j*K+i*K*M]
                            +0.5*dt*D/dd/dd*
                            (-data[(k+1)+j*K+i*K*M]
                            +2.*data[(k)+j*K+i*K*M]
                            -data[(k-1)+j*K+i*K*M]);
                }
                a[0] = 0.0;
                b[0] = 0.0;
                c[0] = 1.0;
                a[K-1] = 0.0;
                b[K-1] = 0.0;
                c[K-1] = 1.0;
                rhs[0] = data[0+j*K+i*K*M];
                rhs[K-1] = data[K-1+j*K+i*K*M];
                solve_matrix(K, a, c, b, rhs, output);
                for (int k = 0; k < K; k++)
                {
                    data_new[k+j*K+i*K*M] = output[k];
                }
            }
        }
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] rhs;
        delete [] output;
        //=========================================================
        return;
    }

    template class Diffusion3d_reg<double>;
    template class Diffusion3d_reg<float>;
}
