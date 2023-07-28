#include "advection_reg.h"

namespace solvers
{
    template<typename T>
    std::pair<T, int> Advection3d_reg<T>::find_inds(T src_pos, int len)
    {
        int ind = -1;
        T distance = 0.0;
        for (int i = 0; i < len; i++)
        {
            if (src_pos < 0.0 || src_pos > (len-1)*dd)
            {
                std::cout << "this is unwanted\n";
                ind = 0;
                distance = -1.0;
                break;
            }
            else if (src_pos < i*dd)
            {
                ind = i;//between i and i-1
                distance = i*dd - src_pos;
                break;
            }
            else if (src_pos == i*dd)
            {
                ind = i;
                distance = 0.0;
                break;
            }
        }
        if (ind == -1) 
        {
            std::cout << "advection::find_inds error\n";
            std::cout << src_pos << " " << len << " " << dd << "\n";
        }
        return std::pair<T, int>{distance, ind};
    }

    template<typename T>
    void Advection3d_reg<T>::set_inds(std::pair<T,int> pair_x, 
                                      std::pair<T,int> pair_y,
                                      std::pair<T,int> pair_z, int idx)
    {
        this->weights[idx].clear();
        this->inds[idx].clear();
        if (pair_x.first<0.0||pair_y.first<0.0||pair_z.first<0.0)
        {
            std::cout << "this is also undwanted\n";
            this->weights[idx].push_back(0.0);
            this->inds[idx].push_back(idx);
        }
        else if (pair_x.first==0&&pair_y.first==0&&pair_z.first==0)
        {
            this->weights[idx].push_back(1.0);
            this->inds[idx].push_back(pair_y.second + pair_z.second*K
                                    + pair_x.second*K*M);
        }
        else
        {
            T tot_dist = 0.0;
            //==================i===j===k==============
            int cur_ind = pair_y.second + pair_z.second*K
                        + pair_x.second*K*M;
            T cur_distance = pow(pair_x.first*pair_x.first 
                               + pair_y.first*pair_y.first 
                               + pair_z.first*pair_z.first, 0.5);
            this->inds[idx].push_back(cur_ind);
            this->weights[idx].push_back(1./cur_distance);
            tot_dist += 1./cur_distance;
            //==================i===j===k-1============
            if (pair_y.first !=0)
            {
                cur_ind = (pair_y.second-1) + pair_z.second*K
                         + pair_x.second*K*M;
                cur_distance = pow((dd-pair_y.first)*(dd-pair_y.first)
                                 + pair_z.first*pair_z.first
                                 + pair_x.first*pair_x.first, 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //==================i===j-1=k==============
            if (pair_z.first !=0)
            {
                cur_ind = pair_y.second + (pair_z.second-1)*K
                        + pair_x.second*K*M;
                cur_distance = pow(pair_y.first*pair_y.first
                                 + (dd-pair_z.first)*(dd-pair_z.first)
                                 + pair_x.first*pair_x.first, 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //==================i-1=j===k==============
            if (pair_x.first !=0)
            {
                cur_ind = pair_y.second + pair_z.second*K
                        + (pair_x.second-1)*K*M;
                cur_distance = pow(pair_y.first*pair_y.first
                           + pair_z.first*pair_z.first
                           + (dd-pair_x.first)*(dd-pair_x.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //==================i===j-1=k-1============
            if (pair_y.first !=0 && pair_z.first !=0)
            {
                cur_ind = (pair_y.second-1) + (pair_z.second-1)*K
                        + pair_x.second*K*M;
                cur_distance = pow((dd-pair_y.first)*(dd-pair_y.first)
                                 + (dd-pair_z.first)*(dd-pair_z.first)
                                 + pair_x.first*pair_x.first, 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //==================i-1=j-1=k==============
            if (pair_z.first !=0 && pair_x.first !=0)
            {
                cur_ind = pair_y.second + (pair_z.second-1)*K
                        + (pair_x.second-1)*K*M;
                cur_distance = pow(pair_y.first*pair_y.first
                           + (dd-pair_z.first)*(dd-pair_z.first)
                           + (dd-pair_x.first)*(dd-pair_x.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //==================i-1=j===k-1============
            if (pair_y.first !=0 && pair_x.first !=0)
            {
                cur_ind = (pair_y.second-1) + pair_z.second*K
                        + (pair_x.second-1)*K*M;
                cur_distance = pow((dd-pair_y.first)*(dd-pair_y.first)
                           + pair_z.first*pair_z.first
                           + (dd-pair_x.first)*(dd-pair_x.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //==================i-1=j-1=k-1============
            if (pair_x.first!=0 && pair_y.first!=0 && pair_z.first!=0)
            {
                cur_ind = (pair_y.second-1) + (pair_z.second-1)*K
                        + (pair_x.second-1)*K*M;
                cur_distance = pow((dd-pair_y.first)*(dd-pair_y.first)
                                 + (dd-pair_z.first)*(dd-pair_z.first)
                                 + (dd-pair_x.first)*(dd-pair_x.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                tot_dist += 1./cur_distance;
            }
            //=========================================
            for (int m=0; m<this->weights[idx].size(); m++)
                this->weights[idx][m] = this->weights[idx][m]/tot_dist;
        }

    }

    template<typename T>
    void Advection3d_reg<T>::set_foam_vels(Vector3d<T> * vel)
    {
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < M-1; j++)
                for (int k = 1; k < K-1; k++)
        {
            int idx = k+j*K+i*K*M;

            vels_foam[idx].x = 0.;
            vels_foam[idx].y = 0.;
            vels_foam[idx].z = 0.;
            for (int l = 0; l < inds_foam[idx].size(); l++)
            {
                vels_foam[idx].x += vel[inds_foam[idx][l]].x * 
                                  weights_foam[idx][l];
                vels_foam[idx].y += vel[inds_foam[idx][l]].y * 
                                  weights_foam[idx][l];
                vels_foam[idx].z += vel[inds_foam[idx][l]].z * 
                                  weights_foam[idx][l];
            }
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::find_foam(Vector3d<T> * c, int n, Vector3d<T> st_point)
    {
        //hardcoded from openfoam (cc.txt file)
        double foam_dy = 0.1;
        double foam_dx = 0.1;
        double foam_dz = 0.1;
        for (int cc = 0; cc < n; cc++)
        {
            //hardcoded from openfoam sim (left points become 0)
            //printf("%f, %f, %f\n", st_point.x, st_point.y, st_point.z);
            float x = c[cc].x-st_point.x;//+dd;
            float y = c[cc].y-st_point.y;//+dd;
            float z = c[cc].z-st_point.z;//+dd;
            if (x<0.0 || y<0.0 || z<0.0) continue;
            if (x>dd*N|| y>dd*K||z>dd*M) continue;

            int ind_N_max = (int) ((x+foam_dx)/dd);
            if (ind_N_max > N-2) ind_N_max = N-2;
            int ind_N_min = (int) ((x-foam_dx)/dd);
            if (ind_N_min < 0) ind_N_min = 0;

            int ind_M_max = (int) ((z+foam_dz)/dd);
            if (ind_M_max > M-2) ind_M_max = M-2;
            int ind_M_min = (int) ((z-foam_dz)/dd);
            if (ind_M_min < 0) ind_M_min = 0;

            int ind_K_max = (int) ((y+foam_dy)/dd);
            if (ind_K_max > K-2) ind_K_max = K-2;
            int ind_K_min = (int) ((y-foam_dy)/dd);
            if (ind_K_min < 0) ind_K_min = 0;

            for (int i = ind_N_min; i < ind_N_max; i++)
                for (int j = ind_M_min; j < ind_M_max; j++)
                    for (int k = ind_K_min; k < ind_K_max; k++)
                    {
                        int ind = (i+1)*M*K+(j+1)*K+(k+1);
                        double dist = pow(pow(i*dd-x,2)+pow(k*dd-y,2)+pow(j*dd-z,2),0.5);

                        if (weights_foam[ind].size()>0 && weights_foam[ind][0] == 1.)
                        {
                            continue;
                        }
                        else if (dist == 0.0)
                        {
                            std::vector<int>().swap(inds_foam[ind]);
                            std::vector<T>().swap(weights_foam[ind]);
                            inds_foam[ind].push_back(cc);
                            weights_foam[ind].push_back(1.);
                        }
                        else
                        {
                            inds_foam[ind].push_back(cc);
                            weights_foam[ind].push_back(1./dist);
                        }
                    }
        }
        for (int i = 0; i < N*M*K; i++)
        {
            T tot = 0.0;
            for (int j = 0; j < this->weights_foam[i].size(); j++)
            {
                tot += this->weights_foam[i][j];
            }
            for (int j = 0; j < this->weights_foam[i].size(); j++)
            {
                this->weights_foam[i][j] /= tot;
            }
        }
    }

    template<typename T>
    void Advection3d_reg<T>::upd_downward_velo(int p)
    {
        int c = 0;
        for (int i = 1; i < N-1; i++)
            for (int j = 1; j < M-1; j++)
                for (int k = 1; k < K-1; k++)
        {
            //if (i > N/2) vels_foam[c].z += -0.3*pow(p, 2./3.);
            vels_foam[c].z += -0.2*pow(p, 2./3.);
            //std::cout << p << " " << -0.1*pow(p, 2./3.) << "\n";
            c++;
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::init_data()
    {
        this->weights_foam = new std::vector<T> [N*M*K];
        this->inds_foam = new std::vector<int> [N*M*K];

        this->weights = new std::vector<T> [N*M*K];
        this->inds = new std::vector<int> [N*M*K];

        this->vels_foam = new Vector3d<T> [N*M*K];
    }

    template<typename T>
    void Advection3d_reg<T>::find_src()
    {
        for (int i = 1; i < N-1; i++)
        {
            for (int j = 1; j < M-1; j++)
            {
                for (int k = 1; k < K-1; k++)
                {
                    int idx = k+j*K+i*K*M;
                    Vector3d<T> cur_pos = {i*dd, k*dd, j*dd};
                    Vector3d<T> src_pos={cur_pos.x-vels_foam[idx].x*dt, 
                                         cur_pos.y-vels_foam[idx].y*dt, 
                                        cur_pos.z-vels_foam[idx].z*dt};
                    auto pair_x = find_inds(src_pos.x, N);
                    auto pair_y = find_inds(src_pos.y, K);
                    auto pair_z = find_inds(src_pos.z, M);
                    //if (i == 0)
                    //{
                    //    std::cout << "i==0\n";
                    //    std::cout << pair_x.first << " " << pair_x.second << "\n";
                    //    std::cout << pair_y.first << " " << pair_y.second << "\n";
                    //    std::cout << pair_z.first << " " << pair_z.second << "\n";
                    //}
                    //if (i == 1)
                    //{
                    //    std::cout << "i==1\n";
                    //    std::cout << pair_x.first << " " << pair_x.second << "\n";
                    //    std::cout << pair_y.first << " " << pair_y.second << "\n";
                    //    std::cout << pair_z.first << " " << pair_z.second << "\n";
                    //}
                    set_inds(pair_x, pair_y, pair_z, idx);
                }
            }
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::iteration(T* data, T* data_new)
    {
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                for (int k = 0; k < K; k++)
        {
            int idx = k+j*K+i*K*M;
            if (i == 0 || j == 0 || k == 0 || i == N-1 || j == M-1 || k == K-1)
            {
                data_new[idx] = data[idx];
                continue;
            }
            data_new[idx] = 0.0;
            for (int l = 0; l < inds[idx].size(); l++)
                data_new[idx] += data[inds[idx][l]]*weights[idx][l];
        }
        return;
    }

    template<typename T>
    Advection3d_reg<T>::~Advection3d_reg()
    {
        delete [] weights;
        delete [] inds;
        delete [] weights_foam;
        delete [] inds_foam;
        delete [] vels_foam;
    }

    template class Advection3d_reg<double>;
    template class Advection3d_reg<float>;
}

