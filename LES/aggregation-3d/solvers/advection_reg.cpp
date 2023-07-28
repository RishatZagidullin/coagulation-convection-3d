#include "advection_reg.h"

namespace solvers
{
    template<typename T>
    std::pair<T, int> Advection3d_reg<T>::find_inds_foam(T src_pos, int len)
    {
        int ind = -1;
        T distance = 0.0;
        for (int i = 0; i < len; i++)
        {
            if (src_pos <= 0.0 || src_pos >= len*dd)
            {
                ind = 0;
                distance = -1.0;
                //std::cout << "this too should probably never happen\n";
                break;
            }
            //between i and i-1
            else if (src_pos < i*dd)
            {
                T d1 = i*dd - src_pos;
                T d2 = src_pos - (i-1)*dd;
                distance = d1 > d2 ? d2 : d1;
                ind = d1 > d2 ? i-1 : i;
                
                break;
            }
            else if (src_pos == i*dd)
            {
                ind = i;
                distance = 0.0;
                break;
            }
        }
        if (ind == -1)//std::cout << "advection::find_inds_foam error\n";
        {
            ind = len-1;
            distance = (src_pos - (len-1)*dd);
            //std::cout << ind << "\n" << distance << "\n";
            //std::cout  << src_pos << " " << len << " " << dd << "\n";
        }
        return std::pair<T, int>{distance, ind};
    }

    template<typename T>
    std::pair<T, int> Advection3d_reg<T>::find_inds(T src_pos, int len)
    {
        int ind = -1;
        T distance = 0.0;
        for (int i = 0; i < len; i++)
        {
            if (src_pos <= 0.0 || src_pos >= (len-1)*dd)
            {
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
        if (ind == -1) std::cout << "advection::find_inds error\n" << src_pos << " " << len << " " << dd << "\n";
        return std::pair<T, int>{distance, ind};
    }

    template<typename T>
    void Advection3d_reg<T>::set_inds(std::pair<T, int> pair_h, std::pair<T, int> pair_w, std::pair<T, int> pair_d, int idx)
    {
        this->weights[idx].clear();
        this->inds[idx].clear();
        if (pair_h.first < 0.0 || pair_w.first < 0.0 || pair_d.first < 0.0)
        {
            this->weights[idx].push_back(0.0);
            this->inds[idx].push_back(0);
        }
        else if (pair_d.first == 0 && pair_w.first == 0 && pair_h.first == 0)
        {
            this->weights[idx].push_back(1.0);
            this->inds[idx].push_back(pair_d.second+pair_w.second*K+pair_h.second*K*M);
        }
        else
        {
            T total_distance = 0.0;
            //==================i===j===k==============
            int cur_ind = pair_d.second+pair_w.second*K+pair_h.second*K*M;
            T cur_distance = pow(pair_d.first*pair_d.first+pair_w.first*pair_w.first+
                                 pair_h.first*pair_h.first, 0.5);
            this->inds[idx].push_back(cur_ind);
            this->weights[idx].push_back(1./cur_distance);
            total_distance += 1./cur_distance;
            //==================i===j===k-1============
            if (pair_d.first !=0)
            {
                cur_ind = (pair_d.second-1)+pair_w.second*K+pair_h.second*K*M;
                cur_distance = pow((dd-pair_d.first)*(dd-pair_d.first)+
                                   pair_w.first*pair_w.first+pair_h.first*pair_h.first, 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //==================i===j-1=k==============
            if (pair_w.first !=0)
            {
                cur_ind = pair_d.second+(pair_w.second-1)*K+pair_h.second*K*M;
                cur_distance = pow(pair_d.first*pair_d.first+(dd-pair_w.first)*(dd-pair_w.first)+
                                   pair_h.first*pair_h.first, 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //==================i-1=j===k==============
            if (pair_h.first !=0)
            {
                cur_ind = pair_d.second+pair_w.second*K+(pair_h.second-1)*K*M;
                cur_distance = pow(pair_d.first*pair_d.first+pair_w.first*pair_w.first+
                                   (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //==================i===j-1=k-1============
            if (pair_d.first !=0 && pair_w.first !=0)
            {
                cur_ind = (pair_d.second-1)+(pair_w.second-1)*K+pair_h.second*K*M;
                cur_distance = pow((dd-pair_d.first)*(dd-pair_d.first)+(dd-pair_w.first)*(dd-pair_w.first)+
                                   pair_h.first*pair_h.first, 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //==================i-1=j-1=k==============
            if (pair_w.first !=0 && pair_h.first !=0)
            {
                cur_ind = pair_d.second+(pair_w.second-1)*K+(pair_h.second-1)*K*M;
                cur_distance = pow(pair_d.first*pair_d.first+(dd-pair_w.first)*(dd-pair_w.first)+
                                   (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //==================i-1=j===k-1============
            if (pair_d.first !=0 && pair_h.first !=0)
            {
                cur_ind = (pair_d.second-1)+pair_w.second*K+(pair_h.second-1)*K*M;
                cur_distance = pow((dd-pair_d.first)*(dd-pair_d.first)+pair_w.first*pair_w.first+
                                   (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //==================i-1=j-1=k-1============
            if (pair_d.first !=0 && pair_w.first !=0 && pair_h.first !=0)
            {
                cur_ind = (pair_d.second-1)+(pair_w.second-1)*K+(pair_h.second-1)*K*M;
                cur_distance = pow((dd-pair_d.first)*(dd-pair_d.first)+(dd-pair_w.first)*(dd-pair_w.first)+
                                   (dd-pair_h.first)*(dd-pair_h.first), 0.5);
                this->inds[idx].push_back(cur_ind);
                this->weights[idx].push_back(1./cur_distance);
                total_distance += 1./cur_distance;
            }
            //=========================================
            for (int m=0; m<this->weights[idx].size(); m++)
                this->weights[idx][m] = this->weights[idx][m]/total_distance;
        }

    }

    template<typename T>
    void Advection3d_reg<T>::set_inds_foam(std::pair<T, int> pair_h, std::pair<T, int> pair_w, std::pair<T, int> pair_d, int idx)
    {
        //boundary source
        if (pair_h.first < 0.0 || pair_w.first < 0.0 || pair_d.first < 0.0)
        {
            //std::cout << "shouldn't happen, i think\n";
        }
        else if (pair_d.first == 0 && pair_w.first == 0 && pair_h.first == 0)
        {
            std::cout << "hope i will not see this\n";
            int cur_ind = pair_d.second+pair_w.second*K+pair_h.second*K*M;
            this->weights_foam[cur_ind].push_back(1.0);
            this->inds_foam[cur_ind].push_back(idx);
        }
        else
        {
            //==================i===j===k==============
            int cur_ind = pair_d.second+pair_w.second*K+pair_h.second*K*M;
            T cur_distance = pow(pair_d.first*pair_d.first+pair_w.first*pair_w.first+
                                 pair_h.first*pair_h.first, 0.5);
            this->inds_foam[cur_ind].push_back(idx);
            this->weights_foam[cur_ind].push_back(1./cur_distance);
        }

    }

    template<typename T>
    void Advection3d_reg<T>::set_foam_vels(Vector3d<T> * vel)
    {
        for (int i = 0; i < N*M*K; i++)
        {
            vels_foam[i].x = 0.;
            vels_foam[i].y = 0.;
            vels_foam[i].z = 0.;
            for (int j = 0; j < inds_foam[i].size(); j++)
            {
                vels_foam[i].x += vel[inds_foam[i][j]].x * 
                                  weights_foam[i][j];
                vels_foam[i].y += vel[inds_foam[i][j]].z * 
                                  weights_foam[i][j];
                vels_foam[i].z += vel[inds_foam[i][j]].y * 
                                  weights_foam[i][j];
            }
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::find_foam(Vector3d<T> * c, int n)
    {
        for (int i = 0; i < n ; i++)
        {
            //hardcoded cube size from openfoam sim
            float x = c[i].x+2.0;
            float y = c[i].y+2.0;
            float z = c[i].z+3.0;
            //w->z h->y but should be temporary
            auto pair_h = find_inds_foam(y, N);
            auto pair_w = find_inds_foam(z, M);
            auto pair_d = find_inds_foam(x, K);
            set_inds_foam(pair_h, pair_w, pair_d, i);
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
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    if (pow((i-K/2)*dd,2)+pow((k-K/2)*dd,2)>=0.16) vels[c].y = -0.005*pow(p, 2./3.);
                    c++;
                }
            }
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

        this->vels = new Vector3d<T> [N*M*K];
        this->vels_foam = new Vector3d<T> [N*M*K];
        int c = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    //data_new[c] = exp(-5.0*((i*dd-2.0)*(i*dd-2.0)+
                    //                    (j*dd-1.0)*(j*dd-1.0)+
                    //                    (k*dd-1.0)*(k*dd-1.0))
                    //          );

                    //vels[c].x = 0.0;
                    //vels[c].y = 0.0;
                    //T r_sq = pow(j*dd-1., 2)+pow(k*dd-1., 2);
                    //if (r_sq<=1.0)
                    //    vels[c].z = 1.0*(1.-r_sq*r_sq);
                    //else
                    //    vels[c].z = 0.0;
                    vels[c].x = 0.0;
                    vels[c].y = ( (j < M/2) && (pow((i-K/2)*dd,2)+pow((k-K/2)*dd,2)<0.09) ) ? 0.3 : 0.0;
                    vels[c].z = ( (j > M/2) ) ? (2.0*(j-M/2))/(M-M/2) : 0.0;
                    //std::cout << vels[c].z << "\n";

                    c++;
                }
            }
        }
    }

    template<typename T>
    void Advection3d_reg<T>::find_src_reg()
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    int idx = k+j*K+i*K*M;
                    Vector3d<T> cur_pos = {k*dd, j*dd, i*dd};
                    Vector3d<T> src_pos={cur_pos.x-vels[idx].x*dt, 
                                         cur_pos.y-vels[idx].y*dt, 
                                        cur_pos.z-vels[idx].z*dt};
                    auto pair_h = find_inds(src_pos.z, N);
                    auto pair_w = find_inds(src_pos.y, M);
                    auto pair_d = find_inds(src_pos.x, K);
                    set_inds(pair_h, pair_w, pair_d, idx);
                }
            }
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::find_src()
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                for (int k = 0; k < K; k++)
                {
                    int idx = k+j*K+i*K*M;
                    Vector3d<T> cur_pos = {k*dd, j*dd, i*dd};
                    Vector3d<T> src_pos={cur_pos.x-vels_foam[idx].x*dt, 
                                         cur_pos.y-vels_foam[idx].y*dt, 
                                        cur_pos.z-vels_foam[idx].z*dt};
                    auto pair_h = find_inds(src_pos.z, N);
                    auto pair_w = find_inds(src_pos.y, M);
                    auto pair_d = find_inds(src_pos.x, K);
                    set_inds(pair_h, pair_w, pair_d, idx);
                }
            }
        }
        return;
    }

    template<typename T>
    void Advection3d_reg<T>::iteration(T* data, T* data_new, int coef, bool boundary)
    {
        //std::swap(data_new, data);
        for (int i = 0; i < N*M*K; i++)
        {
            data_new[i] = 0.0;
            for (int j = 0; j < inds[i].size(); j++)
            {
                if (weights[i][j] == 2.0 && boundary)
                    data_new[i] = 1.0;
                else if (weights[i][j] == 2.0)
                    data_new[i] = data[inds[i][j]];
                else
                    data_new[i] += data[inds[i][j]]*weights[i][j];
            }
        }
        return;
    }

    template<typename T>
    Advection3d_reg<T>::~Advection3d_reg()
    {
        delete [] weights;
        delete [] inds;
        delete [] vels;
    }

    template class Advection3d_reg<double>;
    template class Advection3d_reg<float>;
}

