
%This function runs the optimization algorithm proposed by
%"Salazar et al. Coded Aperture Optimization in Spatial Spectral
%Compressive Spectral Imagers"
%Please cite this artice if used


%%%REQUIRED INPUTS

%L: Number of spectral bands to be recovered
%N: N x N is the size of the coded aperture and sensor
%c1: Varible c1 from Eq. (22) of the paper, related to the on-diagonal
%elements of H
%ctver: Variable that adds importance to the separation of the ON pixels
%and the off-diagonal elements of H (In the paper it is assumed to be 1). 
%shots: Number of captured snapshots

%%Outputs

%t_iter: pattern after each iteration
%t_beg: Inital guess of the coded aperture
%t: Final coded aperture patterns after optimization
%beta: Cost function evolution over the iterations, given by Eq. (22) 
%of the paper.


function [t_iter,t_beg,t,beta]=exe(L,N,ct1,ctver,shots)



s=L/N; %This is valid assuming that beta=1 and Delta_c=Delta_d, in Eq. (8)
%of the paper
load(strcat('Vi_t_128_s',num2str(s),'.mat')); %This matrix contains the percentage of a coded aperture
%that affects a given voxel (pm values in Eq. (7)). For example, if N x  N=128 X 128,
%and L=24, Vi_t would have dimensions of 128*128 x 128*128*24. This is done
%order to avoid a full search over the whole domain (Function Find of the
%proposed algorithm in the paper.

%This Gitub repository only gives access to 5 different cases of such
%function:

%256 x 256 x 24 and s=0.09375
%128 x 128 x 24 and s=0.1875
%128 x 128 x 12 and s=0.09375
%128 x 128 x 6 and s=0.046875

%To have access to other values, please contact the author
[t_iter,t_beg,t,beta]=evaluation_2_blue_fin(N,N,L,shots,ct1,ctver,ctver,Vi_t);

