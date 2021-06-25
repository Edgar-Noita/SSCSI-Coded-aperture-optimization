

%This function calculates the first term and second term of EQ, (22) in
%"Salazar et al. Coded Aperture Optimization in Spatial Spectral
%Compressive Spectral Imagers"
%Please cite this artice if used

%Required Inpunts:

%L: Number of spectral bands to be recovered
%N: N x N is the size of the coded aperture and sensor
%Q:current shot
%Vi_mat2: Variable already defined in evaluation_2_blue_fin.m
%columns: Colums of matrix H that will be used in the current calculations
%t: coded aperture pattern
%ct: Variable c1 in EQ. (22)
%ct_1: Can be used to weight the importance of the off-diagonal elements in the algorithm (Not used in the paper)


%Outputs:

%fnc: Calculated cost function

function [fnc]=cost_fn_2_opt2(N,L,Vi_mat2,Q,columns,t,Q2,ct,ct_1)


indx=(((repmat(0:L-1,[length(columns),1])-floor((repmat(columns',[1,L])-1)/(N^2)))*(N^2)+repmat(columns',[1,L])))';



indx_a=reshape(indx((indx'~=columns')'),[L-1,length(columns)])';

fnc=ct_1*sum(sum(((reshape(Vi_mat2(:,(N^2)*L*(Q-1)+indx_a)'*t,size(indx_a))).*(Vi_mat2(:,(N^2)*L*(Q-1)+columns)'*t)).^2,2)) +ct*sum((sum((reshape(Vi_mat2(:,N^2*L*repmat(0:Q2-1,[length(columns),1])+columns')'*t,[size(indx_a,1),Q2])).^2,2)-1).^2);
