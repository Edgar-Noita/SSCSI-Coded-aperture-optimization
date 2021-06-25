
%This function calculates the first and second term of EQ. (22) in
%"Salazar et al. Coded Aperture Optimization in Spatial Spectral
%Compressive Spectral Imagers"
%Please cite this artice if used


%%%REQUIRED INPUTS

%L: Number of spectral bands to be recovered
%N: N x N is the size of the coded aperture and sensor
%Q:Number of captured snapshots
%Vi_t: Variable already defined in exe.m
%ct: Varible c1 from Eq. (22) of the paper, related to the on-diagonal
%elements of H
%ct_1: Variable that adds importance of the off diagnal elements(Not specified in the paper) 


%%Outputs

%fnc_3: Calculated cost function

function [fnc_3]=cost_fn(N,L,Q,Vi_t,T,ct,ct_1)

div=4; % number of paralell tasks

 
a=speye(Q);
Vi_mat=kron(a,Vi_t); 
outer_sum=Vi_mat'*T;

fnc_t=zeros(div,1);
fnc_t2=zeros(div,1);
fnc_t3=zeros(div,1);
parfor par=1:div

for i=((par-1)/div)*N^2*L+1:((par)/div)*N^2*L
   
    
    %choosing the indexes of the coded apertures used in each iteration
    j=0:L-1;
    indx=(j-floor((i-1)/(N^2)))*(N^2)+i;
    indx_a=indx(indx~=i);  
    func_i=outer_sum(i:N^2*L:(Q-1)*N^2*L+i); 
    fnc_v=zeros(length(indx_a),1); 
   
    %Calculating the first term of  Eq. (22), off diagonal
    %elements
    for sh=1:Q
        func_j=outer_sum(indx_a); 
        fnc_v=fnc_v+((func_i(sh)*func_j)); 
        indx_a=indx_a+(N^2)*L;
    end
   
    fnc_t(par)=fnc_t(par)+sum(fnc_v.^2); %sum of the squared off diagonal elements
    fnc_t2(par)= fnc_t2(par)+(sum(func_i.^2)-1)^2; %Calculation of the second term EQ. (22) (ON diagonal elements)
    fnc_t3(par)= ct_1*fnc_t(par)+ct*fnc_t2(par); %Calculation of the total cost function

end   
        
end

fnc_3=sum(fnc_t3); %sum of the costs function over the parallel tasks
    