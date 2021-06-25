%This function runs the optimization algorithm proposed by
%"Salazar et al. Coded Aperture Optimization in Spatial Spectral
%Compressive Spectral Imagers"
%Please cite this artice if used


%%%REQUIRED INPUTS

%L: Number of spectral bands to be recovered
%N,M: N x M is the size of the coded aperture and sensor
%ct1: Varible c1 from Eq. (22) of the paper, related to the on-diagonal
%elements of H
%Q:Number of cpatured snapshots
%ct_ver: Variable that adds importance to the separation of the ON pixels
%over the Y axis(Not specified in the paper but might be used). 
%ct_1: Variable that adds importance to the off diagonal elements of H, first term, EQ. (22) (Not specified in the paper but might be used). 
%Vi_t: This matrix contains the percentage of a coded aperture
%that affects a given voxel (pm values in Eq. (7)). For more details, see
%exe.m

%%Outputs

%t_beg: Inital guess of the coded aperture
%t: Final coded aperture patterns after optimization
%beta_upd: Cost function evolution over the iterations, given by Eq. (22) 
%of the paper.
function [t_beg,t,beta_upd]=evaluation_2_blue_fin(N,M,L,Q,ct,ct_1,ct_ver,Vi_t)

%Initial guess of coded aperture, must be Boolean
t_beg=(my_bool(N,M,Q));
    for ii=1:Q
        t(:,:,ii)=t_beg{ii};
    end

b=speye(Q);
Vi_mat2=kron(b,Vi_t); %Kroncneker product between the matrix Vi_t
%and a identity matrix of size Q x Q,

cont=1; %variable to control the break of the While loop
it=1; %number of iterations

[beta_ver_2]=cost_fn_ver(L,t,ct_ver); %Function that calculates the cost function
%related to condition: Separation of ON Pixels (third term Eq. (22))
   
[beta_upd_2]=cost_fn(N,L,Q,Vi_t,t(:),ct,ct_1); %Function that calculates the cost function
%related to first and second term of EQ. (22). ct_1 is added to control the
%weight of off diagonal elements in the algorithm. In the paper it is
%asummed 1.
beta_upd(1)=beta_ver_2+beta_upd_2; %Total cost function
epsilon=(0.1/100)*beta_upd(1); %Stop Criteria

%Iterations start
while cont~=-1

it  %print the iteration number

random_access=randperm(N^2); %Random Walk over the 2D code
y=mod(random_access-1,N)+1;
x=floor((random_access-1)/N)+1;


window_size=L; %Size of the window related to third term Eq. (22)
   
    for i=1:N^2 %over the coded aperture elements
      
          
           %Find a value for term 3 of EQ. (22), for a given x(i), y(i)
             blue_hor=x(i)-((window_size)-1):x(i)+((window_size)-1);
             blue_hor=blue_hor(blue_hor>0 & blue_hor<(N+1));
        
             blue_ver=y(i)-((window_size)-1):y(i)+((window_size)-1); 
             blue_ver=blue_ver(blue_ver>0 & blue_ver<(N+1));
             
             [hor,ver]=meshgrid(blue_hor,blue_ver);
    
             weights=1./(sqrt(abs(hor-x(i)).^2+abs(ver-y(i)).^2));
             weights(isinf(weights))=0;
  
             indx=find(Vi_t(random_access(i),:)); 
 
             for sh=1:Q %calculate the cost function for every possible combination
                        %over the shots.. For example, for 3 shots= (1,0,0), (0,1,0), (0,0,1)      
           
                taux(y(i),x(i),:)=zeros(1,1,Q);
                taux(y(i),x(i),sh)=1;
                taux_q=taux(:,:,sh);
        
                [beta_total_p(sh)]=cost_fn_2_opt2(N,L,Vi_mat2,sh,indx,taux(:),Q,ct,ct_1); %calculation of thr first and 
                %second term of EQ. (22)
      
                 beta_ver(sh)=sum(sum(weights.*taux_q(sub2ind(size(taux_q),ver,hor)))); %Calculation of third term of EQ. (22)
                 beta_total(sh)=ct_ver*beta_ver(sh)+beta_total_p(sh); %Summation of the two terms to calculate EQ. (22)
      
             end
       
             beta2=find(beta_total==min(beta_total)); %find all the combinations that throws
             %the minimum value of Eq. (22)
            
             %choose randomly between possible combinations that throw a
             %minimum
             t(y(i),x(i),:)=zeros(1,1,Q); 
             t(y(i),x(i),beta2(randi([1,length(beta2)])))=1;
      
       
     
    end
 
    it=it+1; 
   
    %Calculates the cost function after the present iteration

    [beta_ver_2]=cost_fn_ver(L,t,ct_ver);
   
    [beta_upd_2]=cost_fn(N,L,Q,Vi_t,t(:),ct,ct_1);
     
    beta_upd(it)=beta_ver_2+beta_upd_2;
    
    %checks the stop criteria 
    if (it>2)
      if(abs(beta_upd(it-1)-beta_upd(it))<=epsilon)
        cont=-1;
      end
      
   end
   

end

%emsambles the answer as a cell
for i=1:Q
   t3{i}=t(:,:,i);
end
t=t3;


