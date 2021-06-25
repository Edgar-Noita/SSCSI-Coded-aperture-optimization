
%functions that calculates the third term of EQ. (22)

%%Required Inputs

%L: Number of spectral bands
%t: Coded aperture
%ctver: Weight of the third term in Eq. (22) (In the paper it is assumed to
%be 1)

%Output

%Value of third expression of EQ. (22)
function [fnc]=cost_fn_ver(L,t,ctver)

[N,M,Q]=size(t); %Size of the coded aperture and numer of shots
div=4; %Division of calculation into 4 parallel tasks
window_size=L; %Implemented window size, Window size, if L, size is 2L-1 x 2L-1
fnc_par=zeros(div,1); %Vector that stores the result for every parallel task
y=mod((1:N^2)-1,N)+1; %indexes of the coded aperture over Y
x=floor(((1:N^2)-1)/N)+1; %indexes of the coded aperture over x

%Over the whole coded aperture, using parallel calculation
parfor par=1:div 
     for i=((N^2)/div)*(par-1)+1:((N^2)/div)*par
        
             blue_hor=x(i)-((window_size)-1):x(i)+((window_size)-1); %windows size of 2L-1 x 2L -1
             blue_hor=blue_hor(blue_hor>0 & blue_hor<(N+1));
        
             blue_ver=y(i)-((window_size)-1):y(i)+((window_size)-1); 
             blue_ver=blue_ver(blue_ver>0 & blue_ver<(N+1));
            
             [hor,ver]=meshgrid(blue_hor,blue_ver); %Generation of the indexes grid
   
             weights=1./(sqrt(abs(hor-x(i)).^2+abs(ver-y(i)).^2));
             weights(isinf(weights))=0; %change the middle coefficient from Inf to 0
          
             weights=repmat(weights,[1,Q]); %Weigths repetition for all shots
             fnc_par(par)=fnc_par(par)+ (sum(sum(weights.*t(repmat(sub2ind([N,M],ver,hor),[1,Q])+kron([0:Q-1]*N^2,ones(size(ver))))))); %Calculation of the cost function
     end
end
fnc=ctver*sum(fnc_par); %multiplication of the cost function by a given variable (assumed 1 in the paper)
   