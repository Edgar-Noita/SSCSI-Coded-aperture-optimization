%example to run the optimization process 
clear 
close all
L=24; %number of spectral bands
N=128; %128 x 128 coded aperture and 128 x 128 detector
ct1=10; %variable c1 in Eq. (22)
ctver=1; %this variable should be equal to 1 to follow Eq. (22)
shots=8; %number of captured snapshots

[t_iter,t_beg,t,beta]=exe(L,N,ct1,ctver,shots);

figure
plot(beta)
xlabel('Number of iterations')
ylabel('Cost function (Eq. (22))')

figure

subplot(1,4,1), imagesc(t_beg{1});
title('Initial code for first shot')


t_aux=t_iter{1};

subplot(1,4,2), imagesc(t_aux(:,:,1));
title('Code after iteration 1')


t_aux=t_iter{ceil(length(t_iter)/2)};
subplot(1,4,3), imagesc(t_aux(:,:,1));
title(strcat('Code after iteration ', num2str(ceil(length(t_iter)/2))));


t_aux=t_iter{end};
subplot(1,4,4), imagesc(t_aux(:,:,1));
title('Final optimal code')



