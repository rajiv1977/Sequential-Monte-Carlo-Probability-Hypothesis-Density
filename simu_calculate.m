%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analysis performance of PHD filter based on the simulation
clc
clear
simu_num_start=1;%starting scan of simulation test
simu_num_end=100;%end scan of simulation test
samp_t=1;%sample time, unit second
simu_t=40;%simulation time in one experiment, unit second
tar_num_error_avg=zeros(1,simu_t);%average of absolute error of targets' number's estimation
wass_dist_avg=zeros(1,simu_t);%avverage of wassterstein distance beteween estimation and real targets' state
tar_num_avg=zeros(1,simu_t);%average number of real target
hold on
for simu_k=simu_num_start:simu_num_end
    eval(['load stat_save',num2str(simu_k),'.mat;']);
    eval(['load tar_est_save',num2str(simu_k),'.mat;']);
    eval(['load peak_num_save',num2str(simu_k),'.mat;']);
    eval(['load tar_num_save',num2str(simu_k),'.mat;']);
    eval(['load wass_dist',num2str(simu_k),'.mat;']);
    eval(['load tar_num_error',num2str(simu_k),'.mat;']);
    tar_num_error_avg=tar_num_error_avg+tar_num_error;
    nan_ind=isnan(wass_dist);
    wass_dist(nan_ind)=0;
    wass_dist_avg=wass_dist_avg+wass_dist;
%     plot(wass_dist)
    tar_num_avg=tar_num_avg+tar_num_save;
end
tar_num_error_avg=tar_num_error_avg/(simu_num_end-simu_num_start);
wass_dist_avg=wass_dist_avg/(simu_num_end-simu_num_start);
tar_num_avg=tar_num_avg/(simu_num_end-simu_num_start);