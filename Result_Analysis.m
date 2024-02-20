%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
simu_num=50;
simu_t=40;
average_wass=zeros(1,simu_t);
average_estimate_num=zeros(1,simu_t);
average_abs_error_estimate_num=zeros(1,simu_t);
average_OSAP=zeros(1,simu_t);
current_directory=pwd;
situation_name=['CLE'];
relative_situation_path=['\Simulation_Result\Summary_File'];
situation_path=[current_directory,relative_situation_path];
simulated_input_path=[current_directory,'\Simulated_Input_Data'];
cd(simulated_input_path);
load tar_num_save_run_1.mat;
for simu_k=1:simu_num
    cd(situation_path)
    eval(['load N_k_k_save_',situation_name,'_run',num2str(simu_k),'.mat;']);
    eval(['load wass_dist_',situation_name,'_run',num2str(simu_k),'.mat;']);
    eval(['load OSAP_metric_',situation_name,'_run',num2str(simu_k),'.mat;']);
    average_estimate_num=average_estimate_num+N_k_k_save;
    average_wass(find(isnan(wass_dist)==0))=average_wass(find(isnan(wass_dist)==0))+wass_dist(find(isnan(wass_dist)==0));
    average_abs_error_estimate_num=average_abs_error_estimate_num+abs(N_k_k_save-tar_num_save);
    average_OSAP=average_OSAP+OSAP_metric;
end
average_estimate_num=average_estimate_num/simu_num;
average_wass=average_wass/simu_num;
average_abs_error_estimate_num=average_abs_error_estimate_num/simu_num;
average_OSAP=average_OSAP/simu_num;
eval(['save average_estimate_num_',situation_name,'.mat average_estimate_num;']);
eval(['save average_wass_',situation_name,'.mat average_wass;']);
eval(['save average_OSAP_',situation_name,'.mat average_OSAP;']);
eval(['save average_abs_error_estimate_num_',situation_name,'.mat average_abs_error_estimate_num;']);
cd(current_directory)
% plot(average_estimate_num,'-x')
% hold on
% load tar_num_save_run_1.mat
% plot(tar_num_save,'--')
% figure
% plot(average_wass)
