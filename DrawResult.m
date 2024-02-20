%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clear
current_directory=pwd;
simu_result_path=[current_directory,'\Simulation_Result\Summary_File'];
UnknowC_path=[current_directory,'\Reference_File\UnknownC'];
KnownC_path=[current_directory,'\Reference_File\KnownC'];
simulated_input_path=[current_directory,'\Simulated_Input_Data'];
figure
hold on
cd(UnknowC_path)
load average_estimate_num_UnKnownC.mat;
plot(average_estimate_num,'k-o');
cd(KnownC_path);
load average_estimate_num_KnownC.mat;
plot(average_estimate_num,'-*');
cd(simu_result_path);
load average_estimate_num_CLE.mat;
plot(average_estimate_num,'rs-');
cd(simulated_input_path)
load tar_num_save_run_1.mat;
plot(tar_num_save,'--');

figure
hold on
cd(UnknowC_path)
load average_abs_error_estimate_num_UnKnownC.mat;
plot(average_abs_error_estimate_num,'k-o');
cd(KnownC_path);
load average_abs_error_estimate_num_KnownC.mat;
plot(average_abs_error_estimate_num,'-*');
cd(simu_result_path);
load average_abs_error_estimate_num_CLE.mat;
plot(average_abs_error_estimate_num,'rs-');

figure
hold on
cd(KnownC_path);
load average_wass_KnownC.mat;
plot(average_wass,'-*');
cd(simu_result_path);
load average_wass_CLE.mat;
plot(average_wass,'rs-');

figure
hold on
cd(KnownC_path);
load average_OSAP_KnownC.mat
plot(average_OSAP,'-*');
cd(simu_result_path);
load average_OSAP_CLE.mat
plot(average_OSAP,'rs-');

cd(current_directory);

