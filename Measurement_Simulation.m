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

%%%%%seed for the #1 to #50 run
randn('state',137);
rand('state',761);

simu_num            =1;                     %sequence number of last input path file
samp_t              =1;                     %sample time, unit second
simu_t              =40;                    %simulation time in one experiment, unit second
tar_ini_stan_dev    =sqrt(diag([NaN;1]));   %targets' initial state's variance
stan_dev_n          =1/180*pi;              %the Standard Deviation of measurement noise,
tar_stan_dev        =0.1;                   %standard deviation of target's process noise
prob_d              =1;                     %targets detection probability
F                   =[1 samp_t;0 1];        %state evolve matrix
Tao                 =[samp_t^2/2;samp_t];   %the gain matrix for process noise
plantform_v         =4;
plantform_ini_stat  =[0;5000];
plantform_stan_dev  =[1,1];
measurement_range   =[0 pi];
prob_surv           =0.975;             %target's survive probability
prob_appear         =0.1;               %the parameter for target's appearing poisson process (per sample time)
p_tar_num           =1500;              %particle number/target
p_new_tar           =500;               %particle number/measurement

%generate target state and its corresponding measurement
tar_num_total       =3;
tar_ini             =[800 1400 2500;
                      10  -12  -15];
target_time         =[1  15 20;
                      30 40 40];
%generate clutter 
false_d1            =5;                 %false measurement density, number/sample_time
area_num1           =1;
clut_para1          =[0.9 0.5*pi 12/180*pi];
simu_t1             =40;
cluster_treshold    =9;


current_directory   =pwd;
simulated_data_path =[current_directory,'\Simulated_Input_Data'];
save PHD_parameter.mat;

for kk=1:simu_num
    [stat_save,measure_save,tar_num_save,plan_stat_save,plant_stat_unpertubed_save]=design_path(samp_t,tar_ini,target_time,...
        simu_t,tar_num_total,F,Tao,prob_d,stan_dev_n,tar_stan_dev,plantform_ini_stat,plantform_v,plantform_stan_dev);
    cd(simulated_data_path);
    eval(['save stat_save_run_',num2str(kk),'.mat stat_save;']);
    eval(['save measure_save_run_',num2str(kk),'.mat measure_save;']);
    eval(['save tar_num_save_run_',num2str(kk),'.mat tar_num_save;']);
    eval(['save plan_stat_save_run_',num2str(kk),'.mat plan_stat_save;']);
    eval(['save plant_stat_unpertubed_save_run_',num2str(kk),'.mat plant_stat_unpertubed_save;']);
    cd(current_directory);
end

for kk=1:simu_num
    clutter_measure1=clutter_simu_ellip(area_num1,measurement_range,1,simu_t1,false_d1,clut_para1);
    clutter_measure=[clutter_measure1];
    cd(simulated_data_path);
    eval(['save clutter_measure_run_',num2str(kk),'.mat clutter_measure;']);
    cd(current_directory);
end