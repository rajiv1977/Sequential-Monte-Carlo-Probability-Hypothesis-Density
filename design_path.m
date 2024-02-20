%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stat_save,measure_save,tar_num_save,plan_stat_save,plant_stat_unpertubed_save]=design_path(samp_t,tar_ini,target_time,...
    simu_t,tar_num_total,F,Tao,prob_d,stan_dev_n,tar_stan_dev,plantform_ini_stat,plantform_v,plantform_stan_dev)
stat_save=zeros(2*tar_num_total,simu_t);% form the matrix to save the state for each target at each time index
plan_stat_save=zeros(2,simu_t); %matrix to save the state for sensor plantform
tar_num_save=zeros(1,simu_t);%vector used to save real target number 
measure_save=zeros(tar_num_total,simu_t);%form the matrix to save the measurement for each target at each time index
for k=1:simu_t    
    plan_stat_save(1,k)=plantform_ini_stat(1)+plantform_v*k+randn*plantform_stan_dev(1);
    plan_stat_save(2,k)=plantform_ini_stat(2)+randn*plantform_stan_dev(2);
    plant_stat_unpertubed_save(1,k)=plantform_ini_stat(1)+plantform_v*k;
    plant_stat_unpertubed_save(2,k)=plantform_ini_stat(2);
    tar_num=0;
    for k1=1:tar_num_total
        if k==target_time(1,k1);
            tar_num=tar_num+1;
            stat_save(2*(k1-1)+1:2*k1,k)=tar_ini(:,k1);%initialize target state
            dump_measure=rand;
            if dump_measure<=prob_d
                %get the measurement with the detection probability prob_d
                measure_save(k1,k)=atan2(plan_stat_save(2,k),(stat_save(2*(k1-1)+1,k)-plan_stat_save(1,k)))+stan_dev_n*randn;
            end
            tar_num_save(k)=tar_num;
        elseif (k>target_time(1,k1))&(k<=target_time(2,k1))
            dump_measure=rand;
            tar_num=tar_num+1;
            random_component=tar_stan_dev*randn;
            stat_save(2*(k1-1)+1:2*k1,k)=F*stat_save(2*(k1-1)+1:2*k1,k-1)+Tao*random_component;%update state
            if dump_measure<=prob_d
                %get the measurement with the detection probability prob_d
                measure_save(k1,k)=atan2(plan_stat_save(2,k),(stat_save(2*(k1-1)+1,k)-plan_stat_save(1,k)))+stan_dev_n*randn;
            end
            tar_num_save(k)=tar_num;
        end
    end    
end