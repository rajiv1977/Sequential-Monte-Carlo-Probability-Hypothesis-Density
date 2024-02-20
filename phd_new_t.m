%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_output,p_coeff_output,total_num]=phd_new_t(num_p,prob_appear,prev_measure,...
    tar_ini_stan_dev,F,Tao,tar_stan_dev,stan_dev_n,plant_stat_unpertubed,plantform_stan_dev)
%for each measurement, generate num_p particle
[dump,num_measure]=size(prev_measure);%number of measurements in previous scan
total_num=num_measure*num_p;%number of particles generated for new born target    num_p == p_new_tar=500;%particle number/measurement
p_output=zeros(4,total_num);
p_coeff_output=zeros(1,total_num);

% for k=1:num_measure
%     temp=zeros(4,num_p);
%     temp_coeff=zeros(1,num_p);    
%     v_x=0;
%     stan_dev_v_x=tar_ini_stan_dev(2,2);
%     for kk=1:num_p
%         temp(2,kk)=v_x+stan_dev_v_x*randn;
%         temp(3,kk)=plantform_stan_dev(1)*randn;
%         temp(4,kk)=plantform_stan_dev(2)*randn;
%         temp(1,kk)=plant_stat_unpertubed(1)+temp(3,kk)...
%             +(plant_stat_unpertubed(2)+temp(4,kk))/tan(prev_measure(k));
%     end
%     temp(1:2,:)=F*temp(1:2,:)+Tao*tar_stan_dev*randn(1,num_p);
%     p_output(:,num_p*(k-1)+1:num_p*k)=temp;
% end

for k=1:num_measure
    temp=zeros(4,num_p);
    temp_coeff=zeros(1,num_p);
    x=plant_stat_unpertubed(1)+plant_stat_unpertubed(2)/tan(prev_measure(k));   %sensor location in x and y coordin / tan(previous meas)
    v_x=0;
    stand_dev_x=sqrt(plantform_stan_dev(1)^2+1/tan(prev_measure(k))^2*plantform_stan_dev(2)^2+...
        +plant_stat_unpertubed(2)^2/sin(prev_measure(k))^4*stan_dev_n^2);
    stan_dev_v_x=tar_ini_stan_dev(2,2);
    for kk=1:num_p
        temp(1,kk)=x+stand_dev_x*randn;
        temp(2,kk)=v_x+stan_dev_v_x*randn;
        temp(3,kk)=plantform_stan_dev(1)*randn;
        temp(4,kk)=plantform_stan_dev(2)*randn;
    end
    temp(1:2,:)=F*temp(1:2,:)+Tao*tar_stan_dev*randn(1,num_p);
    p_output(:,num_p*(k-1)+1:num_p*k)=temp;
end
p_coeff_output=prob_appear/total_num*ones(1,total_num);