%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_output,p_coeff_output,total_num]=phd_new_t(num_p,prob_appear,prev_measure,...
    tar_ini_stat,tar_ini_stan_dev,F,Tao,tar_stan_dev,w,g_range)
%for each measurement, generate num_p particle
[dump,num_measure]=size(prev_measure);%number of measurements in previous scan
total_num=num_measure*num_p;%number of particles generated for new born target
p_output=zeros(4,total_num);
for k=1:num_measure
    temp=zeros(4,num_p);
    temp(1,:)=prev_measure(1,k);
    temp(3,:)=prev_measure(2,k);
    temp(2,:)=tar_ini_stat(2)+tar_ini_stan_dev(2,2)*randn(1,num_p);
    temp(4,:)=tar_ini_stat(4)+tar_ini_stan_dev(4,4)*randn(1,num_p);
    temp=F*temp+Tao*tar_stan_dev*randn(2,num_p);
    p_output(:,num_p*(k-1)+1:num_p*k)=temp;
end
p_coeff_output=zeros(1,total_num);
surv_area=(g_range(2)-g_range(1))*(g_range(4)-g_range(3));
%we assume true PHD for new born targets is uniform in whole surveillance region
for k=1:total_num
    dump_pos=p_output([1:2:3],k);
    if (dump_pos(1)>=w(2))&(dump_pos(1)<=w(3))&(dump_pos(2)>=w(4))&(dump_pos(2)<=w(5))
        kappa=w(1)/((w(3)-w(2))*(w(5)-w(4)));
    else
        kappa=(1-w(1))/surv_area;
    end
    p_coeff_output(k)=1/total_num*(prob_appear/surv_area)/...
        (1/surv_area+kappa);
end

