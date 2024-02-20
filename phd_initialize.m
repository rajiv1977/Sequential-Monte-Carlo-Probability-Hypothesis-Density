%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_stat,p_coeff,p_grad]=phd_initialize(measure,prev_measure,w,seed,gauss_num,...
    plant_stat_unpertubed_prev,plant_stat_unpertubed,plantform_stan_dev)

load PHD_parameter.mat

randn('state',seed)
rand('state',seed)

para_len=length(w);

g_range=measurement_range;

[dump,measure_num]=size(measure);
%using finite differential method to finde initial value for
%partial difference of D(x)
%generate particles for new born target based on measurements in previous scan
[p_pred_stat,c_coeff,new_p_total_num]=phd_new_t(p_new_tar,prob_appear,prev_measure,...
    tar_ini_stan_dev,F,Tao,tar_stan_dev,stan_dev_n,plant_stat_unpertubed_prev,plantform_stan_dev);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHD update%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for each measurement,compute weigthed likelihood%%%%%%%%%%%%%%%%%%
%     C_z=zeros(1,measure_num);
[dump,p_pred_num]=size(p_pred_stat);
a_z=zeros(measure_num,p_pred_num);
lh_measure=zeros(measure_num,p_pred_num);
for kk=1:measure_num
    z=measure(:,kk);
    x=atan2((plant_stat_unpertubed(2)+p_pred_stat(4,:)),(p_pred_stat(1,:)-plant_stat_unpertubed(1)-p_pred_stat(3,:)));
    %compute likelihood for measure with index kk
    lh_measure(kk,:)=normpdf(z*ones(1,p_pred_num),x,stan_dev_n*ones(1,p_pred_num));
    a_z(kk,:)=lh_measure(kk,:).*c_coeff;
end
%update particle's PHD weight
update_coeff=zeros(1,p_pred_num);
sum_a_z=(a_z*ones(p_pred_num,1)).';
kappa_k=zeros(1,measure_num);
for mm=1:measure_num%calculate clutter's intensity for each measurement
    for ii=1:gauss_num
        mc=w(3*(ii-1)+1);
        mu1=w(3*(ii-1)+2);
        sigma1=w(3*(ii-1)+3);
        kappa_k(mm)=kappa_k(mm)+cz(measure(:,mm),mc,mu1,sigma1,scaling_factor);
    end
    kappa_k(mm)=kappa_k(mm)+scaling_factor*w(para_len)/(g_range(2)-g_range(1));
end
p_coeff=(1-prob_d)*c_coeff+(1./(kappa_k+prob_d*sum_a_z))*prob_d*a_z;

%based on two different set of parameters to calculate kappa_k and
%using finite difference to initialize partial difference of D(x)
p_grad=zeros(para_len,p_pred_num);
for jj=1:para_len
    kappa_k_plus=zeros(1,measure_num);
    if mod(jj,3)==1
        w_plus=w+0.5/scaling_factor*[zeros(1,jj-1) 1 zeros(1,para_len-jj)];
    else
        w_plus=w+0.005*pi*[zeros(1,jj-1) 1 zeros(1,para_len-jj)];
    end
    for mm=1:measure_num%calculate clutter's intensity for each measurement
        for ii=1:gauss_num
            mc=w_plus(3*(ii-1)+1);
            mu1=w_plus(3*(ii-1)+2);
            sigma1=w_plus(3*(ii-1)+3);
            kappa_k_plus(mm)=kappa_k_plus(mm)+cz(measure(:,mm),mc,mu1,sigma1,scaling_factor);
        end
        kappa_k_plus(mm)=kappa_k_plus(mm)+scaling_factor*w_plus(para_len)/(g_range(2)-g_range(1));
    end
    p_coeff_plus=(1-prob_d)*c_coeff+(1./(kappa_k_plus+prob_d*sum_a_z))*prob_d*a_z;
    kappa_k_neg=zeros(1,measure_num);
    if mod(jj,3)==1
        w_neg=w-0.5/scaling_factor*[zeros(1,jj-1) 1 zeros(1,para_len-jj)];
    else
        w_neg=w-0.005*pi*[zeros(1,jj-1) 1 zeros(1,para_len-jj)];
    end
    for mm=1:measure_num%calculate clutter's intensity for each measurement
        for ii=1:gauss_num
            mc=w_neg(3*(ii-1)+1);
            mu1=w_neg(3*(ii-1)+2);
            sigma1=w_neg(3*(ii-1)+3);
            kappa_k_neg(mm)=kappa_k_neg(mm)+cz(measure(:,mm),mc,mu1,sigma1,scaling_factor);
        end
        kappa_k_neg(mm)=kappa_k_neg(mm)+scaling_factor*w_neg(para_len)/(g_range(2)-g_range(1));
    end
    p_coeff_neg=(1-prob_d)*c_coeff+(1./(kappa_k_neg+prob_d*sum_a_z))*prob_d*a_z;
    if mod(jj,3)==1
        p_grad(jj,:)=(p_coeff_plus-p_coeff_neg)/scaling_factor;    
    else
        p_grad(jj,:)=(p_coeff_plus-p_coeff_neg)/(0.01*pi);    
    end
end

p_stat=p_pred_stat;
