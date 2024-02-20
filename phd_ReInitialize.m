%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_stat,p_coeff1,p_grad]=phd_ReInitialize(p,p_coeff,measure,prev_measure,w,seed,gauss_num)

load PHD_parameter.mat

randn('state',seed)
rand('state',seed)

para_len=length(w);

g_range=[simu_range(1) simu_range(2) simu_range(1) simu_range(2)];
[dump,p_num]=size(p);
[dump,measure_num]=size(measure);
%phd resample+importance sampling for survived particle, single importance distribution is f(x_k|x_{k-1})
N_tar=sum(p_coeff);
p_coeff=p_coeff/N_tar;
p_num=max(round(p_tar_num*N_tar),round(0.5*p_tar_num));
[dump,p_num_max]=size(p);
p_num=min(p_num,p_num_max);
s=cumsum(p_coeff);
rand_dump=rand(1,p_num);
p_surv=zeros(4,p_num);
c_coeff=zeros(1,p_num);
dump_sigma1=samp_t*tar_stan_dev(1,1);
dump_sigma2=samp_t*tar_stan_dev(2,2);
p_coeff=p_coeff*N_tar;
for kk=1:p_num
    rand_kk=rand_dump(kk);
    s_dump=s-rand_kk;
    index_dump=find(s_dump>=0,1,'first');
    p_surv(:,kk)=F*p(:,index_dump)+Tao*tar_stan_dev*randn(2,1);
    c_coeff(kk)=N_tar/p_num;
end
%using finite differential method to finde initial value for
%partial difference of D(x)
%generate particles for new born target based on measurements in previous scan
[new_born_stat,new_born_coeff,new_p_total_num]=phd_new_t(p_new_tar,prob_appear,prev_measure,...
    tar_ini_stat,tar_ini_stan_dev,F,Tao,tar_stan_dev,stan_dev_n);
%combine two part of particles
p_pred_stat=[p_surv new_born_stat];
c_coeff=[c_coeff new_born_coeff];
%PHD update%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for each measurement,compute weigthed likelihood%%%%%%%%%%%%%%%%%%
[dump,p_pred_num]=size(p_pred_stat);
a_z=zeros(measure_num,p_pred_num);
lh_measure=zeros(measure_num,p_pred_num);
for kk=1:measure_num
    z=measure(:,kk);
    x=H*p_pred_stat;
    z_x_lh1=normpdf(z(1)*ones(1,p_pred_num),x(1,:),stan_dev_n*ones(1,p_pred_num));
    %compute likelihood for measure with index kk, x-coordinate
    z_x_lh2=normpdf(z(2)*ones(1,p_pred_num),x(2,:),stan_dev_n*ones(1,p_pred_num));
    %likelehood function of z(with index kk), conditioned on p_pred_state
    lh_measure(kk,:)=z_x_lh1.*z_x_lh2;
    a_z(kk,:)=lh_measure(kk,:).*c_coeff;
end
%update particle's PHD weight
update_coeff=zeros(1,p_pred_num);
sum_a_z=(a_z*ones(p_pred_num,1)).';
kappa_k=zeros(1,measure_num);
for mm=1:measure_num%calculate clutter's intensity for each measurement
    for ii=1:gauss_num
        mc=w(5*(ii-1)+1);
        mu1=w(5*(ii-1)+2);
        sigma1=w(5*(ii-1)+3);
        mu2=w(5*(ii-1)+4);
        sigma2=w(5*(ii-1)+5);
        kappa_k(mm)=kappa_k(mm)+cz(measure(:,mm),mc,mu1,sigma1,mu2,sigma2);
    end
    kappa_k(mm)=kappa_k(mm)+w(para_len)/((g_range(2)-g_range(1))*(g_range(4)-g_range(3)));
end
p_coeff1=(1-prob_d)*c_coeff+(1./(kappa_k+prob_d*sum_a_z))*prob_d*a_z;

%based on two different set of parameters to calculate kappa_k and
%using finite difference to initialize partial difference of D(x)
p_grad=zeros(para_len,p_pred_num);
for jj=1:para_len
    kappa_k_plus=zeros(1,measure_num);
    w_plus=w+0.5*[zeros(1,jj-1) 1 zeros(1,para_len-jj)];
    for mm=1:measure_num%calculate clutter's intensity for each measurement
        for ii=1:gauss_num
            mc=w_plus(5*(ii-1)+1);
            mu1=w_plus(5*(ii-1)+2);
            sigma1=w_plus(5*(ii-1)+3);
            mu2=w_plus(5*(ii-1)+4);
            sigma2=w_plus(5*(ii-1)+5);
            kappa_k_plus(mm)=kappa_k_plus(mm)+cz(measure(:,mm),mc,mu1,sigma1,mu2,sigma2);
        end
        kappa_k_plus(mm)=kappa_k_plus(mm)+w_plus(para_len)/((g_range(2)-g_range(1))*(g_range(4)-g_range(3)));
    end
    p_coeff_plus=(1-prob_d)*c_coeff+(1./(kappa_k_plus+prob_d*sum_a_z))*prob_d*a_z;
    kappa_k_neg=zeros(1,measure_num);
    w_neg=w-0.5*[zeros(1,jj-1) 1 zeros(1,para_len-jj)];
    for mm=1:measure_num%calculate clutter's intensity for each measurement
        for ii=1:gauss_num
            mc=w_neg(5*(ii-1)+1);
            mu1=w_neg(5*(ii-1)+2);
            sigma1=w_neg(5*(ii-1)+3);
            mu2=w_neg(5*(ii-1)+4);
            sigma2=w_neg(5*(ii-1)+5);
            kappa_k_neg(mm)=kappa_k_neg(mm)+cz(measure(:,mm),mc,mu1,sigma1,mu2,sigma2);
        end
        kappa_k_neg(mm)=kappa_k_neg(mm)+w_neg(para_len)/((g_range(2)-g_range(1))*(g_range(4)-g_range(3)));
    end
    p_coeff_neg=(1-prob_d)*c_coeff+(1./(kappa_k_neg+prob_d*sum_a_z))*prob_d*a_z;
    p_grad(jj,:)=(p_coeff_plus-p_coeff_neg);    
end

p_stat=p_pred_stat;