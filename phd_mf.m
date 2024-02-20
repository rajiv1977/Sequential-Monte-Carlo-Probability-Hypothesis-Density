%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_output,p_coeff_output,p_grad_output,log_lh,grad_log_lh]=phd_mf(p,p_coeff,p_grad,frame_num,w,seed,gauss_num)

load PHD_parameter.mat
load measure_buffer.mat
load prev_measure_buffer.mat
load plant_stat_unpertubed_buffer.mat
load prev_plant_stat_unpertubed_buffer.mat

randn('state',seed)
rand('state',seed)

para_len=length(w);

log_lh=zeros(1,frame_num);
grad_log_lh=zeros(para_len,frame_num);
g_range=measurement_range;

for k=1:frame_num
    p_num=length(p_coeff);
    measure=measure_buffer{k};
    prev_measure=prev_measure_buffer{k};
    [dump,measure_num]=size(measure);
    %phd resample+importance sampling for survived particle, single importance distribution is f(x_k|x_{k-1})
    N_tar=sum(p_coeff);
    p_coeff=p_coeff/N_tar;
    p_num=max(round(p_tar_num*N_tar),round(0.5*p_tar_num));
    s=cumsum(p_coeff);
    rand_dump=rand(1,p_num);     
    p_surv=zeros(4,p_num);
    c_coeff=prob_surv*N_tar/p_num*ones(1,p_num);
    d_coeff=zeros(para_len,p_num);     
    dump_sigma1=samp_t*tar_stan_dev(1,1);
    p_coeff=p_coeff*N_tar;
    dump_p_surv=p;
    dump_p_surv(1:2,:)=F*p(1:2,:);    
    dump_part_numerator=prob_surv*N_tar/p_num;
%%%%%%%%%%%%Resampling implementation method 1%%%%%%%%%%%%%%%%%%
%     for kk=1:p_num
%         rand_kk=rand_dump(kk);
%         s_dump=s-rand_kk;
%         index_dump=find(s_dump>=0,1);
%         [dump,index_dump2]=min(abs(s_dump));
%         p_surv(3,kk)=plantform_stan_dev(1)*randn;
%         p_surv(4,kk)=plantform_stan_dev(2)*randn;
%         p_surv(1:2,kk)=dump_p_surv(1:2,index_dump)+Tao*tar_stan_dev*randn;   
%     end
%%%%%%%%%%%Resampling implementation method 2%%%%%%%%%%%%%%%%%%%%
     [Number_New_Particles]=histc(rand_dump,[0 s]);
     ParticleIndex=0;
     for kk=1:length(p_coeff)
         if Number_New_Particles(kk)>0
             p_surv(3,ParticleIndex+1:ParticleIndex+Number_New_Particles(kk))=plantform_stan_dev(1)*randn(1,Number_New_Particles(kk));
             p_surv(4,ParticleIndex+1:ParticleIndex+Number_New_Particles(kk))=plantform_stan_dev(2)*randn(1,Number_New_Particles(kk));
             p_surv(1:2,ParticleIndex+1:ParticleIndex+Number_New_Particles(kk))=dump_p_surv(1:2,kk)*ones(1,Number_New_Particles(kk))+Tao*tar_stan_dev*randn(1,Number_New_Particles(kk)); 
             ParticleIndex=ParticleIndex+Number_New_Particles(kk);
         end
     end
    
    
    %Using IFGT algorithm to calculate d_coeff
    %calculate kernel sum by using predicted PHD and its derivative as weight
    InvCholQ=1/dump_sigma1;
    X_particle=(InvCholQ*dump_p_surv(2,:)).';
    Y_particle=(InvCholQ*p_surv(2,:)).';
    if any(p_grad(:))==1
        NonZeroGrad=1;
    else
        NonZeroGrad=0;
    end
    if NonZeroGrad==1
        %particle coefficients for gradient have non-zero term 
        IFGT_w=[p_coeff;p_grad]/(sqrt(2*pi)*dump_sigma1).';
        IFGT_h=sqrt(2);
        IFGT_epsil=1e-8;
        IFGT_out=computeIFGT(1,X_particle,Y_particle,IFGT_h,IFGT_w,IFGT_epsil);
        for ii=1:para_len
            d_coeff(ii,:)=dump_part_numerator*(IFGT_out(:,ii+1)./IFGT_out(:,1)).';
        end        
        %particle coefficients for gradient are all zero, then do nothing
        %to d_coeff, because under this condition d_coeff should be all zero
    end

    %generate particles for new born target based on measurements in previous scan
    [new_born_stat,new_born_coeff,new_p_total_num]=phd_new_t(p_new_tar,prob_appear,prev_measure,tar_ini_stan_dev,...
        F,Tao,tar_stan_dev,stan_dev_n,prev_plant_stat_unpertubed_buffer{k},plantform_stan_dev);
    %combine two part of particles
    p_pred_stat=[p_surv new_born_stat];
    c_coeff=[c_coeff new_born_coeff];
    %combine two part of particles
    d_coeff=[d_coeff zeros(para_len,new_p_total_num)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PHD update%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for each measurement,compute weigthed likelihood%%%%%%%%%%%%%%%%%%
    [dump,p_pred_num]=size(p_pred_stat);
    a_z=zeros(measure_num,p_pred_num);
    lh_measure=zeros(measure_num,p_pred_num);
    plant_stat_unpertubed=plant_stat_unpertubed_buffer{k};
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
    grad_kappa_k=zeros(para_len,measure_num);
    for mm=1:measure_num%calculate clutter's intensity for each measurement        
        for ii=1:gauss_num
            mc=w(3*(ii-1)+1);
            mu1=w(3*(ii-1)+2);
            sigma1=w(3*(ii-1)+3);
            kappa_k(mm)=kappa_k(mm)+cz(measure(:,mm),mc,mu1,sigma1,scaling_factor);
            grad_kappa_k(3*(ii-1)+1:3*ii,mm)=gradient_cz(measure(:,mm),mc,mu1,sigma1,scaling_factor);
        end
        kappa_k(mm)=kappa_k(mm)+scaling_factor*w(para_len)/(g_range(2)-g_range(1));
        grad_kappa_k(para_len,mm)=scaling_factor/(g_range(2)-g_range(1));        
    end
    p_coeff=(1-prob_d)*c_coeff+(1./(kappa_k+prob_d*sum_a_z))*prob_d*a_z;
    %update particle's gradient weight
    p_grad=d_coeff;
    sum_rho=zeros(para_len,measure_num);
    for kk=1:para_len      
        if NonZeroGrad==1
            dump=d_coeff(kk,:);
            sum_rho(kk,:)=dump*lh_measure.';
            p_grad(kk,:)=prob_d*dump.*((1./(kappa_k+prob_d*sum_a_z))*lh_measure)+(1-prob_d)*dump...
                -((grad_kappa_k(kk,:)+prob_d*sum_rho(kk,:))./((kappa_k+prob_d*sum_a_z).^2))*prob_d*a_z;
        else
            p_grad(kk,:)=-(grad_kappa_k(kk,:)./((kappa_k+prob_d*sum_a_z).^2))*prob_d*a_z;
        end
    end
    %calculate log-likelihood
    for ii=1:gauss_num
        mc=w(3*(ii-1)+1);
        mu1=w(3*(ii-1)+2);
        sigma1=w(3*(ii-1)+3);
        log_lh(k)=log_lh(k)-int_cz(g_range,mc,mu1,sigma1,scaling_factor);
        grad_log_lh(1+3*(ii-1):3*ii,k)=-grad_int_cz(g_range,mc,mu1,sigma1,scaling_factor);
    end
    log_lh(k)=log_lh(k)-w(para_len)*scaling_factor;
    grad_log_lh(3*gauss_num+1,k)=-scaling_factor;
    log_lh(k)=log_lh(k)-sum(c_coeff)+sum(log(kappa_k+prob_d*sum_a_z));
    grad_log_lh(:,k)=grad_log_lh(:,k)-sum(d_coeff,2)+...
        (grad_kappa_k+prob_d*sum_rho)*((1./(kappa_k+prob_d*sum_a_z)).');
    
    p=p_pred_stat;

    if k==1
        %put particle at time 1 as output
        p_output=p;
        p_coeff_output=p_coeff;
        p_grad_output=p_grad;
    end
end