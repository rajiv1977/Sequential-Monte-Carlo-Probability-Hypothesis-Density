%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_output,p_coeff_output]=phd_mf(p,p_coeff,clut_para,false_d,seed,measure,prev_measure,simu_k,clutter_area_type,plant_stat_unpertubed,prev_plant_stat_unpertubed)

load PHD_parameter.mat  %initial infor from the main body

randn('state',seed)  % i need to know what the heck is this
rand('state',seed)

[dump,p_num]=size(p);  %find the state of the target
[dump,measure_num]=size(measure);   %find the size of the measure

if simu_k>2
    %phd resample+importance sampling for survived particle, single importance distribution is f(x_k|x_{k-1})
    N_tar=sum(p_coeff);
    p_coeff=p_coeff/N_tar;
    p_num=max(round(p_tar_num*N_tar),round(0.5*p_tar_num));
    s=cumsum(p_coeff);
    rand_dump=rand(1,p_num);
    p_surv=zeros(4,p_num);
    c_coeff=prob_surv*N_tar/p_num*ones(1,p_num);
    dump_sigma1=samp_t*tar_stan_dev(1,1);
    p_coeff=p_coeff*N_tar;
    dump_p_surv=p;
    dump_p_surv(1:2,:)=F*p(1:2,:);    
    dump_part_numerator=prob_surv*N_tar/p_num;
% %%%%%%%%%%%Resampling implementation method 1%%%%%%%%%%%%%%%%%%
%     for kk=1:p_num
%         rand_kk=rand_dump(kk);
%         s_dump=s-rand_kk;
%         index_dump=find(s_dump>=0,1,'first');
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
    %phd importance sampling for new born target
    %generate particles for new born target based on measurements in previous scan
    [new_born_stat,new_born_coeff,new_p_total_num]=phd_new_t(p_new_tar,prob_appear,prev_measure,tar_ini_stan_dev,...
        F,Tao,tar_stan_dev,stan_dev_n,prev_plant_stat_unpertubed,plantform_stan_dev);
    %combine two part of particles
    p_pred_stat=[p_surv new_born_stat];
    c_coeff=[c_coeff new_born_coeff];
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
    area_num=size(clut_para,1);
    if clutter_area_type==1%elliptic dense clutter area
        for mm=1:measure_num%calculate clutter's intensity for each measurement
            for ww=1:area_num
                kappa_k(mm)=kappa_k(mm)+clut_para(ww,1)*false_d*normpdf(measure(1,mm),clut_para(ww,2),clut_para(ww,3));
            end
            kappa_k(mm)=kappa_k(mm)+false_d*(1-sum(clut_para(:,1)))/(measurement_range(2)-measurement_range(1));
        end
    elseif clutter_area_type==2%rectangular dense clutter area
        for mm=1:measure_num%calculate clutter's intensity for each measurement
            calculated_flag=0;
            for ww=1:area_num
                if (measure(1,mm)>=clut_para(ww,2))&(measure(1,mm)<=clut_para(ww,3))&(calculated_flag==0)
                    kappa_k(mm)=clut_para(ww,1)*false_d/(clut_para(ww,3)-clut_para(ww,2));
                    calculated_flag=1;
                end
            end
            if calculated_flag==0
                kappa_k(mm)=false_d*(1-sum(clut_para(:,1)))/(measurement_range(2)-measurement_range(1));
                calculated_flag=1;
            end
        end
    end
    p_coeff=(1-prob_d)*c_coeff+(1./(kappa_k+prob_d*sum_a_z))*prob_d*a_z;

    p_output=p_pred_stat;
    p_coeff_output=p_coeff;
    
else
    
    %generate particles for new born target based on measurements in previous scan
    [p_pred_stat,c_coeff,new_p_total_num]=phd_new_t(p_new_tar,prob_appear,prev_measure,tar_ini_stan_dev,...
        F,Tao,tar_stan_dev,stan_dev_n,prev_plant_stat_unpertubed,plantform_stan_dev);
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
    area_num=size(clut_para,1);
    if clutter_area_type==1%elliptic dense clutter area
        for mm=1:measure_num%calculate clutter's intensity for each measurement
            for ww=1:area_num
                kappa_k(mm)=kappa_k(mm)+clut_para(ww,1)*false_d*normpdf(measure(1,mm),clut_para(ww,2),clut_para(ww,3));
            end
            kappa_k(mm)=kappa_k(mm)+false_d*(1-sum(clut_para(:,1)))/(measurement_range(2)-measurement_range(1));
        end
    elseif clutter_area_type==2%rectangular dense clutter area
        for mm=1:measure_num%calculate clutter's intensity for each measurement
            calculated_flag=0;
            for ww=1:area_num
                if (measure(1,mm)>=clut_para(ww,2))&(measure(1,mm)<=clut_para(ww,3))&(calculated_flag==0)
                    kappa_k(mm)=clut_para(ww,1)*false_d/(clut_para(ww,3)-clut_para(ww,2));
                    calculated_flag=1;
                end
            end
            if calculated_flag==0
                kappa_k(mm)=false_d*(1-sum(clut_para(:,1)))/(measurement_range(2)-measurement_range(1));
                calculated_flag=1;
            end
        end
    end
    p_coeff=(1-prob_d)*c_coeff+(1./(kappa_k+prob_d*sum_a_z))*prob_d*a_z;

    p_output=p_pred_stat;
    p_coeff_output=p_coeff;
end






