%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                       %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.                  %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                                %
%                                                                                                %
% This program do the simulation for PHD filter (with important targets' path and measurements)  %
% there is only one observer, located at the origin, its output is targets' position             %
% the targets are moving on a 2-D plane, based on same constant velocity model                   %
% the possible target area is [-100,100]-by-[-100,100]                                           %
% all the unit is SI unit                                                                        %
% all the unit is SI unit                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
Measurement_Simulation();
load PHD_parameter.mat;

for simu_k=1:simu_num
    simu_k
    save simu_k.mat simu_k;     %saving the MCR #
    
    clear;
    clc
    
    load simu_k.mat;
    load PHD_parameter.mat;     %loading the above saved info
    cd(simulated_data_path)     %go to the spcified path
    
    eval(['load stat_save_run_',num2str(simu_k),'.mat;']);
    eval(['load measure_save_run_',num2str(simu_k),'.mat;']);
    eval(['load tar_num_save_run_',num2str(simu_k),'.mat;']);
    eval(['load clutter_measure_run_',num2str(simu_k),'.mat;']);
    %   eval(['load plan_stat_save_run_',num2str(kk),'.mat;']);
    eval(['load plant_stat_unpertubed_save_run_',num2str(simu_k),'.mat;']);
    
    cd(current_directory);
    randn('state',37)
    rand('state',13)
    
    N_k_k_save=zeros(1,simu_t);         %vector used to save estimated target number
    peak_num_save=zeros(1,simu_t);      %save the estimated targets' number in every scan
    tar_est_save=zeros(4*20,simu_t);    %matrix used to store the estimation,NOTICE here
    p_stat=zeros(4,1);                  %initialize the particle state matrix
    p_coeff=1;                          %initialize the particle's coefficient
    tar_num_error=zeros(1,simu_t);      %initialize the vector for absolute error of targets' number's estimation
    ass_dist=zeros(1,simu_t);           %wassterstein distance beteween estimation and real targets' state
    
    for k=1:simu_t
        k
        
        if k>1
            prev_measure=measure;
            prev_plant_stat_unpertubed=plant_stat_unpertubed;
        end
        
        non_zero_ind=find(measure_save(:,k)~=0);        %find the non zero value at time index k in measure_save
        tar_measure=measure_save(non_zero_ind,k);       %get measure for the real targets
        measure=[tar_measure' clutter_measure{k}];      %get measurement, include tragets' and clutter's
        [dump,false_num]=size(clutter_measure{k});      %get number of measurement
        [dump,measure_num]=size(measure);               %get number of clutter
        random_index=randperm(measure_num);
        measure=measure(:,random_index);                %do random permutation for measurement
        plant_stat_unpertubed=plant_stat_unpertubed_save(:,k);
        
        if k>1
            
            if k<=simu_t1
                clut_para=clut_para1;
                false_d=false_d1;
                clutter_area_type=1;
            elseif k<=(simu_t1+simu_t2)
                clut_para=clut_para2;
                false_d=false_d2;
                clutter_area_type=1;
            end
            
            [p_stat,p_coeff]=phd_mf_KnowC(p_stat,p_coeff,clut_para,false_d,0,measure,prev_measure,k,clutter_area_type,plant_stat_unpertubed,prev_plant_stat_unpertubed);
            
            N_k_k=sum(p_coeff);
            N_k_k
            N_k_k_save(k)=N_k_k;
            peak_num=round(N_k_k);
            peak_num_save(k)=peak_num;
            
            [stat_dump,stat_dump1,coeff_dump,coeff_dump1,differ_p_num]=phd_resample_V2(p_stat,p_coeff);
            tar_est=state_extract_V2(stat_dump1,coeff_dump1,peak_num,differ_p_num);
            
            dump=tar_est;
            [dump_row,dump_col]=size(dump);
            len_dump=dump_row*dump_col;
            tar_est_save(1:len_dump,k)=dump(:);
            
            
            %calculate the wassterstein distance between the estimation and real target state
            est_i=find(tar_est_save(:,k)~=0);           %get the index for estimation;
            tar_i=find(stat_save(:,k)~=0);              %get the index for real target state
            real_tar_num=round(0.5*length(tar_i));
            est_tar_num=round(0.25*length(est_i));
            if (est_tar_num~=0)&&(real_tar_num~=0)
                
                est_stat=tar_est_save(est_i,k);         %targets' state estimation
                tar_stat=stat_save(tar_i,k);            %real target state
                
                %construct price matrix
                price_mat=zeros(real_tar_num,est_tar_num);
                for tar_k=1:real_tar_num
                    for est_k=1:est_tar_num
                        price_mat(tar_k,est_k)=abs(tar_stat(2*(tar_k-1)+1)-est_stat(4*(est_k-1)+1));
                    end
                end
                
                if est_tar_num==real_tar_num
                    %if estimation of target number is correct,calculaue wasserstein distance based on optimal assignment problem
                    [assignment, cost]=assignmentoptimal(price_mat);
                    wass_dist(k)=cost/real_tar_num;
                else
                    %if estimation of target number is incorrect,calculauewasserstein distance based on simplex method (Linear Convex Optimazation)
                    f=price_mat(:);
                    Aeq=zeros(est_tar_num*real_tar_num,est_tar_num+real_tar_num);
                    for aeq_k=1:(est_tar_num+real_tar_num)
                        if aeq_k<=est_tar_num
                            Aeq((aeq_k-1)*real_tar_num+1:aeq_k*real_tar_num,aeq_k)=ones(real_tar_num,1);
                        else
                            Aeq((aeq_k-est_tar_num):real_tar_num:est_tar_num*real_tar_num,aeq_k)=ones(est_tar_num,1);
                        end
                    end
                    beq=[1/est_tar_num*ones(1,est_tar_num) 1/real_tar_num*ones(1,real_tar_num)];
                    lb=zeros(est_tar_num*real_tar_num,1);
                    [C,fval]=linprog(f,[],[],Aeq.',beq.',lb,[]);
                    wass_dist(k)=fval;
                end
            else
                wass_dist(k)=NaN;
            end
            
        end
    end
    
    current_directory   =pwd;
    simulated_data_path =[current_directory,'\Simulation_Result'];
    cd(simulated_data_path);
    eval(['save tar_est_save_KnownC_run',num2str(simu_k),'.mat tar_est_save;']);
    eval(['save N_k_k_save_KnownC_run',num2str(simu_k),'.mat N_k_k_save;']);
    eval(['save wass_dist_KnownC_run',num2str(simu_k),'.mat wass_dist']);
    cd(current_directory);
    
end


figure(1)
stem(peak_num_save(1:simu_t),'x');
hold on
stem(tar_num_save(1:simu_t),'s');

figure(2)
stem(wass_dist)

figure(3)
hold on
plot([1:30],stat_save(1,1:30),'-x','linewidth',2);
plot([15:40],stat_save(3,15:40),'--*','linewidth',2);
plot([25:40],stat_save(5,25:40),'--o','linewidth',2);
all_zero_flag=0;
k=0;
plot_temp_y1=[];
plot_temp_x1=[];
t_index=1:simu_t;
while all_zero_flag==0
    dump=tar_est_save(4*k+1,:);
    non_zero_ind=find(dump~=0);
    if isempty(non_zero_ind)==0
        plot_temp_y1=[plot_temp_y1,tar_est_save(4*k+1,non_zero_ind)];
        plot_temp_x1=[plot_temp_x1,t_index(non_zero_ind)];
        k=k+1;
    else
        all_zero_flag=1;
    end
end
plot(plot_temp_x1,plot_temp_y1,'ks','MarkerSize',5);

figure(4)
hold on
plot(N_k_k_save(1:simu_t),'-x');
plot(tar_num_save(1:simu_t));
