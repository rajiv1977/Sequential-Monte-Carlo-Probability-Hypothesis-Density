%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this program simulate the target path and its respective measurement used for the simulation of PHD filter
%there is only one observer, located at the origin, its output is targets' position 
%the targets are moving on a 2-D plane, based on same constant velocity model
%the possible target area is [-150,150]-by-[-150,150]
%all the unit is SI unit
clear
clc
randn('state',3)
rand('state',31)
simu_num=100;%total number of simulation test
samp_t=1;%sample time, unit second
simu_t=40;%simulation time in one experiment, unit second
tar_ini_stat=[0;3;0;-3];%the initial state for any target
tar_ini_stan_dev=sqrt(diag([10;1;10;1]));%targets' initial state's variance
tar_stan_dev=diag([1;0.1]);%standard deviation of target's process noise
stan_dev_n=2.5; %the Standard Deviation of measurement noise
prob_surv=0.95;%target's survive probability
prob_appear=0.2;%the parameter for target's appearing poisson process (per sample time)
F=[1 samp_t 0 0;0 1 0 0;0 0 1 samp_t;0 0 0 1];%state evolve matrix
H=[1 0 0 0;0 0 1 0];%measurement matrix
Tao=[samp_t^2/2 0;samp_t,0;0,samp_t^2/2;0,samp_t];%the gain matrix for process noise
max_tar_num=20;%the maximum allowed number for total target in this simulation
prob_d=1;%targets detection probability

for simu_k=1:simu_num
    tar_num=0;%the number of targets' at time K    
    tar_ind=zeros(1,max_tar_num);%indicate the place in stat_save for those target survivied at time k-1
    tar_num_save=zeros(1,simu_t);%vector used to save real target number  
    dump=random('poiss',prob_appear,1,10*simu_t);
    one_begin_ind=find(dump~=0,1,'first');
    new_tar_ind=dump(one_begin_ind:one_begin_ind+simu_t-1);%index for those scans where new target begins, assume at scan1, there is always a target born
    new_tar_ind(find(new_tar_ind~=0))=1;%assume at each scan, only one new taget could bore at most
    tar_num_total=length(find(new_tar_ind~=0));;%the number of all targets that have been appeared
    stat_save=zeros(4*tar_num_total,simu_t);% form the matrix to save the state for each target at each time index
    measure_save=zeros(2*tar_num_total,simu_t);%form the matrix to save the measurement for each target at each time index
    dump=rand(tar_num_total,10*simu_t);
    end_tar_ind=zeros(tar_num_total,simu_t);%index for those scans where target end, here we make sure there is not any targt whose survive length is less than 4 (prob_surv^4=0.8)
    dd_dump=round(log10(0.8)/log10(prob_surv));%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for dump_k=1:tar_num_total
        dump_ind=find(dump(dump_k,:)>=prob_surv);
        dump_lag=dump_ind-[0,dump_ind(1:length(dump_ind)-1)];
        dump_ind=[0,dump_ind(1:length(dump_ind)-1)];
        dump_ind=dump_ind(find(dump_lag>=dd_dump,1,'first'));
        end_tar_ind(dump_k,:)=dump(dump_k,dump_ind+1:dump_ind+simu_t);        
    end
    tar_surv_len=zeros(1,tar_num_total);%vector used to count each target's survive length
    tar_num_total=0;%re-initialize tar_num_total to 0, here it represents the number of all targets that have been appeared sp far
    for k=1:simu_t
        measure=zeros(1,1);%initilaize measurement vector for this time k
        measure_num=0;
        terminate_flag=0;%indicate at time k, whether there are targets have been terminated
        %update state for those already existing target%%%%%%%%%%%%%%%%%%%%%%%%
        if tar_num>0 %if there are some target already existing
            for kk=1:tar_num
                stat_ind=tar_ind(kk);%get target kk's index in stat_save
                dump=tar_surv_len(stat_ind);%get target kk's track length
                temp_dump=end_tar_ind(stat_ind,dump);
                if temp_dump<=prob_surv%test whether the target will survive at this time
                    stat_vec=stat_save(4*(stat_ind-1)+1:4*stat_ind,k-1);%fellowing stat_ind, get the state for target kk from stat_save
                    v=tar_stan_dev*randn(2,1);%get process noise
                    stat_vec=F*stat_vec+Tao*v;%state update
                    stat_save(4*(stat_ind-1)+1:4*stat_ind,k)=stat_vec;%save the state for target kk in stat_save
                    tar_surv_len(stat_ind)=tar_surv_len(stat_ind)+1;%target kk's track length plus one
                else
                    tar_ind(kk)=0;%indicate the target kk has been terminated
                    tar_num=tar_num-1;%change the target number at time k
                    terminate_flag=1;
                end
            end
            if terminate_flag>0%erase the place for thoes has been terminated at k
                temp_dump=zeros(1,max_tar_num);
                dump_k=0;
                for kk=1:max_tar_num
                    if tar_ind(kk)~=0
                        dump_k=dump_k+1;
                        temp_dump(dump_k)=tar_ind(kk);
                    end
                end
                tar_ind=temp_dump;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %generate the spontaneously appearing target at time k%%%%%%%%%%%%%%%%%
        %     if k==1 %assuming there is always a target start at time 1
        %         appear_num=1;
        %     else
        appear_num=new_tar_ind(k);%the number for spontaneously appearing target at time k
        if (appear_num>0)&(k<simu_t-dd_dump)%at end of simulation test, do not generate new target
            if (tar_num_total+appear_num)<=max_tar_num %test whether the total target number has been excedded the maximum allowed number
                tar_num_total=tar_num_total+appear_num;%increase number of all targets that have been appeared
                for kk=1:appear_num
                    tar_num=tar_num+1;%change the target number at time k
                    stat_vec=tar_ini_stat+tar_ini_stan_dev*randn(4,1);
                    stat_save((tar_num_total-1)*4+1:tar_num_total*4,k)=stat_vec;
                    tar_ind(tar_num)=tar_num_total;
                    tar_surv_len(tar_num_total)=tar_surv_len(tar_num_total)+1;%target kk's track length plus one
                end
            end
        end
        tar_num_save(k)=tar_num;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate measurement from target%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kk=1:tar_num
            temp=rand;
            if temp<=prob_d %for the simulation of non-perfect detection
                measure_num=measure_num+1;
                stat_ind=tar_ind(kk);
                stat_vec=stat_save(4*(stat_ind-1)+1:4*stat_ind,k);
                measure_vec=H*stat_vec+stan_dev_n*randn(2,1);%get the measurement with the detection probability prob_d
                measure_save(2*(stat_ind-1)+1:2*stat_ind,k)=measure_vec;%save target kk's measurement at time k
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
    eval(['save stat_save',num2str(simu_k),'.mat stat_save;']);
    eval(['save measure_save',num2str(simu_k),'.mat measure_save;']);
    eval(['save tar_num_save',num2str(simu_k),'.mat tar_num_save;']);
end