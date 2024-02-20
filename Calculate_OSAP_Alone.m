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
simu_num=50;
simu_t=40;
average_wass=zeros(1,simu_t);
average_estimate_num=zeros(1,simu_t);
average_abs_error_estimate_num=zeros(1,simu_t);
relative_situation_path=['\Simulation_Result\Summary_File'];
situation_name=['CLE'];
current_directory=pwd;
situation_path=[current_directory,relative_situation_path];
simulated_input_path=[current_directory,'\Simulated_Input_Data'];
cd(simulated_input_path);
load tar_num_save_run_1.mat;
for simu_k=1:simu_num
    cd(simulated_input_path);
    eval(['load stat_save_run_',num2str(simu_k),'.mat;']);
    cd(situation_path)
    eval(['load tar_est_save_',situation_name,'_run',num2str(simu_k),'.mat;']);
    eval(['load N_k_k_save_',situation_name,'_run',num2str(simu_k),'.mat;']);
    for k=1:simu_t
        est_i=find(tar_est_save(:,k)~=0);%get the index for estimation;
        tar_i=find(stat_save(:,k)~=0);%get the index for real target state
        real_tar_num=round(0.5*length(tar_i));
        est_tar_num=round(0.25*length(est_i));
        %%%%%%%%%%%%%%%%%%%%%Calculate OSAP%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        OSAP_c=10;
        OSAP_order=2;
        if (est_tar_num==0)&(real_tar_num==0)
            OSAP_metric(k)=0;
        elseif ((est_tar_num==0)&(real_tar_num~=0))|((est_tar_num~=0)&(real_tar_num==0))
            OSAP_metric(k)=OSAP_c;
        else
            est_stat=tar_est_save(est_i,k);%targets' state estimation
            tar_stat=stat_save(tar_i,k);%real target state
            est_stat_pos=zeros(1,est_tar_num);%targets' position estimate
            for est_k=1:est_tar_num
                est_stat_pos(:,est_k)=[est_stat(4*(est_k-1)+1)];
            end
            real_stat_pos=zeros(1,real_tar_num);%targets' true position
            for tar_k=1:real_tar_num
                real_stat_pos(:,tar_k)=[tar_stat(2*(tar_k-1)+1)];
            end
            PriceMatrixDim=max(real_tar_num,est_tar_num);
            PriceMatrix=zeros(PriceMatrixDim);
            for tar_k=1:PriceMatrixDim
                for est_k=1:PriceMatrixDim
                    if (tar_k<=real_tar_num)&(est_k<=est_tar_num)
                        PriceMatrix(tar_k,est_k)=min(norm(real_stat_pos(:,tar_k)-est_stat_pos(:,est_k),OSAP_order)^OSAP_order,OSAP_c^OSAP_order);
                    else
                        PriceMatrix(tar_k,est_k)=OSAP_c^OSAP_order;
                    end
                end
            end
            [assignment, cost]=assignmentoptimal(PriceMatrix);
            OSAP_metric(k)=(cost/PriceMatrixDim).^(1/OSAP_order);
        end
    end
    eval(['save OSAP_metric_',situation_name,'_run',num2str(simu_k),'.mat OSAP_metric;']);
end
cd(current_directory)