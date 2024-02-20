%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tar_est=state_extract_V2(stat_dump1,coeff_dump1,peak_num,differ_p_num)
load PHD_parameter.mat
%use improved k-mean algorithm to pick up peak
if peak_num==0
    %in case there is no peak_num needs to be picked up
    tar_est=zeros(4,1);
else
    %in case there is some targets need to be picked up
    tar_est=[];
    tar_already_k=0;%target already peaked
    %clustering p_stat based on each particles position%%%%%%%
    Y = pdist(stat_dump1(1,:)','euclid');
    Z = linkage(Y,'average');
    T = cluster(Z,'cutoff',cluster_treshold,'criterion','distance');
    while tar_already_k<peak_num
        %get expected target number in every cluster
        tar_num_clus=zeros(1,max(T));%sum
        for dump_kk=1:max(T)
            tar_num_clus(dump_kk)=sum(coeff_dump1(find(T==dump_kk)));
        end
        bigger_half=find(tar_num_clus>0.5);
        if isempty(bigger_half)==0            
            max_tries=10;%maximum allowed number for iteration in extraction one peak
            [sort_tar_num_clus,sort_I] = sort(tar_num_clus,'descend');
            for dump_kk=1:max(T)
                tar_est_temp=[];
                if round(sort_tar_num_clus(dump_kk))~=0
                    %use CLEAN algorithm to pickup peak
                    dump_stat=stat_dump1(:,find(T==sort_I(dump_kk)));
                    dump_coeff=coeff_dump1(find(T==sort_I(dump_kk)));
                    sum_dump_coeff=sum(dump_coeff);
                    dump_coeff=dump_coeff/(sum(dump_coeff))*sort_tar_num_clus(dump_kk);
                    target_weight=min([0.9,0.9*sum(dump_coeff)/round(sort_tar_num_clus(dump_kk))]);
                    max_outlier=0;%maximum allowed outlier
                    kk=1;
                    while kk<=round(sort_tar_num_clus(dump_kk))%for kk=1:round(sort_tar_num_clus(dump_kk))
                        [max_coeff max_ind]=max(dump_coeff);%get maximum particle coefficient
                        max_pos=dump_stat(:,max_ind);%state vector correspondes to the maximum particl coefficient
                        if max_coeff<target_weight
                            iter_num=1;
                            dump_flag=0;
                            while (iter_num<max_tries)&(dump_flag==0)
                                rad=sqrt(cluster_treshold)/3*sqrt(iter_num);%get the radius of neighborhood
                                [area_i]=find(abs(dump_stat(1,:)-max_pos(1))<=rad);
                                area_coeff=dump_coeff(area_i);%coefficients for those particles in neighborhood
                                area_stat=dump_stat(:,area_i);%state vectors for those particles in neighborhood
                                neighbor_coeff=sum(area_coeff);
                                if neighbor_coeff>=target_weight
                                    dump_flag=1;
                                end
                                iter_num=iter_num+1;
                            end
                            peak_pos=(area_coeff*area_stat.').'/neighbor_coeff;%get the peak coefficient by the weighted sum
                            if (dump_flag==0)
                                if max_outlier==0
                                    outlier_ind=find(dump_coeff==max_coeff);
                                    dump_coeff(area_i)=0;
                                    tar_est_temp(:,end+1)=[peak_pos;neighbor_coeff];
                                    tar_already_k=tar_already_k+1;
                                    kk=kk+1;
                                    coeff_dump1(find(T==sort_I(dump_kk)))=dump_coeff*sum_dump_coeff/sort_tar_num_clus(dump_kk);
                                else
                                    outlier_ind=find(dump_coeff==max_coeff);
                                    dump_coeff(outlier_ind)=0;
                                    max_outlier=max_outlier-1;
                                    coeff_dump1(find(T==sort_I(dump_kk)))=dump_coeff*sum_dump_coeff/sort_tar_num_clus(dump_kk);
                                end
                            else
                                dump_coeff(area_i)=dump_coeff(area_i)*(1-target_weight/neighbor_coeff);
                                tar_est_temp(:,end+1)=[peak_pos;neighbor_coeff];
                                tar_already_k=tar_already_k+1;
                                kk=kk+1;
                                coeff_dump1(find(T==sort_I(dump_kk)))=dump_coeff*sum_dump_coeff/sort_tar_num_clus(dump_kk);
                            end
                        else
                            tar_est_temp(:,end+1)=[max_pos;target_weight];
                            dump_coeff(max_ind)=dump_coeff(max_ind)-target_weight;
                            tar_already_k=tar_already_k+1;
                            kk=kk+1;
                            coeff_dump1(find(T==sort_I(dump_kk)))=dump_coeff*sum_dump_coeff/sort_tar_num_clus(dump_kk);
                        end
                    end
                    [dump,I]=sort(tar_est_temp(5,:),'descend');
                    TarNum_tar_est_temp=size(tar_est_temp,2);
                    TarNum_tar_est=size(tar_est,2);
                    if ((peak_num-TarNum_tar_est)>TarNum_tar_est_temp)&(peak_num>TarNum_tar_est)
                        tar_est(:,TarNum_tar_est+1:TarNum_tar_est+TarNum_tar_est_temp)=tar_est_temp(1:4,:);
                    elseif peak_num>TarNum_tar_est
                        tar_est(:,TarNum_tar_est+1:peak_num)=tar_est_temp(1:4,I(1:peak_num-TarNum_tar_est));
                    end
                end
            end
        else
            [dump,dump_i]=max(tar_num_clus);
            tar_num_clus(dump_i)=0;
            dump_stat=stat_dump1(:,find(T==dump_i));
            dump_coeff=coeff_dump1(find(T==dump_i));
            tar_already_k=tar_already_k+1;
            TarNum_tar_est=size(tar_est,2);  
            tar_est(:,TarNum_tar_est+1)=dump_coeff*dump_stat.'/sum(dump_coeff);            
        end
    end
end