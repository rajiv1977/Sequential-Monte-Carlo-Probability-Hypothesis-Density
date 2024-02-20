%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tar_est=state_extract(p_stat,stat_dump1,p_coeff,coeff_dump1,peak_num,differ_p_num)
load PHD_parameter.mat
[dump,p_num]=size(p_stat);
stat_dump=p_stat;
%use improved k-mean algorithm to pick up peak
if peak_num==0
    tar_est=zeros(4,1);;
else
    tar_est=zeros(4,peak_num);
    %clustering p_stat based on each particles position%%%%%%%
    Y = pdist(stat_dump1(1:2:3,:)','euclid');
    Z = linkage(Y,'average');
    T1 = cluster(Z,'cutoff',(3*stan_dev_n)^2,'criterion','distance');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %base on clustering result, get the cluster index for p_coeff and p_stat
    T=zeros(1,p_num);
    for kkk=1:differ_p_num
        temp_p=stat_dump1(:,kkk);
        temp_p_mat=temp_p*ones(1,p_num);
        temp_error=sum(abs(stat_dump-temp_p_mat));
        same_stat_ind=find(temp_error==0);
        T(same_stat_ind)=T1(kkk);
    end
    %get expected target number in every cluster
    tar_num_clus=zeros(1,max(T));%sum
    for dump_kk=1:max(T)
        tar_num_clus(dump_kk)=sum(p_coeff(find(T==dump_kk)));
    end
    total_tar_num_temp=round(sum(tar_num_clus));%total target number, NOTE, here it may be bigger than peak_num
    tar_est_temp=zeros(5,total_tar_num_temp);
    max_tries=10;%maximum allowed number for iteration in extraction one peak
    tar_already_k=1;%target already peaked
    for dump_kk=1:max(T)
        if round(tar_num_clus(dump_kk))~=0
            %use CLEAN algorithm to pickup peak
            dump_stat=p_stat(:,find(T==dump_kk));
            dump_coeff=coeff_dump1(find(T==dump_kk));
            dump_coeff=dump_coeff/(sum(dump_coeff))*tar_num_clus(dump_kk);
            target_weight=min([0.9,0.9*sum(dump_coeff)/round(tar_num_clus(dump_kk))]);
            max_outlier=0;%maximum allowed outlier
            kk=1;
            while kk<=round(tar_num_clus(dump_kk))%for kk=1:round(tar_num_clus(dump_kk))
                [max_coeff max_ind]=max(dump_coeff);%get maximum particle coefficient
                max_pos=dump_stat(:,max_ind);%state vector correspondes to the maximum particl coefficient
                if max_coeff<target_weight
                    iter_num=1;
                    dump_flag=0;
                    while (iter_num<max_tries)&(dump_flag==0)
                        rad=stan_dev_n*sqrt(iter_num);%get the radius of neighborhood
                        [area_i]=find((dump_stat(1,:)-max_pos(1)).^2+(dump_stat(3,:)-max_pos(3)).^2<=rad^2);
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
                            dump_coeff(outlier_ind)=0;
                            dump_coeff(area_i)=dump_coeff(area_i)*(1-target_weight/neighbor_coeff);
                            tar_est_temp(:,tar_already_k)=[peak_pos;neighbor_coeff];
                            tar_already_k=tar_already_k+1;
                            kk=kk+1;
                        else
                            outlier_ind=find(dump_coeff==max_coeff);
                            dump_coeff(outlier_ind)=0;
                            max_outlier=max_outlier-1;
                        end
                    else
                        dump_coeff(area_i)=dump_coeff(area_i)*(1-target_weight/neighbor_coeff);
                        tar_est_temp(:,tar_already_k)=[peak_pos;neighbor_coeff];
                        tar_already_k=tar_already_k+1;
                        kk=kk+1;
                    end
                else
                    tar_est_temp(:,tar_already_k)=[max_pos;max_pos];
                    max_coeff(max_ind)=max_coeff(max_ind)-target_weight;
                    tar_already_k=tar_already_k+1;
                    kk=kk+1;
                end
            end
        end
        [dump,I]=sort(tar_est_temp(5,:),'descend');
        tar_est=tar_est_temp(1:4,I(1:peak_num));
    end
end