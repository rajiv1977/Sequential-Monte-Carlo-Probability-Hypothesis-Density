%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stat_dump,stat_dump1,coeff_dump,coeff_dump1,differ_p_num]=phd_resample(p,p_coeff)
load PHD_parameter.mat
%phd resample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_k_k=sum(p_coeff);
p_coeff=p_coeff/N_k_k;
p_num=max(round(p_tar_num*N_k_k),500);
s=cumsum(p_coeff);
rand_dump=rand(1,p_num);
stat_dump=zeros(4,p_num);
stat_dump1=zeros(4,p_num);%particles` state for cluster
coeff_dump=zeros(1,p_num);
coeff_dump1=ones(1,p_num);%particles` coefficient fo cluster
save_index_dump=zeros(1,p_num);
for kk=1:p_num
    rand_kk=rand_dump(kk);
    s_dump=s-rand_kk;
    index_dump=find(s_dump>=0,1,'first');
    save_index_dump(kk)=index_dump;
    coeff_dump(kk)=1/p_num;%p_up_coeff(index_dump);
    coeff_dump1(kk)=p_coeff(index_dump);
    stat_dump(:,kk)=p(:,index_dump);
end
coeff_dump=coeff_dump*N_k_k;
%scanning stat_dump and make sure only those with different state
%vector would be put into stat_dump1
differ_p_num=0;
temp_coeff=coeff_dump1;
while length(find(isnan(temp_coeff)==0))>=1
    differ_p_num=differ_p_num+1;
    non_nan_ind=isnan(temp_coeff);
    dump_p_ind=find(non_nan_ind==0,1,'first');
    temp_p=stat_dump(:,dump_p_ind);
    stat_dump1(:,differ_p_num)=temp_p;
    temp_p_mat=temp_p*ones(1,p_num);
    temp_error=sum(abs(stat_dump-temp_p_mat));
    same_stat_ind=find(temp_error==0);
    temp_coeff(same_stat_ind)=NaN;
end
stat_dump1=stat_dump1(:,1:differ_p_num);