%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stat_dump,stat_dump1,coeff_dump,coeff_dump1,differ_p_num]=phd_resample_V2(p,p_coeff)
load PHD_parameter.mat
%phd resample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_k_k=sum(p_coeff);
p_coeff=p_coeff/N_k_k;
p_num=max(round(p_tar_num*N_k_k),500);
s=cumsum(p_coeff);
rand_dump=rand(1,p_num);
stat_dump=zeros(4,p_num);%save particle after resample
coeff_dump=zeros(1,p_num);%save coefficient after resample
stat_dump1=[];%save particles after resample, whose state is mutually different 
coeff_dump1=[];%save coefficient corresponding those particles after resample, whose state is mutually different
save_index_dump=[];
for kk=1:p_num
    rand_kk=rand_dump(kk);
    s_dump=s-rand_kk;
    index_dump=find(s_dump>=0,1,'first');
    coeff_dump(kk)=1/p_num;
    stat_dump(:,kk)=p(:,index_dump);
    if kk==1
        save_index_dump(end+1)=index_dump;
        coeff_dump1(end+1)=1/p_num;
        stat_dump1(:,end+1)=p(:,index_dump);
    else
        sub_temp=find((save_index_dump-index_dump)==0);
        if isempty(sub_temp)
            %point has not been picked yet
            save_index_dump(end+1)=index_dump;
            coeff_dump1(end+1)=1/p_num;            
            stat_dump1(:,end+1)=p(:,index_dump);
        else
            %point has already been picked
            coeff_dump1(sub_temp)=coeff_dump1(sub_temp)+1/p_num;
        end
    end
end
coeff_dump=coeff_dump*N_k_k;
coeff_dump1=coeff_dump1*N_k_k;
differ_p_num=length(coeff_dump1);