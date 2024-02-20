%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_temp=gradient_cz(z,mc,mu1,sigma1,scaling_factor);%,hm,g_range)
%input z should be an 2*D matrix
[dump,size_z]=size(z);
out_temp=zeros(3,size_z);
dump_cz=cz(z,mc,mu1,sigma1,scaling_factor);%,hm,g_range)-hm/((g_range(2)-g_range(1))*(g_range(4)-g_range(3)));
out_temp(1,:)=dump_cz./mc;
out_temp(2,:)=(z(1,:)-mu1)/sigma1^2.*dump_cz;
% out_temp(4,:)=(z(2,:)-mu2)/sigma2^2.*dump_cz;
out_temp(3,:)=((z(1,:)-mu1).^2-sigma1^2)/sigma1^3.*dump_cz;
% out_temp(5,:)=((z(2,:)-mu2).^2-sigma2^2)/sigma2^3.*dump_cz;
% out_temp(6,:)=ones(1,size_z)/((g_range(2)-g_range(1))*(g_range(4)-g_range(3)));
