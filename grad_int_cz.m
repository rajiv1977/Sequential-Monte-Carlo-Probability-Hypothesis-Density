%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_temp=grad_int_cz(a,mc,mu1,sigma1,scaling_factor)%,hm)
out_temp=zeros(3,1);
dump=int_cz(a,mc,mu1,sigma1,scaling_factor);%,hm)-hm;
norm_x=normpdf([a(1),a(2)],mu1,sigma1);
erf_a1=erf((a(1)-mu1)/sqrt(2*sigma1^2));
erf_a2=erf((a(2)-mu1)/sqrt(2*sigma1^2));
out_temp(1)=dump/mc;
out_temp(2)=mc*scaling_factor*(norm_x(1)-norm_x(2));
out_temp(3)=mc*scaling_factor/sigma1*((mu1-a(2))*norm_x(2)-(mu1-a(1))*norm_x(1));
% out_temp(6)=1;