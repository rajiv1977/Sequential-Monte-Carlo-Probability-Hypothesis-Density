%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_temp=hessian_cz(z,mc,mu1,sigma1,mu2,sigma2)
%input z should be an 2*1 vector
out_temp=zeros(5,5);
out=zeros(15,1);
dump_cz=cz(z,mc,mu1,sigma1,mu2,sigma2);
dump_g_cz=gradient_cz(z,mc,mu1,sigma1,mu2,sigma2);
out_temp(1,1)=0;
out_temp(2:5,1)=dump_g_cz(1)*dump_g_cz(2:5)/dump_cz;
out_temp(1,2:5)=out_temp(2:5,1)';
out_temp(2,2)=(-1/sigma1^2+(z(1)-mu1)^2/sigma1^4)*dump_cz;
out_temp(2,3)=(-3*(z(1)-mu1)/sigma1^3+(z(1)-mu1)^3/sigma1^5)*dump_cz;
out_temp(2,4:5)=dump_g_cz(2)*dump_g_cz(4:5)/dump_cz;
out_temp(3:5,2)=out_temp(2,3:5)';
out_temp(3,3)=((z(1)-mu1)^4/sigma1^6-5*(z(1)-mu1)^2/sigma1^4+2/sigma1^2)*dump_cz;
out_temp(3,4:5)=dump_g_cz(3)*dump_g_cz(4:5)/dump_cz;
out_temp(4:5,3)=out_temp(3,4:5)';
out_temp(4,4)=(-1/sigma2^2+(z(2)-mu2)^2/sigma2^4)*dump_cz;
out_temp(4,5)=(-3*(z(2)-mu2)/sigma2^3+(z(2)-mu2)^3/sigma2^5)*dump_cz;
out_temp(5,4)=out_temp(4,5);
out_temp(5,5)=((z(2)-mu2)^4/sigma2^6-5*(z(2)-mu2)^2/sigma2^4+2/sigma2^2)*dump_cz;
% index=1;
% for kk=1:5
%     for ii=kk:5
%         out(index)=out_temp(kk,ii);
%         index=index+1;
%     end
% end

