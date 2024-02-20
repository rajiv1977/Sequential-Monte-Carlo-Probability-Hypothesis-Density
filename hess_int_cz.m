%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_temp=hess_int_cz(a,mc,mu1,sigma1,mu2,sigma2)
%input a should be an [x1,x2,y1,y2] and represnets a surveillance
%region which confined in [x1,x2]*[y1,y2]
out_temp=zeros(5,5);
out=zeros(15,1);
dump=int_cz(a,mc,mu1,sigma1,mu2,sigma2);
dump_grad=grad_int_cz(a,mc,mu1,sigma1,mu2,sigma2);
norm_x=normpdf([a(1),a(2)],mu1,sigma1);
norm_y=normpdf([a(3),a(4)],mu2,sigma2);
erf_a1=erf((a(1)-mu1)/sqrt(2*sigma1^2));
erf_a2=erf((a(2)-mu1)/sqrt(2*sigma1^2));
erf_a4=erf((a(4)-mu2)/sqrt(2*sigma2^2));
erf_a3=erf((a(3)-mu2)/sqrt(2*sigma2^2));
out_temp(1,1)=0;
out_temp(2:5,1)=dump_grad(1)*dump_grad(2:5)/dump;
out_temp(1,2:5)=out_temp(2:5,1)';
out_temp(2,2)=0.5*mc*(erf_a4-erf_a3)*(-(a(2)-mu1)*norm_x(2)+(a(1)-mu1)*norm_x(1))/sigma1^2;
out_temp(2,3)=mc/(2*sigma1)*(erf_a4-erf_a3)*((1-(a(2)-mu1)^2/sigma1^2)*norm_x(2)-(1-(a(1)-mu1)^2/sigma1^2)*norm_x(1));
out_temp(2,4)=mc*(-norm_y(2)+norm_y(1))*(-norm_x(2)+norm_x(1));
out_temp(2,5)=mc/(sigma2)*(-norm_x(2)+norm_x(1))*((-a(4)+mu2)*norm_y(2)-(-a(3)+mu2)*norm_y(1));
out_temp(3:5,2)=out_temp(2,3:5)';
out_temp(3,3)=-dump_grad(3)/sigma1+mc/(2*sigma1^4)*(erf_a4-erf_a3)*((-a(2)+mu1)*((a(2)-mu1)^2-sigma1^2)*norm_x(2)-...
    (-a(1)+mu1)*((a(1)-mu1)^2-sigma1^2)*norm_x(1));
out_temp(3,4)=mc/(sigma1)*(-norm_y(2)+norm_y(1))*((-a(2)+mu1)*norm_x(2)-(-a(1)+mu1)*norm_x(1));
out_temp(3,5)=mc/(sigma1*sigma2)*((-a(2)+mu1)*norm_x(2)-(-a(1)+mu1)*norm_x(1))*...
    ((-a(4)+mu2)*norm_y(2)-(-a(3)+mu2)*norm_y(1));
out_temp(4:5,3)=out_temp(3,4:5)';
out_temp(4,4)=0.5*mc*(erf_a2-erf_a1)*(-(a(4)-mu2)*norm_y(2)+(a(3)-mu2)*norm_y(1))/sigma2^2;
out_temp(4,5)=mc/(2*sigma2)*(erf_a2-erf_a1)*((1-(a(4)-mu2)^2/sigma2^2)*norm_y(2)-(1-(a(3)-mu2)^2/sigma2^2)*norm_y(1));
out_temp(5,4)=out_temp(4,5);
out_temp(5,5)=-dump_grad(5)/sigma2+mc/(2*sigma2^4)*(erf_a2-erf_a1)*((-a(4)+mu2)*((a(4)-mu2)^2-sigma2^2)*norm_y(2)-...
    (-a(3)+mu2)*((a(3)-mu2)^2-sigma2^2)*norm_y(1));
% index=1;
% for kk=1:5
%     for ii=kk:5
%         out(index)=out_temp(kk,ii);
%         index=index+1;
%     end
% end
