%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f G]=plhgrad(w)
global g_p_stat g_p_coeff g_p_gradient g_frame_num g_gauss_num out_p_stat out_p_coeff out_p_gradient out_log_lh;
[p_stat,p_coeff,p_gradient,log_lh,grad_log_lh]=phd_mf(g_p_stat,g_p_coeff,g_p_gradient,g_frame_num,w,0,g_gauss_num);
f=-sum(log_lh);
G=-sum(grad_log_lh,2);
out_p_stat=p_stat;
out_p_coeff=p_coeff;
out_p_gradient=p_gradient;
out_log_lh=log_lh;
% save p_stat.mat p_stat;
% save p_coeff.mat p_coeff
% save p_gradient.mat p_gradient
% save log_lh.mat log_lh