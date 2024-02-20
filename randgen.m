%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=randgen(m,cov,N)
% mu,means,D-by-1 dimension vector
% cov, covariance matrix, D-by-D, cov must be semi-positive definite
% N, the length of output random variable 
% x, D-dimensions gaussian variable, (in D-by-N matrix)

function x=randgen(mu,cov,N)
[dim_var,dump]=size(cov);% get the dimension for the output random variable
[V,D]= eig(cov);
temp=randn(dim_var,N);
x=V*sqrt(D)*temp+mu*ones(1,N);