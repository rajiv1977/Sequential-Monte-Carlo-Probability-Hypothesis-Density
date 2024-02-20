clc
clear all
clear functions
clear mex
close all


d=2;
N=1000;
M=1000;
d_q=26;
x=randn(N,d);
y=randn(M,d);
q=rand(N,d_q);
h=sqrt(2);
epsil=10^-6;

%naive method
naive_out=zeros(N,d_q);
tic
for kk=1:N
    %calculate transition pdf%%%%%%%%%%%%%%%%%%%%
    x_pdf=normpdf(x(:,1),y(kk,1));
    y_pdf=normpdf(x(:,2),y(kk,2));
    comb_pdf=x_pdf.*y_pdf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj=1:d_q
        naive_out(kk,jj)=comb_pdf.'*q(:,jj);
    end   
end
toc

tic
%IFGT
g1=computeIFGT(d,x,y,h,q,epsil);
IFGT_out=g1/(2*pi);
toc
