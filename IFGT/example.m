clc
clear all
clear functions
clear mex
close all


d=3;
N=1000;
M=1000;
x=randn(N,d);
y=randn(M,d);
q=rand(N,3);
h=10;
epsil=10^-6;

g1=computeIFGT(d,x,y,h,q,epsil);
disp('Evaluated by Automatic Parameter Selection');
disp('Automatic Parameter Selection with Parameters returned');
[g2,p,K,r]=computeIFGT(d,x,y,h,q,epsil);
disp(sprintf('Error wrt first method: %d',sum(sum(g1-g2))))

disp('Semi Automatic Parameter Selection - specifying Truncation Number');
g3=computeIFGT(d,x,y,h,q,epsil,p,'TruncationNumber');
disp(sprintf('Error wrt first method: %d',sum(sum(g1-g3))))

disp('Semi Automatic Parameter Selection - specifying Cluster Size');
g4=computeIFGT(d,x,y,h,q,epsil,K,'ClusterSize');
disp(sprintf('Error wrt first method: %d',sum(sum(g1-g4))))

disp('Semi Automatic Parameter Selection - specifying Truncation Number and Cluster Size ');
g4=computeIFGT(d,x,y,h,q,epsil,p,K);
disp(sprintf('Error wrt first method: %d',sum(sum(g1-g4))))

disp('Explicit parameter specification - specifying Truncation Number, Cluster Size and Cluster Radius');
g4=computeIFGT(d,x,y,h,q,epsil,p,K,r);
disp(sprintf('Error wrt first method: %d',sum(sum(g1-g4))))