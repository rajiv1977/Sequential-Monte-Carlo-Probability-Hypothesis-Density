function [g_ifgt,p_max,K,r]=computeIFGT(d,x,y,h,q,epsil,p1,p2,p3)

% The improved fast gauss transform helps in a fast evaluation of the sums
% of the form Ghat(yj)=sum(qi*e^{|xi-yj|^2/h^2}, j=1...M
% We evaluate G(yj) such that Ghat(yj)-G(yj)|<=Q*eps; where Q=sum(qi)
% 
% Syntax: 
% 
% [g_ifgt,p_max,K,r]=computeIFGT(d,x,y,h,q,epsil,p1,p2,p3)
% 
% Inputs: 
% 
% d        : Data dimension
% x        :  Sources or the set of xi's as above, size N*d or d*N
% y        :  Targets or the set of yi's as above, size M*d or d*M
% h        :  bandwidth, a scalar weight
% q        :  weights or the set of qi's as above, size W*d or d*W, 
%       where W is the number of weights for which the Gauss sum need to be evaluated
% epsil  :  Desired error rate
% p1, p2 and p3 are optional inputs 
% The combination of their values help decide their functionality:
% 
% If you want to specify a limit to the truncation number, 
%         then set p1=desired_truncation and p2='TruncationNumber'
% If you want to specify a value for the cluster size, 
%         then set p1=desired_cluster_size and p2='ClusterSize'
% If you want to specify both the truncation number limit and cluster size,
%         then set p1=desired_truncation and p2=desired_cluster_size 
% If you want to specify truncation number limit, cluster size and cluster radius,
%         then set p1=desired_truncation, p2=desired_cluster_size and 
%         p3 = desired_cluster_radius
% 
% Outputs:
% G_ifgt        : Gauss summation computed by IFGT algorithm
% p_max        : Maximum truncation number
% K             : Cluster size
% r             : Cluster radius
%  
%  References:
% 
% [1] L.F. Greengard and X. Sun, 
% A new version of the fast Gauss transform, 
% Doc. Math. J. DMV Extra Volume ICM (1998) 575–584.
% [2] Yang, C., Duraiswami, R., and Gumerov, N. 
% Improved fast Gauss transform. 
% CS-TR-4495, Dept. of Computer Science, 
% University of Maryland, College Park, 2003
% [3] V.C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov. 
% Fast computation of sums of gaussians in high dimensions. 
% Technical Report CS-TR-4767, Department of Computer Science, 
% University of Maryland, CollegePark, 2005.

%Check parameters for validity
if nargin<6,
    disp('Some mandatory inputs are not specified');
    g_ifgt=NaN;
    return;
else
    if size(x,2)==d && size(x,1)~=d
        x=x';
    end
    if size(y,2)==d && size(y,1)~=d
        y=y';
    end
    if size(x,1)~=d
        g_ifgt=NaN;
        disp('Source dimension not matching given dimensions');
        return;
    elseif size(y,1)~=d
        g_ifgt=NaN;
        disp('Target dimension not matching given dimensions');
        return;
    else
        N=size(x,2);
        M=size(y,2);
    end
    if all(size(h)~=1)
        g_ifgt=NaN;
        disp('h needs to be a scalar');        
        return;
    end
    if any(size(q)==N)        
        if size(q,1)==N
            q=q';
        end
    else
        g_ifgt=NaN;
        disp('Weights q must be of the same size as the sources');
        return;
    end
    W=size(q,1);
end

q=q.';
% Automatic Parameter selection
if nargin==6
    if nargout==1
        g_ifgt = figtree(d,N,M,W,x,h,q,y,epsil);            
    else
        maxRange = max(max([x y],[],2) - min([x y],[],2));
        [K,p_max,r]=figtreeChooseParametersNonUniform(d,N,x,h,epsil,N,maxRange);
        [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, x, K);
        pMax = figtreeChooseTruncationNumber(d, h, epsil, rx,maxRange);
        g_ifgt=figtreeEvaluateIfgt(d, N, M, W, x,h, q, y, pMax, K, clusterIndex, clusterCenter, clusterRadii, r, epsil);
    end
end

% Semi Automatic Parameter selection
if nargin==8
    if ischar(p2)
        if strcmp(p2,'ClusterSize')
            maxRange = max(max([x y],[],2) - min([x y],[],2));
            [K,p_max,r]=figtreeChooseParametersNonUniform(d,N,x,h,epsil,N,maxRange);
            [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, x, p1);
            pMax = figtreeChooseTruncationNumber(d, h, epsil, rx,maxRange);                        
            g_ifgt=figtreeEvaluateIfgt(d, N, M, W, x,h, q, y, pMax, K, clusterIndex, clusterCenter, clusterRadii, r, epsil);
        elseif strcmp(p2,'TruncationNumber')
            maxRange = max(max([x y],[],2) - min([x y],[],2));
            [K,p_max,r]=figtreeChooseParametersNonUniform(d,N,x,h,epsil,N,maxRange);
            [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, x, K);
            pMax=p1;
            g_ifgt=figtreeEvaluateIfgt(d, N, M, W, x,h, q, y, pMax, K, clusterIndex, clusterCenter, clusterRadii, r, epsil);
        else
            disp('Unrecognized string');
            g_ifgt=NaN;
            return;
        end
    else
        maxRange = max(max([x y],[],2) - min([x y],[],2));
        [K,p_max,r]=figtreeChooseParametersNonUniform(d,N,x,h,epsil,N,maxRange);
        [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, x, p2);
        pMax=p1;        
        g_ifgt=figtreeEvaluateIfgt(d, N, M, W, x,h, q, y, pMax, K, clusterIndex, clusterCenter, clusterRadii, r, epsil);
    end
end

if nargin == 9    
    [K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, x, p2);
    pMax=p1; 
    r=p3;
    g_ifgt=figtreeEvaluateIfgt(d, N, M, W, x,h, q, y, pMax, K, clusterIndex, clusterCenter, clusterRadii, r, epsil);
end