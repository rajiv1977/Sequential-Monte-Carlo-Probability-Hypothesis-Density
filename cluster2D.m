%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Sequential Monte Carlo Probability Hypothesis Density (SMC-PHD)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               S.Rajiv, and T.Kirubarajan                                  %
%           ECE Dept., McMaster University, Hamilton, Ontario, L8S 4K1, Canada.             %
%                         sithirr@mcmaster.ca, kiruba@mcmaster.ca                           %
%                                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outIndex centroid] = cluster2D(points, no_of_clusters)

[dummy, no_of_points] = size(points);

% initialise clusters
index = ones(no_of_points,1) + round((no_of_clusters-1)*rand(no_of_points,1));

for i=1:no_of_clusters
    % find the centroid
    clusterindex = find(index==i);    
    centroid(:,i) = mean(points(:,clusterindex)')';
end
oldIndex = [];


while(~isequal(oldIndex, index))
    oldIndex = index;
    for i=1:no_of_points       
        distance = [];
        for j = 1:no_of_clusters
            distance = [distance; norm(centroid(:,j)-points(:,i))];
        end
        % find the closest centriod and change the cluster index
        [value clusterindex] = min(distance);
        index(i) = clusterindex;
	end
	
	for i=1:no_of_clusters
        % find the centroid
        clusterindex = find(index==i);    
        centroid(:,i) = mean(points(:,clusterindex)')';
	end
end

outIndex = index';
