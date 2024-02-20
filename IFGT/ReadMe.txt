Documentation for computeIFGT

The improved fast gauss transform helps in a fast evaluation of the sums
  of the form Ghat(yj)=sum(qi*e^{|xi-yj|^2/h^2}, j=1...M
  We evaluate G(yj) such that Ghat(yj)-G(yj)|<=Q*eps; where Q=sum(qi)
  
  Syntax: 
  
  [g_ifgt,p_max,K,r]=computeIFGT(d,x,y,h,q,epsil,p1,p2,p3)
  
  Inputs: 
  
  d        : Data dimension
  x        :  Sources or the set of xi's as above, size N*d or d*N
  y        :  Targets or the set of yi's as above, size M*d or d*M
  h        :  bandwidth, a scalar weight
  q        :  weights or the set of qi's as above, size W*d or d*W, 
        where W is the number of weights for which the Gauss sum need to be evaluated
  epsil  :  Desired error rate
  p1, p2 and p3 are optional inputs 
  The combination of their values help decide their functionality:
  
  If you want to specify a limit to the truncation number, 
          then set p1=desired_truncation and p2='TruncationNumber'
  If you want to specify a value for the cluster size, 
          then set p1=desired_cluster_size and p2='ClusterSize'
  If you want to specify both the truncation number limit and cluster size,
          then set p1=desired_truncation and p2=desired_cluster_size 
  If you want to specify truncation number limit, cluster size and cluster radius,
          then set p1=desired_truncation, p2=desired_cluster_size and 
          p3 = desired_cluster_radius
  
  Outputs:
  G_ifgt         : Gauss summation computed by IFGT algorithm
  p_max        : Maximum truncation number
  K                 : Cluster size
  r                   : Cluster radius