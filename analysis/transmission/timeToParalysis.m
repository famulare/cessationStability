function [cdf,pdf,t]=timeToParalysis(t)
% Casey 1942 (Figure 2)
% The incubation period in epidemic poliomyelitis
% JAMA
    if nargin==0
        t=1:70;
    end

    d=[1,3.5,9.5,15.5,21.5,27.5,33.5,39.5,max(t)];
    n=cumsum([0,0,4,17,6,1,1,0,0])/29;
    cdf=interp1(d,n,t);
    pdf=[0,diff(cdf)];

end