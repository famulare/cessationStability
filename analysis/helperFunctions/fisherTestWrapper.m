function [h,p,stats]=fisherTestWrapper(pos1,denom1,pos2,denom2,tailStr,alpha)
if nargin<5
    tailStr='both';
end
if nargin<6
    alpha=0.05;
end
[h,p,stats]=fishertest([pos1,denom1-pos1;pos2,denom2-pos2],'tail',tailStr,'alpha',alpha);

end