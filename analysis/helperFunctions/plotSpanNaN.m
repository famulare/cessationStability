function h=plotSpanNaN(x,y,linestr,varargin)
% span nans when plotting
if nargin<3
    linestr='-';
end

for k=1:size(y,1)
    idx=~isnan(y(k,:));
    X=x(idx);
    Y=y(k,idx);
    plot(X,Y,linestr,varargin{:});
    hold on
end
end
