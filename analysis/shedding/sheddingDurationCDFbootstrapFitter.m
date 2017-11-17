function [beta,CI,yCI]=sheddingDurationCDFbootstrapFitter(X,Y,N,reps)

if nargin<4
    reps=100;
end

idx=any(~isnan(Y),2);
X=X(idx);
Y=Y(idx);
N=N(idx);

custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(sheddingDurationCDF(2^b,X,0)) + N.*(1-Y).*log(1-sheddingDurationCDF(2^b,X,0)));
[beta,CI] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',0) ;

betaBoot=zeros(1,reps);
X=unique(X);
N=sum(unique(N))*ones(size(X));
Yboot=zeros(length(X),reps);
custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(sheddingDurationCDF(2^b,X,0)) + N.*(1-Y).*log(1-sheddingDurationCDF(2^b,X,0)));
% figure; hold on;
for k=1:reps
    tmp=nan(length(X),1);
    for n=1:length(X)
        idx=rand<=1-sheddingDurationCDF(2^beta,X,0);
        if any(idx)
            tmp(n)=min(X(idx));
        end
    end
    Yboot(:,k)=(length(tmp)-cumsum(hist(tmp,X)))'/length(tmp);
%     plot(X,Yboot)
    [betaBoot(k)] = mle(Yboot(:,k),'nloglf',custnloglf,'frequency',N,'start',0) ;
end
CI=quantile(betaBoot,[0.025,0.975]);
yCI=quantile(Yboot',[0.025,0.975])';

beta=2^beta;
CI=2.^CI;

end