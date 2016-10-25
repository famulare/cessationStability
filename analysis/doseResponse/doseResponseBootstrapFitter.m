function [beta,CI,yCI,L]=doseResponseBootstrapFitter(X,Y,N,reps,modelDoseRange,NAb,bpAlpha,bpAlphaCI)
if nargin<4
    reps=100;
end
if nargin<5
    modelDoseRange=logspace(-1,6,1e4);
end
if nargin<6
    NAb=1;
end
if nargin<7
    bpAlpha=[];  % takes alpha parameter as given for beta-poisson
    bpAlphaCI=[];
end
if nargin<8 && ~isempty(bpAlpha)
    bpAlphaCI=bpAlpha*ones(reps,1);
end
if nargin==8 && length(bpAlphaCI)>1
    bpAlphaCI=linspace(bpAlphaCI(1),bpAlphaCI(2),reps);
end    

idx=any(~isnan(Y),2);
X=X(idx);
Y=Y(idx);
N=N(idx);

if isempty(bpAlpha)
    custnloglf = @(b,Y,cens,N) -sum(Y.*log(doseResponseModel(X,NAb,1,2.^[b(1),b(2),0])) + (N-Y).*log(1-doseResponseModel(X,NAb,1,2.^[b(1),b(2),0])));
    [beta,CI] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',[7,1]) ;
else
    custnloglf = @(b,Y,cens,N) -sum(Y.*log(doseResponseModel(X,NAb,1,[2.^b(1),bpAlpha,1])) + (N-Y).*log(1-doseResponseModel(X,NAb,1,[2.^b(1),bpAlpha,1])));
    [beta,CI] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',[7]) ;
end
L=-custnloglf(beta,Y,0,N);

% parametric bootstap with resampling around tested doses
betaBoot=nan(length(beta),reps);
[X,~,idx]=unique(X);
totN=zeros(size(X));
for k=1:length(X)
    totN(k)=sum(N(idx==k));
end
Yboot=zeros(length(X),reps);
pBoot=nan(length(modelDoseRange),reps);

for k=1:reps
    if isempty(bpAlpha)
        Yboot(:,k)=binornd(totN,doseResponseModel(X,NAb,1,2.^[beta(1),beta(2),0]));
        custnloglf = @(b,Y,cens,N) -sum(Y.*log(doseResponseModel(X,NAb,1,2.^[b(1),b(2),0])) + (N-Y).*log(1-doseResponseModel(X,NAb,1,2.^[b(1),b(2),0])));
    else
        % normal sample from CI
        ciBin=0;
        while ciBin<1 || ciBin>reps
            ciBin=round(reps/2+reps/4*randn);
        end
        tmpAlpha=bpAlphaCI(ciBin); 
        Yboot(:,k)=binornd(totN,doseResponseModel(X,NAb,1,[2.^beta(1),tmpAlpha,1]));
        custnloglf = @(b,Y,cens,N) -sum(Y.*log(doseResponseModel(X,NAb,1,[2.^b(1),tmpAlpha,1])) + (N-Y).*log(1-doseResponseModel(X,NAb,1,[2.^b(1),tmpAlpha,1])));
    end
    try % inference rarely chokes
        if isempty(bpAlpha)
            [betaBoot(:,k)] = mle(Yboot(:,k),'nloglf',custnloglf,'frequency',totN,'start',[7,1]) ;
            pBoot(:,k)=doseResponseModel(modelDoseRange,NAb,1,2.^[betaBoot(1,k),betaBoot(2,k),0]);
        else
            [betaBoot(:,k)] = mle(Yboot(:,k),'nloglf',custnloglf,'frequency',totN,'start',[7]) ;
            pBoot(:,k)=doseResponseModel(modelDoseRange,NAb,1,[2.^betaBoot(1,k),tmpAlpha,1]);
        end
    catch
        k=k-1;
    end
end
CI=quantile(betaBoot',[0.025,0.975]);
yCI=quantile(pBoot',[0.025,0.975])';

beta=2.^beta;
CI=2.^CI;

% if ~isempty(bpAlpha)
%     beta=[beta(1),bpAlpha];
%     CI=[CI',bpAlphaCI([1,end])'];
% end
end