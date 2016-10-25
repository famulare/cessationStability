function [NAbBeta,CI,yCI,L]=antibodyDependenceDoseResponseBootstrapFitter(NAb,Y,N,naiveBeta,betaCI,reps,modelNAbRange)
if nargin<6
    reps=100;
end
if nargin<7
    modelNAbRange=logspace(0,log10(2^14),1000);
end

tmpCI(:,1)=linspace(betaCI(1,1),betaCI(2,1),reps);
tmpCI(:,2)=linspace(betaCI(1,2),betaCI(2,2),reps);
betaCI=tmpCI;

dose=10^5.7;

idx=any(~isnan(Y),2);
NAb=NAb(idx);
Y=Y(idx);
N=N(idx);

custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(doseResponseModel(dose,NAb,1,[naiveBeta(1),naiveBeta(2),2^b(1)])) + N.*(1-Y).*log(1-doseResponseModel(dose,NAb,1,[naiveBeta(1),naiveBeta(2),2^b(1)])));
[NAbBeta,CI] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',[0]) ;
L=-custnloglf(NAbBeta,Y,0,N);



% parametric bootstap resampling around naive params
betaBoot=nan(length(NAbBeta),reps);
[X,~,idx]=unique(NAb);
totN=zeros(size(X));
for k=1:length(X)
    totN(k)=sum(N(idx==k));
end
Yboot=zeros(length(X),reps);
pBoot=nan(length(modelNAbRange),reps);

for k=1:reps
        % normal sample from CI
        ciBin=[0,0];
        while any(ciBin<1 | ciBin>reps)
            ciBin=round(reps/2+reps/4*randn(1,2));
        end
        tmpBeta(1)=betaCI(ciBin(1),1);
        tmpBeta(2)=betaCI(ciBin(2),2);
        Yboot(:,k)=binornd(totN,doseResponseModel(dose,X,1,[tmpBeta(1),tmpBeta(2),2^NAbBeta]))./totN;
        custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(doseResponseModel(dose,X,1,[tmpBeta(1),tmpBeta(2),2^b])) + N.*(1-Y).*log(1-doseResponseModel(dose,X,1,[tmpBeta(1),tmpBeta(2),2^b])));
    try % inference rarely chokes
            [betaBoot(:,k)] = mle(Yboot(:,k),'nloglf',custnloglf,'frequency',totN,'start',[0]) ;
            pBoot(:,k)=doseResponseModel(dose,modelNAbRange,1,[tmpBeta(1),tmpBeta(2),2^betaBoot(k)]);
    catch
        k=k-1;
    end
end
CI=quantile(betaBoot',[0.025,0.975]);
yCI=quantile(pBoot',[0.025,0.975])';

NAbBeta=2.^NAbBeta;
CI=2.^CI;
%}
end