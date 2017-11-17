% fit WPV in India

data=xls2struct('..\..\Data\India2010\modeledPrevalenceInHealthyContactsOfWPV1Cases.xlsx');

%% fit cumulative incidence


N=data.effectiveSampleSize;
Y=round(data.estimatedPrevalence.*N);
NAb=data.assumedNAb;

custnlogl=@(b,Y,cens,N) -wpvIndiaTransmissionLogLikelihood(b,Y,N,NAb);
[beta] = mle(Y,'nloglf',custnlogl,'frequency',N ...
            ,'start',[-5] ...
            ,'options',statset('maxiter',1e4,'maxfunevals',1e4)) ;
beta

wpvIndiaTransmissionLogLikelihood(beta,Y,N,NAb)

%% parametric bootstrap
reps=1000;
betaBoot=nan(length(beta),reps);
modelBoot=nan(2,reps);
for n=1:reps;
    n
    [~,~,bootY]=wpvIndiaTransmissionLogLikelihood(beta,Y,N,NAb);
    custnlogl=@(b,Y,cens,N) -wpvIndiaTransmissionLogLikelihood(b,Y,N,NAb,n);
    betaBoot(:,n) = mle(bootY,'nloglf',custnlogl,'frequency',N ...
            ,'start',beta ...
            ,'options',statset('maxiter',1e4,'maxfunevals',1e4));
    [~,modelBoot(:,n)]=wpvIndiaTransmissionLogLikelihood(betaBoot(:,n),Y,N,NAb,n);     
end
CI=prctile(betaBoot(betaBoot<0)',[2.5,97.5])
YCI=prctile(modelBoot',[2.5,97.5]);

save('indiaFit.mat','beta','betaBoot','CI','YCI')