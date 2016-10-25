% fit WPV in Louisiana

data=xls2struct('..\..\Data\Louisiana1957\householdSeroconversion.xlsx','seroconversion');

%% fit cumulative incidence


N=data.denominator([1,4,2,5]);
Y=data.seroconversions([1,4,2,5]);
NAb=data.medianPreExposureAntibodyTiter([1,4,2,5]);

custnlogl=@(b,Y,cens,N) -wpvLouisianaTransmissionLogLikelihood(b,Y,N,NAb);
[beta] = mle(Y,'nloglf',custnlogl,'frequency',N ...
            ,'start',[0,0] ...
            ,'options',statset('maxiter',1e4,'maxfunevals',1e4)) ;

10.^beta(1)
beta(2)

wpvLouisianaTransmissionLogLikelihood(beta,Y,N,NAb)

%% parametric bootstrap
reps=1000;
betaBoot=nan(length(beta),reps);
modelBoot=nan(4,reps);
for n=1:reps;
    n
    [~,~,bootY]=wpvLouisianaTransmissionLogLikelihood(beta,Y,N,NAb);
    custnlogl=@(b,Y,cens,N) -wpvLouisianaTransmissionLogLikelihood(b,Y,N,NAb,n);
    betaBoot(:,n) = mle(bootY,'nloglf',custnlogl,'frequency',N ...
            ,'start',beta ...
            ,'options',statset('maxiter',1e4,'maxfunevals',1e4));
    [~,modelBoot(:,n)]=wpvLouisianaTransmissionLogLikelihood(betaBoot(:,n),Y,N,NAb,n);     
end
CI=prctile(betaBoot',[2.5,97.5]);
YCI=prctile(modelBoot',[2.5,97.5]);

10.^CI