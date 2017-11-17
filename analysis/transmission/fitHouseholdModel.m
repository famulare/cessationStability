%% fit transmission model

addpath(genpath('..\Houston1960'))
addpath(genpath('..\helperFunctions'))
addpath(genpath('..\doseResponse'))

[index,sibling,contact]=loadHouston1960();


%% simultaneous fit to all serotypes

timeAxis=repmat([0,7,14,21,28,35]',1,9);
X=timeAxis(:);

s=1;
Y=[0;index(s).SheddingTotal;0;sibling(s).SheddingUnder5;0;contact(s).SheddingUnder5];
N=[max(index(s).DenominatorTotal);index(s).DenominatorTotal;max(sibling(s).DenominatorUnder5);sibling(s).DenominatorUnder5;max(contact(s).DenominatorUnder5);contact(s).DenominatorUnder5];    

for s=2:3
    Y=[Y;0;index(s).SheddingTotal;0;sibling(s).SheddingUnder5;0;contact(s).SheddingUnder5];
    N=[N;max(index(s).DenominatorTotal);index(s).DenominatorTotal;max(sibling(s).DenominatorUnder5);sibling(s).DenominatorUnder5;max(contact(s).DenominatorUnder5);contact(s).DenominatorUnder5];    
end


custnlogl=@(b,Y,cens,N) -householdTransmissionLogLikelihood(b,X,Y,N);
[beta] = mle(Y,'nloglf',custnlogl,'frequency',N ...
            ,'start',[0.1,0.1,0.1,-5,0,0,0] ...
            ,'options',statset('maxiter',1e4,'maxfunevals',1e4)) ;

beta
householdTransmissionLogLikelihood(beta,X,Y,N)

% parametric bootstrap
reps=1000;
betaBoot=nan(length(beta),reps);
modelBoot=nan(101*9,reps);
for n=1:reps;
    n
    [~,~,bootY]=householdTransmissionLogLikelihood(beta,X,Y,N);
    custnlogl=@(b,Y,cens,N) -householdTransmissionLogLikelihood(b,X,Y,N,n);
    betaBoot(:,n) = mle(bootY,'nloglf',custnlogl,'frequency',N ...
            ,'start',beta ...
            ,'options',statset('maxiter',1e4,'maxfunevals',1e4));
    [~,modelBoot(:,n)]=householdTransmissionLogLikelihood(betaBoot(:,n),X,Y,N,n);     
end
CI=prctile(betaBoot',[2.5,97.5]);
YCI=prctile(modelBoot',[2.5,97.5]);

%%

figure(2);
tmax=35;

[~,fit]=householdTransmissionLogLikelihood(beta,0:tmax,1,1);
fit=reshape(fit,101,9);
fit=fit(:,[1,4,7,2,5,8,3,6,9]);
YCI=reshape(YCI,tmax+1,9,2);
YCI=YCI(:,[1,4,7,2,5,8,3,6,9],:);

%%
for s=1:9
    subplot(3,3,s); hold on;
    plotWithCI(0:tmax,fit(1:tmax+1,s),squeeze(YCI(1:tmax+1,s,:)),'k')
end


save('modelFit.mat','beta','betaBoot','CI','YCI')

%% Reff diffs

d{1}=primarySecondaryTertiaryDoseModel('serotype',1);
d{2}=primarySecondaryTertiaryDoseModel('serotype',2);
d{3}=primarySecondaryTertiaryDoseModel('serotype',3);
d{4}=primarySecondaryTertiaryDoseModel('serotype',1,'doseResponseBeta',3);

% Reff12=(sum(d{1}.tertiary.incidence)/sum(d{1}.secondary.incidence))/(sum(d{2}.tertiary.incidence)/sum(d{2}.secondary.incidence))
% Reff32=(sum(d{3}.tertiary.incidence)/sum(d{3}.secondary.incidence))/(sum(d{2}.tertiary.incidence)/sum(d{2}.secondary.incidence))

% can use prevalence instead of incidence since shedding is identical
% across serotypes
Reff12=(sum(d{1}.tertiary.prevalence)/sum(d{1}.secondary.prevalence))/(sum(d{2}.tertiary.prevalence)/sum(d{2}.secondary.prevalence))
Reff32=(sum(d{3}.tertiary.prevalence)/sum(d{3}.secondary.prevalence))/(sum(d{2}.tertiary.prevalence)/sum(d{2}.secondary.prevalence))

sig1=diff(squeeze((sum(YCI(:,7,:))./sum(YCI(:,4,:)))))/4;
sig2=diff(squeeze((sum(YCI(:,8,:))./sum(YCI(:,5,:)))))/4;
sig3=diff(squeeze((sum(YCI(:,9,:))./sum(YCI(:,6,:)))))/4;

Reff12CI(1)=Reff12*(1-2*sqrt(sig1^2/((sum(d{1}.tertiary.prevalence)/sum(d{1}.secondary.prevalence)))+sig2^2/((sum(d{2}.tertiary.prevalence)/sum(d{2}.secondary.prevalence)))));
Reff12CI(2)=Reff12*(1+2*sqrt(sig1^2/((sum(d{1}.tertiary.prevalence)/sum(d{1}.secondary.prevalence)))+sig2^2/((sum(d{2}.tertiary.prevalence)/sum(d{2}.secondary.prevalence)))))

Reff32CI(1)=Reff32*(1-2*sqrt(sig3^2/((sum(d{3}.tertiary.prevalence)/sum(d{3}.secondary.prevalence)))+sig2^2/((sum(d{2}.tertiary.prevalence)/sum(d{2}.secondary.prevalence)))));
Reff32CI(2)=Reff32*(1+2*sqrt(sig3^2/((sum(d{3}.tertiary.prevalence)/sum(d{3}.secondary.prevalence)))+sig2^2/((sum(d{2}.tertiary.prevalence)/sum(d{2}.secondary.prevalence)))))

