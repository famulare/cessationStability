% shedding index

close all; clear all; clc;
addpath(genpath('..\helperFunctions\'))
addpath(genpath('..\doseResponse\'))
addpath(genpath('..\shedding\'))

%% 
NAbRange=(logspace(0,11*log10(2),200));
serotype=2;
age=18;
t=0:80;
t0=0;
reps=10000;

        
doseResponseCI=[logspace(log10(2),log10(30),reps);linspace(0.29,0.83,reps);linspace(0.51,0.57,reps)]';
sheddingViralLoadCI=[linspace(1.26,2.09,reps);linspace(0.01,0.78,reps);linspace(0.01,0.079,reps);linspace(0.08,0.71,reps)]';
sheddingDurationCI=log([linspace(23.6,38.8,reps);linspace(1.13,1.21,reps);linspace(1.57,2.27,reps)])';
for k=1:length(NAbRange)
    sheddingIndex(k)=doseResponseModel(10^5.7,NAbRange(k),serotype)*(sum(sheddingViralLoadGMT(NAbRange(k),age,t,t0).*sheddingDurationCDF(NAbRange(k),t,t0)));
    
    for n=1:reps
        ciIdx=zeros(4);
        while any(any(ciIdx<1 | ciIdx>reps))
            ciIdx=round(reps/2+reps/4*randn(4));
%             % strong anti-correlation between load 2 and load 4
%             ciIdx(2,4)=round(-0.58*(0.063/0.078)*(ciIdx(2,2)-reps/2)+sqrt(1-.58)*reps/4*randn(1));
        end
        for m=1:3
            betaDose(m)=doseResponseCI(ciIdx(1,m),m);
            betaLoad(m)=sheddingViralLoadCI(ciIdx(2,m),m);
            betaDur(m)=sheddingDurationCI(ciIdx(3,m),m);
        end
        betaLoad(4)=sheddingViralLoadCI(ciIdx(2,4),4);
        sheddingIndexBoot(n)=doseResponseModel(10^5.7,NAbRange(k),serotype,betaDose)*(sum(sheddingViralLoadGMT(NAbRange(k),age,t,t0,betaLoad).*sheddingDurationCDF(NAbRange(k),t,t0,betaDur)));
    end
    sheddingIndexCI(k,:)=prctile(sheddingIndexBoot,[2.5,97.5]);
end

%%
figure(2); clf;
plotWithCI(NAbRange,sheddingIndex,sheddingIndexCI,'k')
set(gca,'xscale','log','yscale','log')
xlim([0.8, 2100]); ylim([1e2,10^6.2])
set(gca,'xtick',[1,8,2^5,2^8,2^11])
xlabel('OPV-equivalent humoral antibody titer')
ylabel('mean shedding index [TCID50/g]')

