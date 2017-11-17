function [logL,model,bootY]=wpvIndiaTransmissionLogLikelihood(b,Y,N,NAb,resampleRep)
if nargin<5
    resampleRep=0;
end

t=1:140;
experiment{1} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'serotype',4 ...
                        ,'probDailyPrimarySecondaryContact',1,'t',t ...
                        ,'tertiaryContactAcquire',10^(b(1)+0.97) ...
                        ,'secondaryLog2NAb',log2(NAb(1)) ...
                        ,'secondaryContactAcquire',10^b(1) ...
                        ,'runNetwork',false,'confidenceIntervalSamplerSeed',resampleRep);

experiment{2} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'serotype',4 ...
                        ,'probDailyPrimarySecondaryContact',1,'t',t ...
                        ,'tertiaryContactAcquire',10^(b(1)+0.97) ...
                        ,'secondaryLog2NAb',log2(NAb(2)) ...
                        ,'secondaryContactAcquire',10^b(1) ...
                        ,'runNetwork',false,'confidenceIntervalSamplerSeed',resampleRep);

% weight time from paralysis
[~,pdf]=timeToParalysis(t);
                    
p=nan(length(Y),1);
for k=1:length(Y);
    tmp=conv(experiment{k}.secondary.prevalence,pdf);
    p(k)=mean(tmp(1:88)); % stools collected after onset (weighted by probability, for duration of observation)
end

p=p+10*eps;

model=p;

logL=nansum(Y.*log(p) + (N-Y).*log(1-p));

if nargout==3
    bootY=binornd(N,p);
end

end