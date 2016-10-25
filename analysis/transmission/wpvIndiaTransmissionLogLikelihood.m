function [logL,model,bootY]=wpvIndiaTransmissionLogLikelihood(b,Y,N,NAb,resampleRep)
if nargin<5
    resampleRep=0;
end


experiment{1} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'serotype',4 ...
                        ,'probDailyPrimarySecondaryContact',1,'t',1:70 ...
                        ,'tertiaryContactAcquire',10^(b(1)+0.97) ...
                        ,'secondaryLog2NAb',log2(NAb(1)) ...
                        ,'secondaryContactAcquire',10^b(1) ...
                        ,'runNetwork',false,'confidenceIntervalSamplerSeed',resampleRep);

experiment{2} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'serotype',4 ...
                        ,'probDailyPrimarySecondaryContact',1,'t',1:70 ...
                        ,'tertiaryContactAcquire',10^(b(1)+0.97) ...
                        ,'secondaryLog2NAb',log2(NAb(2)) ...
                        ,'secondaryContactAcquire',10^b(1) ...
                        ,'runNetwork',false,'confidenceIntervalSamplerSeed',resampleRep);

% weight time from paralysis
% Casey 1942 (Figure 2)
% The incubation period in epidemic poliomyelitis
% JAMA
d=[1,3.5,9.5,15.5,21.5,27.5,33.5,70];
n=cumsum([0,0,4,17,6,1,1,0])/29;
cdf=interp1(d,n,1:70);
                    
p=nan(length(Y),1);
for k=1:length(Y);
    p(k)=mean(cdf.*experiment{k}.secondary.prevalence); % stools collected after onset (weighted by probability)
end

p=p+10*eps;

model=p;

logL=nansum(Y.*log(p) + (N-Y).*log(1-p));

if nargout==3
    bootY=binornd(N,p);
end

end