function [logL,model,bootY,model2]=householdTransmissionLogLikelihood(b,X,Y,N,resampleRep)
if nargin<5
    resampleRep=0;
end

% fit with lit-based dose-response for type 1
experiment{1} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',10^(-b(1)) ...
                        ,'secondaryContactAcquire',10^b(4) ...
                        ,'tertiaryContactAcquire',10^b(4) ...
                        ,'numDailySecondaryTertiaryContact',10^b(5) ...
                        ,'runNetwork',false,'serotype',1,'confidenceIntervalSamplerSeed',resampleRep);

% allow dose-response intercept difference for types 2 and 3
experiment{2} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',10^(-b(2)) ...
                        ,'secondaryContactAcquire',10^b(4) ...
                        ,'tertiaryContactAcquire',10^b(4) ...
                        ,'numDailySecondaryTertiaryContact',10^b(5) ...
                        ,'doseResponseBeta',10^b(6) ...
                        ,'runNetwork',false,'serotype',1,'confidenceIntervalSamplerSeed',resampleRep);

experiment{3} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',10^(-b(3)) ...
                        ,'secondaryContactAcquire',10^b(4) ...
                        ,'tertiaryContactAcquire',10^b(4) ...
                        ,'numDailySecondaryTertiaryContact',10^b(5) ...
                        ,'doseResponseBeta',10^b(7) ...
                        ,'runNetwork',false,'serotype',1,'confidenceIntervalSamplerSeed',resampleRep);
 

if any(X==0)
    p=zeros((1+length(experiment{1}.params.t))*9,1);
    idx=ismember(repmat([0,experiment{1}.params.t],1,9),unique(X));
else
    idx=ismember(repmat(experiment{1}.params.t,1,9),unique(X));
    p=zeros(length(experiment{1}.params.t)*9,1);
end

for s=1:3
    if any(X==0)
        tmp=[[0,experiment{s}.primary.prevalence],[0,experiment{s}.secondary.prevalence],[0,experiment{s}.tertiary.prevalence]]';
    else
        tmp=[experiment{s}.primary.prevalence,experiment{s}.secondary.prevalence,experiment{s}.tertiary.prevalence]';
    end
    p(1+(s-1)*length(tmp):s*length(tmp))=tmp;
end
model=p;

p2=p;
for s=1:3
    if any(X==0)
        tmp=[[0,cumsum(experiment{s}.primary.incidence)],[0,cumsum(experiment{s}.secondary.incidence)],[0,cumsum(experiment{s}.tertiary.incidence)]]';
    else
        tmp=[cumsum(experiment{s}.primary.incidence),cumsum(experiment{s}.secondary.incidence),cumsum(experiment{s}.tertiary.incidence)]';
    end
    p2(1+(s-1)*length(tmp):s*length(tmp))=tmp;
end

model2=p2;

p=p(idx);
p=p+10*eps;

logL=sum(Y.*log(p) + (N-Y).*log(1-p));

if nargout>2
    bootY=binornd(N,p);
end

end