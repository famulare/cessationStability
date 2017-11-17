function [experiment]=primarySecondaryTertiaryDoseModel(varargin)

params=struct('t',1:100 ... % duration of trial
            ,'serotype',2 ... % challenge serotype
            ,'vaccinationTime',1 ... % start day of trial
            ,'vaccineDose',10^6 ... % vaccine dose TCID50
            ,'perDoseEfficacy',0.8863... % efficacy scaler to dose response to account for trial arm variations in vaccine take
            ,'secondaryContactAcquire',10^(-5.2943) ... % daily fecal dose (g/day) of secondary from primary contact
            ,'tertiaryContactAcquire',10^(-5.2943) ... % daily fecal dose (g/day) of tertiary from secondary contact
            ,'primaryLog2NAb',0 ... % primary OPV-equivalent immunity
            ,'secondaryLog2NAb',0 ... % secondary OPV-equivalent immunity
            ,'tertiaryLog2NAb',0 ... % tertiary OPV-equivalent immunity
            ,'numDailyPrimarySecondaryContact',1 ... % modify daily contact rates to lower per-day probability
            ,'numDailySecondaryTertiaryContact',8.9685 ...  % uniform probability for how likely a contact is to be visited each day
            ,'primaryAgeMos',12 ... % age of primary in months
            ,'secondaryAgeMos',48 ... % age of tertiary in months
            ,'tertiaryAgeMos',48 ... % age of tertiary in months
            ,'maxNetworkSize',40 ... % maximum size of social network
            ,'runNetwork',true ... % run intimate contacts network? turn off for fitting
            ,'confidenceIntervalSamplerSeed',0 ... % >0: resample parameters from confidence interval for uncertainty propagation?
            ,'doseResponseBeta',[] ... % expose dose response beta param for fitting
            );

params=keyValuePairVararginHandler(params,varargin);

% private parameters internal to model
    if params.serotype~=4
        params.private.probabilityStillInfected.b=[log(30.3),log(1.164),log(1.86)];
        params.private.probabilityStillInfected.bCI=log([23.6,1.13,1.57;38.6,1.21,2.27]);
    else
        params.private.probabilityStillInfected.b=[log(43),log(1.164),log(1.69)];
        params.private.probabilityStillInfected.bCI=log([35.7,1.13,1.53;51.7,1.21,1.94]);
    end
    params.private.sheddingViralLoadGMT.b=[1.64,0.18,0.056,0.32];
    params.private.sheddingViralLoadGMT.bCI=[1.26,0.01,0.01,0.08;2.09,0.78,0.079,0.71];
    params.private.peakSheddingAgeMultiplier.b=[6.67,4.29,9.92];
    params.private.peakSheddingAgeMultiplier.bCI=[5.9,3.5,1;7.5,5.0,33];
    params.private.doseResponseProbability.bCIrho=[0,0.87,0];
    if params.serotype==1
        params.private.doseResponseProbability.b=[14.24,0.444,0.4624];
        params.private.doseResponseProbability.bCI=[3,0.29,0.42;59,0.83,0.50];
    elseif params.serotype==2 
        params.private.doseResponseProbability.b=[7.9611,0.444,0.4624];
        params.private.doseResponseProbability.bCI=[2,0.29,0.42;29.9,0.83,0.50];
    elseif params.serotype==3 
        params.private.doseResponseProbability.b=[17.8232,0.444,0.4624];
        params.private.doseResponseProbability.bCI=[4.66,0.29,0.42;62.79,0.83,0.50];
    elseif params.serotype==4 % WPV 
        params.private.doseResponseProbability.b=[2.31,0.444,0.4624];    
        params.private.doseResponseProbability.bCI=[0.25,0.29,0.42;37,0.83,0.50];
    end
    % resample parameters for uncertainty propagation
    if params.confidenceIntervalSamplerSeed>0;
        rng(params.confidenceIntervalSamplerSeed);
        
        params.private.probabilityStillInfected.b=resampleParametersFromCI(params.private.probabilityStillInfected.bCI);
        params.private.sheddingViralLoadGMT.b=resampleParametersFromCI(params.private.sheddingViralLoadGMT.bCI);
        params.private.peakSheddingAgeMultiplier.b=resampleParametersFromCI(params.private.peakSheddingAgeMultiplier.bCI);
        
        tmp=params.private.doseResponseProbability.bCI;
        tmp(:,1)=log10(tmp(:,1)); % first parameter lives in log space
        params.private.doseResponseProbability.b=resampleParametersFromCI(tmp, params.private.doseResponseProbability.bCIrho);
        params.private.doseResponseProbability.b(1)=10^params.private.doseResponseProbability.b(1);
    end
    % if fitting dose response beta
    if ~isempty(params.doseResponseBeta)
        params.private.doseResponseProbability.b(1)=params.doseResponseBeta;
    end
    
% household transmission model
    % setup
    t = params.t;

    % primary vaccine recipients
        primaryLog2NAb=params.primaryLog2NAb; 

        % primary incidence
        probPrimaryIncidence =  [params.perDoseEfficacy*doseResponseProbability(params.vaccineDose, params.primaryLog2NAb,params),zeros(1,length(t)-1)]; 
        % primary prevalence
        probPrimaryShedding=prevalenceFromIncidence(probPrimaryIncidence,t,primaryLog2NAb,params);

    % secondary contact   
        secondaryLog2NAb=params.secondaryLog2NAb;
        secondaryContactNumPerDay=params.numDailyPrimarySecondaryContact;

        % incidence
        probSecondaryIncidence=incidenceFromContact(probPrimaryIncidence,t,primaryLog2NAb,params.primaryAgeMos,params.secondaryContactAcquire,secondaryLog2NAb,secondaryContactNumPerDay,params);
        % prevalence
        probSecondaryShedding=prevalenceFromIncidence(probSecondaryIncidence,t,secondaryLog2NAb,params);

    % tertiary contact of secondary contact   
        tertiaryLog2NAb=params.tertiaryLog2NAb;
        tertiaryContactNumPerDay=params.numDailySecondaryTertiaryContact;

        % incidence
        probTertiaryIncidence=incidenceFromContact(probSecondaryIncidence,t,secondaryLog2NAb,params.secondaryAgeMos,params.tertiaryContactAcquire,tertiaryLog2NAb,tertiaryContactNumPerDay,params);
        % prevalence
        probTertiaryShedding=prevalenceFromIncidence(probTertiaryIncidence,t,tertiaryLog2NAb,params);

    %output
    experiment.params=params;    
    experiment.primary.incidence = probPrimaryIncidence;
    experiment.primary.prevalence = probPrimaryShedding;

    experiment.secondary.incidence = probSecondaryIncidence;
    experiment.secondary.prevalence = probSecondaryShedding;

    experiment.tertiary.incidence = probTertiaryIncidence;
    experiment.tertiary.prevalence = probTertiaryShedding;

    
% intimate contacts network model
    if params.runNetwork
        % Houston 1960, average family size = 3 children (2 under 5 years),
        % renormalize weights from assuming direct link to accounting for
        % possible indirect links
%         [pss,pis] = adjustContactWeight(sum(experiment.secondary.incidence),sum(experiment.tertiary.incidence));
        pis=sum(experiment.secondary.incidence)/sum(experiment.primary.incidence);
        pss=sum(experiment.tertiary.incidence)/pis;
        
        % expected number infected in family vs family size
        for k=1:params.maxNetworkSize, 
            x(k)=expectedPositiveGivenIndex(k,pis,pss); 
        end
        experiment.household.meanIncidence=x;
        experiment.household.size=2:41;

        % infected number infected in sibling intimate social network 
        % vs network size
        for k=1:params.maxNetworkSize, 
            x(k)=expectedPositiveGivenIndex(k,pss,pss); 
        end
        experiment.siblingFriendNetwork.meanIncidence=x;
        experiment.siblingFriendNetwork.size=2:41;
    end

end

function receivingIncidence=incidenceFromContact(inputIncidence,t,inputLog2NAb,inputAge,contactAcquire,receivingLog2NAb,ContactNumPerDay,params)
    
    inputInfectiousDuration=probabilityStillInfected(t,inputLog2NAb,1,params);
    receivingIncidence=zeros(size(inputIncidence));
    % convolve shedding duration with start probability per day
    for m=t
        % prevalence from secondary infections each day
            stillInfectedHolder=[zeros(1,m-1),inputInfectiousDuration(:,1:length(t)-m+1)];
            inputPrevalence=inputIncidence(m).*stillInfectedHolder;
        % infectious probability from those secondaries
            receivingIncidencePerDailyContact=inputPrevalence.*[zeros(1,m-1),pTransmissionGivenContactShedding(inputLog2NAb,inputAge,contactAcquire,m:t(end),m,receivingLog2NAb,ContactNumPerDay,params)];
        % add up total incidence
            receivingIncidence = receivingIncidence + receivingIncidencePerDailyContact;
    end
end

function prevalence=prevalenceFromIncidence(incidence,t,log2NAb,params)
    prevalence=zeros(size(incidence));
    infectiousDuration=probabilityStillInfected(t,log2NAb,1,params);
    % convolve shedding duration with start probability per day
    for m=t
        stillInfectedHolder=[zeros(1,m-1),infectiousDuration(:,1:length(t)-m+1)];
        prevalence = prevalence + incidence(m).*stillInfectedHolder;
    end
end


function PStillInfected = probabilityStillInfected(t,log2NAb,t0,params)
% shedding duration reverse CDF
    
    b=params.private.probabilityStillInfected.b;
    
    PStillInfected = 1/2*(1-erf(1/sqrt(2*b(3)^2)*(log(abs(t-t0))-(b(1)-b(2)*log2NAb)))).*double(t>=t0);
    PStillInfected(t<(t0+1))=0;
end

function [PTrans] = pTransmissionGivenContactShedding(inputLog2NAb,inputAge,contactAcquire,t,t0,receivingLog2NAb,ContactNumPerDay,params)
% inhomogeneous poisson transmission based on fecal viral load per day

    d=doseResponseProbability(contactAcquire*sheddingViralLoadGMT(inputLog2NAb,t,t0,inputAge,params),receivingLog2NAb,params);
    beta=1-(1-d).^ContactNumPerDay; % sum over number of attempts
    
    compBeta=[1,1-beta(1:end-1)];
    PTrans=beta.*cumprod(compBeta);
    
    PTrans(t<(t0+1))=0;
end

function probabilityInfected = doseResponseProbability(dose,log2NAb,params)
    b=params.private.doseResponseProbability.b;
    probabilityInfected = 1-(1+dose/b(1)).^(-b(2)./(2.^max(0,log2NAb)).^b(3)); 
end

function shedTiter=sheddingViralLoadGMT(log2NAb,t,t0,ageMonths,params)
% shed titer given infected

b=params.private.sheddingViralLoadGMT.b;
 
peakShed = 10.^(peakSheddingAgeMultiplier(ageMonths,params)*(1-b(3)*max(0,log2NAb)));

% lognormal temporal profile in TCID50. Not a probability distrubition
shedTiter= exp(b(1)-b(2)^2/2)./(t-t0+0.01).*exp(-(log(t-t0+0.01)-b(1)).^2./(2*(b(2)+b(4)*log(t-t0)).^2));
shedTiter=peakShed*shedTiter;

% floor
% appears to be a floor. This is an assay artifact in that shedding less that this is a negative, which is accounted for in the sheddingDurationCDF.
shedTiter=max(shedTiter,10^(2.6)); 
shedTiter(t<(t0+1))=0;

end

function peakSheddingAgeFactor = peakSheddingAgeMultiplier(ageMonths,params)
% peak shedding depends on age
b=params.private.peakSheddingAgeMultiplier.b;

peakSheddingAgeFactor =min(b(1),((b(1)-b(2))*exp(-(ageMonths-7)/b(3))+b(2)));

end

function b=resampleParametersFromCI(bCI,rho)
% CI sampler
if nargin==1, 
    rho=zeros(1,size(bCI,2));
end
nBins=100;
b=zeros(1,size(bCI,2));
tmp=zeros(nBins,size(bCI,2));
for n=1:size(bCI,2)
    tmp(:,n)=linspace(bCI(1,n),bCI(2,n),nBins);
end

idx=b;
while any(idx<1 | idx>nBins)
    idx(1)=round(nBins/2+nBins/4*randn);
    for n=2:length(b)
        idx(n)=round(nBins/2+rho(n)*(idx(1)-nBins/2)+sqrt((1-rho(n)^2))*nBins/4*randn);
    end
end
for n=1:length(b)
    b(n)=tmp(idx(n),n);
end

end

function expectedPositive = expectedPositiveGivenIndex(N,pis, pss)
j=1:(N-1); 
expectedPositive= ...
    N* ... % number in network times probability infected per person
     pis*(1 ... % pair interaction: infant to target
     + sum(binom(N-1,j).*pis.^(j-1).*(1-pis).^(N-j).*(1-(1-pss).^j)) ...  % triangle interaction: infant to any sibling, sibling to target
     + sum(binom(N-2,j).*pis.^(j-1).*(1-pis).^(N-1-j).*binom(N-3,j).*pss.^j.*(1-pss).^(N-2-j).*(1-(1-pss).^j) ) )...  % quad interaction: infant to sibling, sibling to sibling, sibling to target
     ; % quad interaction is small compared to triangle for relevant range in this study, so truncated here.
end

function b = binom(N,k)
    b=zeros(size(k));
    idx=k<=N;
    b(idx)=arrayfun(@(x) nchoosek(N,x),k(idx));
end
