% analyze approximate waning data

addpath(genpath('..\helperFunctions'))
addpath(genpath('..\doseResponse'))
addpath(genpath('..\shedding'))

%% build immunity estimates by study

clear all; %close all; clc;


%% tOPVx3 Asturias, tOPVx3 Onorato, tOPVx3 Henry
% natural immunity Abbink2005, (Verlinde1959)
shedData=loadSheddingDurationData;

% add specific trial arms
author={'Henry','Onorato','Asturias','Abbink','Verlinde'};
immPlan={'tOPVx3','natural immunity','seronegative'};
idx=find(ismember(shedData.author,author)&ismember(shedData.primaryImmunization,immPlan) ...
    &(~(ismember(shedData.author,{'Verlinde'})&ismember(shedData.primaryImmunization,{'seronegative'}))));

immunizationPlans=shedData.primaryImmunization(idx);
authors=shedData.author(idx);
waningTime=shedData.challengeAgeMonths(idx)-shedData.ageLastRIMonths(idx);
serotype=shedData.challengeSerotype(idx);
sampleSize=shedData.sampleSize(idx);

NAb=[];
for k=1:length(idx);
    X=shedData.daysSinceChallenge';
    Y=shedData.sheddingCDF(idx(k),:)';
    N=repmat(shedData.sampleSize(idx(k)),length(Y),1);
    [NAb(k),NAbCI(k,:)]=sheddingDurationCDFbootstrapFitter(X,Y,N,1000);
end

waningTime(isnan(waningTime)&ismember(authors,{'Abbink'}))=45*12; % abbink elderly
waningTime(isnan(waningTime)&ismember(authors,{'Verlinde'}))=60; % Verlinde adults 

%% Jafari control

addpath(genpath('..\helperFunctions'));
addpath(genpath('..\..\data\doseResponse'));
addpath(genpath('..\doseResponse'));

doseData=doseResponseData();

waningTime(length(idx)+(1:3))=[2.5,4,60];
authors(length(idx)+(1:3))={'Jafari'};
serotype(length(idx)+(1:3))=1;
sampleSize(length(idx)+(1:3))=doseData.Jafari2014.tOPVxN.N;
count=0;
for k=length(idx)+(1:3)
    count=count+1;
    custnloglf = @(b,Y,cens,N) -sum(Y.*log(doseResponseModel(10^5.7,2^b)) + (N-Y).*log(1-doseResponseModel(10^5.7,2^b)));
    [beta,CI] = mle(doseData.Jafari2014.tOPVxN.pos(count),'nloglf',custnloglf,'frequency',doseData.Jafari2014.tOPVxN.N(count),'start',[7]);
    NAb(k)=round(2^beta);
    NAbCI(k,:)=round(2.^CI);
end


%% plot

figure(1); clf; 
subplot(1,2,1); hold on;

cmap=colormap(lines(3));

for k=1:length(NAb)
    r=10^(0.016*randn);
    if ~any(ismember(authors(k),{'Abbink','Verlinde'}))
        plot(waningTime(k)*r, NAb(k),'.','markersize',10,'markerfacecolor',cmap(serotype(k),:),'markeredgecolor',cmap(serotype(k),:));
    else
        if ~strcmp(immunizationPlans{k},'seronegative')
            plot(waningTime(k)*r, NAb(k),'d','markersize',4,'markerfacecolor',cmap(serotype(k),:),'markeredgecolor',cmap(serotype(k),:));
        else
            plot(waningTime(k)*r, NAb(k),'d','markersize',4,'markeredgecolor',cmap(serotype(k),:));
        end
    end
    plot(waningTime(k)*r*[1,1],NAbCI(k,:),'-','color',cmap(serotype(k),:));
end
set(gca,'yscale','log','xscale','log')
xlim([0.8,840]); ylim([1,10^4])
set(gca,'xtick',[1,6,12,60,600])
ylabel('OPV-equivalent antibody titer')
xlabel('months from last exposure to challenge')

text(120,10^3.3,'Type 1','color',cmap(1,:));
text(120,10^3,'Type 2','color',cmap(2,:));
text(120,10^2.7,'Type 3','color',cmap(3,:));
% colormap(cmap);
% h=colorbar;
% set(h,'ytick',(0.5+[0:length(authors)])/length(authors),'yticklabel',authors)

%% fit

% need proper error propagation
% need way to think about outcome variation
% show variations on included and excluded data. Max uncertainty is 1 order
% of magnitude, equivalent to one dose of tOPV...

t=linspace(1,1200,1200);

fitIdx=true(size(NAb));
% fitIdx=serotype~=3;

% log-linear best-fit
model=@(b,t) b(1)-b(2).*log2(t);
[beta,R,~,COVB]=nlinfit(waningTime(fitIdx),log2(NAb(fitIdx))',model,[11,0],'weights',sampleSize(fitIdx));
YPRED = nlpredci(model,t,beta,R,'covar',COVB);

% boostrap CI
reps=1000;
betaBoot=nan(length(beta),reps);
for k=1:length(NAb)
    tmpCI(:,k)=(linspace(log2(NAbCI(k,1)),log2(NAbCI(k,2)),reps));
end
tmpCI;

for k=1:reps
        % normal sample from CI
        ciBin=zeros(length(NAb),1);
        while any(ciBin<1 | ciBin>reps)
            ciBin=round(reps/2+reps/4*randn(1,length(NAb)));
        end
        for n=1:length(NAb)
            tmpNAb(n)=tmpCI(ciBin(n),n);
        end
        Yboot(:,k)=tmpNAb;
        [betaBoot(:,k),R,~,COVB]=nlinfit(waningTime(fitIdx),(tmpNAb(fitIdx))',model,[11,0],'weights',sampleSize(fitIdx));
        [pBoot(:,k),delta] = nlpredci(model,t,betaBoot(:,k),R,'covar',COVB);
        pBoot(:,k)=pBoot(:,k)+randn/2*delta';
        
end
CI=quantile(betaBoot',[0.025,0.975]);
yCI=quantile(pBoot',[0.025,0.975])';


plotWithCI(t,2.^YPRED',2.^yCI,'k')

%% bOPV

author={'Asturias'};
immPlan={'bOPVx3'};%,'seronegative'};
idx=find(ismember(shedData.author,author)&ismember(shedData.primaryImmunization,immPlan) ...
    &(~(ismember(shedData.author,{'Verlinde'})&ismember(shedData.primaryImmunization,{'seronegative'}))));

immunizationPlans=shedData.primaryImmunization(idx);
authors=shedData.author(idx);
waningTime=shedData.challengeAgeMonths(idx)-shedData.ageLastRIMonths(idx);
serotype=shedData.challengeSerotype(idx);
sampleSize=shedData.sampleSize(idx);

NAb=[];
for k=1:length(idx);
    X=shedData.daysSinceChallenge';
    Y=shedData.sheddingCDF(idx(k),:)';
    N=repmat(shedData.sampleSize(idx(k)),length(Y),1);
    [NAb(k),NAbCI(k,:)]=sheddingDurationCDFbootstrapFitter(X,Y,N,100);
end

for k=1:length(NAb)
    r=10^(0.015*randn);
    plot(waningTime(k)*r, NAb(k),'s','markersize',4,'markerfacecolor',cmap(serotype(k),:),'markeredgecolor',cmap(serotype(k),:));
    plot(waningTime(k)*r*[1,1],NAbCI(k,:),'-','color',cmap(serotype(k),:));
end
fitIdx=true(size(NAb));
% fitIdx=serotype~=3;

% log-linear best-fit
model=@(b,t) b(1)-beta(2).*log2(t);
[betaBOPV]=nlinfit(waningTime(fitIdx),log2(NAb(fitIdx))',model,[2],'weights',sampleSize(fitIdx));

plotWithCI(t,2.^(YPRED'+betaBOPV-YPRED(1)),2.^(yCI+betaBOPV-YPRED(1)),[0.1,0.7,0.1])


%% shedding index

tOPVNAbt=2.^(YPRED');
tOPVNAbtCI=2.^(yCI);

bOPVNAbt=2.^(YPRED'+betaBOPV-YPRED(1));
bOPVNAbtCI=2.^(yCI+betaBOPV-YPRED(1));

tOPVsheddingIndexT=nan(size(tOPVNAbt));
bOPVsheddingIndexT=nan(size(tOPVNAbt));
tOPVsheddingIndexTCI=nan(size(tOPVNAbt,1),2);
bOPVsheddingIndexTCI=nan(size(tOPVNAbt,1),2);

for k=1:length(tOPVsheddingIndexT)
    tOPVsheddingIndexT(k)=doseResponseModel(10^5.7,tOPVNAbt(k))*(sum(sheddingViralLoadGMT(tOPVNAbt(k)).*sheddingDurationCDF(tOPVNAbt(k))));
    bOPVsheddingIndexT(k)=doseResponseModel(10^5.7,bOPVNAbt(k))*(sum(sheddingViralLoadGMT(bOPVNAbt(k)).*sheddingDurationCDF(bOPVNAbt(k))));
    tOPVsheddingIndexTCI(k,1)=doseResponseModel(10^5.7,tOPVNAbtCI(k,1))*(sum(sheddingViralLoadGMT(tOPVNAbtCI(k,1)).*sheddingDurationCDF(tOPVNAbtCI(k,1))));
    bOPVsheddingIndexTCI(k,1)=doseResponseModel(10^5.7,bOPVNAbtCI(k,1))*(sum(sheddingViralLoadGMT(bOPVNAbtCI(k,1)).*sheddingDurationCDF(bOPVNAbtCI(k,1))));
    tOPVsheddingIndexTCI(k,2)=doseResponseModel(10^5.7,tOPVNAbtCI(k,2))*(sum(sheddingViralLoadGMT(tOPVNAbtCI(k,2)).*sheddingDurationCDF(tOPVNAbtCI(k,2))));
    bOPVsheddingIndexTCI(k,2)=doseResponseModel(10^5.7,bOPVNAbtCI(k,2))*(sum(sheddingViralLoadGMT(bOPVNAbtCI(k,2)).*sheddingDurationCDF(bOPVNAbtCI(k,2))));
end
tOPVsheddingIndexT=tOPVsheddingIndexT./(doseResponseModel(10^5.7,1)*(sum(sheddingViralLoadGMT(1).*sheddingDurationCDF(1))));
bOPVsheddingIndexT=bOPVsheddingIndexT./(doseResponseModel(10^5.7,1)*(sum(sheddingViralLoadGMT(1).*sheddingDurationCDF(1))));
tOPVsheddingIndexTCI=tOPVsheddingIndexTCI./(doseResponseModel(10^5.7,1)*(sum(sheddingViralLoadGMT(1).*sheddingDurationCDF(1))));
bOPVsheddingIndexTCI=bOPVsheddingIndexTCI./(doseResponseModel(10^5.7,1)*(sum(sheddingViralLoadGMT(1).*sheddingDurationCDF(1))));

subplot(1,2,2);  hold on;
plotWithCI(t,tOPVsheddingIndexT,tOPVsheddingIndexTCI,'k');
plotWithCI(t,bOPVsheddingIndexT,bOPVsheddingIndexTCI,[0.1,0.7,0.1]);
set(gca,'yscale','log','xscale','log')
xlim([0.8,840]); ylim([1e-4,1])
set(gca,'xtick',[1,6,12,60,600])
ylabel('mean relative shedding index')
xlabel('months from last exposure to challenge')

save('waningWorkspace.mat')