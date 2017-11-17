% explore shedding viral load (concentration in stool) data

close all; clear all;
data=loadSheddingViralLoadData;

%% plot category and serotype
immunizationPlans=unique(data.primaryImmunization);


figure(2); clf;
cmap=colormap(lines(3));
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:3
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeSerotype==n;
        if any(idx)
            plotSpanNaN(data.daysSinceChallenge,10.^data.sheddingGMT(idx,:),'-','color',cmap(n,:))
        end
    end
    ylim(10.^[2.6 6.9])
    xlim([0 30])
    set(gca,'yscale','log','ytick',10.^[3:1:7])
    title(immunizationPlans(k),'fontweight','normal')    
    if k==1,
        text(2,10^5.8,'1','color',cmap(1,:));
        text(2,10^5.5,'2','color',cmap(2,:));
        text(2,10^5.2,'3','color',cmap(3,:));
    end
end

%% plot average by category and serotype
cmap=colormap(lines(3));

figure(3); clf;
for k=1:length(data.immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:3
        idx=find(strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeSerotype==n);
        if any(idx)
            plotSpanNaN(data.daysSinceChallenge,data.meanSheddingGMT(3*(k-1)+n,:),'color',cmap(n,:))
        end
    end
    if k==1,
        text(2,5.8,'1','color',cmap(1,:));
        text(2,5.5,'2','color',cmap(2,:));
        text(2,5.2,'3','color',cmap(3,:));
    end
    ylim([2.6 6.9])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end

%% plot category and age at challenge
immunizationPlans=unique(data.primaryImmunization);
ages=unique(data.challengeAgeMonths);
ageBins=0:36;
cmap=cbrewer('div','RdYlBu',length(ageBins));
for k=1:length(ageBins)
    cmap(k,:)=cmap(k,:)*(0.3*(1-2*k/length(ageBins)).^2+1-0.3);
end
ageIdx=ismember(ageBins,ages);
ageIdx([26,36])=true;
cmap=cmap(ageIdx,:);
cmap(end+1,:)=[0,0,0];

figure(5); clf; colormap(cmap);
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:length(ages)
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeAgeMonths==ages(n);
        if any(idx)
            plotSpanNaN(data.daysSinceChallenge,data.sheddingGMT(idx,:),'color',cmap(n,:))
        end
    end
    ylim([2.6 6.9])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end
h=colorbar;
set(h,'ytick',(0.5+[0:length(ages)])/length(ages),'yticklabel',ages)


%% age exploration: plot mucosal naives grouped by age at challenge
immunizationPlans={'seronegative','unvaccinated','IPVx3'};
ages=unique(data.challengeAgeMonths);
ageBins=0:36;
cmap=cbrewer('div','RdYlBu',length(ageBins));
for k=1:length(ageBins)
    cmap(k,:)=cmap(k,:)*(0.3*(1-2*k/length(ageBins)).^2+1-0.3);
end
ageIdx=ismember(ageBins,ages);
ageIdx([26,36])=true;
cmap=cmap(ageIdx,:);
cmap(end+1,:)=[0,0,0];

figure(5); clf; colormap(cmap);
week1GMT=[];
week1Age=[];
week1Samp=[];
subplot(1,2,1); hold on
for k=1:length(immunizationPlans)
    for n=1:length(ages)
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeAgeMonths==ages(n);
        if any(idx )
            plotSpanNaN(data.daysSinceChallenge,10.^data.sheddingGMT(idx,:),'color',cmap(n,:))
            plot(data.daysSinceChallenge(4),10.^data.sheddingGMT(idx,4),'o','markeredgecolor',cmap(n,:),'markerfacecolor',cmap(n,:));
            
            week1GMT(end+(1:sum(idx)))=data.sheddingGMT(idx,4);
            week1Age(end+(1:sum(idx)))=ages(n);
            week1Samp(end+(1:sum(idx)))=data.sampleSize(idx);
            data.author(idx)
        end
    end
    ylim(10.^[2.6 7.5])
    xlim([0 60])
%     title('no live virus exposure','fontweight','normal')
end
xlabel('days')
ylabel('fecal viral load (TCID50/g)')
set(gca,'yscale','log')

% plot age vs peak titer
[week1Age,idx]=sort(week1Age);
week1GMT=week1GMT(idx);
week1Samp=week1Samp(idx);
subplot(1,2,2); hold on
for n=1:length(ages)
    idx=ages(n)==week1Age;
    if any(idx)
        plot(week1Age(idx)/12,10.^week1GMT(idx),'o','markeredgecolor',cmap(n,:),'markerfacecolor',cmap(n,:));
    end
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('age at challenge (years)')
ylabel('peak viral load (TCID50/g)')
xlim([0,max(week1Age)/12+1])

h=colorbar;
set(h,'ytick',(0.5+[0:length(ages)])/length(ages),'yticklabel',ages)

% fit age-dependent titer model
t=min(week1Age):max(week1Age);
model=@(b,t) ((b(1)-b(2))*exp(-(t-min(week1Age))/b(3)^2)+b(2));
[beta,R,~,COVB]=nlinfit(week1Age,week1GMT,model,[6.7,4.5,1],'weights',week1Samp);
parCI=nlparci(beta,R,'covar',COVB);

AbRange=(exp(linspace(log(t(1)),log(t(end)),20)));
[YPRED, DELTA] = nlpredci(model,AbRange,beta,R,'covar',COVB);

CI=[10.^(YPRED-DELTA);10.^(YPRED+DELTA)]';
% encoded in peakSheddingAgeMultiplier
plotWithCI(AbRange/12,10.^YPRED,CI,[0.4,0.8,0.4],'-',10,false,true)
ylim(10.^[2.6 7.5])

% convert back to proper scale and display
beta(3)=beta(3)^2
parCI(3,:)=parCI(3,:).^2

%% fit temporal profile

ageAdjust=120; %months of age adjustment: saturate age-decline
% find mucosally-naive studies
idx=find(ismember(data.primaryImmunization,{'seronegative','unvaccinated','IPVx2','IPVx3'}));
idx=idx(data.publicationDate(idx)~=2005);% excluding elderly
X=repmat(data.daysSinceChallenge',1,length(idx));
Y=data.sheddingGMT(idx,:)';
N=repmat(data.sampleSize(idx),1,length(Y))';
for k=1:length(idx)
    Y(:,k)=Y(:,k)*peakSheddingAgeMultiplier(ageAdjust)./peakSheddingAgeMultiplier(data.challengeAgeMonths(idx(k)));
end

figure(1); clf
hold on; 
% plot approximate mean shedding duration CDF
meanNaive=nansum(Y.*N,2)./nansum(N.*~isnan(Y),2);
totSamp=nansum(N.*~isnan(Y),2);
meanNaive(meanNaive==0)=nan;
totSamp(isnan(meanNaive))=nan;

plotSpanNaN(data.daysSinceChallenge,10.^meanNaive','-','color',cmap(1,:),'linewidth',2)
hold on 

% sort data for fitting
idx=~isnan(Y) ;
X=reshape(X(idx),1,[])';
Y=reshape(Y(idx),1,[])';
N=reshape(N(idx),1,[])';

% best fit assuming mucosalAb=1 by definition for mucosally-naive
t=1:60;
model=@(b,t) log10(sheddingViralLoadGMT(1,ageAdjust,t,0,[b(1).^2,b(2).^2,0,b(3).^2]));
[beta,R,~,COVB]=nlinfit(X,Y,model,[1.4,0.5,0.02],'weights',N);
parCI=nlparci(beta,R,'covar',COVB);
[YPRED, DELTA] = nlpredci(model,t,beta,R,'covar',COVB);
CI=[10.^(YPRED-2*DELTA);10.^(YPRED+2*DELTA)]';

% covert to natural unit
beta=beta.^2
parCI=parCI.^2

plotWithCI(t,10.^YPRED,CI,[0.4,0.8,0.4],'-',10,false,true)

set(gca,'yscale','log');
ylim([10^2,10^5]);
xlim([0, 50])
xlabel('days')
ylabel('fecal viral load (TCID50/g)')
legend('naive','model')


%% immune exploration (age-adjusted): type 2 only since most complete data

[immunizationPlans,placeIdx]=setdiff(unique(data.primaryImmunization),{'seronegative','IPVx2','IPVx3','natural immunity'});

% load immunity estimates
[~,~,immDat]=xlsread('fitImmunityAllSheddingCDFs.csv');
immDat=immDat(2:end,:);
immDatIdx=ismember(immDat(:,1),immunizationPlans)&cell2mat(immDat(:,2))==2;
bestNAb=cell2mat(immDat(immDatIdx,3))';
CINAb=cell2mat(immDat(immDatIdx,4:5))';

cmap=colormap(lines(length(immunizationPlans)));
cmap(end,:)=[0,0,0];
colormap(cmap);

meanGMT=[];
meanSamp=[];
figure(4); clf; 
subplot(1,2,1), hold all;
hold on 
for k=1:length(immunizationPlans)
    for n=1:3;
        idx=find(strcmp(data.meanSheddingImmunizationPlan,immunizationPlans(k))&data.meanSheddingSerotype==n);
        if any(idx) 
            plotSpanNaN(data.daysSinceChallenge,10.^data.meanSheddingGMT(3*(placeIdx(k)-1)+n,:),'color',cmap(k,:))
            meanGMT(end+(1:length(idx)))=nanmean(data.meanSheddingGMT(idx,:));
            meanSamp(end+(1:length(idx)))=data.meanSheddingSampleSize(idx);
        end
    end
    ylim(10.^[2 6])
    xlim([1 30])
end
xlabel('days')
ylabel('fecal viral load (TCID50/g)')
set(gca,'yscale','log')

subplot(1,2,2); hold on
for k=1:length(meanGMT)
    plot(bestNAb(k),10.^meanGMT(k),'o','markerfacecolor',cmap(k,:),'markeredgecolor',cmap(k,:))
    plot(CINAb(:,k),10.^meanGMT(k)*[1,1],'-','color',cmap(k,:))
end
set(gca,'yscale','log','xscale','log');
xlabel('OPV-equivalent humoral antibody titer')
ylabel('mean viral load (TCID50/g)')
xlim([.5,1e4])

% fit Antibody-dependent titer model and propagate edge uncertainty
AbRange=(exp(linspace(log(min(bestNAb)),log(max(bestNAb)),100)));
model=@(b,NAb) 10.^(b(1)*(1-b(2)*log2(NAb)));

% best-fit
[beta,R,~,COVB]=nlinfit(bestNAb,10.^meanGMT,model,[5,0.05],'weights',meanSamp);
[YPRED,DELTA] = nlpredci(model,AbRange,beta,R,'covar',COVB);

% CI sampler
reps=1000;
betaBoot=zeros(2,reps);
yCI=zeros(length(AbRange),reps);
for k=1:reps
    for n=1:length(CINAb)
        tmp=logspace(log10(CINAb(1,n)),log10(CINAb(2,n)),reps);
        idx=0;
        while idx<1 || idx>reps
            idx=round(reps/2+randn*reps/4);
        end
        bootNAb(n)=tmp(idx);
    end
    [betaBoot(:,k),R,~,COVB]=nlinfit(bootNAb,10.^meanGMT,model,beta,'weights',meanSamp);
    yCI(:,k)=nlpredci(model,AbRange,betaBoot(:,k),R,'covar',COVB);
end

parCI=prctile(betaBoot',[2.5,97.5])
CI=prctile(yCI',[2.5,97.5]);

plotWithCI(AbRange,YPRED',CI',[0.4,0.8,0.4],'-',10,false,true)
% encoded in shedViralLoadGMT
ylim(10.^[1,6])


h=colorbar;
set(h,'ytick',(0.5+[0:length(immunizationPlans)])/length(immunizationPlans),'yticklabel',immunizationPlans)

subplot(1,2,1); hold on
plot(t,sheddingViralLoadGMT(2^0,12,t),'--r')
plot(t,sheddingViralLoadGMT(2^3,12,t),'--g')
plot(t,sheddingViralLoadGMT(2^11,12,t),'--b')

