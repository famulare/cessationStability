% fit shedding durations

clear all; close all;
data=loadSheddingDurationData;

%% fit mucosally-naive shedding duration data WPV
cmap=colormap(lines(4));
cmap=cmap([1,2,4],:);

% find mucosally-naive studies
idx=ismember(data.meanSheddingImmunizationPlan,{'naive-wpv'});
X=repmat(data.daysSinceChallenge',1,sum(idx));
Y=data.meanSheddingCDF(idx,:)';
N=repmat(data.meanSheddingSampleSize(idx),1,length(Y))';

figure(1); %clf
hold all; 
% plot approximate mean shedding duration CDF
meanY=nansum(Y.*N,2)./nansum(N.*~isnan(Y),2);
meanY(meanY==0)=nan;
plot(data.daysSinceChallenge,meanY,'color',cmap(1,:),'linewidth',2)

% sort data for fitting
idx=~isnan(Y);
X=reshape(X(idx),1,[])';
Y=reshape(Y(idx),1,[])';
N=reshape(N(idx),1,[])';

% best fit assuming mucosalAb=1 by definition for mucosally-naive
custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(sheddingDurationCDF(1,X,0,[b(1),0,b(2)])) + N.*(1-Y).*log(1-sheddingDurationCDF(1,X,0,[b(1),0,b(2)])));
[beta] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',[log(30),log(2)]) ;

% bootstrap CI since binomial likelihood isn't proper survival curve
% estimator, but Kaplan-Meier and its uncertainty estimators aren't
% really available because I don't have the individual data
reps=1000;
betaBoot=nan(2,reps);
X=unique(X);
Yboot=zeros(length(X),reps);
N=sum(unique(N))*ones(size(X));
custnloglf = @(b,D,cens,N) -sum(D.*N.*log(sheddingDurationCDF(1,X,0,[b(1),0,b(2)])) + N.*(1-D).*log(1-sheddingDurationCDF(1,X,0,[b(1),0,b(2)])));
for k=1:reps
    tmp=nan(length(X),1);
    for n=1:length(X)
        idx=rand<=1-sheddingDurationCDF(1,X,0,[beta(1),0,beta(2)]);
        if any(idx)
            tmp(n)=min(X(idx));
        end
    end
    Yboot(:,k)=(length(tmp)-cumsum(hist(tmp,X)))'/length(tmp);
    try % rarely, inference chokes
        [betaBoot(:,k)] = mle(Yboot(:,k),'nloglf',custnloglf,'frequency',N,'start',beta) ;
    end
end
CI=quantile(betaBoot',[0.025,0.975]);

yCI=quantile(Yboot',[0.025,0.975])';


exp(beta)
exp(CI)

tmpMean=sheddingDurationCDF(1,X,0,[beta(1),0,beta(2)]);
plotWithCI(X,tmpMean,yCI,cmap(1,:))
xlabel('days since OPV challenge')
ylabel('probability shedding')
% title('Shedding duration ','fontweight','normal');


%% fit mucosally-naive shedding duration data OPV

% find mucosally-naive studies
idx=ismember(data.meanSheddingImmunizationPlan,{'seronegative','unvaccinated'});
X=repmat(data.daysSinceChallenge',1,sum(idx));
Y=data.meanSheddingCDF(idx,:)';
N=repmat(data.meanSheddingSampleSize(idx),1,length(Y))';

% plot approximate mean shedding duration CDF
meanY=nansum(Y.*N,2)./sum(N.*~isnan(Y),2);
meanY(meanY==0)=nan;
plot(data.daysSinceChallenge,meanY,'color',cmap(2,:),'linewidth',2)

% sort data for fitting
idx=~isnan(Y);
X=reshape(X(idx),1,[])';
Y=reshape(Y(idx),1,[])';
N=reshape(N(idx),1,[])';

% best fit assuming mucosalAb=1 by definition for mucosally-naive
custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(sheddingDurationCDF(1,X,0,[b(1),0,b(2)])) + N.*(1-Y).*log(1-sheddingDurationCDF(1,X,0,[b(1),0,b(2)])));
[beta] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',[log(30),log(2)]) ;

% bootstrap CI since binomial likelihood isn't proper survival curve
% estimator, but Kaplan-Meier and its uncertainty estimators aren't
% really available because I don't have the individual data
reps=1000;
betaBoot=zeros(2,reps);
X=unique(X);
Yboot=zeros(length(X),reps);
N=sum(unique(N))*ones(size(X));
custnloglf = @(b,D,cens,N) -sum(D.*N.*log(sheddingDurationCDF(1,X,0,[b(1),0,b(2)])) + N.*(1-D).*log(1-sheddingDurationCDF(1,X,0,[b(1),0,b(2)])));
for k=1:reps
    tmp=nan(length(X),1);
    for n=1:length(X)
        idx=rand<=1-sheddingDurationCDF(1,X,0,[beta(1),0,beta(2)]);
        if any(idx)
            tmp(n)=min(X(idx));
        end
    end
    Yboot(:,k)=(length(tmp)-cumsum(hist(tmp,X)))'/length(tmp);
    try % rarely, inference chokes
        [betaBoot(:,k)] = mle(Yboot(:,k),'nloglf',custnloglf,'frequency',N,'start',beta) ;
    end
end
CI=quantile(betaBoot',[0.025,0.975]);

yCI=quantile(Yboot',[0.025,0.975])';


exp(beta)
exp(CI)

tmpMean=sheddingDurationCDF(1,X,0,[beta(1),0,beta(2)]);
plotWithCI(X,tmpMean,yCI,cmap(2,:))
% xlabel('days since OPV challenge')
% ylabel('probability shedding')
% title('Shedding duration ','fontweight','normal');



%% Maximum mucosal immunity: fit type 2 after 3 tOPV
idx=ismember(data.author,{'Asturias'})&ismember(data.primaryImmunization,{'tOPVx3'});
X=repmat(data.daysSinceChallenge',1,sum(idx));
Y=data.sheddingCDF(idx,:)';
N=repmat(data.sampleSize(idx),1,length(Y))';

figure(1); 
% plot approximate mean shedding duration CDF
meanY=nansum(Y.*N,2)./sum(N.*~isnan(Y),2);
meanY(meanY==0)=nan;
plot(data.daysSinceChallenge,meanY,'color',cmap(3,:),'linewidth',2)

idx=~isnan(Y(:,1));
X=reshape(X(idx,:),1,[])';
Y=reshape(Y(idx,:),1,[])';
N=reshape(N(idx,:),1,[])';

% best fit at maximal immunity: 1 month after tOPVx3 at mucosalAb=2^11
custnloglf = @(b,Y,cens,N) -sum(Y.*N.*log(sheddingDurationCDF(2^11,X,0,[beta(1),b,beta(2)])) + N.*(1-Y).*log(1-sheddingDurationCDF(2^11,X,0,[beta(1),b,beta(2)])));
[beta2,CI2] = mle(Y,'nloglf',custnloglf,'frequency',N,'start',0) ;

% bootstrap CI 
reps=1000;
beta2Boot=zeros(1,reps);
X=unique(X);
Yboot=zeros(length(X),reps);
N=sum(unique(N))*ones(size(X));
custnloglf = @(b,D,cens,N) -sum(D.*N.*log(sheddingDurationCDF(2^11,X,0,[beta(1),b,beta(2)])) + N.*(1-D).*log(1-sheddingDurationCDF(2^11,X,0,[beta(1),b,beta(2)])));
for k=1:reps
    tmp=nan(length(X),1);
    for n=1:length(X)
        idx=rand<=1-sheddingDurationCDF(2^11,X,0,[beta(1),beta2,beta(2)]);
        if any(idx)
            tmp(n)=min(X(idx));
        end
    end
    Yboot(:,k)=(length(tmp)-cumsum(hist(tmp,X)))'/length(tmp);
    try % rarely, inference chokes
        [beta2Boot(:,k)] = mle(Yboot(:,k),'nloglf',custnloglf,'start',beta2);
    end
end
CI2=quantile(beta2Boot',[0.025,0.975]);

yCI2=quantile(Yboot',[0.025,0.975])';

exp(beta2)
exp(CI2)

tmpMean=sheddingDurationCDF(2^11,X,0,[beta(1),beta2,beta(2)]);
plotWithCI(X,tmpMean,yCI2,cmap(3,:))
xlabel('days since OPV challenge')
ylabel('probability shedding')
% title('Shedding duration ','fontweight','normal');

text(40,1,'naive, WPV','color',cmap(1,:))
text(40,0.9,'naive, OPV','color',cmap(2,:))
text(40,0.8,'max immunity, OPV','color',cmap(3,:))
xlim([0 60]);

%% best fit to all vaccine schedules

cmap=colormap(lines(3));

beta=nan(length(data.immunizationPlans),3);
CI=nan(length(data.immunizationPlans),2,3);

medianShedDur=beta;
ciShedDur=CI;
sampleSize=beta;

[~,idx]=sort(nanmean(beta,2));

figure(14); clf; 
for k=5;1:length(data.immunizationPlans)    
    idx=find(ismember(data.meanSheddingImmunizationPlan,data.immunizationPlans(k)));
    subplot(5,4,k)
    
    hold all
    for n=2;1:length(idx)
        if ~all(isnan(data.meanSheddingCDF(idx(n),:)))
            plot(data.daysSinceChallenge,data.meanSheddingCDF(idx(n),:),'color',cmap(n,:))
           
            X=repmat(data.daysSinceChallenge',1,length(idx(n)));
            Y=data.meanSheddingCDF(idx(n),:)';
            N=repmat(data.meanSheddingSampleSize(idx(n)),1,length(Y))';

            try
                sampleSize(k,n)=data.meanSheddingSampleSize(idx(n));
                [beta(k,n),CI(k,:,n),yCI]=sheddingDurationCDFbootstrapFitter(X,Y,N,1000);
                plot(data.daysSinceChallenge,sheddingDurationCDF(beta(k,n),data.daysSinceChallenge),'--','color',cmap(n,:));
                text(30,0.9-0.1*(n-1),['N_{Ab} = ',num2str(round(beta(k,n))),' (',num2str(round(CI(k,1,n))),',',num2str(round(CI(k,2,n))),')'],'color',cmap(n,:),'fontsize',8);
                
                % median shedding duration
                medianShedDur(k,n)=data.daysSinceChallenge(find(0.5>sheddingDurationCDF(beta(k,n),data.daysSinceChallenge),1));
                ciShedDur(k,1,n)=data.daysSinceChallenge(find(0.5>sheddingDurationCDF(CI(k,1,n),data.daysSinceChallenge),1));
                ciShedDur(k,2,n)=data.daysSinceChallenge(find(0.5>sheddingDurationCDF(CI(k,2,n),data.daysSinceChallenge),1));
            end
        end
    end
    xlim([0,70])
    if k==1
        text(2,0.5,'Type 1','color',cmap(1,:),'fontsize',8);
        text(2,0.35,'Type 2','color',cmap(2,:),'fontsize',8);
        text(2,0.2,'Type 3','color',cmap(3,:),'fontsize',8);
    end
    title(data.immunizationPlans{k},'fontweight','normal','fontsize',10);

end

%% write out immunity fits
fid=fopen('fitImmunityAllSheddingCDFs.csv','w');
fprintf(fid,'vaccineSchedule,challengeSerotype,bestNab,lowerCINAb,upperCINAb,medianDur,lowerCIDur,upperCIDur,sampleSize\r\n');
for k=1:length(data.immunizationPlans) 
    for n=1:3
        if ~isnan(beta(k,n))
            fprintf(fid,'%s,%f,%f,%f,%f,%f,%f,%f,%f\r\n',data.immunizationPlans{k} ...
                    ,n,beta(k,n),CI(k,1,n),CI(k,2,n),medianShedDur(k,n) ...
                    ,ciShedDur(k,2,n),ciShedDur(k,1,n),sampleSize(k,n));
        end
    end
end
fclose(fid);

%% reload and plot all
[shedDur,immSched]=xlsread('fitImmunityAllSheddingCDFs.csv');
figure(17); clf; 

% sort by mean mucosalAb
[schedules,~,idx]=unique(immSched(2:end,1));
sortMean=zeros(size(schedules));
for k=1:length(schedules)
    sortMean(k)=mean(shedDur(idx==k,2));
end
[~,sortIdx]=sort(sortMean);

% plot each result
cmap=colormap(lines(3));
for n=1:length(shedDur) 

    pos=find(idx(n)==sortIdx);
    sero=shedDur(n,1);
    
    subplot(1,2,2)
    hold on;
    plot(log2(shedDur(n,2)),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(log2(shedDur(n,3:4)),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))

    subplot(1,2,1)
    hold on;
    plot(shedDur(n,5),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(shedDur(n,6:7),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))
    
end
subplot(1,2,2)
xlim([-4 14])
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
set(gca,'xtick',-2:2:12,'xticklabel',2.^(-2:2:13))
xlabel('OPV-equivalent mucosal antibody titer')
set(gca,'ygrid','on')

subplot(1,2,1)
xlim([0 50]);
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
xlabel('median shedding duration (days)')
text(2,3.5,'Type 1','color',cmap(1,:));
text(2,2.5,'Type 2','color',cmap(2,:));
text(2,1.5,'Type 3','color',cmap(3,:));

set(gca,'ygrid','on')


%% reload and plot OPV
[shedDur,immSched]=xlsread('fitImmunityAllSheddingCDFs.csv');

idx=~cellfun(@isempty,regexp(immSched(:,1),'tOPV')) ...
    | ~cellfun(@isempty,regexp(immSched(:,1),'mOPV')) ...
    | ~cellfun(@isempty,regexp(immSched(:,1),'unvaccinated'));
shedDur=shedDur(idx(2:end),:);
immSched=immSched(idx,:);
figure(18); clf; 

% sort by mean mucosalAb
[schedules,~,idx]=unique(immSched(:,1));
sortMean=zeros(size(schedules));
for k=1:length(schedules)
    sortMean(k)=mean(shedDur(idx==k,2));
end
[~,sortIdx]=sort(sortMean);

subplot(1,2,2)
plot([0,0],[0.5,length(schedules)+0.5],'color',[0,0,0,0.15])

% plot each result
cmap=colormap(lines(3));
for n=1:length(shedDur) 

    pos=find(idx(n)==sortIdx);
    sero=shedDur(n,1);
    
    subplot(1,2,2)
    hold on;
    plot(log2(shedDur(n,2)),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(log2(shedDur(n,3:4)),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))

    subplot(1,2,1)
    hold on;
    plot(shedDur(n,5),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(shedDur(n,6:7),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))
    
end
subplot(1,2,2)
xlim([-4 14])
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
set(gca,'xtick',-2:2:12,'xticklabel',2.^(-2:2:13))
xlabel('OPV-equivalent mucosal antibody titer')
set(gca,'ygrid','on')
box off

subplot(1,2,1)
xlim([0 50]);
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
xlabel('median shedding duration (days)')
text(2,3.5,'Type 1','color',cmap(1,:));
text(2,2.5,'Type 2','color',cmap(2,:));
text(2,1.5,'Type 3','color',cmap(3,:));

set(gca,'ygrid','on')


%% reload and plot IPV
[shedDur,immSched]=xlsread('fitImmunityAllSheddingCDFs.csv');

idx=(~cellfun(@isempty,regexp(immSched(:,1),'IPV')) &  cellfun(@isempty,regexp(immSched(:,1),'OPV'))) ...
    | ~cellfun(@isempty,regexp(immSched(:,1),'unvaccinated')) ;
shedDur=shedDur(idx(2:end),:);
immSched=immSched(idx,:);
figure(19); clf; 

% sort by number of IPV
[schedules,~,idx]=unique(immSched(:,1));
sortIdx=[4,1,2,3];

subplot(1,2,2)
plot([0,0],[0.5,length(schedules)+0.5],'color',[0,0,0,0.15])

% plot each result
cmap=colormap(lines(3));
for n=1:length(shedDur) 

    pos=find(idx(n)==sortIdx);
    sero=shedDur(n,1);
    
    subplot(1,2,2)
    hold on;
    plot(log2(shedDur(n,2)),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(log2(shedDur(n,3:4)),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))

    subplot(1,2,1)
    hold on;
    plot(shedDur(n,5),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(shedDur(n,6:7),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))
    
end
subplot(1,2,2)
xlim([-4 14])
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
set(gca,'xtick',-2:2:12,'xticklabel',2.^(-2:2:13))
xlabel('OPV-equivalent mucosal antibody titer')
set(gca,'ygrid','on')
box off

subplot(1,2,1)
xlim([0 50]);
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
xlabel('median shedding duration (days)')
text(2,3.5,'Type 1','color',cmap(1,:));
text(2,2.5,'Type 2','color',cmap(2,:));
text(2,1.5,'Type 3','color',cmap(3,:));

set(gca,'ygrid','on')


%% reload and plot bOPV cross-protection
[shedDur,immSched]=xlsread('fitImmunityAllSheddingCDFs.csv');

idx=~cellfun(@isempty,regexp(immSched(:,1),'bOPV')) ...
    | (~cellfun(@isempty,regexp(immSched(:,1),'tOPVx1')) & [0;shedDur(:,1)]==2)...
    | (~cellfun(@isempty,regexp(immSched(:,1),'unvaccinated'))& [0;shedDur(:,1)]==2) ;
shedDur=shedDur(idx(2:end),:);
immSched=immSched(idx,:);
figure(20); clf; 

% sort by mean mucosalAb
[schedules,~,idx]=unique(immSched(:,1));
sortMean=zeros(size(schedules));
for k=1:length(schedules)
    sortMean(k)=mean(shedDur(idx==k,2));
end
[~,sortIdx]=sort(sortMean);
% sortIdx=[7,3,1,2,4,5,6]; % tOPVx1 on top

subplot(1,2,2)
plot([0,0],[0.5,length(schedules)+0.5],'color',[0,0,0,0.15])

% plot each result
cmap=colormap(lines(3));
for n=1:length(shedDur) 

    pos=find(idx(n)==sortIdx);
    sero=shedDur(n,1);
    
    subplot(1,2,2)
    hold on;
    plot(log2(shedDur(n,2)),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(log2(shedDur(n,3:4)),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))

    subplot(1,2,1)
    hold on;
    plot(shedDur(n,5),pos+(sero-2)*.15,'.','color',cmap(sero,:));
    plot(shedDur(n,6:7),(pos+(sero-2)*.15)*[1,1],'-','color',cmap(sero,:))
    
end
subplot(1,2,2)
xlim([-4 14])
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
set(gca,'xtick',-2:2:12,'xticklabel',2.^(-2:2:13))
xlabel('OPV-equivalent mucosal antibody titer')
set(gca,'ygrid','on')
box off

subplot(1,2,1)
xlim([0 50]);
ylim([0.5,length(schedules)+0.5]);
set(gca,'ytick',1:length(schedules))
set(gca,'yticklabel',schedules(sortIdx));
xlabel('median shedding duration (days)')
text(2,2.5,'Type 2','color',cmap(2,:));

set(gca,'ygrid','on')



%% agregates by serotype

cmap=colormap(lines(3));

[~,idx]=sort(nanmean(beta,2));

immPlanGroups.name={'mucosal-naive','OPVx1','bOPV','bOPV+IPV','OPVx2','OPVx3','natural immunity'};
immPlanGroups.set={{'IPV/seronegative','naive','IPVx3','IPVx4'} ...
                    ,{'mOPV1x1','tOPVx1'} ...
                    ,{'bOPV13x3','IPVx2&bOPV13x1','IPVx1&bOPV13x2'} ...
                    ,{'bOPV13x3&IPVx1','bOPV13x3&IPVx2'} ...
                    ,{'IPVx2&tOPVx2','tOPVx2'} ...
                    ,{'tOPVx3','IPVx3&tOPVx3'} ...
                    ,{'natural immunity'}};
immPlanGroups.exclude={{'IPVx3',1},{},{},{},{},{},{}};                


beta=nan(length(immPlanGroups.name),3);
CI=nan(length(immPlanGroups.name),2,3);

for k=1:length(immPlanGroups.name)    
    for n=1:3
        idx=find(ismember(data.meanSheddingImmunizationPlan,immPlanGroups.set{k}) & data.meanSheddingSerotype==n);
        if ~isempty(immPlanGroups.exclude{k})
            for m=1:size(immPlanGroups.exclude{k},1)
                exIdx=find(ismember(data.meanSheddingImmunizationPlan,immPlanGroups.exclude{k}{m,1}) & data.meanSheddingSerotype==immPlanGroups.exclude{k}{m,2});
            end
            idx=setdiff(idx,exIdx);
        end


        X=repmat(data.daysSinceChallenge',length(idx),1);
        N=reshape(repmat(data.meanSheddingSampleSize(idx),1,length(data.daysSinceChallenge))',[],1);
        Y=reshape(data.meanSheddingCDF(idx,:)',[],1);

        X=X(~isnan(Y));
        N=N(~isnan(Y));
        Y=Y(~isnan(Y));
        if ~isempty(Y)
            [beta(k,n),CI(k,:,n)]=sheddingDurationCDFbootstrapFitter(X,Y,N);    
        end
    end

end

figure(18); clf; 
hold on;
[~,idx]=sort(nanmean(beta,2));
cmap=colormap(lines(3));
for k=1:length(immPlanGroups.name) 
    for n=1:3
%         subplot(1,3,n); hold on;
        plot(log2(beta(idx(k),n)),k+(n-2)*.15,'.','color',cmap(n,:));
        plot(log2(CI(idx(k),:,n)),(k+(n-2)*.15)*[1,1],'-','color',cmap(n,:))
    
        xlim([-2 14])
        ylim([0.5,length(immPlanGroups.name(idx))+0.5]);
        set(gca,'ytick',1:length(immPlanGroups.name(idx)))
        set(gca,'yticklabel',immPlanGroups.name(idx));
        set(gca,'xtick',0:2:12,'xticklabel',2.^(0:2:12))
        xlabel('antibody titer')
        
        text(log2(2000),4.5-n,['Type ',num2str(n)],'color',cmap(n,:));
        set(gca,'ygrid','on')
    end
end