function data = loadSheddingViralLoadData()


% load
addpath(genpath('..\helperFunctions'));
addpath(genpath('..\..\data\shedding'));

[numData,~,allData]=xlsread('sheddingTitersGivenPositive.xlsx');

%% structure data
allData=allData(1:28,:); % for some reason, data is loading lots of nan rows
excludeIdx=strcmp(allData(:,38),'yes'); % remove studies excluded from analysis
allData=allData(~excludeIdx,:);
numData=numData(~excludeIdx,:);

% extract useful arrays
data.daysSinceChallenge=numData(1,2:26);
data.sampleSize=numData(2:end,1);
data.sheddingGMT=numData(2:end,2:26);


data.ageLastRIMonths=round(numData(2:end,27)*12/365);
data.challengeAgeMonths=round(numData(2:end,28)*12/365);
data.challengeSerotype=numData(2:end,29);
data.challengeDose=allData(2:end,37);

data.primaryImmunization=allData(2:end,35);
data.immunizationSchedule=allData(2:end,36);
data.publicationDate=numData(2:end,32);
data.location=allData(2:end,33);
data.author=allData(2:end,31);

% average by category

data.immunizationPlans=unique(data.primaryImmunization);
data.meanSheddingGMT=nan(length(data.immunizationPlans)*3,length(data.daysSinceChallenge));
data.meanSheddingSampleSize=nan(length(data.immunizationPlans)*3,1);
data.meanSheddingSerotype=nan(length(data.immunizationPlans)*3,1);
data.meanSheddingImmunizationPlan=cell(length(data.immunizationPlans)*3,1);

for k=1:length(data.immunizationPlans)
    for n=1:3
        data.meanSheddingSerotype(3*(k-1)+n)=n;
        data.meanSheddingImmunizationPlan(3*(k-1)+n)=data.immunizationPlans(k);
        idx=find(strcmp(data.primaryImmunization,data.immunizationPlans(k))&data.challengeSerotype==n);
        if ~isempty(idx)
            tmp=data.sheddingGMT(idx,:);
            % age adjustment to 12 months
            for m=1:length(idx)
                tmp(m,:)=tmp(m,:)*peakSheddingAgeMultiplier(12)/peakSheddingAgeMultiplier(data.challengeAgeMonths(idx(m)));
            end
            
            for m=1:length(idx)
                tmp(m,:)=tmp(m,:)*data.sampleSize(idx(m));
            end
            sampSizeTmp=nan(size(tmp));
            for m=1:length(sampSizeTmp)
                sampSizeTmp(:,m)=sum(data.sampleSize(idx(~isnan(tmp(:,m)))));
            end
                        
            tmp2=nan(size(data.daysSinceChallenge));
            for m=1:length(tmp2)
                tmp2(m)=nansum(tmp(:,m))/sum(data.sampleSize(idx(~isnan(tmp(:,m)))));
            end
            tmp2(tmp2==0)=nan;
            sampSizeTmp(sampSizeTmp==0)=nan;
            data.meanSheddingGMT(3*(k-1)+n,:)=tmp2;
            data.meanSheddingSampleSize(3*(k-1)+n,:)=round(nanmean(sampSizeTmp(1,:)));
        end
    end
end



end