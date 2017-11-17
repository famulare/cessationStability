% explore shedding duration data

data=loadSheddingDurationData;

%% plot all data by category
immunizationPlans=unique(data.primaryImmunization);
cmap=cbrewer('qual','Paired',length(immunizationPlans));

figure(1); clf;
for k=1:length(immunizationPlans)
    subplot(5,4,k)
    idx=strcmp(data.primaryImmunization,immunizationPlans(k));
    plot(data.daysSinceChallenge,data.sheddingCDF(idx,:),'color',cmap(k,:))
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end


%% plot category and serotype
immunizationPlans=unique(data.primaryImmunization);


figure(2); clf;
cmap=colormap(lines(3));
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:3
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeSerotype==n;
        if any(idx)
            plot(data.daysSinceChallenge,data.sheddingCDF(idx,:),'color',cmap(n,:))
        end
    end
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end

%% plot average by category and serotype
cmap=colormap(lines(3));

figure(3); clf;
for k=1:length(data.immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:3
        idx=find(strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeSerotype==n);
        if any(idx)
            plot(data.daysSinceChallenge,data.meanSheddingCDF(3*(k-1)+n,:),'color',cmap(n,:))
        end
    end
    if k==1,
        text(2,0.5,'1','color',cmap(1,:));
        text(2,0.35,'2','color',cmap(2,:));
        text(2,0.2,'3','color',cmap(3,:));
    end
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end


%% plot category and study date
immunizationPlans=unique(data.primaryImmunization);
dates=unique(data.publicationDate);
dateBins=min(dates):max(dates);
cmap=cbrewer('div','RdYlBu',length(dateBins));
for k=1:length(dateBins)
    cmap(k,:)=cmap(k,:)*(0.3*(1-2*k/length(dateBins)).^2+1-0.3);
end
dateIdx=ismember(dateBins,dates);
cmap=cmap(dateIdx,:);


figure(4); clf; colormap(cmap);
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:length(dates)
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&data.publicationDate==dates(n);
        if any(idx)
            plot(data.daysSinceChallenge,data.sheddingCDF(idx,:),'color',cmap(n,:))
        end
    end
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
%     if k==6
        legend(data.location(strcmp(data.primaryImmunization,immunizationPlans(k))))
%     end
end
h=colorbar;
set(h,'ytick',(0.5+[0:length(dates)])/length(dates),'yticklabel',dates)

%% plot category and age at challenge
immunizationPlans=unique(data.primaryImmunization);
ages=unique(data.challengeAgeMonths);
ageBins=0:26;
cmap=cbrewer('div','RdYlBu',length(ageBins));
for k=1:length(ageBins)
    cmap(k,:)=cmap(k,:)*(0.3*(1-2*k/length(ageBins)).^2+1-0.3);
end
ageIdx=ismember(ageBins,ages);
cmap=cmap(ageIdx,:);
cmap(end+1,:)=[0,0,0];

figure(5); clf; colormap(cmap);
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:length(ages)
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&data.challengeAgeMonths==ages(n);
        if any(idx)
            plot(data.daysSinceChallenge,data.sheddingCDF(idx,:),'color',cmap(n,:))
        end
    end
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end
h=colorbar;
set(h,'ytick',(0.5+[0:length(ages)])/length(ages),'yticklabel',ages)



%% plot category and challenge type
immunizationPlans=unique(data.primaryImmunization);
doses=unique(data.challengeDose);
doseBins=1:6;
cmap=colormap(lines(length(doseBins)));
cmap=cmap([1,3:6,2],:);

figure(6); clf; colormap(cmap);
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:length(doses)
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&ismember(data.challengeDose,doses(n));
        if any(idx)
            plot(data.daysSinceChallenge,data.sheddingCDF(idx,:),'color',cmap(n,:))
        end
    end
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
%     if k==6
%         legend(data.location(strcmp(data.primaryImmunization,immunizationPlans(k))))
%     end
end
h=colorbar;
set(h,'ytick',(0.5+[0:length(doses)])/length(doses),'yticklabel',doses)


%% plot category and time between last dose and challenge
immunizationPlans=unique(data.primaryImmunization);
ages=unique(data.challengeAgeMonths-data.ageLastRIMonths);
ageBins=0:7;
cmap=cbrewer('div','RdYlBu',length(ageBins));
for k=1:length(ageBins)
    cmap(k,:)=cmap(k,:)*(0.3*(1-2*k/length(ageBins)).^2+1-0.3);
end
ageIdx=ismember(ageBins,ages);
cmap=cmap(ageIdx,:);
cmap(end+1,:)=[0,0,0];

figure(5); clf; colormap(cmap);
for k=1:length(immunizationPlans)
    subplot(5,4,k); hold all;
    for n=1:length(ages)
        idx=strcmp(data.primaryImmunization,immunizationPlans(k))&(data.challengeAgeMonths-data.ageLastRIMonths)==ages(n);
        if any(idx)
            plot(data.daysSinceChallenge,data.sheddingCDF(idx,:),'color',cmap(n,:))
        end
    end
    ylim([0 1])
    xlim([0 45])
    title(immunizationPlans(k),'fontweight','normal')
end
h=colorbar;
set(h,'ytick',(0.5+[0:sum(~isnan(ages))])/sum(~isnan(ages)),'yticklabel',ages)

