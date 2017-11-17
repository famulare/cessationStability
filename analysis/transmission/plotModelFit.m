% plotModelFit


addpath(genpath('..\Houston1960'))
addpath(genpath('..\helperFunctions'))
addpath(genpath('..\doseResponse'))

[index,sibling,contact]=loadHouston1960();

cmap(1,:)=[229,38,40]/255;
cmap(2,:)=[62,130,186]/255;
cmap(3,:)=[71,172,68]/255;

timeAxis=repmat([0,7,14,21,28,35]',1,9);
X=timeAxis(:);

s=1;
Y=[0;index(s).SheddingTotal;0;sibling(s).SheddingUnder5;0;contact(s).SheddingUnder5];
N=[max(index(s).DenominatorTotal);index(s).DenominatorTotal;max(sibling(s).DenominatorUnder5);sibling(s).DenominatorUnder5;max(contact(s).DenominatorUnder5);contact(s).DenominatorUnder5];    

for s=2:3
    Y=[Y;0;index(s).SheddingTotal;0;sibling(s).SheddingUnder5;0;contact(s).SheddingUnder5];
    N=[N;max(index(s).DenominatorTotal);index(s).DenominatorTotal;max(sibling(s).DenominatorUnder5);sibling(s).DenominatorUnder5;max(contact(s).DenominatorUnder5);contact(s).DenominatorUnder5];    
end


load('modelFit.mat')
tmax=100;
reps=1000;
modelBoot=nan((tmax+1)*9,reps);
for n=1:reps;
    [~,modelBoot(:,n)]=householdTransmissionLogLikelihood(betaBoot(:,n),X,Y,N,n);     
end
YCI=prctile(modelBoot',[2.5,97.5]);


[~,fit,~,fit2]=householdTransmissionLogLikelihood(beta,0:tmax,1,1);
fit=reshape(fit,tmax+1,9);
fit2=reshape(fit2,tmax+1,9);
for k=1:9
    idx=(k-1)*(tmax+1)+(1:tmax+1);
    yCI{k}=YCI(:,idx)';
end

%% sabin
figure(1); clf; 

for k=1:3
    
    % prevalence
    subplot(2,3,k);cla; hold on
    
    [p,ci]=binofit(index(k).SheddingTotal,index(k).DenominatorTotal);
    plotWithCI(-0.5+index(k).Day,p,ci,cmap(1,:),'.',10,true,false)
    idx=1+3*(k-1);
    plotWithCI(0:tmax,fit(:,idx),yCI{idx},cmap(1,:))
    
    [p,ci]=binofit(sibling(k).SheddingUnder5,sibling(k).DenominatorUnder5);
    plotWithCI(sibling(k).Day,p,ci,cmap(2,:),'.',10,true,false)
    idx=2+3*(k-1);
    plotWithCI(0:tmax,fit(:,idx),yCI{idx},cmap(2,:))
    
    [p,ci]=binofit(contact(k).SheddingUnder5,contact(k).DenominatorUnder5);
    plotWithCI(0.5+contact(k).Day,p,ci,cmap(3,:),'.',10,true,false)
    idx=3+3*(k-1);
    plotWithCI(0:tmax,fit(:,idx),yCI{idx},cmap(3,:))
    
    ylim([0,1.0])
    xlim([0,38])
    ylabel('fraction shedding')
    xlabel('days since index infection')
    set(gca,'xtick',0:7:35)
    
    % incidence
    subplot(2,3,3+k);cla; hold on
    idx=1+3*(k-1);
    plot(0:tmax,fit2(:,idx),'color',cmap(1,:))
    
    idx=2+3*(k-1);
    plot(0:tmax,fit2(:,idx),'color',cmap(2,:))
    
    idx=3+3*(k-1);
    plot(0:tmax,fit2(:,idx),'color',cmap(3,:))
    
    
    ylim([0,1.0])
    xlim([0,38])
    ylabel('incidence')
    xlabel('days since index infection')
    set(gca,'xtick',0:7:35)
end

%% louisiana

data=xls2struct('..\..\Data\Louisiana1957\householdSeroconversion.xlsx','seroconversion');

load('LouisianaFit.mat')

N=data.denominator([1,4,3,2,5]);
Y=data.seroconversions([1,4,3,2,5]);
NAb=data.medianPreExposureAntibodyTiter([1,4,3,2,5]);


% fit naive OlderSibling
experiment{1} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1 ...
                        ,'secondaryLog2NAb',log2(NAb(1)) ...
                        ,'doseResponseBeta',10^beta(1) ...
                        ,'runNetwork',false,'serotype',4);

% fit immune OlderSibling
experiment{2} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1 ...
                        ,'secondaryLog2NAb',log2(NAb(2)) ...
                        ,'doseResponseBeta',10^beta(1) ...
                        ,'runNetwork',false,'serotype',4);

% fit naive YoungerSibling
experiment{3} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1 ...
                        ,'secondaryContactAcquire',10^(-5.2943+beta(3)) ...
                        ,'secondaryLog2NAb',log2(NAb(3)) ...
                        ,'doseResponseBeta',10^beta(1) ...
                        ,'runNetwork',false,'serotype',4);
                    
                 
% fit naive adult
experiment{4} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+beta(2)) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+beta(2)) ...
                    ,'secondaryLog2NAb',log2(NAb(4)) ...
                    ,'doseResponseBeta',10^beta(1) ...
                    ,'runNetwork',false,'serotype',4);

% fit immune adult
experiment{5} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+beta(2)) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+beta(2)) ...
                    ,'secondaryLog2NAb',log2(NAb(5)) ...
                    ,'doseResponseBeta',10^beta(1) ...
                    ,'runNetwork',false,'serotype',4);

figure(4); clf; 
figure(2); clf;

for k=1:5

    figure(2);
    tmax=100;
    reps=1000;
    modelBoot=nan((tmax)*2,reps);
    for n=1:reps;
        if k<3
            ex = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'secondaryLog2NAb',log2(NAb(k)),'doseResponseBeta',10^betaBoot(1,n),'runNetwork',false,'serotype',4,'confidenceIntervalSamplerSeed',n);
        elseif k==3
            ex = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'secondaryContactAcquire',10^(-5.2943+betaBoot(3,n)),'secondaryLog2NAb',log2(NAb(k)),'doseResponseBeta',10^(betaBoot(1,n)),'runNetwork',false,'serotype',4,'confidenceIntervalSamplerSeed',n);
        else
            ex = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'secondaryContactAcquire',10^(-5.2943+betaBoot(2,n)),'tertiaryContactAcquire',10^(-4.3299+betaBoot(2,n)),'secondaryLog2NAb',log2(NAb(k)),'doseResponseBeta',10^(betaBoot(1,n)),'runNetwork',false,'serotype',4,'confidenceIntervalSamplerSeed',n);
        end
            
        modelBoot(:,n)=[cumsum(ex.primary.incidence(1:tmax)),cumsum(ex.secondary.incidence(1:tmax))];     
    end
    YCI=prctile(modelBoot',[2.5,97.5]);

    clear yCI
    for n=1:2
        idx=(n-1)*(tmax)+(1:tmax);
        yCI{n}=YCI(:,idx)';
    end

    % prevalence
    subplot(2,5,k);cla; hold on
    plot(1:tmax,experiment{k}.primary.prevalence(1:tmax),'color',cmap(1,:))
    plot(1:tmax,experiment{k}.secondary.prevalence(1:tmax),'color',cmap(2,:))
    plot(1:tmax,experiment{k}.tertiary.prevalence(1:tmax),'color',cmap(3,:))
    
    
    ylim([0,1.0])
    xlim([0,38])
    ylabel('fraction shedding')
    xlabel('days since index infection')
    set(gca,'xtick',0:7:35)
    
    % incidence
    subplot(2,5,5+k);cla; hold on
    
    [p,ci]=binofit(Y,N);
    plotWithCI(30,p(k),ci(k,:),cmap(2,:),'.',10,true,false)
    
    plot(1:tmax,cumsum(experiment{k}.primary.incidence(1:tmax)),'color',cmap(1,:))
    plotWithCI(1:tmax,cumsum(experiment{k}.secondary.incidence(1:tmax)),yCI{2},cmap(2,:))
    plot(1:tmax,cumsum(experiment{k}.tertiary.incidence(1:tmax)),'color',cmap(3,:))
    
    ylim([0,1.0])
    xlim([0,38])
    ylabel('incidence')
    xlabel('days since index infection')
    set(gca,'xtick',0:7:35)
    
    drawnow
    
    
    figure(4); 

    subplot(1,2,1); hold on;
    plotWithCI(k-.2,p(k),ci(k,:),cmap(2,:),'.',10,true,false)
    plotWithCI(k+.2,sum(experiment{k}.secondary.incidence(1:30)),yCI{2}(30,:),cmap(2,:),'s',10,true,false)
    ylim([0,1])
    
end




%% india

data=xls2struct('..\..\Data\India2010\modeledPrevalenceInHealthyContactsOfWPV1Cases.xlsx');
load('indiaFit.mat')

N=data.effectiveSampleSize;
Y=round(data.estimatedPrevalence.*N);
NAb=data.assumedNAb;

tmax=140;
experiment{1} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'serotype',4 ...
                        ,'probDailyPrimarySecondaryContact',1,'t',1:tmax ...
                        ,'tertiaryContactAcquire',10^(beta(1)+0.97) ...
                        ,'secondaryLog2NAb',log2(NAb(1)) ...
                        ,'tertiaryLog2NAb',log2(NAb(1)) ...
                        ,'secondaryContactAcquire',10^beta(1) ...
                        ,'runNetwork',false);

experiment{2} = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'serotype',4 ...
                        ,'probDailyPrimarySecondaryContact',1,'t',1:tmax ...
                        ,'tertiaryContactAcquire',10^(beta(1)+0.97) ...
                        ,'secondaryLog2NAb',log2(NAb(2)) ...
                        ,'tertiaryLog2NAb',log2(NAb(2)) ...
                        ,'secondaryContactAcquire',10^beta(1) ...
                        ,'runNetwork',false);


[~,pdf]=timeToParalysis(1:tmax);
                    
tmax=88; % post-paralysis observation window
q=nan(length(Y),1);
for k=1:length(Y);
    tmp=conv(experiment{k}.secondary.prevalence,pdf);
    q(k)=mean(tmp(1:tmax)); % stools collected after onset (weighted by probability)
end

figure(3); clf;

for k=1:2

    figure(3);
    reps=1000;
    modelBoot=nan(1,reps);
    for n=1:reps;
        ex = primarySecondaryTertiaryDoseModel('perDoseEfficacy',1,'secondaryLog2NAb',log2(NAb(k)),'tertiaryLog2NAb',log2(NAb(k)),'secondaryContactAcquire',10^betaBoot(n),'tertiaryContactAcquire',10^(betaBoot(n)+0.97),'runNetwork',false,'serotype',4,'confidenceIntervalSamplerSeed',n,'t',1:tmax);
        tmp=conv(ex.secondary.prevalence,pdf);
        modelBoot(:,n)=mean(tmp(1:tmax));
    end
    YCI=prctile(modelBoot',[2.5,97.5]);

    % prevalence
    subplot(2,2,k);cla; hold on
    plot(1:tmax,experiment{k}.primary.prevalence(1:tmax),'color',cmap(1,:))
    plot(1:tmax,experiment{k}.secondary.prevalence(1:tmax),'color',cmap(2,:))
    plot(1:tmax,experiment{k}.tertiary.prevalence(1:tmax),'color',cmap(3,:))
    
    [p,ci]=binofit(Y,N);
    
%     plotWithCI(1:tmax,cumsum(experiment{k}.secondary.incidence(1:tmax)),yCI{2},cmap(2,:))
    t=find(cdf>0.5,1):tmax;
    plotWithCI(t,p(k)*ones(size(t))',repmat(ci(k,:),length(t),1),cmap(2,:))%,'.',10,true,false)
    plotWithCI(t,q(k)*ones(size(t))',repmat(YCI,length(t),1),cmap(3,:))%,'.',10,true,false)
    
    ylim([0,1.0])
    xlim([0,tmax])
    ylabel('fraction shedding')
    xlabel('days since index infection')
    set(gca,'xtick',0:7:tmax)
    
    % incidence
    subplot(2,2,2+k);cla; hold on
    
    plot(1:tmax,cumsum(experiment{k}.primary.incidence(1:tmax)),'color',cmap(1,:))
    plot(1:tmax,cumsum(experiment{k}.secondary.incidence(1:tmax)),'color',cmap(2,:))
    plot(1:tmax,cumsum(experiment{k}.tertiary.incidence(1:tmax)),'color',cmap(3,:))
    
    ylim([0,1.0])
    xlim([0,tmax])
    ylabel('incidence')
    xlabel('days since index infection')
    set(gca,'xtick',0:7:tmax)
    
    drawnow
    
    figure(4); 

    subplot(1,2,2); hold on;
    plotWithCI(k-.2,p(k),ci(k,:),cmap(2,:),'.',10,true,false)
    plotWithCI(k+.2,q(k),YCI,cmap(2,:),'s',10,true,false)
	ylim([0,1])
end



