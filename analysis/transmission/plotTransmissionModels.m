% plot transmission models


addpath(genpath('..\Houston1960'))
addpath(genpath('..\helperFunctions'))
addpath(genpath('..\doseResponse'))
addpath(genpath('..\shedding'))
[index,sibling,contact]=loadHouston1960();

%%
typeStr={'primary','secondary','tertiary'};
infantLog2NAb=[0,log2(512),log2(8),log2(8)];
siblingLog2NAb=[0,log2(512),log2(256),log2(2)];

titleStr={'naive & IPV-only','tOPVx3',{'1 year after cessation','(mixed bOPV & tOPV)'},'bOPVx3'};

% figure 1 
cmap=cbrewer('qual','Set1',3);
cmap2=cbrewer('seq','YlOrRd',4);

fecalOffset=2.8; % india

figure(10); clf;
for n=1:4
    
    experiment{n} = primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',infantLog2NAb(n) ...
                    ,'secondaryLog2NAb',siblingLog2NAb(n),'tertiaryLog2NAb',siblingLog2NAb(n) ...
                    ,'perDoseEfficacy',1);
    
    subplot(2,4,n); hold all
    for m=1:3        
        plot(experiment{n}.params.t,experiment{n}.(typeStr{m}).prevalence,'color',cmap(m,:),'linewidth',1.5)
%         plot(experiment{n}.params.t,cumsum(experiment{n}.(typeStr{m}).incidence),'color',cmap(m,:),'linewidth',1.5)
    end
    
    experiment{n+3} = primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',infantLog2NAb(n) ...
                    ,'secondaryLog2NAb',siblingLog2NAb(n),'tertiaryLog2NAb',siblingLog2NAb(n) ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset) ...
                    ,'perDoseEfficacy',1);
    
    subplot(2,4,4+n); hold all
    for m=1:3        
        plot(experiment{n+3}.params.t,experiment{n+3}.(typeStr{m}).prevalence,'color',cmap(m,:),'linewidth',1.5)
%         plot(experiment{n+3}.params.t,cumsum(experiment{n+3}.(typeStr{m}).incidence),'color',cmap(m,:),'linewidth',1.5)
    end
end

for n=1:4
    for k=1:2
        subplot(2,4,n+4*(k-1))
        xlim([1 60])
        ylim([0,1])
        xlabel('days')
        ylabel('prevalence')
        title(titleStr{n},'fontweight','normal')
    end
    if n==1
        subplot(2,4,1)
        text(25,0.95,'index case','color',cmap(1,:))
        text(25,0.85,'household contact','color',cmap(2,:))
        text(25,0.75,'extrafamilial contact','color',cmap(3,:))
    end 
    
end

%% shedding index vs pTrans
serotype=2;

fecalOffset=2.8;

clear sheddingIndex pis pss idx
NAbRange=(logspace(0,11*log10(2),200));
age=18;
t=1:80;
t0=1;
sheddingIndex=nan(size(NAbRange,1),3);
pss=nan(length(NAbRange),3);
for k=1:length(NAbRange)
    sheddingIndex(k,1)=doseResponseModel(10^5.7,NAbRange(k),serotype)*(sum(sheddingViralLoadGMT(NAbRange(k),age,t,t0).*sheddingDurationCDF(NAbRange(k),t,t0)));
    sheddingIndex(k,2)=sqrt(sheddingIndex(k,1)*sheddingIndex(1,1));
    sheddingIndex(k,3)=sqrt(sheddingIndex(k,1)*sheddingIndex(1,1));
    
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',0 ...
                    ,'secondaryLog2NAb',log2(NAbRange(k)),'tertiaryLog2NAb',log2(NAbRange(k)),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset) );
    pss(k,1)=sum(d.tertiary.incidence)/sum(d.secondary.incidence);
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(NAbRange(k)) ...
                    ,'secondaryLog2NAb',0,'tertiaryLog2NAb',log2(NAbRange(k)),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset) );
    pss(k,2)=sum(d.tertiary.incidence)/sum(d.secondary.incidence);
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',0 ...
                    ,'secondaryLog2NAb',log2(NAbRange(k)),'tertiaryLog2NAb',0,'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset) );
    pss(k,3)=sum(d.tertiary.incidence)/sum(d.secondary.incidence);
end

cmap=cbrewer('qual','Set1',3);
% cmap2=[[1,0,0];[0,0.8,0];[0,0,1]];
cmap2=colormap(lines(5));
cmap2=cmap2([1,2,4,5],:);

figure(11);
clf; 
subplot(1,2,1)
hold all;
for k=1:3
    plot(sheddingIndex(:,1),pss(:,k)/pss(1,1),'-','color',cmap(k,:))
end
axis tight
set(gca,'xscale','log','yscale','log')
ylabel('relative R_{eff}')
xlabel('mean shedding index')
set(gca,'xtick',[1e2,1e3, 1e4,1e5,1e6])
text(0.5e4,10^(-.3),'both immunized','color',cmap(1,:),'fontsize',8)
text(0.5e4,10^(-.6),'index naive, contact immunized','color',cmap(2,:),'fontsize',8)
text(0.5e4,10^(-.9),'index immunized, contact naive','color',cmap(3,:),'fontsize',8)



%% figure 3: IPV boosting
cmap=cbrewer('qual','Set1',3);
cmap2=cbrewer('seq','YlOrRd',4);

fecalOffset=2.8; % india
siblingLog2NAb=[log2(512),log2(680),log2(2048)];
middle=[0,2,0];
pde=5.62;

% fecalOffset=0; % houston
% siblingLog2NAb=[log2(8),log2(8*8),log2(1028)];
% middle=[0,8,0];
% pde=1.16;

figure(14); clf;



for n=1:3
    
    experiment{n} = primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',siblingLog2NAb(n) ...
                    ,'secondaryLog2NAb',middle(n),'tertiaryLog2NAb',siblingLog2NAb(n) ...
                    ,'probDailyPrimarySecondaryContact',1,'perDoseEfficacy',pde,'serotype',4 ...
                    ,'secondaryContactAcquire',10^(-4.3299+fecalOffset) ...
                    ,'tertiaryContactAcquire',10^(-5.2943+fecalOffset) );
    
    subplot(2,3,n); hold all
    for m=1:3        
%         plot(experiment{n}.params.t,cumsum(experiment{n}.(typeStr{m}).incidence),'color',cmap(m,:),'linewidth',1.5)
        plot(experiment{n}.params.t,experiment{n}.(typeStr{m}).prevalence,'color',cmap(m,:),'linewidth',1.5)
    end
    subplot(2,3,n+3); hold all
    for m=1:3        
        plot(experiment{n}.params.t,cumsum(experiment{n}.(typeStr{m}).incidence),'color',cmap(m,:),'linewidth',1.5)
%         plot(experiment{n}.params.t,experiment{n}.(typeStr{m}).prevalence,'color',cmap(m,:),'linewidth',1.5)
    end
%     plot(experiment{n}.params.t,cumsum(experiment{n}.(typeStr{3}).incidence)./sum(experiment{n}.(typeStr{2}).incidence),'color',cmap(m,:),'linewidth',1.5)
end

ts={'6+ tOPV',{'6+ tOPV &','mOPV boost'},{'6+ tOPV &','IPV boost'}};
for n=1:3
        subplot(2,3,n)
        xlim([1 60])
        ylim([0,1])
        xlabel('days')
        ylabel('prevalence')
        title(ts{n},'fontweight','normal')
        
        subplot(2,3,n+3)
        xlim([1 60])
        ylim([0,1])
        xlabel('days')
        ylabel('incidence')
        title(ts{n},'fontweight','normal')
    if n==1
        subplot(2,3,1)
        text(25,0.95,'index case','color',cmap(1,:))
        text(25,0.85,'household contact','color',cmap(2,:))
        text(25,0.75,'extrafamilial contact','color',cmap(3,:))
    end 
    
end



%% color plots
clear R

serotype=4;

% Reff vs fecal contamination and family size
fecalOffset=linspace(-2,3.5,100);
groupSize=linspace(3,12,100);


% Sabin 2, naive
R{1}=nan(length(fecalOffset),length(groupSize));
for k=1:length(R{1})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',0,'serotype',serotype ...
                    ,'secondaryLog2NAb',0,'tertiaryLog2NAb',0,'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset(k)) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset(k)) );
    R{1}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end

% Sabin 2, plenty of tOPV
R{2}=nan(length(fecalOffset),length(groupSize));
for k=1:length(R{2})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(512),'serotype',serotype ...
                    ,'secondaryLog2NAb',log2(512),'tertiaryLog2NAb',log2(512),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset(k)) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset(k)) );
    R{2}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end


% Sabin 2, 5 years after
R{4}=nan(length(fecalOffset),length(groupSize));
for k=1:length(R{4})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(8),'serotype',serotype ...
                    ,'secondaryLog2NAb',log2(8),'tertiaryLog2NAb',log2(8),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset(k)) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset(k)) );
    R{4}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end


% Sabin 2, 1 year after
R{3}=nan(length(fecalOffset),length(groupSize));
for k=1:length(R{3})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(8),'serotype',serotype ...
                    ,'secondaryLog2NAb',log2(256),'tertiaryLog2NAb',log2(256),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(-5.2943+fecalOffset(k)) ...
                    ,'tertiaryContactAcquire',10^(-4.3299+fecalOffset(k)) );
    R{3}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end

Rmin=min(min(R{2}));
Rmax=max(max(R{1}));

Rmin=1e-4;
Rmax=15;


ca=linspace(log10(Rmin),log10(Rmax),256);
idx1=find(ca>=-1,1);
idx2=find(ca>=0,1);
c1=cbrewer('div','RdYlBu',(idx2-idx1)*2+13);
c1=c1(end:-1:end/2,:);
% c0=[interp1([1,idx1],[0.1,c1(1,1)],1:idx1)',interp1([1,idx1],[0.1,c1(1,2)],1:idx1)',interp1([1,idx1],[0.2,c1(1,3)],1:idx1)'];
c0=[interp1([1,idx1],[c1(1,1),c1(13,1)],1:idx1)',interp1([1,idx1],[c1(1,2),c1(13,2)],1:idx1)',interp1([1,idx1],[c1(1,3),c1(13,3)],1:idx1)'];
% c0=c0(end:-1:1,:);
c1=c1(13:end,:);
c2=cbrewer('seq','YlOrRd',256-idx2);
cmap=[c0;c1;c2];

figure(12+serotype); clf;
colormap(cmap)
subplot(1,4,4)
imagesc(-4.33+fecalOffset,groupSize,log10(R{1}'))
title('naive & IPV-only','fontweight','normal')

subplot(1,4,2)
imagesc(-4.33+fecalOffset,groupSize,log10(R{3}'))
title('mixed bOPV and tOPV','fontweight','normal')

subplot(1,4,3)
imagesc(-4.33+fecalOffset,groupSize,log10(R{4}'))
title('bOPV','fontweight','normal')

subplot(1,4,4)
imagesc(-4.33+fecalOffset,groupSize,log10(R{1}'))
title('naive & IPV-only','fontweight','normal')

subplot(1,4,1)
imagesc(-4.33+fecalOffset,groupSize,log10(R{2}'))
title('>3 tOPV','fontweight','normal')

for k=1:4
    subplot(1,4,k)
    axis xy
    xlabel('fecal-oral exposure (mg/day)')
    ylabel('intimate contact group size')
    caxis(log10([Rmin,Rmax]))
    set(gca,'xtick',[-5,-4,-3,-2,-1],'xticklabel',{'0.01','0.1','1','10','100'})
    if k==1
        c1=[1,1,1];
        c2=[0,0,0];
    else
        c1=[0,0,0];
        c2=[1,1,1];
    end
    text(-4.33,4,{'Houston','1960'},'HorizontalAlignment','center','fontsize',10,'color',c1);
%     text(-4.33+1,5,{'Matlab'},'HorizontalAlignment','center','fontsize',10,'color',c1);
%     text(-4.33+1.7,8,{'Mirpur'},'HorizontalAlignment','center','fontsize',10,'color',c2);
    text(-4.33+2.8,10,{'Uttar','Pradesh'},'HorizontalAlignment','center','fontsize',10,'color',c2);
end

h=colorbar;
g=get(h,'ytick');
set(h,'ytick',g(1:2:end),'yticklabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10'})
ylabel(h,'local R_{eff}')

