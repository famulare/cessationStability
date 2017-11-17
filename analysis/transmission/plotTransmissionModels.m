% plot transmission models


addpath(genpath('..\Houston1960'))
addpath(genpath('..\helperFunctions'))
addpath(genpath('..\doseResponse'))
addpath(genpath('..\shedding'))
[index,sibling,contact]=loadHouston1960();


t=1:100;


% surface colormap
    Rmin=1e-4;
    Rmax=15;

    ca=linspace(log10(Rmin),log10(Rmax),256);
    idx1=find(ca>=-1,1);
    idx2=find(ca>=0,1);
    c1=cbrewer('div','RdYlBu',(idx2-idx1)*2+13);
    c1=c1(end:-1:end/2,:);
    c0=[interp1([1,idx1],[c1(1,1),c1(13,1)],1:idx1)',interp1([1,idx1],[c1(1,2),c1(13,2)],1:idx1)',interp1([1,idx1],[c1(1,3),c1(13,3)],1:idx1)'];
    c1=c1(13:end,:);
    c2=cbrewer('seq','YlOrRd',256-idx2);
    cmap=[c0;c1;c2];

% line colormaps
    cmap3=cbrewer('qual','Set1',3);
    cmap2=cbrewer('seq','YlOrRd',4);
    %% dose and immunity
clear R

serotype=4;

% Reff vs fecal contamination and family size
fecalOral=linspace(-7.5,-2.7,100);
log2NAb=linspace(0,11,100);
groupSize=12;

clear R

for k=1:length(fecalOral)
    for n=1:length(log2NAb)

        d = primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',log2NAb(n) ...
                        ,'secondaryLog2NAb',log2NAb(n),'tertiaryLog2NAb',log2NAb(n) ...
                        ,'secondaryContactAcquire',10^fecalOral(k) ...
                        ,'tertiaryContactAcquire',10^(fecalOral(k)) ...
                        ,'perDoseEfficacy',1,'serotype',serotype);
        R{1}(n,k)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;

    end
end



figure(2+serotype); clf;
colormap(cmap)

imagesc(log2NAb,fecalOral,log10(R{1}'))
title('dose vs immunity','fontweight','normal','fontsize',9)

    axis xy
    ylabel({'fecal-oral dose [\mug/day]'},'fontsize',10)
    xlabel('OPV-equivalent antibody titer','fontsize',10)
    caxis(log10([Rmin,Rmax]))
    set(gca,'ytick',[-7:1:-1],'yticklabel',{'10^{-1}','','10^1','','10^3',''})
    set(gca,'fontsize',9)
    xlim([0,11])
    set(gcf,'Units','in','Pos',[1,1,2,1.5],'papersize',[7.5,2],'paperunits','in','paperPosition',[0,0,7.5,2])
%     export_fig(['fecalVsImmunity_',num2str(serotype),'.pdf'],'-pdf','-transparent')

    
%% color plots
clear R

serotype=3;

% Reff vs fecal contamination and family size
fecalOral=linspace(-7.5,-2.7,100);
groupSize=linspace(3,11,100);


% Sabin 2, naive
R{1}=nan(length(fecalOral),length(groupSize));
for k=1:length(R{1})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',0,'serotype',serotype ...
                    ,'secondaryLog2NAb',0,'tertiaryLog2NAb',0,'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(fecalOral(k)) ...
                    ,'tertiaryContactAcquire',10^(fecalOral(k)) );
    R{1}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end

% define low transmission as R{1}<1;
lowMicrogams=10^(6+fecalOral(find(R{1}(:,1)<1,1,'last')))

% Sabin 2, plenty of tOPV
R{2}=nan(length(fecalOral),length(groupSize));
for k=1:length(R{2})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(512),'serotype',serotype ...
                    ,'secondaryLog2NAb',log2(512),'tertiaryLog2NAb',log2(512),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(fecalOral(k)) ...
                    ,'tertiaryContactAcquire',10^(fecalOral(k)) );
    R{2}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end


% Sabin 2, 5 years after
R{4}=nan(length(fecalOral),length(groupSize));
for k=1:length(R{4})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(8),'serotype',serotype ...
                    ,'secondaryLog2NAb',log2(8),'tertiaryLog2NAb',log2(8),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(fecalOral(k)) ...
                    ,'tertiaryContactAcquire',10^(fecalOral(k)) );
    R{4}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end

% define high transmission as R{4}>1;
highMicrogams=(10^(6+fecalOral(find(R{4}(:,1)>=1,1,'first'))))


% Sabin 2, 1 year after
R{3}=nan(length(fecalOral),length(groupSize));
for k=1:length(R{3})
    d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2(8),'serotype',serotype ...
                    ,'secondaryLog2NAb',log2(256),'tertiaryLog2NAb',log2(256),'runNetwork',false,'perDoseEfficacy',1 ...
                    ,'secondaryContactAcquire',10^(fecalOral(k)) ...
                    ,'tertiaryContactAcquire',10^(fecalOral(k)) );
    R{3}(k,:)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize;
end


figure(12+serotype); clf;
colormap(cmap)

subplot(1,4,2)
imagesc(fecalOral,groupSize,log10(R{3}'))
title('mixed bOPV and tOPV','fontweight','normal','fontsize',9)

subplot(1,4,3)
imagesc(fecalOral,groupSize,log10(R{4}'))
title('bOPV era','fontweight','normal','fontsize',9)

subplot(1,4,4)
imagesc(fecalOral,groupSize,log10(R{1}'))
title('naive','fontweight','normal','fontsize',9)

subplot(1,4,1)
imagesc(fecalOral,groupSize,log10(R{2}'))
title('tOPV era','fontweight','normal','fontsize',9)

for k=1:4
    subplot(1,4,k)
    axis xy
    xlabel({'fecal-oral exposure','[\mug/day]'},'fontsize',10)
    ylabel('contact network size','fontsize',10)
    caxis(log10([Rmin,Rmax]))
    set(gca,'xtick',[-7:1:-1],'xticklabel',{'10^{-1}','','10^1','','10^3',''})
    set(gca,'fontsize',9)
    if k<3
        c1=[1,1,1];
        c2=[0,0,0];
    else
        c1=[0,0,0];
        c2=[1,1,1];
    end
    text(-5.3,4,{'Houston','1960'},'HorizontalAlignment','center','fontsize',6,'color',c1);
    text(-5.3+2.8,10,{'UP&Bihar','2003-8'},'HorizontalAlignment','center','fontsize',6,'color',c2);
    AX{k}=get(gca,'pos');
end

for k=1:4
    subplot(1,4,k)
    AX{k}(2)=AX{k}(2)+0.1;
    AX{k}(3:4)=AX{k}(3:4).*[0.9,0.8];
    set(gca,'pos',AX{k})
end

h=colorbar;
set(h,'ytick',[-4:1],'yticklabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10'},'fontsize',8)
% ylabel(h,'local R_{eff}')
set(h,'pos',[0.91,AX{4}(2),0.01,AX{4}(4)])

set(gcf,'Units','in','Pos',[1,1,7.5,1.5],'papersize',[7.5,2],'paperunits','in','paperPosition',[0,0,7.5,2])
export_fig(['localReffSweep_type',num2str(serotype),'.pdf'],'-pdf','-transparent')

%%
typeStr={'primary','secondary','tertiary'};
infantLog2NAb=[log2(512),log2(8),log2(8),0];
siblingLog2NAb=[log2(512),log2(256),log2(2),0];

titleStr={'tOPVera','bOPV & tOPV era','bOPV era','naive'};

% figure 1 

fecalOral=[0.2,6,230]*1e-6;
clear YCI

serotype=2;
reps=100;
figure(10); clf;
for k=1:3
    for n=1:4

        experiment{k,n} = primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',infantLog2NAb(n) ...
                        ,'secondaryLog2NAb',siblingLog2NAb(n),'tertiaryLog2NAb',siblingLog2NAb(n) ...
                        ,'secondaryContactAcquire',fecalOral(k) ...
                        ,'tertiaryContactAcquire',fecalOral(k) ...
                        ,'numDailySecondaryTertiaryContact',8.96 ...
                        ,'perDoseEfficacy',1,'serotype',serotype);
        
        modelBoot=nan(length(experiment{k,n}.params.t),3,reps);
        for q=1:reps;
            tmp=primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',infantLog2NAb(n) ...
                        ,'secondaryLog2NAb',siblingLog2NAb(n),'tertiaryLog2NAb',siblingLog2NAb(n) ...
                        ,'secondaryContactAcquire',fecalOral(k) ...
                        ,'tertiaryContactAcquire',fecalOral(k) ...
                        ,'numDailySecondaryTertiaryContact',8.96 ...
                        ,'perDoseEfficacy',1,'serotype',serotype,'confidenceIntervalSamplerSeed',q);  
            for m=1:3
                modelBoot(:,m,q)=tmp.(typeStr{m}).prevalence;
            end
        end

        subplot(3,4,4*(k-1)+n); hold all
        for m=1:3        
            YCI{m}=prctile(squeeze(modelBoot(:,m,:))',[2.5,97.5])';
%             plot(experiment{k,n}.params.t,experiment{k,n}.(typeStr{m}).prevalence,'color',cmap3(m,:),'linewidth',1.5)
            plotWithCI(experiment{k,n}.params.t,experiment{k,n}.(typeStr{m}).prevalence,YCI{m},cmap3(m,:))
        end
    end
end

for n=1:4
    for k=1:3
        subplot(3,4,4*(k-1)+n);
        xlim([1 60])
        ylim([0,1])
        xlabel('days')
        ylabel({'probability of shedding'})
        title(titleStr{n},'fontweight','normal')
    end
    if n==1
        subplot(3,4,1)
        text(25,0.95,'index case','color',cmap3(1,:))
        text(25,0.85,'household contact','color',cmap3(2,:))
        text(25,0.75,'extrafamilial contact','color',cmap3(3,:))
    end 
    
end

set(gcf,'Units','in','Pos',[1,1,7.5,4.2],'papersize',[7.5,5],'paperunits','in','paperPosition',[0,0,7.5,5])
export_fig(['scenarios_type',num2str(serotype),'.pdf'],'-pdf','-transparent')

%% HID50

clear R AX

serotype=2;

% Reff vs HID50 and immunity
log2NAb=linspace(0,11,100);
hid=linspace(5,75,100);
beta=hid;
for k=1:length(beta)
    beta(k)=fzero(@(x) 0.5-doseResponseModel(hid(k),1,1,[x^2,0.444,0.545]), 2)^2;
end

% fecalOral=[0.5,6,230]*1e-6;
groupSize=[4,5,10];

figure(22+serotype); clf;
colormap(cmap)
for m=1:3
    R{m}=nan(length(log2NAb),length(hid));
    for k=1:size(R{m},1)
        for n=1:size(R{m},2)
            % shedding duration depends on reversion
            % small effect relative to HID50
%             if beta(n)<7.9611
%                 stype=4;
%             else
                stype=2;
%             end
            d=primarySecondaryTertiaryDoseModel('t',t,'primaryLog2NAb',log2NAb(k),'serotype',serotype ...
                        ,'secondaryLog2NAb',log2NAb(k),'tertiaryLog2NAb',log2NAb(k),'runNetwork',false,'perDoseEfficacy',1 ...
                        ,'secondaryContactAcquire',fecalOral(m) ...
                        ,'tertiaryContactAcquire',fecalOral(m) ...
                        ,'numDailySecondaryTertiaryContact',8.96 ...
                        ,'doseResponseBeta',beta(n),'serotype',stype);
            R{m}(k,n)=groupSize(m)*sum(d.tertiary.incidence)/sum(d.primary.incidence);
        end
    end
    subplot(1,3,m)
    imagesc(log2NAb,hid,log10(R{m}'))
    drawnow;
end

for k=1:3
    subplot(1,3,k)
    axis xy
    ylim([min(hid),max(hid)])
    xlabel({'OPV-equivalent','antibody titer'},'fontsize',10)
    ylabel('HID50','fontsize',10)
    caxis(log10([Rmin,Rmax]))
    set(gca,'xtick',[0,3,6,9,11],'xticklabel',{'1','8','64','512','2048'})
    set(gca,'ytick',(5:5:75),'yticklabel',{'','10','','','','30','','','','50','','','','70'})
    set(gca,'fontsize',9)
    AX{k}=get(gca,'pos');
    AX{k}(2)=AX{k}(2)+0.1;
    AX{k}(3:4)=AX{k}(3:4).*[0.9,0.8];
    set(gca,'pos',AX{k})
end

h=colorbar;
set(h,'ytick',[-4:1],'yticklabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10'},'fontsize',8)
% ylabel(h,'local R_{eff}')
set(h,'pos',[0.91,AX{3}(2),0.01,AX{3}(4)])

rs=7.5*.72*3/4+7.5*0.28;
set(gcf,'Units','in','Pos',[1,1,rs,1.75],'papersize',[7.5,2],'paperunits','in','paperPosition',[0,0,rs,2])
% export_fig('HID50sweep.pdf','-pdf','-transparent')

%% directionality of immunization

clear R AX

serotype=2;

fecalOral=[5,230]*1e-6;
groupSize=[40,100]; 
log2NAb=linspace(0,11,100);

figure(36+serotype); clf;
colormap(cmap)
for k=1:2
    R{k}=nan(length(log2NAb));
    for m=1:length(log2NAb)
        for n=1:length(log2NAb) 
            d=primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',log2NAb(m),'serotype',serotype ...
                        ,'secondaryLog2NAb',log2NAb(m),'tertiaryLog2NAb',log2NAb(n),'runNetwork',false,'perDoseEfficacy',1 ...
                        ,'secondaryContactAcquire',fecalOral(k) ...
                        ,'tertiaryContactAcquire',fecalOral(k)/10 ...
                        ,'numDailyPrimarySecondaryContact',1 ... % modify daily contact rates to lower per-day probability
                        ,'numDailySecondaryTertiaryContact',1/30);
            R{k}(m,n)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize(k);
        end
    end
    subplot(1,3,k)
    imagesc(log2NAb,log2NAb,log10(R{k}'))
    drawnow;
end

max(max(R{2}))

for k=1:2
    subplot(1,3,k)
    axis xy
    xlabel('index immunity')
    ylabel('contact immunity')
    caxis(log10([Rmin,Rmax]))
    set(gca,'xtick',[0,3,6,9,11],'xticklabel',{'1','8','64','512','2048'})
    set(gca,'ytick',[0,3,6,9,11],'yticklabel',{'1','8','64','512','2048'})
    set(gca,'fontsize',9)
    AX{k}=get(gca,'pos');
    AX{k}(2)=AX{k}(2)+0.1;
    AX{k}(3:4)=AX{k}(3:4).*[0.9,0.8];
    set(gca,'pos',AX{k})
end

h=colorbar;
set(h,'ytick',[-4:1],'yticklabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10'},'fontsize',8)
% ylabel(h,'local R_{eff}')
set(h,'pos',[0.91,AX{2}(2),0.01,AX{2}(4)])

rs=7.5*.72*3/4+7.5*0.28;
set(gcf,'Units','in','Pos',[1,1,rs,1.75],'papersize',[7.5,2],'paperunits','in','paperPosition',[0,0,rs,2])
export_fig('directionality3.pdf','-pdf','-transparent')

%% dose response

clear R

serotype=2;

log2NAb=linspace(0,11,101);
dose=linspace(-1,6,100);

figure(42+serotype); clf;
cmap4=cbrewer('div','RdYlBu',100);
cmap4=cmap4(end:-1:1,:);
% cmap4=cmap4(:,[2,1,3])

% cmap4=cbrewer('seq','Reds',27);
% cmap4=cmap4(end:-1:1,:);
% for k=28:36
%     cmap4(k,:)=cmap4(27,:);
% end
% cmap4=(cmap4+cbrewer('seq','Blues',36))/2;
% tmp=cbrewer('seq','Greens',100);
% cmap4(37:100,:)=tmp(37:100,:);
% cmap4(1,:)=mean(cmap4(1,:))*[1,1,1];

colormap(cmap4)

R=nan(length(log2NAb),length(dose));
for k=1:size(R,1)
    for n=1:size(R,2)
         R(k,n)=doseResponseModel(10^dose(n),2^log2NAb(k),serotype);
    end
end
imagesc(log2NAb,dose,R')

axis xy
ylim([min(dose),max(dose)])
xlabel({'OPV-equivalent antibody titer'},'fontsize',10)
ylabel('dose [TCID50]','fontsize',10)
caxis([0,1])
set(gca,'xtick',[0,3,5,8,11],'xticklabel',{'1','8','32','256','2048'})
set(gca,'ytick',(-1:1:6),'yticklabel',{'10^{-1}','','10^1','','10^3','','10^5',''})
set(gca,'fontsize',9)
AX=get(gca,'pos');
AX(2)=AX(2)+0.1;
AX(3:4)=AX(3:4).*[0.9,0.8];
set(gca,'pos',AX)

h=colorbar;
set(h,'ytick',[0:0.2:1],'fontsize',8)
set(h,'pos',[0.91,AX(2),0.03,AX(4)])
ylabel(h,'probability of shedding')

set(gcf,'Units','in','Pos',[1,1,3.5,3],'papersize',[3.5,3],'paperunits','in','paperPosition',[0,0,3.5,3])
export_fig('doseResponse.pdf','-pdf','-transparent')



%% directionality of immunization

clear R AX

serotype=2;

fecalOral=logspace(-3,0,100);
contactRate=logspace(-3,0,100);

% groupSize=linspace(3,11,100);

f=[5,230]*1e-6;
groupSize=[4,10]; 
log2NAb=0;

figure(36+serotype); clf;
colormap(cmap)
for k=1:2
    R{k}=nan(length(fecalOral));
    for m=1:length(fecalOral)
        for n=1:length(contactRate) 
            d=primarySecondaryTertiaryDoseModel('t',1:80,'primaryLog2NAb',log2NAb,'serotype',serotype ...
                        ,'secondaryLog2NAb',log2NAb,'tertiaryLog2NAb',log2NAb,'runNetwork',false,'perDoseEfficacy',1 ...
                        ,'secondaryContactAcquire',f(k) ...
                        ,'tertiaryContactAcquire',fecalOral(m)*f(k) ...
                        ,'numDailyPrimarySecondaryContact',1 ... % modify daily contact rates to lower per-day probability
                        ,'numDailySecondaryTertiaryContact',8.9685*contactRate(n));
            R{k}(n,m)=sum(d.tertiary.incidence)/sum(d.primary.incidence)*groupSize(k);
        end
    end
    subplot(1,3,k)
    imagesc(log10(contactRate*8.9685),log10(fecalOral*f(k)),log10(R{k}'))
    drawnow;
end

max(max(R{2}))

for k=1:2
    subplot(1,3,k)
    axis xy
    ylabel('fecal-oral dose')
    xlabel('extrafamilial contact rate')
    caxis(log10([Rmin,Rmax]))
%     set(gca,'ytick',[-3:1:0],'yticklabel',{'10^{-3}','10^{-2}','10^{-1}','1'})          
    set(gca,'xtick',[-2:1:1],'xticklabel',{'10^{-2}','10^{-1}','1','10'})          
    set(gca,'fontsize',9)
    AX{k}=get(gca,'pos');
    AX{k}(2)=AX{k}(2)+0.1;
    AX{k}(3:4)=AX{k}(3:4).*[0.9,0.8];
    set(gca,'pos',AX{k})
end

h=colorbar;
set(h,'ytick',[-4:1],'yticklabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','1','10'},'fontsize',8)
% ylabel(h,'local R_{eff}')
set(h,'pos',[0.91,AX{2}(2),0.01,AX{2}(4)])

rs=7.5*.72*3/4+7.5*0.28;
set(gcf,'Units','in','Pos',[1,1,rs,1.75],'papersize',[7.5,2],'paperunits','in','paperPosition',[0,0,rs,2])
export_fig('longDistance.pdf','-pdf','-transparent')
