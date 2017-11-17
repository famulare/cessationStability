% script to visualize and provide basic analysis of Houston shedding data


addpath(genpath('..\helperFunctions'));
%% load data
[index,sibling,contact]=loadHouston1960();

%% plot shedding by subject type and age
% columns: serotype
% rows: index, sibling, primary contact
% index age bins: under 7mos, 7+
% Sibling/contact age bins: under 5, 5-9, 10+
cmap=colormap(lines(6));

% supplemental figure (data as reported)
figure(1); clf;
for k=1:3
    % index
    subplot(3,3,k)
    hold all
    [p,ci]=binofit(index(k).SheddingUnder7mo,index(k).DenominatorUnder7mo);
    plotWithCI(index(k).Day+0.25,p,ci,cmap(1,:));
    [p,ci]=binofit(index(k).Shedding7plus,index(k).Denominator7plus);
    plotWithCI(index(k).Day-0.25,p,ci,cmap(2,:));
    
    xlim([0 37])
    ylim([0,1])
    ylabel('fraction shedding')
    xlabel('days')
    set(gca,'xtick',[0:7:35]);
    title('index','fontweight','normal')
    if k==3
        text(1,0.25,'under 7 months','color',cmap(1,:))
        text(1,0.1,'7 and over','color',cmap(2,:))        
    end
    
    % sibling
    subplot(3,3,k+3)
    hold all
    [p,ci]=binofit(sibling(k).SheddingUnder12mo,sibling(k).DenominatorUnder12mo);
    plotWithCI(sibling(k).Day+0.5,p,ci,cmap(1,:));
    [p,ci]=binofit(sibling(k).Shedding12to23,sibling(k).Denominator12to23);
    plotWithCI(sibling(k).Day+0.3,p,ci,cmap(2,:));
    [p,ci]=binofit(sibling(k).Shedding23to35,sibling(k).Denominator23to35);
    plotWithCI(sibling(k).Day+0.1,p,ci,cmap(3,:));
    [p,ci]=binofit(sibling(k).Shedding36to59,sibling(k).Denominator36to59);
    plotWithCI(sibling(k).Day-0.1,p,ci,cmap(4,:));
    [p,ci]=binofit(sibling(k).Shedding60to107,sibling(k).Denominator60to107);
    plotWithCI(sibling(k).Day-0.3,p,ci,cmap(5,:));
    [p,ci]=binofit(sibling(k).Shedding108plus,sibling(k).Denominator108plus);
    plotWithCI(sibling(k).Day-0.5,p,ci,cmap(6,:));
    xlim([0 37])
    ylim([0,1])
    ylabel('fraction shedding')
    xlabel('days')
    set(gca,'xtick',[0:7:35]);
    title('sibling','fontweight','normal')
    if k==3
        text(1,0.85,'under 12 months','color',cmap(1,:))
        text(1,0.7,'12 to 23','color',cmap(2,:))        
        text(1,0.55,'24 to 25','color',cmap(3,:))        
        text(1,0.4,'36 to 59','color',cmap(4,:))        
        text(1,0.25,'60 to 107','color',cmap(5,:))        
        text(1,0.1,'108+','color',cmap(6,:))        
    end
    
    subplot(3,3,k+6)
    hold all
    % aggregated as reported
    [p,ci]=binofit(contact(k).SheddingTotal,contact(k).DenominatorTotal);
    plotWithCI(sibling(k).Day,p,ci,cmap(1,:));
    
    xlim([0 37])
    ylim([0,0.5])
    ylabel('fraction shedding')
    xlabel('days')
    set(gca,'xtick',[0:7:35]);
    title('extrafamilial contact','fontweight','normal')
    if k==3
        text(1,0.35,'all ages','color',cmap(1,:))
    end
         
end

%% main figure (active analysis groups only)
figure(2); %clf;

for k=1:3
    % index
    subplot(3,3,k)
    hold all
    [p,ci]=binofit(index(k).SheddingTotal,index(k).DenominatorTotal);
    plotWithCI(index(k).Day,p,ci,cmap(1,:),'-',10,false,false);
    
    xlim([0 37])
    ylim([0,1])
    ylabel('fraction shedding')
    xlabel('days')
    set(gca,'xtick',[0:7:35]);
    title('index','fontweight','normal')
    if k==3
        text(1,0.1,'age 2 to 18 months','color',cmap(1,:))
    end
    
    % sibling
    subplot(3,3,k+3)
    hold all
    [p,ci]=binofit(sibling(k).SheddingUnder5,sibling(k).DenominatorUnder5);
    plotWithCI(sibling(k).Day+0.25,p,ci,cmap(1,:),'-',10,false,false);
    [p,ci]=binofit(sibling(k).Shedding60to107,sibling(k).Denominator60to107);
    plotWithCI(sibling(k).Day-0.25,p,ci,cmap(2,:),'-',10,false,false);
    xlim([0 37])
    ylim([0,1])
    ylabel('fraction shedding')
    xlabel('days')
    set(gca,'xtick',[0:7:35]);
    title('sibling','fontweight','normal')
    if k==3
        text(1,0.85,'age <5 years','color',cmap(1,:))
        text(1,0.7,'5 to 9 years','color',cmap(2,:))        
    end
    
    % contacts
        % age correction factor for contacts
        % contacts have similar demographics to siblings, but demographic
        % breakdown is not presented in Benyesh-Melnick1967 paper
         estimatedDenominatorFractionUnder5 = (max([sibling.DenominatorUnder5])./max([sibling.DenominatorTotal]));
         estimatedSheddingFractionUnder5 = (max([sibling.SheddingUnder5])./max([sibling.SheddingTotal]));
    
    subplot(3,3,k+6)
    hold all
    % adjusted under 5y
    [p,ci]=binofit(contact(k).SheddingUnder5,contact(k).DenominatorUnder5);
    plotWithCI(contact(k).Day+0.25,p,ci,cmap(1,:),'-',10,false,false);
    % adjusted over 5
    [p,ci]=binofit(contact(k).Shedding5plus,contact(k).Denominator5plus);
    plotWithCI(contact(k).Day-0.25,p,ci,cmap(2,:),'-',10,false,false);
    xlim([0 37])
    ylim([0,0.5])
    ylabel('fraction shedding')
    xlabel('days')
    set(gca,'xtick',[0:7:35]);
    title('extrafamilial contact','fontweight','normal')
    if k==3
        text(1,0.375,'age <5 years (est.)','color',cmap(1,:))
        text(1,0.3,'5 to 9 years (est.)','color',cmap(2,:))        
    end
         
end
