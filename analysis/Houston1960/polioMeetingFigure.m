% polioMeetingFigure

BMaxis=[0,7,14,21,28,35];
BMdata= [[0,0,0]; [0.88, 0.18, 0.02];... % digitized figure 4
         [0.84, 0.59, 0.1]; ...
         [0.74, 0.52, 0.1]; ...
         [0.64, 0.46, 0.17]; ...
         [0.46, 0.29, 0.13]];


% figure
hold on
cmap=colormap(lines(3));
plot(1:35,(experiment.primary.probInfectedGivenVaccine(1:35)),'color',cmap(1,:),'lineWidth',2)
plot(1:35,sum(experiment.secondary.probGenotypeShed(:,1:35),1),'lineWidth',2,'color',cmap(2,:))
plot(1:35,sum(experiment.tertiary.probGenotypeShed(:,1:35),1),'lineWidth',2,'color',cmap(3,:))
for k=1:3
    plot(BMaxis,BMdata(:,k),'--','linewidth',2,'color',cmap(k,:))
end

% need uncertainty

h=legend('index','sibling','contact');
set(h,'box','off','fontsize',7)
A=get(h,'pos');
A(1)=0.35;
A(2)=0.8;
set(h,'pos',A);

% plot([26, 32]+.2, [0.75 0.75]+.04,'k--','linewidth',2); text(33,0.79,'data','fontsize',7)
% plot([26, 32]+.2, 0.69*[1,1]+.04,'k-','linewidth',2); text(33,0.73,'fit','fontsize',7)
axis tight
ylim([0,1])
ylabel('probability shedding')
xlabel('days since challenge')
% print('-dpdf','Benyesh2WithFit.pdf')