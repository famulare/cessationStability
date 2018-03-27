% explore and fit boosting after infection


close all; clear all;

addpath(genpath('..\helperFunctions'));
addpath(genpath('..\..\data\Louisiana1957'));

[data,hdrStr]=xlsread('householdSeroconversion.xlsx','prePostAntibodyLineList');
hdrStr=hdrStr(1,1:2);


%log diff
hdrStr{3}='postPreRatio';
data(:,3)=data(:,2)./data(:,1);

df=mat2dataset(data,'varnames',hdrStr);

df.log2PreExposureTiter = log2(df.preExposureTiter);
df.log2PostExposureTiter = log2(df.postExposureTiter);
df.logePostPreRatio = log(df.postPreRatio);

%% explore non-infected data
seroconversionThreshold=4; % 
idx=df.postPreRatio>=seroconversionThreshold;

LM=fitlm(df(idx,:),'logePostPreRatio ~ log2PreExposureTiter')

figure(1); clf; hold all;

jitter=0.1*randn(length(df.postPreRatio(idx)),2);

subplot(1,2,1); hold all;
h=LM.plot();
set(h(1),'visible','off')
plot(df.log2PreExposureTiter(idx)+jitter(:,1),df.logePostPreRatio(idx)+jitter(:,2),'x')

subplot(1,2,2); hold all;
test=mat2dataset([0:15]','varnames','log2PreExposureTiter');
[test.theta,test.CI]=LM.predict(test);
test.theta=exp(test.theta);
test.CI=exp(test.CI);
plotWithCI(test.log2PreExposureTiter,test.theta,test.CI,'r')
ylabel('theta')
xlabel('log2PreExposureTiter')
title({'WPV mean response model','from Gelfand 1957b'})

% jafari
plot(test.log2PreExposureTiter,exp(4.28-.2*test.log2PreExposureTiter),'k')

% inset
axes('position',[.7,.5,.25,.3]); hold on
plotWithCI(test.log2PreExposureTiter,test.theta,test.CI,'r')
plot(test.log2PreExposureTiter,exp(4.28-.2*test.log2PreExposureTiter),'k')
ylabel('theta')
xlabel('log2PreExposureTiter')
xlim([6,14])

%% is log titer for non-coverters roughly gaussian?
figure(2); hold all;
[y,x]=hist(data(~idx,3),-2:2);
plot(x,y/sum(y),'k');
% for k=unique(data(~idx,1))'
%     [y,x]=hist(data(~idx & data(:,1)==k,3),-2:2);
%     plot(x,y/sum(y));
% end
xlabel('log2(post - pre)')
ylabel('frequency')
title('non-seroconverters','fontweight','normal')

% sure...