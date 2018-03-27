% explore and fit boosting after infection


close all; clear all;

addpath(genpath('..\helperFunctions'));
addpath(genpath('..\..\data\Louisiana1957'));

[data,hdrStr]=xlsread('householdSeroconversion.xlsx','prePostAntibodyLineList');
hdrStr=hdrStr(1,1:2);
data=log2(data);
hdrStr={'preExposureLog2NAb','postExposureLog2NAb'};

%log diff
hdrStr{3}='log2Diff';
data(:,3)=data(:,2)-data(:,1);

df=mat2dataset(data,'varnames',hdrStr);

%% explore non-infected data
seroconversionThreshold=4; % 
idx=df.log2Diff>=log2(seroconversionThreshold);

LM=fitlm(df(idx,:),'log2Diff~preExposureLog2NAb')

figure(1); clf; hold all;

jitter=0.1*randn(length(df.log2Diff(idx)),2);

subplot(1,2,1); hold all;
h=LM.plot();
set(h(1),'visible','off')
plot(df.preExposureLog2NAb(idx)+jitter(:,1),df.log2Diff(idx)+jitter(:,2),'x')

subplot(1,2,2); hold all;
test=mat2dataset([0:14]','varnames','preExposureLog2NAb');
[test.theta,test.CI]=LM.predict(test);
test.theta=2.^test.theta;
test.CI=2.^test.CI;
plotWithCI(test.preExposureLog2NAb,test.theta,test.CI,'r')
ylabel('theta')
xlabel('preInfection log2NAb')
title({'WPV mean response model','from Gelfand 1957b'})



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