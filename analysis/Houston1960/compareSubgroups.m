% script for various subgroup comparisons

addpath(genpath('..\helperFunctions'));
%% load data
[index,sibling,contact]=loadHouston1960();

%% index: compare under 7 with over 7
clc
for k=1:3
    [~,p]=fisherTestWrapper(sum(index(k).SheddingUnder7mo),sum(index(k).DenominatorUnder7mo),sum(index(k).Shedding7plus),sum(index(k).Denominator7plus))
end
% punchline: diffence is significant for types 2 and 3, but effect size is small

%% index: serotype comparison

clc
% total shedding
% type 1 vs type 2
[~,p]=fisherTestWrapper(sum(index(1).SheddingTotal),sum(index(1).DenominatorTotal),sum(index(2).SheddingTotal),sum(index(2).DenominatorTotal))
% type 2 vs type 3
[~,p]=fisherTestWrapper(sum(index(3).SheddingTotal),sum(index(3).DenominatorTotal),sum(index(2).SheddingTotal),sum(index(2).DenominatorTotal))
% type 1 vs type 3
[~,p]=fisherTestWrapper(sum(index(1).SheddingTotal),sum(index(1).DenominatorTotal),sum(index(3).SheddingTotal),sum(index(3).DenominatorTotal))

% first week
% type 1 vs type 2
[~,p]=fisherTestWrapper(index(1).SheddingTotal(1),index(1).DenominatorTotal(1),index(2).SheddingTotal(1),index(2).DenominatorTotal(1))
% type 2 vs type 3
[~,p]=fisherTestWrapper(index(3).SheddingTotal(1),index(3).DenominatorTotal(1),index(2).SheddingTotal(1),index(2).DenominatorTotal(1))
% type 1 vs type 3
[~,p]=fisherTestWrapper(index(1).SheddingTotal(1),index(1).DenominatorTotal(1),index(3).SheddingTotal(1),index(3).DenominatorTotal(1))

% punchline: type 2 shedding by index children is significantly higher than
% types 1 or 3. No significant differences between type 1 and type 3


%% sibling: compare under 5y with 5-9

clc
for k=1:3
    % compare 5+ to under 5
    [~,p]=fisherTestWrapper(sum(sibling(k).SheddingUnder5),sum(sibling(k).DenominatorUnder5),sum(sibling(k).Shedding60to107),sum(sibling(k).Denominator60to107))    
    
%     [~,p]=fisherTestWrapper(sum(sibling(k).SheddingUnder12mo),sum(sibling(k).DenominatorUnder12mo),sum(sibling(k).Shedding12to23),sum(sibling(k).Denominator12to23))    
%     [~,p]=fisherTestWrapper(sum(sibling(k).Shedding12to23),sum(sibling(k).Denominator12to23),sum(sibling(k).Shedding23to35),sum(sibling(k).Denominator23to35))
%     [~,p]=fisherTestWrapper(sum(sibling(k).Shedding23to35),sum(sibling(k).Denominator23to35),sum(sibling(k).Shedding36to59),sum(sibling(k).Denominator36to59))
%     [~,p]=fisherTestWrapper(sum(sibling(k).Shedding36to59),sum(sibling(k).Denominator36to59),sum(sibling(k).Shedding60to107),sum(sibling(k).Denominator60to107))
%     [~,p]=fisherTestWrapper(sum(sibling(k).Shedding60to107),sum(sibling(k).Denominator60to107),sum(sibling(k).Shedding108plus),sum(sibling(k).Denominator108plus))    
%     pause
end
% punchline: under 5s shed more than over 5s. All under-5's equivalent.


%% sibling: serotype comparison

clc
% total shedding
% type 1 vs type 2
[~,p]=fisherTestWrapper(sum(sibling(1).SheddingUnder5),sum(sibling(1).DenominatorUnder5),sum(sibling(2).SheddingUnder5),sum(sibling(2).DenominatorUnder5))
% type 2 vs type 3
[~,p]=fisherTestWrapper(sum(sibling(3).SheddingUnder5),sum(sibling(3).DenominatorUnder5),sum(sibling(2).SheddingUnder5),sum(sibling(2).DenominatorUnder5))
% type 1 vs type 3
[~,p]=fisherTestWrapper(sum(sibling(1).SheddingUnder5),sum(sibling(1).DenominatorUnder5),sum(sibling(3).SheddingUnder5),sum(sibling(3).DenominatorUnder5))

% punchline: type 2 shedding by sibling children is significantly higher than
% types 1 or 3. No significant differences between type 1 and type 3

%% contact: serotype comparison

clc
% total shedding
% type 1 vs type 2
[~,p]=fisherTestWrapper(sum(contact(1).SheddingUnder5),sum(contact(1).DenominatorUnder5),sum(contact(2).SheddingUnder5),sum(contact(2).DenominatorUnder5))
% type 2 vs type 3
[~,p]=fisherTestWrapper(sum(contact(3).SheddingUnder5),sum(contact(3).DenominatorUnder5),sum(contact(2).SheddingUnder5),sum(contact(2).DenominatorUnder5))
% type 1 vs type 3
[~,p]=fisherTestWrapper(sum(contact(1).SheddingUnder5),sum(contact(1).DenominatorUnder5),sum(contact(3).SheddingUnder5),sum(contact(3).DenominatorUnder5))

% punchline: type 2 shedding by contact children is significantly higher than
% types 1 or 3. No significant differences between type 1 and type 3



















