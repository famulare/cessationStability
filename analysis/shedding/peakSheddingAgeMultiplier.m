function peakSheddingAgeFactor = peakSheddingAgeMultiplier(ageMonths)

b=[6.67,4.29,9.92];
% one-factor model:  estimated timescale matches with move to solid food and immunity maturation
peakSheddingAgeFactor =min(b(1),((b(1)-b(2))*exp(-(ageMonths-7)/b(3))+b(2)));

end