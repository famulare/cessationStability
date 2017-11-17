function shedTiter=sheddingViralLoadGMT(NAb,ageMonths,t,t0,b)
if nargin<1
    NAb=1;
end
if nargin<2
    ageMonths=18;
end
if nargin<3
    t=0:80;
end
if nargin<4
    t0=0;
end
if nargin<5
    b=[1.65,0.17,0.056,0.32];
end


peakShed = 10.^(peakSheddingAgeMultiplier(ageMonths)*(1-b(3)*log2(NAb)));

% lognormal temporal profile in TCID50. Not a probability distrubition
shedTiter= exp(b(1)-b(2)^2/2)./(t-t0+0.01).*exp(-(log(t-t0+0.01)-b(1)).^2./(2.*(b(2)+b(4)*log(t-t0+0.01)).^2));
shedTiter=peakShed*shedTiter;

% floor
% appears to be a floor at tOPVx3. This is an assay artifact in that shedding less that this is a negative.
shedTiter=max(shedTiter,10^(2.6)); 
shedTiter(t<(t0+2))=0;

end