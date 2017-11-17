function PStillInfected=sheddingDurationCDF(NAb,t,t0,b)
if nargin<1
    NAb=1;
end
if nargin<2
    t=0:80;
end
if nargin<3
    t0=0;
end
if nargin<4
%     b=[log(33),log(1.2205),log(1.8965)]; % Behrend unpublished, but
%     actually is derived in supplemental code to Behrend et al 2014.
    b=[log(30.3),log(1.164),log(1.86)];
end
M=b(1)-b(2)*log2(NAb);
sigma=b(3);

PStillInfected = 1/2*(1-erf(1/sqrt(2*sigma^2)*(log(abs(t-t0))-M))).*double(t>=t0);
PStillInfected(t<(t0+2))=0;

end