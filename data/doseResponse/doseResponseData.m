function data = doseResponseData()
% summary of studies that investigate OPV dose response. If you're reading
% this, I (MGF) apologize for coding data into matlab instead of a
% spreadsheet. 

addpath(genpath('..\helperFunctions'));
addpath(genpath('..\..\data\doseResponse'));

% TYPE 1 challenge %%%%%%%%%%%%%%%%%%%%
% Minor1981, OPV1
% Minor, T E, Allen, C I, Tsiatis, A A, Nelson, D B, D'Alessio, D J
% Human infective dose determinations for oral poliovirus type 1 vaccine in infants.
% Journal of clinical microbiology 13(2), 1981.
%   data from Table 1, combine across pre-challenge serum antibody
% excluded from study because sample size per dose is too small.
data.Minor1981.serotype=1;
data.Minor1981.ageMos=2;
data.Minor1981.vaccSchedule='first dose challenge';
% reported data
    data.Minor1981.dose=[280,210,160,90,80,65,55,50,42,27,16,7];
    data.Minor1981.pos=[1,2,3,3,1,0,1,3,0,0,0,0];
    data.Minor1981.N=[1,2,3,3,1,6,3,3,1,2,2,1];
% smoothed in bins of roughly equal log width
    binEdgeIdx=[0,cumsum(hist(log(data.Minor1981.dose),4))];
    for k=2:length(binEdgeIdx)
        tmpIdx=(binEdgeIdx(k-1)+1):binEdgeIdx(k);
        data.Minor1981.smooth.dose(k-1)=geomean(data.Minor1981.dose(tmpIdx));
        data.Minor1981.smooth.pos(k-1)=sum(data.Minor1981.pos(tmpIdx));
        data.Minor1981.smooth.N(k-1)=sum(data.Minor1981.N(tmpIdx));
    end

% onorato1991, OPV1
% Ida M. Onorato, John F. Modlin, A. Marshall McBean, Mary Lou Thoms, Genevieve A. Losonsky, Roger H. Bernier
% Mucosal Immunity Induced by Enhanced-Potency Inactivated and Oral Polio Vaccines 
% Journal of Infectious Diseases. 163(1), 1991.
%   data from table 1    
data.Onorato1991.serotype=1;
data.Onorato1991.dose=[5.8e5,6.5e2]; % median challenge dose

data.Onorato1991.IPVx3.ageMos=26.9;
data.Onorato1991.IPVx3.vaccSchedule='eIPV 2,4,18 months';
data.Onorato1991.IPVx3.pos=[37,22];
data.Onorato1991.IPVx3.N=[45, 48];

data.Onorato1991.tOPVx3.ageMos=23.6;
data.Onorato1991.tOPVx3.vaccSchedule='tOPV 2,4,18 months';
data.Onorato1991.tOPVx3.pos=[14,6];
data.Onorato1991.tOPVx3.N=[45, 34];

% Henry1966
% Henry, J L, Jaikaran, E S, Davies, J R, Tomlinson, A J, Mason, P J, Barnes, J M, Beale, A J
% A study of poliovaccination in infancy: excretion following challenge with
% live virus by children given killed or living poliovaccine.
% The Journal of Hygiene, 64(1) 1966. 
%   data from Table 6

data.Henry1966.serotype=1;
data.Henry1966.dose=[10^5.7,10^4.7,10^3.7,10^2.7,10^1.7];

data.Henry1966.IPVx3.ageMos=6;
data.Henry1966.IPVx3.vaccSchedule='DTP+IPV 2,3,4 months';
data.Henry1966.IPVx3.pos=[9,9,10,7,7];
data.Henry1966.IPVx3.N=[10,9,11,9,10];

data.Henry1966.unvaccinated.vaccSchedule='no polio (DTP 2,3,4 months)';
data.Henry1966.unvaccinated.ageMos=6;
data.Henry1966.unvaccinated.pos=[9,9,11,6,5];
data.Henry1966.unvaccinated.N=[9,10,11,8,10];

data.Henry1966.IPVx4.vaccSchedule='DTP+IPV 2,3,4,15 months';
data.Henry1966.IPVx4.ageMos=16;
data.Henry1966.IPVx4.pos=[8,9,4,4,3];
data.Henry1966.IPVx4.N=[8,9,8,10,8];

data.Henry1966.tOPVx3.vaccSchedule='tOPV 7,8,9 months';
data.Henry1966.tOPVx3.ageMos=16;
data.Henry1966.tOPVx3.pos=[8,7,1,0,0];
data.Henry1966.tOPVx3.N=[11,12,10,8,9];


% this isn't a dose response experiment, but it puts important modern
% data about IPV boosting on the plot at high dose
% Jafari2014
% Jafari, Hamid, Deshpande, Jagadish M, Sutter, Roland W, Bahl, Sunil, Verma, Harish
% Ahmad, Mohammad, Kunwar, Abhishek, Vishwakarma, Rakesh, Agarwal, Ashutosh,
% Jain, Shilpi, Estivariz, Concepcion, Sethi, Raman, Molodecky, Natalie A,
% Grassly, Nicholas C, Pallansch, Mark A, Chatterjee, Arani, Aylward, R Bruce
% Science, 245(6199) 2014.
%   data from Table S1, ever excreting, Type 1 only 
data.Jafari2014.serotype=1;
data.Jafari2014.dose=10^5.7*[1 1 1];
data.Jafari2014.ageMos=[10,65,125]; % data for each age bin

data.Jafari2014.IPVboost.vaccSchedule='tOPV background, IPV boost, then challenge 4 weeks later';
data.Jafari2014.IPVboost.pos=[9,10,14];
data.Jafari2014.IPVboost.N=[102,110,104];

data.Jafari2014.IPVboost.vaccSchedule='tOPV background, bOPV boost, then challenge 4 weeks later';
data.Jafari2014.bOPVboost.pos=[25,23,28];
data.Jafari2014.bOPVboost.N=[102,108,111];

data.Jafari2014.tOPVxN.vaccSchedule='tOPV background, enroll but no boost, then challenge 4 weeks later';
data.Jafari2014.tOPVxN.pos=[15,26,55];
data.Jafari2014.tOPVxN.N=[104,108,105];


% TYPE 2 %%%%%%%%%%%%%%%%%%%

% Dane1961, OPV2
% Dane, D S, Dick, G W, Briggs, Moya, Nelson, R, McAlister, J, Connolly, J H, 
% Haire, Margaret, Mckeown, Florence, Field, C M B
% Vaccination against poliomyelitis with live virus vaccine: 
% 6. Changes in Sabin Type II oral vaccine virus after multiplication in the intestinal tract
% British Medical Journal, 2(5247) 1961.
%   data from Table X
data.Dane1961.serotype=2;

data.Dane1961.unvaccinated.Sabin2.dose=[10^4.6,10^3.6,10^2.6];
data.Dane1961.unvaccinated.Sabin2.vaccSchedule='naive';
data.Dane1961.unvaccinated.Sabin2.ageMos=11; % mean (5-17 range)
data.Dane1961.unvaccinated.Sabin2.pos=[7,9,6];
data.Dane1961.unvaccinated.Sabin2.N=[7,9,9];

% Notes (MGF): fascinating study, ethically interesting.  They did a dose
% response study with type 2 virus isolated from stool 5 days after Sabin 2
% challenge.  The monkey neurovirulence assays in the paper show that
% little genetic reversion had happened by this point, and, correpsondingly
% they find no difference in dose response in the overlapping region. From
% a cVDPV modeling point of view, it would have been fascinating had they
% used virus from stool a few weeks after challenge, when neurovirulence
% was high (and modern understanding is the two key attenuating nucleotides
% are likely to be reverted).  The idea of 
% challenging after any human passage at any time is ethically dubious to me, 
% but as the paper shows, in this setting transmission guaranteed exposure,
% and so the dose response study is no less ethical than any other study design
% that only partially vaccinated the ward.  

data.Dane1961.unvaccinated.HumanPassaged5DayStool.dose=[10^3,10^2,10^1];
data.Dane1961.unvaccinated.HumanPassaged5DayStool.vaccSchedule='naive';
data.Dane1961.unvaccinated.HumanPassaged5DayStool.ageMos=11; % mean (5-17 range)
data.Dane1961.unvaccinated.HumanPassaged5DayStool.pos=[9,4,2];
data.Dane1961.unvaccinated.HumanPassaged5DayStool.N=[9,9,9];

% This isn't a dose response study, but it puts important modern points
% on the figure

% O'Ryan2015, OPV2
% O'Ryan, Miguel, Bandyopadhyay, Ananda S, Villena, Rodolfo, Espinoza, Mónica
% Novoa, José, Weldon, William C, Oberste, M Steven, Self, Steve, Borate, Bhavesh R
% Asturias, Edwin J, Clemens, Ralf, Orenstein, Walter, Jimeno, José, Rüttimann, Ricardo
% Costa Clemens, Sue Ann
% Inactivated poliovirus vaccine given alone or in a sequential schedule with bivalent oral poliovirus vaccine in Chilean infants: a randomised, controlled, open-label, phase 4, non-inferiority study
% The Lancet Infectious Diseases. 3099(15), 2015.

data.ORyan2015.serotype=2;
data.ORyan2015.dose=10^5.7;
data.ORyan2015.ageMos=6.5;

data.ORyan2015.IPVx3.vaccSchedule='IPV 2,3,4 months';
data.ORyan2015.IPVx3.pos=158;
data.ORyan2015.IPVx3.N=171;

data.ORyan2015.IPVx2bOPVx1.vaccSchedule='IPV 2,3,4 months';
data.ORyan2015.IPVx2bOPVx1.pos=139;
data.ORyan2015.IPVx2bOPVx1.N=179;

data.ORyan2015.IPVx1bOPVx2.vaccSchedule='IPV 2,3,4 months';
data.ORyan2015.IPVx1bOPVx2.pos=132;
data.ORyan2015.IPVx1bOPVx2.N=164;


% Asturias2016, OPV2
%Asturias, Edwin J, Bandyopadhyay, Ananda S, Self, Steve, Rivera, Luis, 
% Saez-llorens, Xavier, Lopez, Eduardo, Melgar, Mario, Gaensbauer, James T,
% Weldon, William C, Oberste, M Steven, Borate, Bhavesh R, Gast, Chris, 
% Clemens, Ralf, Orenstein, Walter, G, Miguel O Ryan, Jimeno, José
% Humoral and intestinal immunity induced by new schedules of bivalent oral
% poliovirus vaccine and one or two doses of inactivated poliovirus vaccine 
% in Latin American infants : an open-label randomised controlled trial
data.Asturias2016.serotype=2;
data.Asturias2016.dose=10^5.7;
data.Asturias2016.ageMos=4;

% positive for shedding 1 week after challenge
data.Asturias2016.bOPVx3.vaccSchedule='bOPV 6,10,14 weeks'; % Group 1
data.Asturias2016.bOPVx3.pos=[150,146];
data.Asturias2016.bOPVx3.N=[197,187];

data.Asturias2016.bOPVx3IPVx1.vaccSchedule='bOPV + 1 IPV 6,10,14 weeks'; % Group 4
data.Asturias2016.bOPVx3IPVx1.pos=140;
data.Asturias2016.bOPVx3IPVx1.N=193;

data.Asturias2016.bOPVx3IPVx2.vaccSchedule='bOPV + 1 IPV 6,10,14,36 weeks'; % Group 5
data.Asturias2016.bOPVx3IPVx2.pos=141;
data.Asturias2016.bOPVx3IPVx2.N=189;

data.Asturias2016.tOPVx3.vaccSchedule='tOPV IPV 6,10,14 weeks'; % Group 3
data.Asturias2016.tOPVx3.pos=7;
data.Asturias2016.tOPVx3.N=186;

end