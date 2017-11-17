% explore secondary contacts in BM1967

% secondary contacts pos
%page 121

S2posCount=15;
N=280;

%tOPV pos frac
tOPVposFrac=0.13;

% weight N by arm size
weightsN=[666,323,198,475];
weightsPos=[1.5,11,1,14];

% attribute 15 pos to mOPV2 arm vs tOPV arm

tOPVcontactN=(N*weightsN(4)/sum(weightsN));
tOPVposCount=(0.13*tOPVcontactN);

mOPV2contactN=round(N*weightsN(2)/sum(weightsN));
mOPV2posCount=S2posCount-tOPVposCount;

mOPV2posFrac=mOPV2posCount/mOPV2contactN