function probabilityInfected = doseResponseModel(dose,NAb,serotype,b)
if nargin<1
    dose=10^5.7; %TCID50
end
if nargin<2
    NAb=1;
end
if nargin<3
    serotype=2;
    b=[14.24,0.444,0.545];
end
if nargin<4 
    if serotype==1
        b=[14.24,0.444,0.545];
    elseif serotype==2 
        b=[7.9611,0.444,0.545];
    elseif serotype==3 
        b=[17.8232,0.444,0.545];
    else
        b=[1.83,0.444,0.545];
    end
end

% beta-poisson dose response model
probabilityInfected = 1-(1+dose/b(1)).^(-b(2)./NAb.^b(3)); 

end