% calculate HID50 from dose response model
reps=1000;
for k=1:4

    if k==1
        b=[14.24,0.444,0.4624];
        bCI=[3,0.29,0.42;59,0.83,0.50];
    elseif k==2 
        b=[7.9611,0.444,0.4624];
        bCI=[2,0.29,0.42;29.9,0.83,0.50];
    elseif k==3 
        b=[17.8232,0.444,0.4624];
        bCI=[4.66,0.29,0.42;62.79,0.83,0.50];
    elseif k==4 % WPV 
        b=[1.83,0.444,0.4624];    
        bCI=[0.20,0.29,0.42;29,0.83,0.50];
    end

    HID50(k)=fzero(@(x) 0.5-doseResponseModel(x.^2,1,k,b), 1).^2;

    
    CI=nan(reps,1);
    boot=nan(3,reps);
    % strong parameter correlation
        rho=0.87;
        boot(1,:)=logspace(log10(bCI(1,1)),log10(bCI(2,1)),reps);
        boot(2,:)=logspace(log10(bCI(1,2)),log10(bCI(2,2)),reps);
        L=reps/2*[1.12,0.99];
    for n=1:reps
        rng(n);
        idx=0;
        while any(idx<1 | idx>reps)
            idx(1)=round(L(1)+reps/4*randn);
            idx(2)=round(L(2)+rho*(idx(1)-L(1))+sqrt((1-rho^2))*reps/4*randn);
        end

        CI(n)=fzero(@(x) 0.5-doseResponseModel(x.^2,1,k,[boot(1,idx(1)),boot(2,idx(2)),b(3)]), 1).^2;

    end
    HIDCI(k,:)=quantile(CI,[0.025,0.975]);
end