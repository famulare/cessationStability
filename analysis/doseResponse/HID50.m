% calculate HID50 from dose response model
reps=1000;
for k=1:4

    if k==1
        b=[14.24,0.444,0.545];
        bCI=[3,0.29,0.51;59,0.83,0.57];
    elseif k==2 
        b=[7.9611,0.444,0.545];
        bCI=[2,0.29,0.51;29.9,0.83,0.57];
    elseif k==3 
        b=[17.8232,0.444,0.545];
        bCI=[4.66,0.29,0.51;62.79,0.83,0.57];
    elseif k==4 % WPV 
        b=[1.83,0.444,0.545];    
        bCI=[0.20,0.29,0.51;29,0.83,0.57];
    end

    HID50(k)=fsolve(@(x) 0.5-doseResponseModel(x,1,k,b), 1);

    CI=nan(reps,1);
    boot=nan(3,reps);
    boot(1,:)=logspace(log10(bCI(1,1)),log10(bCI(2,1)),reps);
    boot(2,:)=linspace(bCI(1,2),bCI(2,2),reps);
    boot(3,:)=linspace(bCI(1,3),bCI(2,3),reps);
    for n=1:reps
        idx=0;
        while any(idx<1 | idx>reps)
            idx=round([reps/2+reps/4*randn(3,1)]);
        end

        CI(n)=fsolve(@(x) 0.5-doseResponseModel(x,1,k,[boot(1,idx(1)),boot(2,idx(2)),boot(3,idx(3))]), 1);

    end
    HIDCI(k,:)=quantile(CI,[0.025,0.975]);
end