%% loading data -> 4 batches
b1= load("batch1_10sim.mat");
b2= load("batch2_10sim.mat");
b3= load("batch3_20sim.mat");
b4= load("batch4_20sim.mat");

%% Check for nan values (b3 and b4 has nan values)

hasNaN= zeros(1,20);

% b3
for i = [1:20];
    temp_arg3 = any(isnan(b3.mainmatrix(:,:,:,i)));  % it gives an array of logicals-> true-isnan, false- not nan
    hasNaN3(i)= any(temp_arg3(:)==true); % gives a single value if any of the element is true
end

% b4
for i = [1:20];
    temp_arg4 = any(isnan(b4.mainmatrix(:,:,:,i)));  % it gives an array of logicals-> true-isnan, false- not nan
    hasNaN4(i)= any(temp_arg4(:)==true); % gives a single value if any of the element is true
end

%% array for no nan values
no_nan_array3= [1:20];
no_nan_array4= [1:20];

no_nan_array3 = no_nan_array3(~hasNaN3);
no_nan_array4 = no_nan_array4(~hasNaN4);

%%
mainmatrix= cat(4,cat(4,cat(4,b1.mainmatrix,b2.mainmatrix),b3.mainmatrix(:,:,:,no_nan_array3)),b4.mainmatrix(:,:,:,no_nan_array4));
allperUP= mean(mainmatrix,4);
allthetaes= b1.allthetaes;
allbetas=b1.allbetas;
allthetaas= b1.allthetaas;

%% check if nothing in all per up is NaN
any(isnan(allperUP(:)))

%% 
save("allbetas.mat","allbetas")
save("allthetaas.mat", "allthetaas")
save("allthetaes.mat", "allthetaes")
save("allperUP.mat", "allperUP")

