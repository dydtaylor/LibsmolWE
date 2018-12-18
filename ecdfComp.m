load('brianCode/Output/WE_flux.mat');
fluxes = NaN(7500,10);
for j = [1:7,9,10]
    fluxes(:,j) = load("ECDFs/NR500/ecdfFlux"+j+".txt");
end
fluxes(:,8) = [];
figure()
hold on
for j = 1:9
    [f,x] = ecdf(fluxes(5000:end,j));
    plot(x,f)
end
title('libsmol ECDF')
figure()
ecdf(flux_est(1900:end))
title('Brian ECDF')
