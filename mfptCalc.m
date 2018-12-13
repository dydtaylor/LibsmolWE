fluxes = NaN(18, 7500);
for j = 2:19
    filename = "Data/flux"+j+".txt";
    data = load(filename);
    fluxes((j-1),:) = data;
end

meanFluxes = mean(fluxes,2);
MFPT = ones(size(meanFluxes))./(meanFluxes * 10^5);
figure(2)
hold on
semilogy(2:19, MFPT,'ro','linewidth',2);