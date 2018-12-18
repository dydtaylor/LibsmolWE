fluxes = NaN(18, 7500);
bruteMeans = NaN(18,1);
bruteErr = NaN(18,1);
for j = 2:18
    filename = "Data/flux"+j+".txt";
    data = load(filename);
    fluxes((j-1),:) = data;
    bruteData = load("oldData/nSweep/"+j+".txt");
    bruteMeans(j) = mean(bruteData);
    bruteErr(j) = std(bruteData/sqrt(numel(bruteData)));
end

meanFluxes = mean(fluxes,2);
MFPT = ones(size(meanFluxes))./(meanFluxes * 10^5);
figure()
hold on
plot(2:19, MFPT,'ro','linewidth',2);
errorbar(2:18,bruteMeans(2:18),bruteErr(2:18),'bo','linewidth',2);
set(gca,'yscale','log')
title('MFPT')

figure()
plot(2:18,log(bruteMeans(2:18))-log(MFPT(2:18)))
title('Log Brute Force - Log WE')