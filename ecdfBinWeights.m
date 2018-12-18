weights = zeros(15,10);
for j = 1:10
    data = load("ECDFs/NR500/ecdfSim"+j+".txt");
    for k = 1:length(data)
        bin = data(k,1);
        weights(bin,j) = weights(bin,j) + data(k,2);
    end
end
load("brianCode/output/BinWeights.mat");
meanOld = mean(BinWeights,2);