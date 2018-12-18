for j = [1 2 5 10 20 50 100]
    data = load("tauSweep/fluxTau"+j+"N15.txt");
    figure()
    plot(1:7500,data)
    title("tau = "+j)
end