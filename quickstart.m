FluxFilename = "simFluxes.txt";

DynParamsID = fopen('dynamicsParams.txt');
DynParams = textscan(DynParamsID, '%q %f');
DynParams = DynParams{2};
fclose(DynParamsID);
WEParamsID = fopen('WEParams.txt');
WEParams = textscan(WEParamsID, '%q %f');
WEParams = WEParams{2};
fclose(WEParamsID);
dt = DynParams(1);
tau = WEParams(1);

TimescaleFactor = 1/(dt*tau);

fluxData = load(FluxFilename)*TimescaleFactor;

meanFluxes = mean(fluxData(end/2:end));

MFPT = 1/meanFluxes;

fprintf('The MFPT measured from the dataset is: %f \n The domain length is: %f \n The number of molecules in the domain is: %f \n The domain density is: %f \n The dynamics timestep (dt) is: %f \n The WE timestep (tau) is: %f (units of dt) \n', MFPT, DynParams(2), DynParams(8), DynParams(8)/(DynParams(2)^2), DynParams(1), WEParams(1))

figure
plot(fluxData)
set(gca,'yscale','log')
title('Flux measured at timestep')
ylabel('WE Step')
xlabel('Flux')