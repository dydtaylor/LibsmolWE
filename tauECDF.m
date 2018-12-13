fluxes = NaN(7500, 10);
CDFS = cell(10,1);
x = [];
for j = 1:10
    data = load("ECDFs/ecdfFlux"+j+".txt");
    fluxes(:,j) = data;
    [CDFS{j},x]= ecdf(fluxes(:,j));
    %figure(j)
    %ecdf(fluxes(:,j))
    %set(gca,'Yscale','log')
end
clc