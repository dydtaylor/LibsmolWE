function generateStartingTrajectoriesAndWeights()

%numParticles = 13; %the total number of particles in the box
%L=300; % box length of side
%L0=100; % roi radius 

    propertiesFromText = dlmread('./input/parameters.txt')
    L0 = propertiesFromText(3);
    L = propertiesFromText(4);
    numParticles = propertiesFromText(5);

    xvals = 0:1:numParticles;
    yval = binopdf(xvals,numParticles,(pi/9));
    [i,j] = max(yval);
    bintransition = j-1;
    save bintransition bintransition;
    CGBins = [0:1:numParticles numParticles+1]
    reps = 500; %initial reps, WE Parameter, to define later, taken from input
    yval = binopdf(xvals,numParticles,(pi/9));
    yvs = yval;
    weights = repmat(yval/reps, reps, 1);
    weights = reshape(weights,size(weights,1)*size(weights,2),1);
    save('./input/weights.mat','weights')

    trajectories = [];

    posInConfigs = [];
    posOutConfigs = [];


    while size(posInConfigs,1) < 1000 & size(posOutConfigs,1) < 1000

            NumRand=10000;
            Xconfigs=(rand(NumRand,1)-0.5)*L;
            Yconfigs=(rand(NumRand,1)-0.5)*L;
            Rconfigs=sqrt(Xconfigs.^2+Yconfigs.^2);
            IN = find(Rconfigs<L0);
            OUT = find(Rconfigs>=L0);



            posInConfigs = [Xconfigs(IN) Yconfigs(IN)];
            posOutConfigs = [Xconfigs(OUT) Yconfigs(OUT)];


            if size(posInConfigs,1) >=1000


                posInConfigs = posInConfigs(1:1000,:);

            end



            if size(posOutConfigs,1) >=1000


                posOutConfigs = posOutConfigs(1:1000,:);

            end


    end

    sizetrrows = length(yvs)*reps;


    %numParticles = findNumParticles(numparticles);

    %trajectoriesNP = repmat(numParticles,sizetrrows,1);

    %allNP = trajectoriesNP;
    nin=[]
    for i=1:length(yvs)
        for j=1:(reps)
                    N_In = randi([CGBins(i) CGBins(i+1)-1],1,1)
                    nin = [nin N_In];
                    trajC = [datasample(posInConfigs,N_In,1);datasample(posOutConfigs,numParticles-N_In,1 )];
                    tr = reshape(trajC',[1 numParticles*2]);
                    trajectories = [trajectories; tr];
                    if numberInside(tr,numParticles*2,L0) ~= N_In
                        disp('incorrect')
                    end

        end
    end


    sizetrrows = length(trajectories);

    numMonomers =numParticles;
    numDimers = 0;

    trajectoriesNP = repmat([numMonomers numDimers],sizetrrows ,1);
    save('./input/trajectories.mat','trajectories')
    save('./input/trajectoriesNP.mat','trajectoriesNP')
    save('./input/numParticles.mat','numParticles')

end




