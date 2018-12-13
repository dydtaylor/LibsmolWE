function WE_main(strTau,strNumreplicas,type,inputflag,stopflag,strIter)
    % timestep , the weighted ensemble timestep
    % NP, the total number of particles
    % Kon, dimerization on rate
    % Koff, dimerization off rate
    % DT, the smoldyn timestep


    function [flux_est]=  WE_equil(WENitersMax,weParams,Ggrid,GBinDef,GSimDef,nzconfigs,nznp,nzweights)
        %%%% helper functions begin %%%%
        function propagate(replicasCopyNumbers,replicasXYPositions,weight,tau)
            
            % call simulation program, returns probability entering B per unit time
            trajsize = size(replicasXYPositions,1);
            seeds = randi(intmax('int16'),1 ,trajsize);
            stopflag = 0; %off
            for ii = 1:trajsize % propagate each replica
              trajrow = replicasXYPositions(ii,:);
              trajNProw = replicasCopyNumbers(ii,:);
              seed = seeds(ii);
              
              [newTrajRow,trajNPRow,endtime] = GSimDef.doSimulation( trajrow,trajNProw,weParams,seed,stopflag); % calls smoldyn cpp evacuation simulation

              replicasXYPositions(ii,:)= newTrajRow;
              replicasCopyNumbers(ii,:) = trajNPRow;

            end


            weightend = size(weight,1);
            for i =1:weightend % add new configurations and weights back to the Grid
                config = replicasXYPositions(i,:); %x y coordinates
                configNP = replicasCopyNumbers(i,:); 
                [binIndex] = GBinDef.assignToBin(config,configNP);
                Ggrid.addToGrid(config,configNP,weight(i),binIndex);
            end
        end

     %%%% helper functions end %%%%
        

        %%%% main code start %%%%
        bw = []; % container for total weight in each bin at each iteration

 

        
        for j=1:numel(nzweights) % initialize the grid with starting trajectories and weightsj
            config = nzconfigs(j,:);
            configNP = nznp(j,:);
            [binIndex] = GBinDef.assignToBin(config,configNP);
            Ggrid.addToGrid(config,configNP,nzweights(j,:),binIndex);
        end

        Ggrid.balance();
        taus_elapsed = 0
        BinWeights = [];
        ElapsedTime = []
        iter=0
        totalProbIntoB = 0;
        tauavgs = [];
        while (size(ElapsedTime,1) <WENitersMax)
            
            tic


        
            [weights,replicasXYPositions,replicasCopyNumbers] = Ggrid.getTrajectoriesAndWeightsClear();
        

            

            
            

            propagate(replicasCopyNumbers,replicasXYPositions,weights,tau) % simulate all trajectories and update grid, reinitialize in A if reached

            Ggrid.balance();

            [binweights,tr,trajNP] = Ggrid.getTrajectoriesAndWeights(); 
            BinWeights = [BinWeights binweights];
            save('./output/BinWeights.mat', 'BinWeights')


            elapsedtime = toc;
            ElapsedTime = [ElapsedTime elapsedtime];
            save('./output/ElapsedTime.mat', 'ElapsedTime')
            

            
            
        end

    end 

    function [flux_est]=  WE_flux(WENitersMax,weParams,Ggrid,GBinDef,GSimDef,nzconfigs,nznp,nzweights,stopflag)
        %%%% helper functions begin %%%%
        
        
        function [ProbIntoB] =  propagateAndReintroduce(replicasCopyNumbers,replicasXYPositions,weight,tau,replicasXYPositionsA,weightA,replicasCopyNumbersA)
            
            % call simulation program, returns probability entering B per unit time
            trajsize = size(replicasXYPositions,1);
            ProbIntoB = 0;
            seeds = randi(intmax('int16'),1 ,trajsize);
            for ii = 1:trajsize % propagate each replica
              trajrow = replicasXYPositions(ii,:);
              trajNProw = replicasCopyNumbers(ii,:);
              seed = seeds(ii);
              [newTrajRow,trajNPRow,endtime] = GSimDef.doSimulation( trajrow,trajNProw,weParams,seed,stopflag); % calls smoldyn cpp evacuation simulation

              replicasXYPositions(ii,:)= newTrajRow;
              replicasCopyNumbers(ii,:) = trajNPRow;

              
              if (ismember(GBinDef.assignToBin(replicasXYPositions(ii,:),replicasCopyNumbers(ii,:)),GBinDef.basinB) ) % if replica ends up in basin B
                    
                    [aConfig,aNP] = GBinDef.chooseConfigFromBasinA(replicasXYPositionsA,replicasCopyNumbersA,weightA); %choose new configuration to insert from basin A
                    replicasXYPositions(ii,:) = aConfig;
                    replicasCopyNumbers(ii,:)  = aNP;
                    ProbIntoB = ProbIntoB +weight(ii);  %add weight from this particle to probablity entering B
              end



            end


            weightend = size(weight,1);
            for i =1:weightend % add new configurations and weights back to the Grid
                config = replicasXYPositions(i,:); %x y coordinates
                configNP = replicasCopyNumbers(i,:); 
                [binIndex] = GBinDef.assignToBin(config,configNP);
                Ggrid.addToGrid(config,configNP,weight(i),binIndex);
            end
        end

        %%%% helper functions end %%%%
        

        %%%% main code start %%%%
        flux_est=[]; % container for flux estimate at each iteration
        bw = []; % container for total weight in each bin at each iteration

 

        
        for j=1:numel(nzweights) % initialize the grid with starting trajectories and weightsj
            config = nzconfigs(j,:);
            configNP = nznp(j,:);
            [binIndex] = GBinDef.assignToBin(config,configNP);
            Ggrid.addToGrid(config,configNP,nzweights(j,:),binIndex);
        end

        Ggrid.balance();
        taus_elapsed = 0
        BinWeights = [];
        ElapsedTime = []
        iter=0
        totalProbIntoB = 0;
        tauavgs = [];
        while (size(flux_est,1) <WENitersMax)
            
            tic


        
            [weightsA,replicasXYPositionsA,replicasCopyNumbersA] = Ggrid.getTrajectoriesAndWeightsFromA(GBinDef.basinA);
            [weights,replicasXYPositions,replicasCopyNumbers] = Ggrid.getTrajectoriesAndWeightsClear();
        

            

            
            

            [ProbIntoB] = propagateAndReintroduce(replicasCopyNumbers,replicasXYPositions,weights,tau,replicasXYPositionsA,weightsA,replicasCopyNumbersA) % simulate all trajectories and update grid, reinitialize in A if reached

            Ggrid.balance();

            [binweights,tr,trajNP] = Ggrid.getTrajectoriesAndWeights(); 
            BinWeights = [BinWeights binweights];
            save('./output/BinWeights.mat', 'BinWeights')

            flux_est = [flux_est ; ProbIntoB/tau];
            save('./output/WE_flux.mat', 'flux_est')

            elapsedtime = toc;
            ElapsedTime = [ElapsedTime elapsedtime];
            save('./output/ElapsedTime.mat', 'ElapsedTime')
            

            
            
        end




    end

%%%% parameters 
%addpath('cpp')
WENitersMax = str2num(strIter);  
tau =str2num(strTau)
numreplicas=str2num(strNumreplicas);
weParams = struct('tau',tau,'numreplicas',numreplicas);
GSimDef = SimulationParameterRetrieval();
GBinDef =BinDefinition();
if strcmp(inputflag,'false')
    generateStartingTrajectoriesAndWeights()
end
load('input/weights.mat', 'weights')
load('input/trajectories.mat', 'trajectories')
load('input/trajectoriesNP.mat', 'trajectoriesNP');
if strcmp(stopflag,'false')
    stopflag = 0
else
    stopflag = 1
end
nzconfigs = trajectories;
nzweights = weights;
nznp = trajectoriesNP;
nBins = length(GBinDef.binEdges) -1;
Ggrid = Grid(nBins,numreplicas,GBinDef.basinA);
if strcmp(type,'flux')
    flux_est = WE_flux(WENitersMax,weParams,Ggrid,GBinDef,GSimDef,nzconfigs,nznp,nzweights,stopflag);
elseif strcmp(type,'equil')
    WE_equil(WENitersMax,weParams,Ggrid,GBinDef,GSimDef,nzconfigs,nznp,nzweights);
end

end
