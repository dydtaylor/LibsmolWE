classdef SimulationParameterRetrieval<handle
   properties
        % begin modify 
        kon;
        koff;
        difcm;
        difcd;
        radiusROI;
        boxLength ;
        numberParticles;
        smoldynTimeStep;
        % end modify
    end

    methods
            function obj = SimulationParameterRetrieval()
                % begin modify
                propertiesFromText = dlmread('./input/parameters.txt')
                obj.kon = propertiesFromText(1);
                obj.koff = propertiesFromText(2);
                obj.radiusROI = propertiesFromText(3);
                obj.boxLength = propertiesFromText(4);
                obj.numberParticles = propertiesFromText(5);
                obj.smoldynTimeStep = propertiesFromText(6);
                obj.difcm = propertiesFromText(7);
                obj.difcd = propertiesFromText(8);
                % end modify
            end


            function [dataOut,dataCounts,endtime] = doSimulation(obj,traj,trajCounts,weParam,seed,stopflag)
                %begin modify
                simulationParams = struct('kon',obj.kon,'koff',obj.koff,'difcd',obj.difcd,'difcm',obj.difcm,'stopflag',stopflag,'radiusROI', obj.radiusROI, 'boxLength', obj.boxLength);
                modelParams = struct('numberParticles',obj.numberParticles,'dt',obj.smoldynTimeStep);
                %end modify
                [dataOut,dataCounts,endtime] = DynamicsEngine(traj,trajCounts,weParam,simulationParams,modelParams,seed);

            end
        


    end

end
