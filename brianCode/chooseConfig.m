        function [inputConfig,numDimers] =  chooseConfig(traj,weights)



            config = chooseRandom(weights,1);
            inputConfig =traj(config,:);
            numDimers = config-1;

        end
