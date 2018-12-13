classdef BinDefinition<handle
   properties
        binEdges
        basinA
        basinB
        % begin modify
        radius
        nParticles 
        % end modify
    end

    methods

        % the constructor define properties to be used in functions in below 
        % first 3 arguments define the state space 
        % 4th argument are model specific parameters in a 1d array
        function obj = BinDefinition()
            % begin modify
            propertiesFromText = dlmread('./input/parameters.txt')
            obj.radius = propertiesFromText(3);
            obj.nParticles = propertiesFromText(5);
            xvals = 0:1:obj.nParticles;
            yval = binopdf(xvals,obj.nParticles,(pi/9));
            [i,j] = max(yval);
            bintransition = j-1;
            obj.binEdges = [0:obj.nParticles 1000] % bins defined by edges
            obj.basinA = [2:length(obj.binEdges)-1]; % defined with respect to the bin edges
            obj.basinB = [1]; % defined with respect to the bin edges

            %obj.binEdges = [0 linspace(1,bintransition+1,bintransition+2-1) (bintransition + floor((obj.nParticles-bintransition)/2)) 1000] % bins defined by edges
            %obj.basinA = [2:length(obj.binEdges)-1]; % defined with respect to the bin edges
            %obj.basinB = [1]; % defined with respect to the bin edges
            % end modify

        end
       
        function [inputConfig,dNP] =  chooseConfigFromBasinA(obj,traj,trajNP,weights)
        % chooses a configuration from basin A using the weights     
        % traj, input configurations which are from basin A
        % trajNP, input particle type counts
        % weights, weight of each input
            
            config = obj.selectConfigFromWeights(weights,1);
            inputConfig =traj(config,:);
            dNP = trajNP(config,:);
             
        end

        function x = selectConfigFromWeights(obj,p,n)
        % selecting a random config according to weights
            if nargin == 1
                n = 1;
            end
            k = length(p);
            size(p)
            
            p = reshape(p,k,1);
            
            x = sum(repmat(rand(1,n),k,1)> repmat(cumsum(p)/sum(p),1,n),1)+1;
        end
        % general functions

        function [binIndex] = assignToBin(obj,traj,trajNP)
                orderParameter = obj.findOrderParameter(traj,trajNP);
                binIndex = find(histc(orderParameter,obj.binEdges)); % 
        end

        function [orderParameter] = findOrderParameter(obj,traj,trajNP)
        % input: a data trajectory
        % output: an order parameter 
            % begin modify
            particleLength = 2*sum(trajNP); 
            orderParameter = numberInside(traj,particleLength,obj.radius); % for my specific implementation, i compute the number inside the ROI as my order parameter using the properties defined in the constructor
            % end modify

        end


        % my specific implementations

        function x = numberInside(coords,particleLength,radius)
        % compute the number inside the ROI
          x = 0 ;
            
          for iii = 1:2:particleLength
            
            dist = sqrt(coords(iii)^2 + coords(iii+1)^2);
            
            if dist < radius
            
              x = x +1;

            end

          end 
        end

    end
end
