classdef Grid<handle
    % a data structure that collects configurations from the simulation
    properties
        weights % cell array containing the configurations
        bins % cell array containing coordinates
        binsNP;
        numberOfBins1;
        target;
        %np;
        basinA;
    end
    
    
    
    methods
        
        
            
            
        function obj = Grid(nBins1,Mtarg,basina)
            % contructor that creates a cell array bins and weights to store trajectories and weights
            % nBins1, total number  of bins
            % Mtarg, target number of replicas
            % np, the total number of particles
            obj.numberOfBins1 = nBins1;
            %obj.np = np;

            obj.weights = cell(nBins1,1);
            obj.bins = cell(nBins1,1);
            obj.binsNP = cell(nBins1,1);
            obj.basinA = basina;


            obj.target = Mtarg;
            for i=1:numel(obj.weights)
                obj.weights{i} = [];
                obj.bins{i} = [];
                obj.binsNP{i}=[];
            end
            
        end


      % function [linearIndex] = convertToLinearIndex(obj,inputMatrix)
      %      ind1 = inputMatrix(1)+1;
      %      ind2 = inputMatrix(2)+1;
      %      matrix = obj.numberOfBins1.*ones(obj.np,1)
      %      linearIndex = sub2ind(matrix,ind1,ind2)
      % end
        
            function balance(obj)
                for i=1:numel(obj.weights)
                    if size(obj.weights{i},1) == 0
                        continue
                    else
                        obj.idealSplit(i);
                        obj.idealMerge(i);
                          if size(obj.weights{i},1) <  obj.target
                              obj.split(i);
                          elseif size(obj.weights{i},1) >  obj.target
                              obj.merge(i);
                          end
                    end
                end
            end
        
        
            function idealSplit(obj,indexBin)
            % splits replicas according to ideal weight in bin 
                idealWeight = sum(obj.weights{indexBin})/obj.target;
                [B,I]  = sort(obj.weights{indexBin},'descend');
                obj.weights{indexBin} = B;
                N = size(obj.weights{indexBin},1);
                z = sum(obj.weights{indexBin});
                obj.bins{indexBin} = obj.bins{indexBin}(I,:);
                obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);
                

                
                insertConfigs = [];
                indsNotSplit = [];
                insertIndices = [];
                for j=1:N
                    weight = obj.weights{indexBin}(j,:);
                    split = weight;
                    ind = 1;
                    while split > 2*idealWeight
                        
                        
                        
                        split = weight/(ind+1);
                        ind = ind +1;
                       
                        
                    end
                    if ind==1
                        indsNotSplit=[indsNotSplit;j];
                        ind = 0;
                    end
                    for i=1:ind
                        insertConfigs = [insertConfigs;split];
                        insertIndices = [insertIndices;j];
                    end
                                    
                end
                if isempty(insertConfigs) == 0
                    obj.weights{indexBin} = vertcat(obj.weights{indexBin}(indsNotSplit,:),insertConfigs);
                    obj.bins{indexBin} = vertcat(obj.bins{indexBin}(indsNotSplit,:),obj.bins{indexBin}(insertIndices,:));
                    obj.binsNP{indexBin} = vertcat(obj.binsNP{indexBin}(indsNotSplit,:),obj.binsNP{indexBin}(insertIndices,:));

                end
                
            end
        
         function idealMerge(obj,indexBin)
         %merge replicas according to ideal weight
             [B,I]  = sort(obj.weights{indexBin},'ascend');
             obj.weights{indexBin} = B;
             obj.bins{indexBin} = obj.bins{indexBin}(I,:);
             obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);

             tw = sum(B);
             idealWeight = sum(obj.weights{indexBin})/obj.target;
             
             while obj.weights{indexBin}(1) < 0.5*idealWeight && size(obj.weights{indexBin},1) > 1
                 totalWeight = 0;
                
                 N = size(obj.weights{indexBin},1);
                 stopIndex = 1;
                 for i=1:N
                     totalWeight = totalWeight +obj.weights{indexBin}(i);
                     if (totalWeight < idealWeight) && (obj.weights{indexBin}(i) < 0.5*idealWeight)
                         
                         stopIndex = stopIndex +1;
                     elseif (totalWeight < idealWeight) && (obj.weights{indexBin}(i) > 0.5*idealWeight) && (totalWeight+obj.weights{indexBin}(i) < 1.5*idealWeight)
                         
                         stopIndex = stopIndex +1;
                     else
                         break;
                     end
                 end
                 
                 indexKeep = discreternd(obj.weights{indexBin}(1:stopIndex,:));
                 
                 obj.weights{indexBin} = obj.weights{indexBin}(stopIndex:end,:);
                 obj.weights{indexBin}(1) = totalWeight;
                 obj.bins{indexBin}(stopIndex,:) = obj.bins{indexBin}(indexKeep,:); 
                 obj.bins{indexBin} = obj.bins{indexBin}(stopIndex:end,:);
                 obj.binsNP{indexBin}(stopIndex,:) = obj.binsNP{indexBin}(indexKeep,:); 
                 obj.binsNP{indexBin} = obj.binsNP{indexBin}(stopIndex:end,:);

                 [B,I]  = sort(obj.weights{indexBin},'ascend');
                 obj.bins{indexBin} = obj.bins{indexBin}(I,:);
                 obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);

                 tw2 =sum(B);
                 obj.weights{indexBin} = B;
                 idealWeight = sum(obj.weights{indexBin})/obj.target;
             end
             function x = discreternd(p,n)
                 if nargin == 1
                     n = 1;
                 end
                 p = p/sum(p);
                 k = length(p);
                 p = reshape(p,k,1);
                 x = sum(repmat(rand(1,n),k,1)> repmat(cumsum(p)/sum(p),1,n),1)+1;
             end
         end
       
        function merge(obj,indexBin)
        % merge replicas to meet target number
           [B,I]  = sort(obj.weights{indexBin},'descend');
           obj.weights{indexBin} = B;
           
           
           while size(obj.weights{indexBin},1) > obj.target
               if obj.weights{indexBin}(end-1,:)/(obj.weights{indexBin}(end-1,:) + obj.weights{indexBin}(end,:)) > rand
                   obj.weights{indexBin}(end-1,:) = obj.weights{indexBin}(end-1,:) + obj.weights{indexBin}(end,:);
                   obj.weights{indexBin} = obj.weights{indexBin}(1:end-1,:);
                   I = I(1:end-1);
                    obj.bins{indexBin} = obj.bins{indexBin}(I,:);
                    obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);

                   [B,I]  = sort(obj.weights{indexBin},'descend');
                    obj.weights{indexBin} = B;
                    obj.bins{indexBin} = obj.bins{indexBin}(I,:);
                    obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);
                   

                   
               else
                   obj.weights{indexBin}(end-1,:) = obj.weights{indexBin}(end-1,:) + obj.weights{indexBin}(end,:);
                   obj.weights{indexBin} = obj.weights{indexBin}(1:end-1,:);
                   I = [I(1:end-2); I(end)];
                   obj.bins{indexBin} = obj.bins{indexBin}(I,:);
                   obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);

                   [B,I]  = sort(obj.weights{indexBin},'descend');
                    obj.weights{indexBin} = B;
                    obj.bins{indexBin} = obj.bins{indexBin}(I,:);
                    obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);
                   

               end
               
               
           end
           
               
        end
        
        
        
        
        function split(obj,indexBin)
        % split to meet replica target number
           [B,I]  = sort(obj.weights{indexBin},'ascend');
           obj.weights{indexBin} = B;
           disp('size');
          
           while size(obj.weights{indexBin},1) < obj.target
               
         
              
               obj.weights{indexBin}(end,:) = obj.weights{indexBin}(end,:)/2;
               obj.weights{indexBin}(end+1,:) = obj.weights{indexBin}(end,:);
                [B,I]  = sort(obj.weights{indexBin},'ascend');
                obj.weights{indexBin} = B;
                 obj.bins{indexBin} = vertcat(obj.bins{indexBin},obj.bins{indexBin}(end,:));
               obj.bins{indexBin} = obj.bins{indexBin}(I,:);
               obj.binsNP{indexBin} = vertcat(obj.binsNP{indexBin},obj.binsNP{indexBin}(end,:));
               obj.binsNP{indexBin} = obj.binsNP{indexBin}(I,:);
               

           end
          
               
        end
        
   
        
        function addToGrid(obj,config,NPconfig,weight,indexBin)
            %add the trajectory data and weight into a bin at indexBin
            %indexBin = obj.convertToLinearIndex(indexBin);

            % add a config at an interface
            bs = obj.bins{indexBin};
            if (size(bs,1)==0) %isempty(obj.bins{indexBin})%if I get error I think bc it is[]and I should try command equivalent to isempty 
                obj.bins{indexBin} = config;
                obj.binsNP{indexBin} = NPconfig;

                obj.weights{indexBin} = weight;
            elseif ~(size(obj.bins{indexBin},1)==0)
                obj.bins{indexBin} = vertcat(obj.bins{indexBin},config);
                obj.binsNP{indexBin} = vertcat(obj.binsNP{indexBin},NPconfig);

                obj.weights{indexBin} = vertcat(obj.weights{indexBin},weight);
            end
        end

        function [tweights,tr,tNP]  = getTrajectoriesAndWeightsFromA(obj,basinA)
        % retrieve all weights and trajectories (both have same order) and clear all bins 
            tweights=[];
            tr = [];
            tNP = [];
            for i=1:numel(obj.basinA)
                indexA = basinA(i);
                tweights = [tweights;obj.weights{indexA}];
                tr= [tr; obj.bins{indexA}];
                tNP= [tNP; obj.binsNP{indexA}];

            end
               
        end

        
        function [tweights,tr,tNP]  = getTrajectoriesAndWeightsClear(obj)
        % retrieve all weights and trajectories (both have same order) and clear all bins 
            tweights=[];
            tr = [];
            tNP = [];
            for i=1:numel(obj.weights)
                tweights = [tweights;obj.weights{i}];
                tr= [tr; obj.bins{i}];
                tNP= [tNP; obj.binsNP{i}];

                obj.bins{i} = [];
                obj.binsNP{i} = [];
                obj.weights{i} = [];

            end
               
        end

        
    function [tweights,tr,tNP]  = getTrajectoriesAndWeights(obj)
    % retrieve all weights and trajectories (both have same order without clearing bins)
            tweights=[];
            tNP = [];
            tr = [];
            for i=1:numel(obj.weights)
                
                    tweights = [tweights;sum(obj.weights{i})];
                    tr= [tr; obj.bins{i}];
                    tNP= [tNP; obj.binsNP{i}];
                

                %obj.bins{i} = [];
                %obj.weights{i} = [];
            end
               
        end
        
        
    end
    
end
