classdef gene < handle
    %GENE
    %   Genome class used to generate population within genetic algorithm.
    %   Properties:
    %       alleles    = Array of important identifying information associated
    %                    with the gene. Alleles determine genome performance
    %                    and alleles that optimize the fitness function are 
    %                    selected for by the algorithm.
    %
    %       states     = The states associated with the alleles of the gene.
    %
    %       constCheck = Boolean array of what allele constraints are met.

    properties
        alleles
        states
        constCheck
    end

    methods
        function obj = gene(allele_seed)
        % CONSTRUCTOR
        %   Populates alleles with allele seed

        obj.alleles = allele_seed;

        end

        function crossover(obj1, obj2)
        %CROSSOVER Summary of this function goes here
        %   Performs crossover between two gene objects
        
            p1 = randi(size(obj1.alleles,1));
            p2 = randi(size(obj1.alleles,1));
            
            temp = obj1.alleles(p1:p2);
            obj1.alleles(p1:p2) = obj2.alleles(p1:p2);
            obj2.alleles(p1:p2) = temp;
            clear temp
        end

        function mutate(obj, mut_rate, seed)
        % MUTATE applies random mutation to solution set of genetic algorithm
        %   Performs mutation on alleles of gene. 
        %   Note: seed correpsonds to the seed value of the gene pool
        %         (initially this is just Bryson's Rule but could change as
        %         persistent data becomes avaliable).
            
            maxVal = seed * 1e+04;
            minVal = ones(size(seed)) * 1.0e-08;
            sigma = (maxVal - minVal)/6;
            
            for i = 1:length(obj.alleles)
                if rand < mut_rate
                    obj.alleles(i) = obj.alleles(i) + sigma(i) * randn;
                end
            end
            
            obj.alleles(obj.alleles > maxVal) = maxVal(obj.alleles > maxVal);
            obj.alleles(obj.alleles < minVal) = minVal(obj.alleles < minVal);
        end

        function constraintCheck(obj, stateLimits)
        %CONSTRAINT CHECK
        %   Checks against constraint array to ensure solution meats pre-defined 
        %   constraints.
        
            %check against state limits
            minCheck = obj.states > stateLimits(:, 1);
            maxCheck = obj.states < stateLimits(:, 2);
        
            % account for states with no limits (isnan is a logical array with a 1 
            % for any NaN in the passed array so adding it to constCheck will make 
            % all states with no limit true
            obj.constCheck = [minCheck, maxCheck] + isnan(stateLimits);
        end

        function fitness = cost_function(xsegment, segArray)
        %COST_FUNCTION
        %   Evaluates fitness of genetic algorithm solution based on state
        %   information from simulation and reference critical value. Also cross
        %   checks state information with state min and state max vectors from
        %   parameter cell array.
        
            xcrit2 = segArray{5}(:, 2);
            weights = segArray{21};
            
            %evaluate fitness of solution
            zsegment = xcrit2 - xsegment;
            OS = max(abs(zsegment), [], 2);
            FE = abs(zsegment(:,end));
            fitness = 1/(weights(1) * (OS(1) + OS(2) + OS(3)) + weights(2)*(FE(1) + FE(2) + FE(3)));
        end
    end
end