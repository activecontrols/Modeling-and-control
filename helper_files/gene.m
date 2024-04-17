classdef gene
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
        alleles % array of gene specific data that is optimized for by GA
        fitness % fitness value of gene based on a fitness function defined at the population level
        lqrFail % boolean recording if the genome can be successfully used in the LQR function
        odeStopped % boolean recording if the genome results in successful integration using ODE45
        states % states vs. time of simulated trajectory based on gene alleles (might want to remove to make gene more generally applicable)
        constCheck % logical array of which state constraints are met (might want to remove to make gene more generally applicable)
        parameters % All important parameters for our systems (segArray)
    end

    methods
        % CONSTRUCTOR
        %   Populates alleles with allele seed
        function obj = gene(allele_seed, parameters)
            obj.alleles = allele_seed;
            obj.parameters = parameters;

        end

        %CROSSOVER Summary of this function goes here
        %   Performs crossover between two gene objects
        %   
        %   RETURNS: array of alleles for child1 and child2
        function [child1, child2] = crossover(obj1, obj2)
            p1 = randi(size(obj1.alleles,1));
            p2 = randi(size(obj1.alleles,1));
            
            child1 = obj1.alleles;
            child2 = obj2.alleles;
            
            temp = child1(p1:p2);
            child1(p1:p2) = child2(p1:p2);
            child2(p1:p2) = temp;
            clear temp
        end

        % MUTATE applies random mutation to solution set of genetic algorithm
        %   Performs mutation on alleles of gene. 
        %   Note: seed correpsonds to the seed value of the gene pool
        %         (initially this is just Bryson's Rule but could change as
        %         persistent data becomes avaliable).
        function child = mutate(obj, mut_rate, seed, mut_factor)
            maxVal = seed * mut_factor;
            minVal = ones(size(seed)) * 1.0e-08; % I didn't use zero to avoid making Q not positive-definite
            sigma = (maxVal - minVal)/3;
            child = obj.alleles;
            
            for i = 1:length(child)
                if rand < mut_rate
                    child(i) = child(i) + sigma(i) * randn;
                end
            end
            
            child(child > maxVal) = maxVal(child > maxVal);
            child(child < minVal) = minVal(child < minVal);
        end

        %CONSTRAINT CHECK
        %   Checks against constraint array to ensure solution meats pre-defined 
        %   constraints.
        function check = constraintCheck(obj, stateLimits)
            %check against state limits
            minCheck = obj.states > stateLimits(:, 1);
            maxCheck = obj.states < stateLimits(:, 2);
        
            % account for states with no limits (isnan is a logical array with a 1 
            % for any NaN in the passed array so adding it to constCheck will make 
            % all states with no limit true
            minCheck = minCheck + isnan(stateLimits(:, 1)) .* ones(size(minCheck));
            maxCheck = maxCheck + isnan(stateLimits(:, 2)) .* ones(size(maxCheck));
            check = [minCheck, maxCheck];
        end
    end
end