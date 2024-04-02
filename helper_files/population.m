classdef population < handle & matlab.mixin.Copyable
    %POPULATION
    %   Collection of genes used in genetic algorithm

    properties
        allele_seed % starting seed for alleles of each gene in population
        popSize % population size
        mut_rate1 % initial mutation rate
        mut_rate2 % final mutation rate
        mut_func % function that evaluates mutation rate b/w mut_rate1 and mut_rate2 based on current pop size
        gen_cut % general cutoff location (lower percentile killed)
        elite_cut = 0 % elite cutoff location (upper percentile perserved without radnom changes)
        fit_func % function to evaluate fitness of a gene within the population
        nodes % cell array of nodes that make up the population
        nodeRanking % array of node indexes sorted from most fit to least fit
        numAlive % number of nodes alive in population
    end

    methods
        % CONSTRUCTOR
        %   Constructs initial population
        function obj = population(allele_seed, popSize, mut_rate1, mut_rate2, mut_func, gen_cut, elite_cut, fit_func)
            % Add user defined params
            obj.allele_seed = allele_seed;
            obj.popSize = popSize;
            obj.mut_rate1 = mut_rate1;
            obj.mut_rate2 = mut_rate2;
            obj.mut_func = mut_func;
            obj.gen_cut = gen_cut;
            obj.elite_cut = elite_cut;
            obj.fit_func = fit_func;
            obj.nodeRanking = 1:popSize; % default to unordered list of nodes
            obj.numAlive = popSize;

            % Generate nodes
            for i = 1:obj.popSize
                % Create node
                obj.nodes{i} = node([], [], 1, 0, true, gene(allele_seed), []);

            end

            % Initial mutation and crossover (not reproduction! this just
            % introduces some initial randomness to the population)
            for i = 1:length(obj.nodes)
                child = obj.nodes{i}.genome.mutate(mut_func(mut_rate1, mut_rate2, popSize), allele_seed);
                obj.nodes{i}.genome.alleles = child;
            end

            for i = 1:2:length(obj.nodes)-2
                [child1, child2] = obj.nodes{i}.genome.crossover(obj.nodes{i+2}.genome);
                obj.nodes{i}.genome.alleles = child1;
                obj.nodes{i+2}.genome.alleles = child2;
            end
        end

        % FIT EVAL
        %   Evaluates fitness of a population using given fitness function
        %   handle. Returns ordered list of node indexes in the population
        %   for use in the KILL method
        function fitEval(pop)
            ns = pop.nodes; %makes code cleaner
            fitArray = zeros(size(ns)); %initialize fitArray
            
            for i = 1:length(pop.nodes)
               fit = pop.fit_func(ns{i}.genome);
               ns{i}.genome.fitness = fit;
               fitArray(i) = fit;
            end
            
            pop.nodeRanking = conditional_quicksort(1:length(ns), fitArray);

            % CONDITIONAL_QUICKSORT
            %   Uses quick sort to sort vector elements corresponding to the sorted
            %   condition array from greatest to least
            function [vector] = conditional_quicksort(vector, condition)
                if numel(vector) <= 1 % vectors with one or less elements are sorted
                    return
                else
                    pivot=condition(1);  % take first value as pivot element
                                      % randomization would help avoiding worst case
                                      % runtime
                    % We need three partitions in order to make use of Matlabs
                    % in-place processing feature.
                    vector = [ conditional_quicksort(vector(condition > pivot), condition(condition > pivot))...
                               vector(condition == pivot)...
                               conditional_quicksort(vector(condition < pivot), condition(condition < pivot)) ]; 
                end
            end
        end

        % KILL
        %   Kills nodes based on general cuttoffs and marks nodes as elite
        %   based on elite cuttoffs
        function kill(pop)
            % Make sure we are only evaluating currently alive genes
            r = pop.nodeRanking(1:pop.numAlive);

            % Mark elite genes
            i_elite = r(1:floor(length(r)*pop.elite_cut)); 
            if isempty(i_elite); i_elite = 1; end % need to make sure i_elite always at least 1
            for i = i_elite
                pop.nodes{i}.elite = true;
            end
            
            % Mark dead genes
            if (ceil(length(r)*pop.gen_cut) + 1) > length(r) % Minimum of last gene always dies
                i_kill = r(end);
            else
                i_kill = r((ceil(length(r)*pop.gen_cut) + 1):end);
            end

            for i = i_kill
                pop.nodes{i}.alive = false;
            end

            % Total number of alive genes
            pop.numAlive = length(r) - length(i_kill);
        end

        % REPRODUCE
        %   Creates new population through mutation and crossover. Tracks
        %   lineage.
        function newPop = reproduce(oldPop)
            % Copy old pop (Can't just use copy() without node not being a
            % handle class which is needed to track parentage as far as I
            % can tell
            newPop = population(oldPop.allele_seed, oldPop.numAlive, oldPop.mut_rate1, oldPop.mut_rate2, oldPop.mut_func, oldPop.gen_cut, oldPop.elite_cut, oldPop.fit_func);
            r = oldPop.nodeRanking;
            for i = 1:length(r)
                if oldPop.nodes{r(i)}.alive == true
                    newPop.nodes{i} = copy(oldPop.nodes{r(i)}); % newPop nodes' default order should match the order specified by oldPop.nodeRanking
                end
            end

            for i = 1:length(newPop.popSize)
                if newPop.nodes{i}.alive == true && newPop.nodes{i}.elite == false
                    % Perform mutation
                    child = newPop.nodes{i}.genome.mutate(newPop.mut_func(newPop.mut_rate1, newPop.mut_rate2, newPop.popSize), newPop.allele_seed);
                    newPop.nodes{i}.genome.alleles = child;

                    % Track parentage
                    oldPop.nodes{r(i)}.children = [oldPop.nodes{r(i)}.children, newPop.nodes{i}];
                    newPop.nodes{i}.parent1 = oldPop.nodes{r(i)};
                end
            end

            for i = ceil(length(newPop.popSize)*newPop.elite_cut):2:length(newPop.popSize)-2
                if newPop.nodes{r(i)}.alive == true && newPop.nodes{r(i+1)}.alive == true
                    % Perform crossover
                    [child1, child2] = newPop.nodes{r(i)}.genome.crossover(newPop.nodes{r(i+1)}.genome);
                    newPop.nodes{r(i)}.genome.alleles = child1;
                    newPop.nodes{r(i+1)}.genome.alleles = child2;

                    % Track parentage again
                    oldPop.nodes{r(i)}.children = [oldPop.nodes{r(i)}.children, newPop.nodes{r(i+1)}];
                    newPop.nodes{r(i)}.parent2 = oldPop.nodes{r(i+1)};
                    
                    oldPop.nodes{r(i+1)}.children = [oldPop.nodes{r(i+1)}.children, newPop.nodes{r(i)}];
                    newPop.nodes{r(i+1)}.parent2 = oldPop.nodes{r(i)};
                end
            end
        end
    end
end