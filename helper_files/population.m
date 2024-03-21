classdef population < handle
    %POPULATION
    %   Collection of genes used in genetic algorithm

    properties
        allele_seed
        popSize
        mut_rate1
        mut_rate2
        mut_func
        gen_cut
        elite_cut = 0
        nodes
    end

    methods
        function obj = population(allele_seed, popSize, mut_rate1, mut_rate2, mut_func, gen_cut, elite_cut)
        % CONSTRUCTOR
        %   Constructs 

            % Add user defined params
            obj.allele_seed = allele_seed;
            obj.popSize = popSize;
            obj.mut_rate1 = mut_rate1;
            obj.mut_rate2 = mut_rate2;
            obj.mut_func = mut_func;
            obj.gen_cut = gen_cut;
            obj.elite_cut = elite_cut;

            % Generate nodes
            for i = 1:obj.popSize
                % Create node
                obj.nodes{i} = node([], [], 1, 0, true, gene(allele_seed), []);

            end

            % Initial mutation and crossover
            for i = 1:length(obj.nodes)
                obj.nodes{i}.genome.mutate(mut_func(mut_rate1, mut_rate2, popSize), allele_seed)
                
            end

            for i = 1:length(obj.nodes)-1
                obj.nodes{i}.genome.crossover(obj.nodes{i+1}.genome)
            end
        end

        function reproduce(obj)
        % REPRODUCE
        %   Creates new population through mutation and crossover. Tracks
        %   lineage.
            for i = 1:length(obj.nodes)
                obj.nodes(i).genome.mutate(obj.mut_func(obj.mut_rate1, obj.mut_rate2, obj.popSize), obj.allele_seed)
                
            end

            for i = 1:length(obj.nodes)-1
                obj.nodes(i).genome.crossover(obj.nodes(i+1).genome)
            end
        end
    end
end