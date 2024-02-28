classdef population < handle
    %POPULATION
    %   Collection of genes used in genetic algorithm

    properties
        allele_seed
        popSize = 100
        mut_rate1 = 0.5
        mut_rate2 = 0.5
        mut_func
        gen_cut = 0.5
        elite_cut = 0
        genes
    end

    methods
        function generate_genes(obj)
        %GENERATE_GENES
        %   Populates genes property with values based on an initial seed.
            for i = 1:obj.popSize
                obj.genes{i} = gene;
                obj.genes{i}.alleles = obj.allele_seed;
            end
        end

        function reproduce(obj)
            % work with nodes instead so child-parent relationships established continuously
            for i = 1:length(obj.genes)
                obj.genes{i}.mutate(obj.mut_func(obj.mut_rate1, obj.mut_rate2, obj.popSize), obj.allele_seed)
            end

            for i = 1:length(obj.genes)-1
                obj.genes{i}.crossover(obj.genes{i+1})
            end
        end
    end
end