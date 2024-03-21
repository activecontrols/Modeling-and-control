% push to helperfiles under traj qr opt branch in modeling and control
classdef node
    properties
        parent1 % node type
        parent2 % node type
        generation % int
        batch % numeric
        alive = true; % boolean
        genome % gene type
        children % array of nodes        
        % if root node:
        % generation = 0
        % genome = most fit in the tree
        % children = all of the starting population
    end
    methods
        % constructor
        function obj = node(parent1, parent2, generation, batch, alive, genome, children)
            arguments
                parent1 
                parent2 
                generation int32 {mustBeInteger}
                batch double {mustBeNumeric}
                alive logical {mustBeNumericOrLogical}
                genome
                children
            end
            obj.parent1 = parent1;
            obj.parent2 = parent2; 
            obj.generation = generation; 
            obj.batch = batch;
            obj.alive = alive;
            obj.genome = genome; 
            obj.children = children;
        end
        % initializes a tree
        function root = init_tree(batch_num, starting_pop)
            arguments
                batch_num class {mustBeNumeric}
                starting_pop % array of genes
            end
            %conv staring pop to array of nodes
            root = node([], [], 0, batch_num, true, [], starting_pop);
        end
        % adds a new generation to the tree
        % assumes last_gen is already in tree
        function new = new_gen(last_gen, new_gen)
            arguments
                last_gen % array of nodes
                new_gen % array of genes
            end
            batch_num = last_gen(0).batch_num;
            generation = 1 + last_gen(0).generation;
            parent_num1 = 0;
            parent_num2 = 1;
            new = zeros(length(new_gen));
            for i = 1:length(new_gen)
                while (~last_gen(parent_num1).alive)
                    parent_num1 = parent_num1 + 1;
                    parent_num2 = parent_num2 + 1;
                end
                while (~last_gen(parent_num2).alive)
                    parent_num2 = parent_num2 + 1;
                end
                new(i) = node(last_gen(parent_num1), last_gen(parent_num2), generation, batch_num, true, new_gen(i), []);
            end
        end
        % set final genome? to root
        % find lca (last common ancestor) of two given nodes
        function lca = LCA(node1, node2)
            arguments
                node1 
                node2
            end
            if (node1.generation > node2.generation)
                temp = node1;
                node1 = node2;
                node2 = temp;
            end
            ancestry1 = {node1};
            ancestry2 = {};
            while (node1.generation ~= node2.generation)
                temp = {};
                for i = 1:length(ancestry2)
                    if ~ismember(ancestry2(i).parent1, temp)
                        temp.append(ancestry2(i).parent1)
                    end
                    if ~ismember(ancestry2(i).parent2, temp)
                        temp.append(ancestry2(i).parent2)
                    end    
                end
                ancestry2 = temp;
                node2 = temp(i);
            end
            while (node1.generation ~= 0)
                %edit
                for i = 1:length(ancestry2)
                    if ~ismember(ancestry2(i).parent1, temp)
                        temp.append(ancestry2(i).parent1)
                    end
                    if ~ismember(ancestry2(i).parent2, temp)
                        temp.append(ancestry2(i).parent2)
                    end    
                end
            end
            if (node1.generation == 0)
                lca = node1;
            end
        end
        % traverse array + apply func
        function traverse(node, func)
            arguments
                node class {mustBeA(node)}
                func class {mustBeFunction}
            end
            disp(node.genome);
            func(node);
            for i = 1:length(node.children)
                traverse(node.children(i));
            end
        end
    end
end