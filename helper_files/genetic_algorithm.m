function [Q, R, root] = genetic_algorithm(x, A, B, segArray)
% GENETIC_ALGORITHM
%   Unique Inputs: popSize = initial population size
%                  mut_rate1 = initial mutation rate
%                  mut_rate2 = final mutation rate
%                      - mutation rate decreases on a power scale from 
%                        mut_rate1 to mut_rate2 as the population decreases
%                  gen_cut = cuttoff point for general population
%                  elite_cut = cuttoff point for elite population (top 1 
%                              solution will always be preserved regardless
%                              of the cuttoff rate)
%
%   General Info: The genetic algorithm has X stages. First, an initial 
%                 population is generated using Bryson's Rule as a basis
%                 then performing mutation and crossover on all but one 
%                 solution to create variation. Second, the dynamics of the
%                 system are modeled in order to evaluate fitness. Third,
%                 the population undergoes the culling/reproduction stage.
%                 The top gen_cut solutions are kept and undergo mutation
%                 and crossover with eachother. The top elite_cut solutions
%                 do not undergo mutation or crossover with any other
%                 solutions.
    MOI = segArray{1};
    constants = [segArray{3}; segArray{2}; segArray{4}];
    xcrit1 = segArray{5}(:, 1);
    xcrit2 = segArray{5}(:, 2);
    limits = segArray{7};
    numPoints = segArray{10};
    ti = segArray{11};
    tf = segArray{12};
    Qbry = segArray{13};
    Rbry = segArray{14};
    popSize = segArray{16};
    mut_rate1 = segArray{17};
    mut_rate2 = segArray{18};
    gen_cut = segArray{19};
    elite_cut = segArray{20};
    mut_func = segArray{21};
    fit_func = segArray{22};
    mut_factor = segArray{23};
    
    pop = population(1, [Qbry; Rbry], popSize, mut_rate1, mut_rate2, mut_factor, mut_func, gen_cut, elite_cut, fit_func, segArray);
    root = node.init_tree(int64(0), pop.nodes);
    
    while pop.popSize > 1
        fprintf("\nPOPULATION: %d\n", pop.generation)
        
        nodes = pop.nodes; % We can pass pop.nodes but not the entire population into the parfor loop for some reason
        % Simulate dynamics for each solution
        %% REPLACE LINEAR SIM WITH 2PBVP
        parfor i = 1:pop.popSize
            lqrFail = false;
            Q = diag(nodes{i}.genome.alleles(1:size(x,1)));
            R = diag(nodes{i}.genome.alleles(size(x,1) + 1:end));
            
            %Optimal control gain matrix K, solution S, and poles P
            try
                [Ksegment, ~, ~] = lqr(A, B, Q, R);
            catch
                fprintf("LQR Fail   ")
                lqrFail = true;
            end
            nodes{i}.genome.lqrFail = lqrFail;
    
            if ~lqrFail
                %Simulate using input data
                %Utilizes dynamics' translational symmetry to approach critical points
                [xsegment, ~, ~, odeStopped] = simulate(ti, tf, numPoints, Ksegment, constants, MOI, xcrit1-xcrit2, limits); 
                nodes{i}.genome.states = xsegment;
                nodes{i}.genome.odeStopped = odeStopped

                if (odeStopped == true) 
                    fprintf("ODE STOPPED    "); 
                else
                    fprintf("ODE SUCCESS    ");
                end
            end
        end
        pop.nodes = nodes;
        
        % Remaining GA processes
        pop.fitEval;
        pop.kill;
        pop = pop.reproduce;
    end
    
    Q = diag(pop.nodes{1}.genome.alleles(1:size(x,1)));
    R = diag(pop.nodes{1}.genome.alleles(size(x,1) + 1:end));
    root.genome = pop.nodes{1}.genome;
end