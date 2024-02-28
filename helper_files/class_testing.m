dbstop in population at 25
clear

pop = population;
pop.allele_seed = [10,100,1000];
pop.mut_rate1 = 1;
pop.mut_rate2 = 1;
pop.mut_func = @mfunc;
pop.popSize = 3;

pop.generate_genes;
pop.reproduce

function mrate = mfunc(mrate1, mrate2, popSize)
    mrate = mrate2 + ((mrate1 - mrate2)/2^popSize)*2^popSize;

end