clear

allele_seed = [10,100,1000];
mut_rate1 = 1;
mut_rate2 = 1;
mut_func = @mfunc;
gen_cut = 0.5;
elite_cut = 0;
popSize = 3;

pop = population(allele_seed, popSize, mut_rate1, mut_rate2, mut_func, gen_cut, elite_cut);

function mrate = mfunc(mrate1, mrate2, popSize)
    mrate = mrate2 + ((mrate1 - mrate2)/2^popSize)*2^popSize;

end