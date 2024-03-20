function writeTrajectory(filename, K_set, x_set, u_set, t_set)
    % TODO: add quick stabilization matrix
    fp = fopen(filename,'w');
    [n, k] = size(x_set);
    K_size = size(K_set);
    m = K_size(1);
    N = K_size(2) - n;
    p = K_size(3);
    % k, p, m, n and then N on one line
    fprintf(fp, "%d,%d,%d,%d,%d\n", k, p, m, n, N);
    
    % unfortunately the indexing is different than the one in C
    % so we have to do this
    K_set_p = permute(K_set, [2 1 3]); % N+n, then m, then p, fast to slow
    x_set_p = permute(x_set, [2 1]);
    u_set_p = permute(u_set, [2 1]);
    fprintf(fp, "%f,", K_set_p);
    % add q matrices here eventually
    fprintf(fp, "%f,", x_set_p);
    fprintf(fp, "%f,", u_set_p);
    fprintf(fp, "%f,", t_set);

    fclose(fp);
end

