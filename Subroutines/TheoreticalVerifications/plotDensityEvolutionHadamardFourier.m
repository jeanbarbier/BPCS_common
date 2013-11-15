% plot of the full or seeded density evolution compared with the algorithm MSE(t)

if (numBlockC == 1)
    MSEtheoric(1) = rho;
    
    for (t = 2 : max(size(MSEtAlgo) ) )
        if (strcmp('real', type) ); [E_, V_, MSEtheoric(t), V(t) ] = DensEvoSparseGauss('full', 0, 1, rho, sum(Mblock) / N, My.var_noise, MSEtheoric(t - 1), MSEtheoric(t - 1) );
        else [MSEtheoric(t) ] = DensEvoSparseGaussComplex('full', rho, sum(Mblock) ./ N, My.var_noise, MSEtheoric(t - 1) ); end
    end
    
    if (strcmp('real', type) ); plot_(1) = semilogy(MSEtAlgo, 'or'); hold on; plot_(2) = semilogy(MSEtheoric, '-k'); hold off;
    else plot_(1) = semilogy(MSEtAlgo, 'ob'); hold on; plot_(2) = semilogy(MSEtheoric, '-k'); hold off; end
    
else
    MSEtheoric{1}(1 : numBlockC) = rho;
    
    for c = 1 : numBlockC; blocAlgo(1, c) = MSEblockAlgo{1}(c); blocTh(1, c) = MSEtheoric{1}(c); end
    
    for (t = 2 : max(size(MSEtAlgo) ) )
        if (strcmp('real', type) ) [E_, V_, MSEtheoric{t}, V(t) ] = DensEvoSparseGauss('seeded', 0, 1, rho, sum(Mblock) / N, My.var_noise, MSEtheoric{t - 1}, MSEtheoric{t - 1}, numBlockC, numBlockL, Nblock, Mblock, N, J);
        else [MSEtheoric{t} ] = DensEvoSparseGaussComplex('seeded', rho, sum(Mblock) / N, My.var_noise, MSEtheoric{t - 1}, numBlockC, numBlockL, Nblock, Mblock, N, J); end
        
        for c = 1 : numBlockC; blocAlgo(t, c) = MSEblockAlgo{t}(c); blocTh(t, c) = MSEtheoric{t}(c); end
    end
    
    for c = 1 : numBlockC; plot_(1) = semilogy(blocAlgo(:, c), 'or'); hold on; plot_(2) = semilogy(blocTh(:, c), '-k'); end
    
    hold off;
end