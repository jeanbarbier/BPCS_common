% initialisation of all the quantities

if (strcmpi(opt.method,'AMP') || (strcmpi(opt.method,'AMPhb') ) || (strcmpi(opt.method,'AMPcomplex') ) || (strcmpi(opt.method,'AMPhadamard') ) || (strcmpi(opt.method,'AMPseededHadamard') ) || (strcmpi(opt.method,'AMPseededHadamardTranspose') ) || (strcmpi(opt.method,'AMPseededFourier') ) )
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = ones(1,M);
    R_init = zeros(1,N); S2_init = zeros(1,N); av_mess_init = zeros(1,N); var_mess_init = opt.signal_rho .* ones(1,N);
end

if (strcmpi(opt.method,'AMPseededHadamardTransposeA') )
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = ones(1,N);
    R_init = zeros(1,M); S2_init = zeros(1,M); av_mess_init = zeros(1,M); var_mess_init = opt.signal_rho .* ones(1,M);
end

if (strcmpi(opt.method,'AMPh') )
    if (opt.remove_mean == 0); W = Y; else W = Yeff; end; V = 1;
    R_init = zeros(1,N); S2_init = zeros(1,N); av_mess_init = zeros(1,N); var_mess_init = opt.signal_rho .* ones(1,N);
end