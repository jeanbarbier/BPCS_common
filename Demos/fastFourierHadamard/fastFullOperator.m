%% readme
% Demo of fast, memory efficient and quasi-optimal CS reconstruction of a sparse-gauss real
% or complex signal by use of full structured operators.
% For the signal, each component is i.i.d with p(x) = (1 - rho) *
% delta(x) + rho * ComplexNormal(x | mGauss, varGauss).
% The imaginary part is zero in the real case.
% The final estimated signal is results.av_mess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM AND PARAMETERS DEFINITION

close all; clear;

%% problem parameters
type = 'real'; % type of the signal : 'real' (will use an Hadamard based operator) or 'complex' (will use a Fourier based operator)
N = 2^14; % size of the signal (POWER OF 2 if type = 'real')
rho = 0.1; % ..and it's density/sparsity
mGauss = 0; % mean of the real signal or of the real and imaginary parts of the complex one
varGauss = 1; % variance of the real signal or of the real and imaginary parts of the complex one
alpha = 0.3; % measurement rate
trueMat = 1; % want to compare with true random matrix (1/0)? (for size N <= 2^15)
var_noise = 0; % noise variance

%% algorithm parameters
My = CSBP_Solver_Opt();
My.nb_iter = 300; % max number of iterations
My.print = 1; % frequency of printing results
My.conv = 1e-10; % convergence accuracy
My.dump_mes = 0.5; % dumping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO MODIFICATIONS AFTER THIS LINE ARE REQUIRED

%% full matrix required parameters
numBlockC = 1; w = 1; J = 1 ./ N; numBlockL = 1; My.MSEbyBlock = -1; Nblock = N; Mblock = floor(alpha * Nblock);

%% creation of the signal
if (strcmp(type, 'complex') ); S = S_SparseGaussComplex(N, rho, mGauss, varGauss);
else S = S_2Gauss(N, rho, 0, mGauss, 0, varGauss); end

%% lines and signs randomization of the operators
[rp, flipedSigns] = createRandomLinesAndSignsPermutationForOperators(numBlockC, numBlockL, J, Mblock, Nblock);

%% algorithm inputs
if (strcmp(type, 'complex') ); My.prior = 'Complex'; My.method = 'AMPseededFourier'; else My.prior = 'SparseGauss'; My.method = 'AMPseededHadamard'; end
My.M = sum(Mblock);
My.N = N;
My.J = J;
My.numBlockL = numBlockL;
My.numBlockC = numBlockC;
My.Mblock = Mblock;
My.Nblock = Nblock;
My.rp = rp;
My.noBlockError = flipedSigns;
My.signal = S;
My.signal_rho = rho;
My.m_gauss = mGauss;
My.var_gauss = varGauss;
My.var_noise = var_noise;

%% measure
if (strcmp(type, 'complex') ); Y = MultSeededFourier(S, J, numBlockL, numBlockC, Mblock, Nblock, rp, flipedSigns);
else Y = MultSeededHadamard(S, J, numBlockL, numBlockC, Mblock, Nblock, rp, flipedSigns); end
if (My.var_noise > 0); Y = Y + randn(size(Y) ) .* sqrt(My.var_noise); end

%% reconstruction
tic; [results, n_and_e, MSEtAlgo, MSEblockAlgo] = CSBP_Solver(Y, [], My); toc;

%% density evolution
plotDensityEvolutionHadamardFourier(); drawnow; hold on;

%% comparison with random matrix
if ((N <= 2^15) && (trueMat == 1) )
    My.save_memory = 1; My.save_speed = 0;
    if (strcmp(type, 'complex') ); My.method = 'AMPcomplex'; else My.method = 'AMP'; end
    G = createSeededRandomMatrix(type, J, Mblock, Nblock); % creation of the seeded matrix
    Y = G * S.'; % measure
    if (My.var_noise > 0); Y = Y + randn(size(Y) ) .* sqrt(My.var_noise); end
    tic; [results, n_and_eR, MSEtAlgoR, MSEblockAlgoR] = CSBP_Solver(Y, G, My); toc; % reconstruction
    plot_(3) = semilogy(MSEtAlgoR(:), 'om'); % plot
    
    if (strcmp('real', type) ); legend(plot_([1,3,2]), 'Hadamard', 'Random', 'Density evolution');
    else legend(plot_([1,3,2]), 'Fourier', 'Random', 'Density evolution'); end
end