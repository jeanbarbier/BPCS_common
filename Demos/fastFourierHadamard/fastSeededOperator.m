%% readme
% Demo of fast, memory efficient and quasi-optimal CS reconstruction of a sparse-gauss real
% or complex signal by use of seeded (spatially coupled) structured operators.
% For the signal, each component is i.i.d with p(x) = (1 - rho) *
% delta(x) + rho * ComplexNormal(x | mGauss, varGauss).
% The imaginary part is zero in the real case.
% The final estimated signal is results.av_mess


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM AND PARAMETERS DEFINITION

if(numBlockC > 1); if (alphaCs(1) < alphaCs(2) ); error('The first block must have a measurement rate bigger or equal to the one of the bulk blocks'); end; end

close all; clear; load line_BP.mat;

%% problem parameters
type = 'real'; % type of the signal : 'real' (will use an Hadamard based operator) or 'complex' (will use a Fourier based operator)
N = 2^15; % size of the signal (POWER OF 2 if type = 'real')
rho = 0.1; % ..and it's density/sparsity
[a, closest] = min(abs(rho - line_BP(1,:) ) );
mGauss = 0; % mean of the real signal or of the real and imaginary parts of the complex one
varGauss = 1; % variance of the real signal or of the real and imaginary parts of the complex one
alphaGlobal = 0.2; % global measurement rate (the true rate will be a little bit different)
numBlockC = 8; % number of blocks (must be >= 2) for the columns of the seeding matrix, it must divide N (POWER OF 2 if type = 'real')
w = 1; % coupling window (number of sub-diagonal blocks)
JJ = 0.1.^2; % coupling variance (not rescaled by N, this will be done automatically)
numBlockL = numBlockC + w - 1; % number of blocks for the rows of the seeding matrix
alphaCs(1) = 0.1 + line_BP(2, closest); % measurement rate 1st block/seed, taken into acount if numBlockC > 1 : can be modified
alphaCs(2 : numBlockL) = (numBlockC .* alphaGlobal - alphaCs(1) ) ./ (numBlockL - 1); % measurement rate of the bulk blocks : not to be modified
trueMat = 1; % want to compare with true random matrix (1/0)? (for size N <= 2^15)
var_noise = 0; % noise variance

%% algorithm parameters
My = CSBP_Solver_Opt();
My.nb_iter = 300; % max number of iterations
My.print = 1; % frequency of printing results
My.conv = 1e-6; % convergence accuracy
My.dump_mes = 0.5; % dumping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO MODIFICATIONS AFTER THIS LINE ARE REQUIRED

%% creation of the signal
if (strcmp(type, 'complex') ); S = S_SparseGaussComplex(N, rho, mGauss, varGauss);
else S = S_2Gauss(N, rho, 0, mGauss, 0, varGauss); end

%% matrix/operator block sizes
Nblock = N / numBlockC; % number of variables per block in the seeded matrix, must divide N (POWER OF 2 if type = 'real')
for (i = 1 : numBlockL) Mblock(i) = floor(alphaCs(i) * Nblock); end % number of lines of the blocks in the seeded matrix
disp('true measurement rate'); disp(sum(Mblock) / N);

%% block variance matrix
J = createSeededJ(numBlockL, numBlockC, JJ, w, N);

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
My.MSEbyBlock = 1;
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
    for c = 1 : numBlockC; blocAlgoR(1, c) = MSEblockAlgoR{1}(c); end % plot
    for (t = 1 : max(size(MSEtAlgoR) ) )
        for c = 1 : numBlockC; blocAlgoR(t, c) = MSEblockAlgoR{t}(c); end
    end
    for c = 1 : numBlockC; plot_(3) = semilogy(blocAlgoR(:, c), 'om'); hold on; end
    
    if (strcmp('real', type) ); legend(plot_([1,3,2]), 'Hadamard', 'Random', 'Density evolution');
    else legend(plot_([1,3,2]), 'Fourier', 'Random', 'Density evolution'); end
end