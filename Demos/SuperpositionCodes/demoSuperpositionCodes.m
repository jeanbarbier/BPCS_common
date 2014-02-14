%% readme
% Demo of the AMP decoder for the superposition codes : a binary message S
% of size N = L * B made of L sections of size B, each of them having only one 1
% inside it is sent through the Additive White Gaussian Noise channel.
% The section error rate is the fraction of incorrectly decoded sections.
%
% In this demo (if ampSingleInstance = 1), a signal with given B and L is randomly generated,
% coded by multiplication with a random Gaussian matrix made of i.i.d with zero mean
% and variance 1 / L (to ensure a signal power ~ 1), corrupted by gaussian
% noise of variance 1 / snr (snr is the signal to noise ratio) and then decoded by AMP.
% You can track the evolution of the ser (error).
%
% In addition (if ampDensityEvolution = 1), the asymptotic ser (L -> infinity) of AMP can be
% predicted as well by Density Evolution.

close all; clear all;

%% parameters
ampSingleInstance = 1; % do you want to decode one single random instance using AMP ? 1 if yes, 0 for no
ampDensityEvolution = 1; % do you want to get the asymptotic section error rate of AMP using the Density Evolution technic ? 1 if yes, 0 for no
rate = 1.4; % rate in bits/channel use
B = 2^2; % number of variables per section
L = 5000; % number of sections
snr = 15; % signal to noise ratio
nbIterMax = 50; % maximum number of iterations for the algorithm/Density Evolution
sizeMonteCarlo = 1e5; % sample size for the estimation of the integrals by monte carlo in the Density Evolution
accuracyDecoding = 1e-3; % threshold of the section error rate for convergence
printFrequency = 1; % printing frequency of the section error rate

%% print informations
alpha = log2(B) / (rate * B); % measurement rate
Delta = 1 / snr; % noise variance
capacity = log2(1 + snr) / 2; % Shannon capacity of the channel
pr = [capacity  rate]; disp('Capacity  Rate'); disp(pr);

%%%%%%%%%%%%%%%%%% no modifications required after this line %%%%%%%%%%%%%%%%%%

%% single instance decoding by AMP
if (ampSingleInstance)
    N = B * L; % total signal size
    if (N > 2^15); error('too big N, the memory of the computer can saturate'); end % error if too big size, comment this line for bigger sizes
    word = ceil(B * rand(1, L) ); % select the support of the message
    S = zeros(N, 1); for i = 1 : L; pos = word(i); S((i - 1) * B + pos) = 1; end % create the message
    G = randn(ceil(alpha * N), N) / sqrt(L); % the measurement matrix (scaled such that the power of the signal is ~ 1)
    Y = G * S; % the encoding step
    Y = Y + sqrt(Delta) * randn(size(Y) ); % additive white gaussian noise
    
    % properties of the algorithm
    My = CSBP_Solver_Opt(); % default values of the option field
    My.nb_iter = nbIterMax;
    My.conv = accuracyDecoding;
    My.signal = S;
    My.dump_mes = 0.;
    My.print = printFrequency;
    My.var_noise = Delta;
    My.signal_rho = 1 / B;
    My.method = 'AMP';
    My.prior = 'SuperpositionCode';
    My.NbSC = B;
    My.nonZeroValues = ones(1, L);
    
    tic; [X, dumb, dum2, serSingleInstance] = CSBP_Solver(Y', G, My); toc; % reconstruction
end


%% asymptotic behavior of AMP by Density Evolution
if (ampDensityEvolution)
    nbIterMaxDE = nbIterMax; if (ampSingleInstance); nbIterMaxDE = max(size(serSingleInstance) ); end
    Ebiased = 1 / B; serDensityEvolution(1) = 1;
    for t = 2 : nbIterMaxDE;
        [Ebiased, serDE] = DensEvoSuperpositionCode(B, alpha, Delta, Ebiased, sizeMonteCarlo, 'sectionErrorRate');
        serDensityEvolution(t) = serDE;
    end
end

%% plots
if (ampDensityEvolution && ampSingleInstance == 0); semilogy(max(serDensityEvolution, 1e-50),'ro-'); legend('Asymptotic behavior');
elseif (ampDensityEvolution == 0 && ampSingleInstance); semilogy(max(serSingleInstance, 1e-50),'bo-'); legend('Single instance');
elseif (ampDensityEvolution && ampSingleInstance); semilogy(max(serSingleInstance, 1e-50),'bo-'); hold on;  semilogy(max(serDensityEvolution, 1e-50),'ro-'); legend('Single instance', 'Asymptotic behavior'); end
title(['rate = ' num2str(rate) ', snr = ' num2str(snr) ]);