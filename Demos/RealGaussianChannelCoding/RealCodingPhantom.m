%% readme
% Demo related to the article about robust error-correction by J.Barbier, F.Krzakala, P.Zhang and L.Zdeborova : http://arxiv.org/pdf/1304.6599.pdf
% It demonstrates how a picture (Phantom) can be sent through a noisy channel and
% decoded by reconstruction of the noise by use of AMP. The noise model is
% a background noise (that affects all the components) with a big Gaussian that corrupts a fraction rho of the components : p(z) = Gaussian(0, var_small) + rho * Gaussian(0, 1).
% AMP decoding is compared with what one would obtain without the big
% gaussian in the noise but only the small background, it's the ideal
% decoding, the L1 decoding and the naive one, which is
% directly estimating the signal without taking care of the noise.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM AND PARAMETERS DEFINITION

clear all; close all;

%% problem parameters
N = 32.^2; % size of the phantom (power of 2, N >= 64^2 becomes very demanding in memory)
hadamard = 1; % do you want to use hadamard matrix (way more fast) or gaussian one?
seeded = 1; % seeded or homogeneous matrix (if seeded = 0)? for both there are random and hadamard matrices
rho_corrupt = 0.1; % density of the big noise
alphaGlobal = 0.25; % measurement rate
var_small = 1e-6; % variance of the small additive gaussian noise
if (seeded == 1)
    JJ = 0.1.^2; % variance parameter between the blocks
    w = 1; % number of sub-diagonal blocks
    numBlockC = 8; % number of blocks for the columns of the seeding matrix (POWER OF 2)
    numBlockL = numBlockC + w - 1; % number of blocks for the rows of the seeding matrix
    alphaCs(1) = 3 .* rho_corrupt; % measurement rate seed
    alphaCs(2 : numBlockL) = (numBlockC .* alphaGlobal - alphaCs(1) ) ./ (numBlockL - 1); % measurement rate bulk blocks
else
    alphaCs(1) = alphaGlobal; JJ = []; numBlockC = 1; numBlockL = 1; w = 0; % parameters for non-seeded matrix
end

%% algorithm properties
nb_iter = 300; % maximum number of iterations
print = 1; % printing frequency
dump_mes = 0.5; % dumping
save_ = 0; % if you want to save the results, save_ = 1

%% do you want to use l1magic pack (that must be added to the path) or the inluded l1 AMP solver?
% default is zero (AMP l1 solver), be sure to download l1_magic and put it in the path if you use it: http://users.ece.gatech.edu/~justin/l1magic/
l1_magic = 0; l1_AMP = 1 - l1_magic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% NO MODIFICATIONS AFTER THIS LINE ARE REQUIRED

if(numBlockC > 1); if (alphaCs(1) < alphaCs(2) ); error('The first block must have a measurement rate bigger or equal to the one of the bulk blocks'); end; end

%% Initialisation
M = ceil(N ./ (1 - alphaGlobal) ); % number of sent components through the noisy channel
if (hadamard == 1)
    if (seeded == 0)
        NN = 2^(ceil(log(M) / log(2) ) );  % number of variables extended to the first next power of two
        Nblock = NN; % number of variables per block in the seeded matrix (POWER OF 2)
        Mblock = floor(alphaGlobal * Nblock); % number of measurements
    else
        NN = 2^(ceil(log(M) / log(2) ) ); % number of variables extended to the first next power of two
        Nblock = 2^(ceil(log(NN / numBlockC) / log(2) ) ); % number of variables per block in the seeded matrix (POWER OF 2)
        Mblock = floor(alphaCs * Nblock); % number of lines per block
    end
else
    if (seeded == 0)
        NN = M;
        Nblock = NN;
        Mblock = floor(alphaGlobal * Nblock); % number of lines per block
    else
        NN = ceil(M ./ numBlockC); % number of variables per block
        Nblock = NN; % number of variables per block in the seeded matrix
        Mblock = floor(alphaCs * Nblock); % number of lines per block
    end
    
end
J = createSeededJ(numBlockL, numBlockC, JJ, w, M); % block variance matrix
[randomPermutations, flipedSigns] = createRandomLinesAndSignsPermutationForOperators(numBlockC, numBlockL, J, Mblock, Nblock); % lines and signs randomization of the operators
Phantom = phantom(sqrt(N) ); % image loading
I_L1 = zeros(size(Phantom) ); I_CS = zeros(size(Phantom) ); I_id = zeros(size(Phantom) );
if (hadamard == 0); F = createSeededRandomMatrix('real', J, Mblock, Nblock); else F = createSeededHadamardMat(Mblock, Nblock, numBlockL, numBlockC, J, randomPermutations, flipedSigns); end % Generate the measurement matrix that will measure the noise only (= null space of the coding one)
% if (max(size(F) ) - min(size(F) ) < max(size(Phantom) ) );

disp('Lena decoding');
disp('true gamma'); disp(1 / (1 - sum(Mblock) / (Nblock .* numBlockC) ) );
disp('true alpha'); disp(sum(Mblock) / (Nblock .* numBlockC) );

%% Generate the signal and the perfect measure
S = Phantom; % original picture
rp = randperm(N); S = reshape(S, N, 1); S_mixed = S(rp); % randomization of the signal
A = null(F).'; % coding matrix
comp = 0; while (min(size(A) ) > max(size(S_mixed) ) ); S_mixed = [S_mixed; randn]; comp = comp + 1; end % add components to the signal to fulfill size constraints
Y_perfect = A.' * S_mixed; MM = max(size(Y_perfect) ); % perfect measure for comparison (non-noisy)

%% Error due to noisy channel transmission
k = floor(rho_corrupt .* MM);
rp2 = randperm(MM); big_error = [[randn(k, 1); zeros(MM - k, 1)] ]; error_ = big_error(rp2); big_error = error_;
small_error = randn(MM, 1) .* sqrt(var_small);
Y_err = Y_perfect + small_error + big_error;
Y_null = F * Y_err; % noisy measure

%% Ideal coding and decoding
Y_id = Y_perfect + small_error;
X_id = A.' \ Y_id; X_id_(rp) = X_id(1 : max(size(X_id) ) - comp, 1); X_id = X_id_.';

%% Naive decoding (with no error correction)
X_naive = A.' \ Y_err; X_naive_(rp) = X_naive(1:max(size(X_naive) ) - comp, 1); X_naive = X_naive_.';

%% BPCS coding and decoding by reconstruction of the sparse error vector
% algorithm properties
My = CSBP_Solver_Opt();
if (hadamard == 1); My.method = 'AMPseededHadamard'; else My.method = 'AMP'; end
MSEbyBlock = 1;
My.conv = var_small .* 10;
My.signal_rho = rho_corrupt;
My.signal = big_error + small_error;
My.MSEbyBlock = MSEbyBlock;
My.noBlockError = flipedSigns;
My.N = NN;
My.M = sum(Mblock);
My.J = J;
My.numBlockL = numBlockL;
My.numBlockC = numBlockC;
My.Mblock = Mblock;
My.Nblock = Nblock;
My.rp = randomPermutations;
My.alphaBig = 0;
My.nb_iter = nb_iter;
My.print = print;
My.dump_mes = dump_mes;
% 2Gauss
My.m_2_gauss = 0;
My.var_1_gauss = var_small;
My.var_2_gauss = 1 + var_small;

%% AMP decoding with 2Gauss prior
My.prior = '2Gauss';
if (hadamard == 0); [results2G, n_and_e] = CSBP_Solver(Y_null, F, My); else [results2G, n_and_e] = CSBP_Solver(Y_null, [], My); end;
Y_CS = Y_err - results2G.av_mess';
X_CS = A.' \ Y_CS; X_CS_(rp) = X_CS(1:max(size(X_CS) ) - comp, 1); X_CS = X_CS_.';

%% L1
if (l1_AMP == 1)
    My.prior = 'L1';
    if (hadamard == 0); [resultsL1, n_and_e] = CSBP_Solver(Y_null, F, My); else [resultsL1, n_and_e] = CSBP_Solver(Y_null, [], My); end;
    Y_L1 = Y_err - resultsL1.av_mess';
    X_L1 = A.' \ Y_L1; X_L1_(rp) = X_L1(1 : max(size(X_L1) ) - comp, 1); X_L1 = X_L1_.';
else
    err_L1 = l1dantzig_pd(F * Y_null, F', [], Y_null, 3e-3, 1e-2, 1000); % Dantzig selector
    err_index = abs(err_L1) > sqrt(var_small .* M); F_t = F'; % regression step
    F_index_t = F_t(:,err_index);
    err_L1 = zeros(size(Y_err) ); err_L1(err_index) = F_index_t \ Y_null;
    Y_L1 = Y_err - err_L1;
    X_L1 = A \ Y_L1; X_L1_(rp) = X_L1(1 : max(size(X_L1) ) - comp, 1); X_L1 = X_L1_.';
end

%% Results
S(rp) = S(1 : max(size(S_mixed) ) - comp, 1);
rho_id_CS = norm(X_CS - S, 2) ./ norm(X_id - S, 2);
rho_id_L1 = norm(X_L1 - S, 2) ./ norm(X_id - S, 2);
mse_id = mean((X_id - S).^2);
mse_CS = mean((X_CS - S).^2);
mse_L1 = mean((X_L1 - S).^2);

results_(1) = {reshape(X_id, sqrt(N), sqrt(N) )};
results_(2) = {reshape(X_CS, sqrt(N), sqrt(N) )};
results_(3) = {reshape(X_L1, sqrt(N), sqrt(N) )};
results_(4) = {[mse_id, mse_CS, mse_L1, rho_id_CS, rho_id_L1]};

I_CS = reshape(X_CS, sqrt(N), sqrt(N) );
I_L1 = reshape(X_L1, sqrt(N), sqrt(N) );
I_id = reshape(X_id, sqrt(N), sqrt(N) );
I_naive = reshape(X_naive, sqrt(N), sqrt(N) );

%% Plots
clf
subplot(3, 2, 1); imshow(Phantom); xlabel('Original Image');
subplot(3, 2, 2); imshow(I_id); xlabel('Ideal');
subplot(3, 2, 3); imshow(I_CS); xlabel('AMP 2Gauss');
subplot(3, 2, 4); imshow(I_L1); xlabel('L1');
subplot(3, 2, 5); imshow(I_naive); xlabel('No error correction');

%% Save
if (save_ == 1); save('Phantom_decoding','I_id','I_CS','I_L1','I_naive','results_'); end