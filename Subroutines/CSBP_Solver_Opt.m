% EMBP_partial_Opt :  Function to set EMBP_Solver_Opt to their default values
%
%   Details of the option :
%
%   .method             method of resolution of the inference problem (see
%                       the article by Krzakala and al.). It can be 'AMP', 'AMPtap' for statistically
%                       homogeneous measurements matrices, 'AMPcomplex' for complex numbers, 'AMPseededFourier' for Fourier seeded operator
%                       or 'AMPseededHadamard' for Hadamard seeded operator. ['AMP']
%
%   .save_memory        option that does not allow to create new objects of the
%                       size comparable with the measurement matrix during
%                       the procedure, such that the memory required by the
%                       routine is lowered. [1]
%
%   .save_speed         option that allows pre-computation of new objects (just presents for some methods) of the
%                       size comparable with the measurement matrix, such that the memory required by the
%                       routine is multiplied by a factor 3 to 4 times the one required to store the measurement matrix
%                       but the algorithm is speeded-up. It is usefull for not too big matrices or high memory capacities. [0]
%                       WARNING: save_memory and save_speed are not
%                       compatible: One must be equal to 1, the other to 0.
%
%   .nb_iter            max number of iterations. [default : 1000]
%
%   .print              print results every opt.print iterations (0 -> never). [10]
%
%   .conv               convergence criterion. [1e-8]
%
%   .learn              learning parameters or not. [0]
%
%   .signal_rho         first estimate of the density of the signal. [M ./ (10 .* N)]
%
%   .var_noise          first estimate of the variance vector of the noise. [1e-10]
%
%   .signal             a solution to compare to while running.
%
%   .dump_learn         dumping coefficient of the learning. [0]
%
%   .dump_mes           dumping of the messages. [0.5]
%
%   .option_noise       Value is 1/~=1 [0, i.e not activated]. If activated,
%                       it deals with matrix uncertainty.
%
%   .remove_mean        When activited (1 or 2 instead of 0), this
%                       allows the algorithm to deal with non zero
%                       mean matrices that can have a different mean on every column (1).
%                       With value 2 it assumes the same average value for each of them.
%                       If provided (see below), it uses the one given by the user. [0]
%
%   .Gmean              see .remove_mean
%
%   .Ymean              see .remove_mean
%
%   .Nvec               This can be used to specify a structure
%                       to the signal, and change the printing output. [1]
%
%   .Mvec               This can be used to specify a structure
%                       to the solution vector Y (this is used in
%                       Active option.noise for instance).
%
%   .varG               Used if AMPtap is activated (for bloc matrices). It
%                       can be given in two different forms : Or it is a full sparse [M N] matrix where the element varG(i,j) is the variance of
%                       the bloc to which G(i,j) belong to, or it a small [L C]
%                       matrix where the element varG(i,j) is the variance
%                       of the bloc (i,j). In the second case, the user
%                       must also give the opt.Mvec_bm and opt.Nve_bm.
%
%   .Nvec_bm            Used when AMPtap is activated. It is a vector containing the number of columns in each block.
%
%   .Mvec_bm            Used when AMPtap is activated. It is a vector containing the number of lines in each block.
%
%   .alphaBig           Must be equal to one if alpha (measurement rate >=1). [0]
%
%   .adaptDump          Adaptative dumping (for the moment, only for 'SparseGauss' or '2Gauss' prior with 'AMP' method). [0]
%
%   .prior              prior on the data. ['SparseGauss']
%                       'SparseGauss' : Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var_gauss) * exp(-(x - m_gauss)^2 / (2 * var_gauss) ) : param_1 = m_gauss; param_2 = var_gauss;
%                       'SparseGaussCut' : Gaussian sparse prior enforcing value inside a symetric interval : p(x) ~ [(1 - rho) * delta(x) + rho / sqrt(2 * pi * var_gauss) * exp(-(x - m_gauss)^2 / (2 * var_gauss) )] * I(|x| < cut) : param_1 = m_gauss; param_2 = var_gauss; param_3 = cut;
%                       'SparseGaussPositive' : Positive Gaussian sparse prior : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var_gauss) * exp(-(x - m_gauss)^2 / (2 * var) ) * I(x > 0) : param_1 = m_gauss; param_2 = var_gauss;
%                       '2Gauss' : Mixture of two gaussians : p(x) ~ (1 - rho) * exp(-(x - m_1)^2 / (2 * var_1) ) / sqrt(2 * pi * var_1) + rho * exp(-(x - m_2)^2 / (2 * var_2) ) / (sqrt(2 * pi * var_2) ) : param_1 = m_1; param_2 = m_2; param_3 = var_1; param_4 = var_2;
%                       'SparseBinary' : Binary prior : p(x) ~ (1 - rho) * delta(x) + rho * delta(x - 1);
%                       'SparseExponential' : Exponential sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(x > 0) * exp(-expo * x), expo > 0 : param_1 = expo;
%                       'SparseConstant' : Unity inside a finite interval sparse prior : p(x) ~ (1 - rho) * delta(x) + rho * I(c_down < x < c_up) : param_1 = c_down; param_2 = c_up;
%                       'Laplace' : Laplace prior : p(x) ~ 2 / beta * exp{-beta * |x|} : param_1 = beta;
%                       'L1' : L1 optimization (soft tresholding) : p(x) ~ lim_{beta -> infinity} exp{-beta * |x|}, where the x values are bounded by [min, max] : param_1 = min; param_2 = max;
%                       'Binary1' : Plus or minus one prior with fraction half of each : p(x) ~ 0.5 * delta(x - 1) + 0.5 * delta(x + 1);
%                       'Complex' : Complex prior (x complex number) : p(x) ~ (1 - rho) * delta(x) + rho / sqrt(2 * pi * var_gauss) * exp(-(|x| - m_gauss)^2 / (2 * var_gauss) ) : param_1 = m_gauss; param_2 = var_gauss;
%
%   .m_gauss            Mean of the gaussian for the prior 'SparseGauss' [0].
%   .var_gauss          Variance .. [1.]
%
%   .m_gaussP           Mean of the gaussian for the prior 'SparseGaussPositive' [0].
%   .var_gaussP         Variance .. [1.]
%
%   .expo               Exposent for the prior 'SparseExponential' [1];
%
%   .c_down             Lower bound of the signal for the prior 'SparseConstant'. [0]
%   .c_up               Upper bound .. [1]
%
%   .m_gaussC           Mean of the gaussian for the prior 'SparseGaussCut'. [0]
%   .var_gaussC         Variance .. [1]
%   .Cut                Cut/bound for the values of the x. [3]
%
%   .m_2_gauss          Mean of the second gaussian in the prior '2Gauss', so the one of big components of density rho. [0]
%   .var_2_gauss        Variance .. [1]
%   .var_1_gauss        Variance of the small gaussian of mean 0 that models the small components of density (1 - rho). [1e-5]
%
%   .min                Minimum bound for the L1 reconstruction. [-3]
%   .max                Max bound for the L1 reconstruction. [3]
%
%   .beta               Inverse temperature for the prior 'Laplace'. [1]
%
%   .m_gaussComp        Mean of the gaussian for the prior 'Complex'. [0]
%
%   .var_gaussComp      Variance .. [1]

function opt = CSBP_Solver_Opt()
opt.nb_iter = 1000;
opt.print = 10;
opt.conv = 1e-8;
opt.learn = 0;
opt.signal_rho = -1;
opt.var_noise = 1e-10;
opt.signal = [];
opt.dump_learn = 0.;
opt.dump_mes = 0.5;
opt.prior = 'SparseGauss';
opt.save_memory = 1;
opt.save_speed = 0;
opt.option_noise = 0;
opt.remove_mean = 0;
opt.Gmean = [];
opt.Ymean = [];
opt.Nvec = [1];
opt.Mvec = [1];
opt.varG = [];
opt.method = 'AMP';
opt.MSEbyBlock = 0;
opt.alphaBig = 0;
opt.adaptDump = 0;

% SparseGauss
opt.m_gauss = 0;
opt.var_gauss = 1;

% SparseGaussCut
opt.m_gaussC = 0;
opt.var_gaussC = 1;
opt.Cut = 3;

% SparseGaussPositive
opt.m_gaussP = 0;
opt.var_gaussP = 1;

% SparseExponential
opt.expo = 1;

% SparseConstant
opt.c_down = 0;
opt.c_up = 1;

% 2Gauss
opt.m_1_gauss = 0;
opt.m_2_gauss = 0;
opt.var_1_gauss = 1e-5;
opt.var_2_gauss = 1;

% L1
opt.min = -3;
opt.max = 3;

% Laplace
opt.beta = 1;

% Complex
opt.m_gaussComp = 0;
opt.var_gaussComp = 1;

end