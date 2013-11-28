% option field
My = CSBP_Solver_Opt();

% problem definition
size = 1000;
rho = 0.2;
mean = 0;
var = 1;
measure_rate = 0.4;

% signal generation
S = S_2Gauss(size, rho, 0, mean, 0, var);
My.signal = S;

% gaussian matrix generation
G = randn(ceil(measure_rate .* size), size);

% measure
Y = G * S';

% algorithm
X = CSBP_Solver(Y,G,My);