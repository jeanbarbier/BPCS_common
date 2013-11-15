function [E_new_cav] = DensEvoSparseGaussComplex(type, rho, alpha, Delta, E, varargin)
% Density evolution iteration for a complex gaussian signal with zero mean and unit
% variance, in the full or seeded case (type) under the NISHIMORI condition

minn = -20; maxx = -minn;
E_new_cav = 0;

alpha_ = @(S2, R, rho) abs(R).^2 ./ S2 - abs(R ./ (S2 + 1) ).^2 ./ (S2./ (1 + S2) );
Z = @(S2, R, rho) (1 - rho) .* exp(-0.5 .* abs(R).^2 ./ S2) + rho .* S2 ./ (S2 + 1) .* exp(-alpha_(S2, R, rho) ./ 2);
fa = @(S2, R, rho) (rho .* S2 ./ (S2 + 1) .* R ./ (S2 + 1) .* exp(-alpha_(S2, R, rho) ./ 2) ) ./ Z(S2, R, rho);
fc = @(S2, R, rho) (rho ./ Z(S2, R, rho) .* S2 ./ (S2 + 1) .* (2 .* S2 ./ (1 + S2) + abs(R ./ (S2 + 1) ).^2) .* exp(-alpha_(S2, R, rho) ./ 2) - abs(fa(S2, R, rho) ).^2) ./ 2;
G = @(z1, z2) exp(-0.5 .* (z1.^2 + z2.^2) ) ./ (2 .* pi);

if (strcmp(type, 'full') )
    
    int1 = integral2(@(z1, z2) G(z1, z2) .* fc((Delta + E) ./ alpha, (z1 + 1i .* z2) .* sqrt((E + Delta) ./ alpha), rho), minn, maxx, minn, maxx, 'AbsTol', 1e-10);
    int2 = integral2(@(z1, z2) G(z1, z2) .* fc((Delta + E) ./ alpha, (z1 + 1i .* z2) .* sqrt((E + Delta) ./ alpha + 1), rho), minn, maxx, minn, maxx, 'AbsTol', 1e-10);
    E_new_cav = (1 - rho) .* int1 + rho .* int2;
    
elseif (strcmp(type, 'seeded') )
    
    numBlockC = varargin{1}; numBlockL = varargin{2}; Nblock = varargin{3}; Mblock = varargin{4}; N = varargin{5}; J = varargin{6} .* N;
    
    for c = 1 : numBlockC;
        m(c) = 0;
        
        for l = 1 : numBlockL; m(c) = m(c) + Mblock(l) ./ Nblock .* J(l, c) ./ numBlockC ./ (Delta + sum(J(l, :) .* E) ./ numBlockC); end
                     
        int_ = integral(@(z) exp(-0.5 .* z.^2) ./ sqrt(2 .* pi) .* z.^2 ./ (rho + (1 - rho) .* exp(-0.5 .* z.^2 .* m(c) ) .* sqrt(m(c) + 1) ), minn, maxx, 'AbsTol', 1e-10);
        int1 = integral2(@(z1, z2) G(z1, z2) .* fc(1 ./ m(c), (z1 + 1i .* z2) .* sqrt(1 ./ m(c) ), rho), minn, maxx, minn, maxx, 'AbsTol', 1e-15);
        int2 = integral2(@(z1, z2) G(z1, z2) .* fc(1 ./ m(c), (z1 + 1i .* z2) .* sqrt(1 ./ m(c) + 1), rho), minn, maxx, minn, maxx, 'AbsTol', 1e-15);
        E_new_cav(c) = (1 - rho) .* int1 + rho .* int2;
    end
    
end

end