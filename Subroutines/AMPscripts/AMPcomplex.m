% AMP for complex number. To do : fast version and remove mean

V_new = (abs(G).^2 * prior.var_mess.').';
W_new = dumping(W,(G * prior.av_mess.').' - ((Y - W) ./ (n_and_e.var_noise + V) ) .* V_new,opt.dump_mes);
V_new = dumping(V,V_new,opt.dump_mes);
var_1 = (1 ./ (n_and_e.var_noise + V_new) ) * abs(G).^2;
var_2 = ((Y - W_new) ./ (n_and_e.var_noise + V_new) ) * conj(G);

prior.R = var_2 ./ var_1 + prior.av_mess; prior.S2 = 1 ./ var_1; prior = F(prior); V = V_new; W = W_new;
% Iteration done ----