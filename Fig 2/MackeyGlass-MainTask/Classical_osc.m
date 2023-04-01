function da = Classical_osc(t, a, ut, u, param)

K = param(1); kappa = param(2); alpha = param(3);

u = interp1(ut, u, t);

da = -1j * K * (1+2*abs(a)^2)*a - kappa/2 * a - 1j* alpha * u;

end