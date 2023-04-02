function [rho] = coherent_state(alpha,d)

psi = zeros(d,1);

for i = 1:d
    psi(i) = alpha^(i-1)/sqrt(factorial(i-1));
end

rho = psi * psi';
rho = rho/real(trace(rho));
end