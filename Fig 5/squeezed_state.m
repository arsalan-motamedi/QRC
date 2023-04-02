function [rho] = squeezed_state(alpha, zeta, d)
a = diag(sqrt(1:d-1),1);
S = expm( 1/2* (zeta' * a*a - zeta * a' * a'));
D = expm( alpha * a' - alpha' * a);

psi = zeros(d,1);
psi(1) = 1;

psi = D * S * psi;
rho = psi * psi';

rho = rho/trace(rho);

end