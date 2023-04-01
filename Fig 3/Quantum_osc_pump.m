
function dy = Quantum_osc_pump(t, y, ut, u, a,param )
rho = reshape(y, size(a));
u = interp1(ut, u, t);

Num = a' * a;
%Oscillator parameters
K = param(1); kappa = param(2);kappa2 = param(3) ; alpha = param(4);

%evolution
H = K * Num^2 + alpha * u * (a + a'); %Hamiltonian
D_rho = a * rho * a' - 1/2 * ( Num*rho + rho * Num) ;
D2 =  a' * rho * a  - 1/2 * ((a * a') * rho + rho * (a * a')); %Noise term
drho = -1j * (H*rho - rho * H) + kappa(1) * D_rho + kappa2 * D2;

dy = drho(:);

end