% test for ode45

function dy = Quantum_osc(t, y, ut, u, a,param )
rho = reshape(y, size(a));
u = interp1(ut, u, t);

Num = a' * a;
%Oscillator parameters
K = param(1); kappa = param(2); alpha = param(3);

%evolution
H = K * Num^2 + alpha * u * (a + a'); %Hamiltonian
D_rho = a * rho * a' - 1/2 * ( Num*rho + rho * Num); %Noise term
drho = -1j * (H*rho - rho * H) + kappa * D_rho;

dy = drho(:);

end