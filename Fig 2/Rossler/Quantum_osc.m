% test for ode45

function dy = Quantum_osc(t, y, ut, u, a,param)
rho = reshape(y, size(a));
u = interp1(ut, u, t);

Num = a' * a;
%Oscillator parameters
K = param(1); kappa = param(2); alpha = param(3);% lambda = param(4);

%evolution
H = K * Num^2 + alpha * u * (a + a'); %Hamiltonian
D_rho = a * rho * a' - 1/2 * ( Num*rho + rho * Num); %Noise term
%dephasing = -lambda * exp(-lambda * t) * rho + lambda * exp(-lambda * t) * diag(diag(rho));
%dephasing = zeros(size(Num));
%for i = 1:length(dephasing)
%    ket_n = zeros(length(dephasing),1);
%    ket_n(i) = 1;
%    dephasing = dephasing + rho(i,i) *(ket_n * ket_n') - ((ket_n * ket_n') * rho + rho * (ket_n * ket_n'))/2;
%end
%dephasing = dephasing * lambda;
drho = -1j * (H*rho - rho * H) + kappa * D_rho;% + dephasing;

dy = drho(:);

end