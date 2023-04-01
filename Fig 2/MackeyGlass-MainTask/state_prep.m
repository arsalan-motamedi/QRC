function rho = state_prep(d, alpha, M, c)
%This function prepares the desired initial state rho, to be either the
%coherent(c=1) or incohernt (c=0) superposition of two ket{alpha}'s.

V = zeros(d,M);
omega = exp(1j*2*pi/M);

for i = 1:d
    temp = alpha^(i-1)/sqrt(factorial(i-1));
    for j = 1:M
        V(i, j) = temp * (omega^((i-1)*(j-1)));
    end
end

if( c == 1 )
    temp = zeros(d,1);
    for i = 1:M
        temp = temp + V(:,i);
    end
    rho = temp*temp';
    rho = rho / trace(rho);
else
    rho = zeros(d,d);
    for i = 1:M
        rho = rho + V(:,i) * V(:,i)';
    end
    rho = rho/trace(rho);
end

end