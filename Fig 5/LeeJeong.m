function I = LeeJeong(rho)

s = size(rho);
d = sqrt(s(2)); %dimensionality of rho
nt = s(1); %number of time steps
I = zeros(1,nt);
a = diag(sqrt(1:d-1),1);

for i = 1:nt
    D = reshape(rho(i,:), [d,d]);
    L = a * D * a' - 1/2 * D * (a') * a - 1/2 * (a') * a * D;
    I(i) = real(-trace(D*L));
end

end