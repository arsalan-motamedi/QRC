function out = phase_space(X, n)

N = length(X);
out = zeros(2,N-n);
for i = 1:N-n
    out(1,i) = X(i); out(2,i) = X(i+n);
end

end