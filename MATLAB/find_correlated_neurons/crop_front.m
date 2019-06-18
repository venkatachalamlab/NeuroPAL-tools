function y = crop_front(x, N)
x(1:N) = x(N+1);
y = x;