function y = highpass(x, N)

x = [repmat(x(1), [1 N-1]) x];
m = median(x, 'omitnan');

z = lowpass(x, N);
y = x - z + m;

y = y(N:end);

