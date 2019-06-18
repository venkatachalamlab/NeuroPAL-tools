function y = lowpass(x, N)

x = [repmat(x(1), [1 N-1]) x];
y = conv(x, ones(1,N)/N, 'valid');

