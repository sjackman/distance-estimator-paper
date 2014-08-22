%L = 6950;
%L = 6200;
L = 8000;
delta = 1;
k = 0;
gap = delta + 2*k;

% Population PMF
mu = 5000;
X = normpdf(1:L, mu, 0.1 * mu);
xstar = min(L, 2 * mu);

C = 0;
for i = gap : L
%	C += (i - gap) * X(i);
	C += (i - gap) / (L - i + 1) * X(i);
end
clear w
for i = gap : xstar
%	w(i) = (i - gap) / C;
	w(i) = (i - gap) / (L - i + 1) / C;
end
%plot(1:xstar, w);

S = X(1:xstar) .* w;
%plot(1:xstar, S);

plot(1:xstar, X(1:xstar), '-', 1:xstar, S, '-');
