l1 = 100;
l2 = 200;
gap = 350;

% Population PMF
mu = 500;
X = normpdf(1 : 500 + l1 + l2, mu, 0.1 * mu);

%theta = 200;
figure
hold on
xstar = min(gap + l1 + l2, 2 * mu);
plot(1:xstar, X(1:xstar), 'g-');
for theta = 200 : 50 : 500
	C = 0;
	clear l p;
	for i = 1 : l1 + l2
%		p(i + theta) = X(i + theta) * l(i);
%		C += p(i + theta);
		l(i) = max(0, min([i, l1, l2, l1 + l2 - i]));
		p(i) = X(i + theta) * l(i);
%		p(i) = X(i + theta);
		C += p(i);
	end
	p = p / C;

	%xstar = min(theta + l1 + l2, 2 * mu);
	xstar = l1 + l2;
	plot(1:xstar, p(1:xstar));
	%plot(theta + (1:l1+l2), l(1:l1+l2) / max(l) * max(X));
end
hold off
%	plot(1:xstar, X(1:xstar), '-', 1:xstar, p(1:xstar), '-',
%		theta + (1:l1+l2), l(1:l1+l2) / max(l) * max(X));
