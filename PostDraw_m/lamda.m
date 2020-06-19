function g = lamda(A)
R = qr(A);
W = diag(R);
w1 = min(abs(W));
w2 = max(abs(W));
B = 0.00025;
g = B*(w1^2/w2^2)*((w2^2+1)/2);
end