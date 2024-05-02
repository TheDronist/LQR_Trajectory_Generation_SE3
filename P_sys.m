function F = mRiccati(t, X, A, B, Q, R)
[n, m] = size(A);
Y = reshape(X , n, m+1);
P = Y(1:n, 1:m); %Convert from "n^2"-by-1 to "n"-by-"n"
N = Y(1:n, m+1); 
dPdt = -(A.'*P + P*A - P*B/R*B.'*P + Q);%Determine derivative
dNdt = (-A.' + P*B/R*B.')*N;
F = [dPdt(:); dNdt(:)];
