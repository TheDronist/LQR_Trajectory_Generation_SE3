function X_sol = dynamics(t,x,A,B,R,P,N)
[n m]= size(A);
P_ = reshape(P,size(A));
N_ = reshape(N,n,1);
u = -inv(R)*B'*(P_*x+N_) ;
dxdt = A*x + B*u;
X_sol = dxdt(:);
