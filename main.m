clear
close all
clc

% A: The state matrix
A = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% B: The input matrix 
B = [0 0 0;
     0 0 0;
     0 0 0; 
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];

% C: The output matrix
C = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
 
g = 9.81;                     % gravity in kgm/s²
m = 4.34;                     % mass of the quadrotor in kg
Sp = [-5; 0; 4];              % center of the sphere
e3 = [0; 0; 1];               % the direction in z-axis
f = m*g;                      % the nominal force of the quad
n = 4;                        % the number of waypoints
y = cell(1,5);                % initializing the output vector containing the position and the acceleration of the quad
y{1} = [0; 0; 0; 0; 0; 0];

for l = 1:n
    %o = i-1;
    b{l} = [4*((l+2)/2)*cos((pi*4*l)/20); 6*((l+2)/2)*sin((pi*4*l)/20); 0.6*4*l];    % The positions of the waypoints in 3D (Eq. 34)
    s{l} = (Sp - b{l})/norm((Sp - b{l}));                                            % The pointing direction at the waypoints (Eq. 35)
    q{l} = (e3 - (e3'*s{l}*s{l}))/norm((e3 - (e3'*s{l}*s{l})));                      % The thrust direction at the waypoints (Eq. 6)
    a{l} = (f/m)*q{l} - g*e3;                                                        % The acceleration of the quad at the waypoints (Eq. 10)
    y{l+1} = [b{l}(1); b{l}(2); b{l}(3); a{l}(1); a{l}(2); a{l}(3)];                 % Update the output vector at the waypoints (Eq. 9(b))
    p{l} = cross(s{l},q{l});
    Ri{l} = [p{l} s{l} q{l}];
    v{l} = [b{l}'; Sp'];
end

%% Generation of the position trajectory of the quadrotor
% The weighting matrices of the LQR generator
S = 200*eye(6);
R = 20*eye(3);
Q = 1*eye(15);

% Simulation time variables
t_f = 16;
dt = 0.01;
t = 0:dt:t_f;

% The solution of the riccati equation at the last waypoint (4th waypoint)
% in (Eq. 23)
P_f1= C'*S*C;
N_f1 = -C.'*S*y{5};
P0 = P_f1;
N0 = N_f1;

% Reshape the solution of the riccati equation at the last waypoint
P_all_res1(:,1)= P_f1(:);
N_all_res1(:,1)= N_f1(:);

% The ricatti equation solution for the time period [16s, 12s] using Euler method

t_res1 = t_f:-dt:12;
for i = 2:1:length(t_res1)
    
    P1 =reshape(P_all_res1(:,i-1),size(A));
    N1 =reshape(N_all_res1(:,i-1),size(N_f1));
    
    dPdt1 = (A.'*P1 + P1*A - P1*B*inv(R)*B.'*P1 + Q);      %(Eq. 19)
    dNdt1 = -(-A' + P1*B*inv(R)*B.')*N1;                   %(Eq. 20)
    
    P1 = P1 + dt*(dPdt1);
    N1 = N1 + dt*(dNdt1);
    
    P_all_res1(:,i)= P1(:);
    N_all_res1(:,i)= N1(:);
end

P_all_res1 = P_all_res1';
N_all_res1 = N_all_res1';
P_f2 = P_all_res1(400,:);
N_f2 = N_all_res1(400,:);
P_f2 = reshape(P_f2, size(A));
N_f2 = N_f2';

% Update the solution of the riccati equation at the 3rd waypoint in (Eq. 23)
P_f2 = P_f2 + C'*S*C;
N_f2 = N_f2 - C'*S*y{4};

u{1} = y{3};
u{2} = y{2};
u{3} = y{1};

% The ricatti equation solution for the time period [12s, 8s], [8s, 4s], [4s, 0s]
for j = 1:1:n-1
    t_res = 4*j+dt:dt:4*(j+1);
    P_all_res{1,j}(:,1)= P_f2(:);
    N_all_res{1,j}(:,1)= N_f2(:);
    
    for k = 2:1:length(t_res)
        P =reshape(P_all_res{1,j}(:,k-1),size(A));
        N =reshape(N_all_res{1,j}(:,k-1),size(N_f1));    %(Eq. 19)
        dPdt = (A.'*P + P*A - P*B*inv(R)*B.'*P + Q);     %(Eq. 20)
        dNdt = -(-A' + P*B*inv(R)*B.')*N;
        P = P + dt*(dPdt);
        N = N + dt*(dNdt);
        P_all_res{1,j}(:,k)= P(:);
        N_all_res{1,j}(:,k)= N(:);
    end
    P_all_res{1,j} = P_all_res{1,j}';
    N_all_res{1,j} = N_all_res{1,j}';
    P_f = P_all_res{1,j}(399,:);
    N_f = N_all_res{1,j}(399,:);
    P_f = reshape(P_f, size(A));
    N_f = N_f';
    P_f2 = P_f + C'*S*C;
    N_f2 = N_f - C'*S*u{j};
end

%Concatenate the obtained solutions and flip them
P_sol = [ flipud(P_all_res{1,3}); flipud(P_all_res{1,2}); flipud(P_all_res{1,1}); flipud(P_all_res1)];
N_sol = [ flipud(N_all_res{1,3}); flipud(N_all_res{1,2}); flipud(N_all_res{1,1}); flipud(N_all_res1)];

%solve the system dynamics
X0 = zeros(15,1);
X_eul(:,1) = X0;

% add t_dis array 
t_dis = 0:dt:t_f;

for i = 2:length(t_dis)
    P_eul = reshape(P_sol(i-1,:),size(A));
    N_eul = reshape(N_sol(i-1,:),size(N_f1));
    U_eul(:,i-1) = -inv(R)*B.'*(P_eul*X_eul(:,i-1)+ N_eul);              %(Eq. 22)
    %X_eul(:,i) = X_eul(:,i-1) + dt*(((A-B*inv(R)*B'*P_eul)*X_eul(:,i-1))-B*inv(R)*B'*N_eul);
    X_eul(:,i) = X_eul(:,i-1) + dt* ((A*X_eul(:,i-1))+ B*U_eul(:,i-1) ); %(Eq. 21)
end

% Plot the Figure 2(b)
figure
plot3(X_eul(1,:), X_eul(2,:), X_eul(3,:), 'LineWidth', 2);
xlabel('x (m)') 
ylabel('y (m)') 
zlabel('z (m)')
hold on
b0 = [0; 0; 0];  %Starting waypoint
plot3(b0(1), b0(2), b0(3),'o'); 
hold on
%plot the sphere centered at Sp
[x,y,z] = sphere;  
x = x*1;
y = y*1;
z = z*1;
surf(x-5,y,z+4)
hold on
for i = 1:n
    plot3(b{i}(1),b{i}(2),b{i}(3), 'o')  % Plot the waypoints at bi's
    plot3(v{i}(:,1),v{i}(:,2),v{i}(:,3),'--');
    quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1]);         % plot the coordinate at the origine
    quiver3(b{i}(1),b{i}(2),b{i}(3),s{i}(1),s{i}(2),s{i}(3),1.5,'r','filled'); % plot the vectors si's
    quiver3(b{i}(1),b{i}(2),b{i}(3),q{i}(1),q{i}(2),q{i}(3),1.5,'k','filled'); % plot the vectors qi's
    quiver3(b{i}(1),b{i}(2),b{i}(3),p{i}(1),p{i}(2),p{i}(3),1.5,'g','filled'); % plot the vectors si's cross qi's
end

%% Generation of the attitude trajectory of the quadrotor

for o = 1:length(t_dis)
    a_t = [X_eul(7,:); X_eul(8,:); X_eul(9,:)];                                                % The acceleration vector obtained from the state
    a_t_dot = [X_eul(10,:); X_eul(11,:); X_eul(12,:)];                                         % The derivative of the acceleration
    q_t(:,o)= (g*e3 + a_t(:,o))/norm(g*e3 + a_t(:,o));                                         % The thrust vector in (Eq. 24)
    q_t_dot(:,o) = (a_t_dot(:,o) - (a_t_dot(:,o)'*q_t(:,o))*q_t(:,o))/norm(g*e3 + a_t(:,o));   % the derivative of the thrust vector in (Eq. 27(b))
    U_t(:,o) = cross(q_t(:,o),q_t_dot(:,o));                                                   % when alpha_t = 0
end

for o = 1:length(t_dis)
    U_t(:,o) = cross(q_t(:,o),q_t_dot(:,o));% + alpha_t(:,o)*q_t(:,o);
    S_U_t(o,:) = [0 U_t(3,o) -U_t(2,o) -U_t(3,o) 0 U_t(1,o) U_t(2,o) -U_t(1,o) 0];
end

%Solving the rotational dynamics R_dot = S(U)*R in (Eq. 28(a))
R0 = eye(3);    
R_sol(:,1) = R0(:);

for h = 2:1:length(t_dis)
    R_s = reshape(R_sol(:,h-1),size(R0));
    R_s = Rodrigues(dt*U_t(:,h))*R_s;
    R_sol(:,h)= R_s(:);
end

R_c_i{1} = eye(3);

for l = 1:n
    s_c_i{l} = (s{1,l}-(s{1,l}'*q_t(:,400*l+1))*q_t(:,400*l+1))/norm(s{1,l}-...
    (s{1,l}'*q_t(:,400*l+1))*q_t(:,400*l+1));                                     % The corrected pointing direction vectors at every waypoint eq. (26)
    R_c_i{l+1} = [cross(s_c_i{l},q_t(:,400*l+1)) s_c_i{l} q_t(:,400*l+1)];        % THe corrected waypoints rotations eq.(25)
    Rd{l+1} = R_c_i{l+1}*reshape(R_sol(:,400*l+1),[3,3])';                        % Obtain the difference between the propagated atitude and the desired corrected one eq.(28)
    Log_SO3_R{l+1} = logm(Rd{l+1});
    Theta{l+1} = [Log_SO3_R{l+1}(3,2);Log_SO3_R{l+1}(1,3);Log_SO3_R{l+1}(2,1)];   % The vee mapping eq.(29) 
    c_i{l+1} = ((30*Theta{l+1}(3,1))/(4^5));                                      % The coefficients c_i of alpha_t in eq.(33)
end

% Obtain the function alpha(t) in [0s, 4s]
temps{1} = 0:dt:4-dt;
alpha{1,1} = c_i{1,2}*((temps{1,1}(:,1:400)).^2).*((temps{1,1}(:,1:400)-4).^2);

% Obtain the function alpha(t) in [4s, 12s]
for w = 2:3
    temps{w} = 4*(w-1):dt:4*w-dt;
    alpha{1,w} = c_i{1,w+1}*((temps{1,w}(:,1:400)-4*(w-1)).^2).*((temps{1,w}(:,1:400)-4*w).^2);
end

% Obtain the function alpha(t) in [12s, 16s]
temps{4} = 4*3:dt:16;
alpha{1,4} = c_i{1,4+1}*((temps{1,4}(:,1:401)-4*(4-1)).^2).*((temps{1,4}(:,1:401)-4*4).^2);

alpha_t = [alpha{1,1} alpha{1,2} alpha{1,3} alpha{1,4}];         % Obtain the equation of alpha_t in eq.(32)

% Update all of the vector U = q*q_dot + alpha_t*q in (Eq. 27(b))
for o = 1:length(t_dis)
    U_t(:,o) = cross(q_t(:,o),q_t_dot(:,o)) + alpha_t(:,o)*q_t(:,o);
end

%Solve again the rotational dynamics but with new vector U
R0 = eye(3);
R_sol(:,1) = R0(:);

for h = 2:1:length(t_dis)
    R_s = reshape(R_sol(:,h-1),size(R0));
    R_s = Rodrigues(dt*U_t(:,h))*R_s;
    R_sol(:,h)= R_s(:);
end

% The rotational matrix of quad
J = [0.820 0 0; 0 0.0845 0; 0 0 0.1377];

% Ploting the Figure 2(b)
s_t(:,1) = [0;1;0];
for h = 2:1:length(t_dis)    
    R_ss = reshape(R_sol(:,h),size(R0));
    s_t(:,h)  =R_ss*[0;1;0];
    q_t2(:,h) = R_ss*[0;0;1];
end
figure
plot(t, s_t', 'LineWidth', 2)
xlabel('t(s)')
ylabel('s(t)')

% Ploting the Figure 2(c)
q_t2(:,1) = [0;0;1];
for i = 1:length(t_dis)
    R_d = reshape(R_sol(:,i), 3,3);
    f_d = (g*e3 + a_t)';
    R_dd = R_d*e3;
    f_dd(:,i) = m*(f_d(i,:)*R_dd);
end

figure
plot(t,q_t2', 'LineWidth', 2)
xlabel('t(s)')
ylabel('q(t)')

% Ploting the Figure 2(d)
figure
plot(t, f_dd, 'LineWidth', 2)
xlabel('t(s)')
ylabel('f(N)')