clc
close all
clear all 

 A = [zeros(12,3), eye(12,12);
     zeros(3,3), zeros(3,12)];
 B = [zeros(12,3); eye(3,3)];

 C = [eye(3,3), zeros(3,3), zeros(3,3), zeros(3,6);
      zeros(3,3), zeros(3,3), eye(3,3), zeros(3,6)];

g = 9.81;    %gravity
m = 4.34;    %mass of the drone
Sp = [-5; 0; 4];
e3 = [0; 0; 1];
f_x = m*g;
n = 4;
yi{1} = [0; 0; 0; 0; 0; 0];
for i = 1:n
    b{i} = [4*((i+2)/2)*cos((pi*4*i)/20); 6*((i+2)/2)*sin((pi*4*i)/20); 0.6*4*i];
    s{i} = (Sp - b{i})/norm(Sp - b{i});
    q{i} = (e3 - (e3'*s{i})*s{i})/norm((e3 - (e3'*s{i})*s{i}));
    a{i} = (f_x/m)*q{i} - g*e3;
    yi{i+1} = [b{i}(1); b{i}(2); b{i}(3); a{i}(1); a{i}(2); a{i}(3)];
    p{i} = cross(s{i},q{i});
    Ri{i} = [p{i} s{i} q{i}];
    v{i} = [b{i}'; Sp'];
end

%% Generation of position trajectory of the quadrotor
S = 10000*eye(6);
R = 1*eye(3);
Q = 0.01*eye(15);

P_f = C'*S*C;
N_f = -C.'*S*yi{5};
P0 = P_f;
N0 = N_f;
t_f = 16;
dt = 0.01;

%Riccati equation solution from 16s to 12s
t_res_1 = 16:-dt:12;
[t_1, P_1] = ode45(@(t,X)P_sys(t,X,A,B,Q,R),t_res_1,[P0 N0]); %solve the first riccati equation 
P1 = P_1((4/dt),1:225); %solution of riccati equation at ti+ = 12s + dt.
N1 = P_1((4/dt),226:240);
P1 = reshape(P1, 15, 15); %make P1 a 15*15 matrix
P1 = P1 + C'*S*C; %updating solution (equation 23 in the paper).
N1 = N1' - C'*S*yi{4};

%Riccati equation solution from 12s to 8s
t_res_2 = 12-dt:-dt:8;
[t_2, P_2] = ode45(@(t,X)P_sys(t,X,A,B,Q,R),t_res_2,[P1 N1]); %solve the first riccati equation 
P2 = P_2((4/dt),1:225); %solution of riccati equation at ti+ = 12s + dt.
N2 = P_2((4/dt),226:240);
P2 = reshape(P2, 15, 15); %make P1 a 15*15 matrix
P2 = P2 + C'*S*C; %updating solution (equation 23 in the paper).
N2 = N2' - C'*S*yi{3};

%Riccati equation solution from 8s to 4s
t_res_3 = 8-dt:-dt:4;
[t_3, P_3] = ode45(@(t,X)P_sys(t,X,A,B,Q,R),t_res_3,[P2 N2]); %solve the first riccati equation 
P3 = P_3((4/dt),1:225); %solution of riccati equation at ti+ = 12s + dt.
N3 = P_3((4/dt),226:240);
P3 = reshape(P3, 15, 15); %make P1 a 15*15 matrix
P3 = P3 + C'*S*C; %updating solution (equation 23 in the paper).
N3 = N3' - C'*S*yi{2};

%Riccati equation solution from 4s to 0s
t_res_4 = 4-dt:-dt:0;
[t_4, P_4] = ode45(@(t,X)P_sys(t,X,A,B,Q,R),t_res_4,[P3 N3]); %solve the first riccati equation 
P4 = P_4((4/dt),1:225); %solution of riccati equation at ti+ = 12s + dt.
N4 = P_4((4/dt),226:240);
P4 = reshape(P4, 15, 15); %make P1 a 15*15 matrix
P4 = P4 + C'*S*C; %updating solution (equation 23 in the paper).
N4 = N4' - C'*S*yi{1};

%Ricatti equation solution from 16s to 0s
P = [P_1(:,1:225); P_2(:,1:225); P_3(:,1:225); P_4(:,1:225)];
N = [P_1(:,226:240); P_2(:,226:240); P_3(:,226:240); P_4(:,226:240)];

%Initial solution of the system dynamics
X0 = zeros(15,1);
X_eul(:,1) = X0;

t_dis = [t_1; t_2; t_3; t_4];
t_dis1 = 0:dt:t_f;

%Solution of the system dynamics
[t,X_eul]=ode45(@(t,x)dynamics(t,x,A,B,R...
    ,linearInter(t,flipud(t_dis),flipud(P)),linearInter(t,flipud(t_dis),flipud(N))),t_dis1, X0);
X_eul = X_eul';

%Plot Figure 2(a) of the paper
figure
plot3(X_eul(1,:), X_eul(2,:), X_eul(3,:));
xlabel('x (m)') 
ylabel('y (m)') 
zlabel('z (m)')
set(gca,'xLim',[-20 20],'yLim',[-5 15],'zLim',[0 15]);
view(62,-3)
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
