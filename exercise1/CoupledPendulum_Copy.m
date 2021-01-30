% Simulation of coupled pendulum by Lagrangian mechanics
close all; clear; clc
% generalized coordinates
syms t dum_
th1 = str2sym('th1(t)');
th2 = str2sym('th2(t)');
q = str2sym('q(t)');
ic = str2sym('ic(t)');
% constants, length, mass, g, geometry
L1 = 5;
L2 = 6;
M = 20;
m_1 = 2;
m_2 = 2;
m_3 = 4;
g = 9.81;
d_0 = 7.5; % rest length spring
Kr = 0.7; % Motor constant
r = 0.1; % Wheels radii
R = 0.5; %Resistance
Li = 2.7e-3; %Inductance
% positions and velocities as function of the generalized coordinates

% First mass
x1 = q + L1 * sin(th1);
y1 = -L1 * cos(th1);

% Second mass
x2 = q + L1 * sin(th2); 
y2 = -L1 * cos(th2); 

% Distance from Cart to third mass

d1 = L1*cos((th2-th1)/2);
d2 = sqrt(L2^2 - L1^2*(sin((th2-th1)/2)^2));
d=d1+d2;

% Third mass

x3 = d*sin((th2+th1)/2) + q;
y3 = -d*cos((th2+th1)/2);

q_dot = diff(q, t);
x1_dot = diff(x1, t);
x2_dot = diff(x2, t);
x3_dot = diff(x3, t);

y1_dot = diff(y1, t);
y2_dot = diff(y2, t);
y3_dot = diff(y3, t);

% Bars Inertias


% Angles

alph1 = (acos(d2/L2))+(th2 + th1)/2;
alph2 = - (acos(d2/L2))+(th2 + th1)/2;
alph1_dot = diff(alph1,t);
alph2_dot = diff(alph2,t);

% kinetic and potential energy
T = M/2 *(diff(q,t))^2 + m_1/2 * ((diff(x1,t))^2 + (diff(y1,t))^2) +...
    m_2/2 * ((diff(x2,t))^2 + (diff(y2,t))^2) + m_3/2 * ((diff(x3,t))^2 + (diff(y3,t))^2);
k = 30;
V = m_1 * g * y1 + m_2 * g * y2 +  m_3 * g * y3 +...
    1/2 * k * (sqrt((x2 - x1)^2 + (y2 - y1)^2) - d_0)^2;

% Lagrangian
L = T - V;

% dL/d(qdot)
% dL_dqdot = diff(L, diff(q,t));
dL_dqdot =subs(diff(subs(L, diff(q, t), dum_), dum_), dum_, diff(q, t));
% dL_dth1dot = diff(L, diff(th1,t));
dL_dth1dot =subs(diff(subs(L, diff(th1, t), dum_), dum_), dum_, diff(th1, t));
% dL_dth2dot = diff(L, diff(th2,t));
dL_dth2dot =subs(diff(subs(L, diff(th2, t), dum_), dum_), dum_, diff(th2, t));


% dL/dq
% dL_dq = diff(L, q);
dL_dq = subs(diff(subs(L, q, dum_), dum_), dum_, q);
% dL_dth1 = diff(L, th1);
dL_dth1 = subs(diff(subs(L, th1, dum_), dum_), dum_, th1);
% dL_dth2 = diff(L,th2);
dL_dth2 = subs(diff(subs(L, th2, dum_), dum_), dum_, th2);


% dFdq
b = 5; % Joints dissipation constant
k_M = 5; % Wheels dissipation constant

F = 1/2 * b * (diff(th1,t)^2+diff(th2,t)^2+(alph1_dot-diff(th1,t))^2+...
    (alph2_dot - diff(th2,t))^2);
F_M = 1/2 * k_M * (diff(q,t))^2;
% dF_dqdot = diff(F, diff(q,t));
dF_dqdot = subs(diff(subs(F_M, diff(q, t), dum_), dum_), dum_, diff(q, t));
% dF_dth1dot = diff(F, diff(th1,t));
dF_dth1dot = subs(diff(subs(F, diff(th1, t), dum_), dum_), dum_, diff(th1, t));
% dF_dth2dot = diff(F, diff(th2,t));
 dF_dth2dot = subs(diff(subs(F, diff(th2, t), dum_), dum_), dum_, diff(th2, t));
 
% F external
U = 0;
% U = 5*cos(t/2);
% U = 5*sign(2-t)+5;
% F_ext = 100*cos(t);
F_ext= Kr*ic/r;
% generalized equations of motion
deq_1 = diff(dL_dth1dot, t) - dL_dth1 + dF_dth1dot;
deq_2 = diff(dL_dth2dot, t) - dL_dth2 + dF_dth2dot;
% deq_3 = diff(dL_dqdot, t) - dL_dq + dF_dqdot + F_ext ;
deq_3 = diff(dL_dqdot, t) - dL_dq + dF_dqdot;
deq_4 = Li*diff(ic,t) + R*ic + (Kr*q_dot)/r;

deqs = struct('th1', deq_1, 'th2', deq_2, 'q', deq_3, 'ic', deq_4); 
var = [th1;th2;q;ic];

leqs = [deqs.th1 == 0; deqs.th2 == 0; deqs.q == F_ext; deqs.ic == U];

[eqs,vars,m] = reduceDifferentialOrder(leqs,var);

[MassMatrix,ForceMatrix] = massMatrixForm(eqs,vars);
MassMatrix = simplify(MassMatrix);
ForceMatrix = simplify(ForceMatrix);

MM = odeFunction(MassMatrix, vars);
FF = odeFunction(ForceMatrix, vars);
%% Solve non linear ode system with losses
time = linspace(0, 60, 4000);
% initial conditions [th1, th2, q, ic, th1dot, th2dot, qdot]
x_0 = [-pi/3 pi/3 0 0 0 0 0];
opts=odeset('Mass', MM, 'Stats','on');
[~, x] = ode45(FF, time, x_0, opts);
% Calculate positions as function of generalized coordinates

TH1 = x(:,1);
TH2 = x(:,2);
X1 = x(:, 3) + L1 * sin(x(:, 1));
Y1 = -L1 * cos(x(:, 1));
X2 = x(:, 3) + L1 * sin(x(:, 2));
Y2 = -L1 * cos(x(:, 2));
Q = x(:, 3);
Q_dot = x(:,7);
Ic = x(:,4);
TH1_DOT = x(:,5);
TH2_DOT = x(:,6);

TH3= (x(:,2)-x(:,1))/2;
TH4= (x(:,1)+x(:,2))/2;

D1 = L1*cos(TH3);
D2=sqrt(L2^2-L1^2*sin(TH3).^2);
D= D1+ D2;

X3 = Q + D.*sin(TH4);
Y3 = -D.*cos(TH4);


%% Set Equations of motion without losses

% External Inputs

U = 0;
% U = 5*cos(t/2);
% U = 5*sign(2-t)+5;

F_ext= Kr*ic/r;

% generalized equations of motion

deq_1 = diff(dL_dth1dot, t) - dL_dth1;
deq_2 = diff(dL_dth2dot, t) - dL_dth2;
% deq_3 = diff(dL_dqdot, t) - dL_dq + dF_dqdot + F_ext ;
deq_3 = diff(dL_dqdot, t) - dL_dq;
deq_4 = L*diff(ic,t) + (Kr*diff(q,t))/r;

deqs = struct('th1', deq_1, 'th2', deq_2, 'q', deq_3, 'ic', deq_4); 
var = [th1;th2;q;ic];

leqs = [deqs.th1 == 0; deqs.th2 == 0; deqs.q == F_ext; deqs.ic == U];

[eqs,vars] = reduceDifferentialOrder(leqs,var);

[MassMatrix,ForceMatrix] = massMatrixForm(eqs,vars);
MassMatrix = simplify(MassMatrix);
ForceMatrix = simplify(ForceMatrix);

MM = odeFunction(MassMatrix, vars);
FF = odeFunction(ForceMatrix, vars);

%% Solve non linear ode system without losses
time = linspace(0, 60, 4000);
% initial conditions [th1, th2, q, ic, th1dot, th2dot, qdot]
x_0 = [-pi/3 pi/3 0 0 0 0 0]; 
opts=odeset('Mass', MM, 'Stats','on');
[~, x_nl] = ode45(FF, time, x_0, opts);
% Calculate positions as function of generalized coordinates
X1_nl = x_nl(:, 3) + L1 * sin(x_nl(:, 1));
Y1_nl = -L1 * cos(x_nl(:, 1));
X2_nl = x_nl(:, 3) + L1 * sin(x_nl(:, 2));
Y2_nl = -L1 * cos(x_nl(:, 2));
Q_nl = x_nl(:, 3);
Q_dot_nl = x_nl(:,7);
Ic_nl = x_nl(:,4);

TH3_nl= (x_nl(:,2)-x_nl(:,1))/2;
TH4_nl= (x_nl(:,1)+x_nl(:,2))/2;

D1_nl = L1*cos(TH3_nl);
D2_nl=sqrt(L2^2-L1^2*sin(TH3_nl).^2);
D_nl= D1_nl+ D2_nl;

X3_nl = Q_nl + D_nl.*sin(TH4_nl);
Y3_nl = -D_nl.*cos(TH4_nl);


%% Calculate Electrical power loss

%Input voltage in absolute value

% Uin=5*cos(time/2);

 Uin = 5*sign(2-time)+5;

% Calculate Ic_dot

pp = bsxfun(@plus,R*Ic,Kr*Q_dot/r);

Ic_dot = bsxfun(@minus,Uin',pp)/Li;

% Calculate power loss in the inductor

P_i = bsxfun(@times,Ic,Ic_dot)*Li;

% Input electrical power with losses

P_e = bsxfun(@times,Uin',Ic);

% Input electrical power without losses

P_e_nl = bsxfun(@times,Uin',Ic_nl);

%Output mechanical power with losses
 
P_m = bsxfun(@times,Q_dot,Ic)*Kr/r;

%Output mechanical power without losses
 
P_m_nl = bsxfun(@times,Q_dot_nl/r,Ic_nl*Kr);

% Electrical losses with losses

L_elec = bsxfun(@times,Ic,Ic)*R;

L_elec_th = P_e - P_m;

Dif = bsxfun(@minus,L_elec_th,L_elec);

Compare = [Dif P_i];

% Electrical losses without friction losses

L_elec_nl = bsxfun(@times,Ic_nl,Ic_nl)*R;

L_elec_th_nl = P_e_nl - P_m_nl;

%% Plot graphs
 figure(1)
 plot(X1)
 hold on
 plot(X1_nl)
 hold off
 
 figure(2)
 plot(500*abs(cos(time/2))*x(:,7))

%% Calculate Mechanical Power Loss

X1_dot_=L1.*cos(x(:,1)).*x(:,5) + x(:,7);
Y1_dot_=L1.*sin(x(:,1)).*x(:,5);

X2_dot_=L1.*cos(x(:,2)).*x(:,6) + x(:,7);
Y2_dot_=L1.*sin(x(:,2)).*x(:,6);


X3_dot_= Q_dot(:,1) - sin(x(:,1)/2 + x(:,2)/2).*(5*sin(x(:,1)/2 - x(:,2)/2).*...
    (x(:,5)/2 - x(:,6)/2) + (25*cos(x(:,1)/2 - x(:,2)/2).*...
    sin(x(:,1)/2 - x(:,2)/2).*(x(:,5)/2 - x(:,6)/2))./...
    (36 - 25*sin(x(:,1)/2 - x(:,2)/2).^2).^(1/2)) + cos(x(:,1)/2 + x(:,2)/2).*...
    (x(:,5)/2 + x(:,6)/2).*...
    (5*cos(x(:,1)/2 - x(:,2)/2) + (36 - 25*sin(x(:,1)/2 - x(:,2)/2).^2).^(1/2));

Y3_dot_ = cos(x(:,1)/2 + x(:,2)/2).*(5*sin(x(:,1)/2 - x(:,2)/2).*...
    (x(:,5)/2 - x(:,6)/2) + (25*cos(x(:,1)/2 - x(:,2)/2).*...
    sin(x(:,1)/2 - x(:,2)/2).*(x(:,5)/2 - x(:,6)/2))./...
    (36 - 25*sin(x(:,1)/2 - x(:,2)/2).^2).^(1/2)) + sin(x(:,1)/2 + x(:,2)/2).*...
    (x(:,5)/2 + x(:,6)/2).*...
    (5*cos(x(:,1)/2 - x(:,2)/2) + (36 - 25*sin(x(:,1)/2 - x(:,2)/2).^2).^(1/2));

% kinetic and potential energy

T_ = M/2 * Q_dot.^2 + m_1/2 * (X1_dot_.^2 + Y1_dot_.^2) +...
    m_2/2 * (X2_dot_.^2 + Y2_dot_.^2) + m_3/2 * (X3_dot_.^2 + Y3_dot_.^2);

V_ = m_1 * g * Y1 + m_2 * g * Y2 +  m_3 * g * Y3 +...
    1/2 * k * (sqrt((X2 - X1).^2 + (Y2 - Y1).^2) - d_0).^2;

% Lagrangian
L_ = T_ + V_;

% Mechanical Losses

% ALPH1_DOT_= TH1_DOT/2 + TH2_DOT/2 +...
%     (25*cos(TH1/2 - th2(t)/2)*sin(TH1/2 - th2(t)/2)*(TH1_DOT/...
%     2 - TH2_DOT/2))/(6*((25*sin(TH1/2 - th2(t)/2)^2)/36)^(1/2)*...
%     (36 - 25*sin(TH1/2 - th2(t)/2)^2)^(1/2));
% 
% ALPH2_DOT_= TH1_DOT/2 + TH2_DOT/2 -...
%     (25*cos(TH1/2 - th2(t)/2)*sin(TH1/2 - th2(t)/2)*(TH1_DOT/...
%     2 - TH2_DOT/2))/(6*((25*sin(TH1/2 - th2(t)/2)^2)/36)^(1/2)*...
%     (36 - 25*sin(TH1/2 - th2(t)/2)^2)^(1/2));
% 
% F = 1/2 * b * (diff(th1,t))^2+diff(th2,t)^2+(alph1_dot-diff(th1,t))^2+...
%     (alph2_dot - diff(th2,t))^2);
% F_M = 1/2 * k_M * (diff(q,t))^2;


%% Check Mechanical Power Loss=0 in no losses case 


X1_dot_nl=L1.*cos(x_nl(:,1)).*x_nl(:,5) + x_nl(:,7);
Y1_dot_nl=L1.*sin(x_nl(:,1)).*x_nl(:,5);

X2_dot_nl=L1.*cos(x_nl(:,2)).*x_nl(:,6) + x_nl(:,7);
Y2_dot_nl=L1.*sin(x_nl(:,2)).*x_nl(:,6);


X3_dot_nl= Q_dot_nl(:,1) - sin(x_nl(:,1)/2 + x_nl(:,2)/2).*(5*sin(x_nl(:,1)/2 - x_nl(:,2)/2).*...
    (x_nl(:,5)/2 - x_nl(:,6)/2) + (25*cos(x_nl(:,1)/2 - x_nl(:,2)/2).*...
    sin(x_nl(:,1)/2 - x_nl(:,2)/2).*(x_nl(:,5)/2 - x_nl(:,6)/2))./...
    (36 - 25*sin(x_nl(:,1)/2 - x_nl(:,2)/2).^2).^(1/2)) + cos(x_nl(:,1)/2 + x_nl(:,2)/2).*...
    (x_nl(:,5)/2 + x_nl(:,6)/2).*...
    (5*cos(x_nl(:,1)/2 - x_nl(:,2)/2) + (36 - 25*sin(x_nl(:,1)/2 - x_nl(:,2)/2).^2).^(1/2));

Y3_dot_nl= cos(x_nl(:,1)/2 + x_nl(:,2)/2).*(5*sin(x_nl(:,1)/2 - x_nl(:,2)/2).*...
    (x_nl(:,5)/2 - x_nl(:,6)/2) + (25*cos(x_nl(:,1)/2 - x_nl(:,2)/2).*...
    sin(x_nl(:,1)/2 - x_nl(:,2)/2).*(x_nl(:,5)/2 - x_nl(:,6)/2))./...
    (36 - 25*sin(x_nl(:,1)/2 - x_nl(:,2)/2).^2).^(1/2)) + sin(x_nl(:,1)/2 + x_nl(:,2)/2).*...
    (x_nl(:,5)/2 + x_nl(:,6)/2).*...
    (5*cos(x_nl(:,1)/2 - x_nl(:,2)/2) + (36 - 25*sin(x_nl(:,1)/2 - x_nl(:,2)/2).^2).^(1/2));

% kinetic and potential energy

T_nl = M/2 * Q_dot_nl(:,1).^2 + m_1/2 * (X1_dot_nl(:,1).^2 + Y1_dot_nl(:,1).^2) +...
    m_2/2 * (X2_dot_nl(:,1).^2 + Y2_dot_nl(:,1).^2) + m_3/2 * (X3_dot_nl.^2 + Y3_dot_nl.^2);

V_nl = m_1 * g * Y1_nl(:,1) + m_2 * g * Y2_nl(:,1) +  m_3 * g * Y3_nl(:,1) +...
    1/2 * k * (sqrt((X2_nl(:,1) - X1_nl(:,1)).^2 + (Y2_nl(:,1) - Y1_nl(:,1)).^2) - d_0).^2;

% Lagrangian
L_nl = T_nl - V_nl;




%% plot simulation
set(gcf, 'color', 'w')
set(gcf, 'position', [10, 100, 750, 750])
figure(3)
h = plot([]);
hold on
box on
axis equal
for i = 1 : numel(time)
    if ~ishghandle(h)
        break
    end
    cla
    plot([Q(i), X1(i)], [0, Y1(i)], 'b', 'Linewidth', 2);
    plot(X1(i), Y1(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4 * m_1);
    plot([Q(i), X2(i)], [0, Y2(i)], 'k', 'Linewidth', 2);
    plot(X2(i), Y2(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4 * m_2);
    plot([X2(i), X3(i)], [Y2(i),Y3(i)], 'r', 'Linewidth', 2);
    plot(X3(i), Y3(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4 * m_3);
    plot([X1(i), X3(i)], [Y1(i),Y3(i)], 'g', 'Linewidth', 2);
    axis([-20, 20, -15, 15]);
    h = draw_spring_2D([X1(i); Y1(i)], [X2(i); Y2(i)], 10, 0.5);
    drawnow
end
function h = draw_spring_2D(A, B, number_of_coils, y_amplitude)    
    persistent t
    
    normalvector_AB = (B - A) / norm(B - A);
    offset_A = A + 1.25 * normalvector_AB;
    offset_B = B - 1.25 * normalvector_AB;
    distance_between_offsets = norm(offset_B - offset_A);
    
    t = linspace(-pi, number_of_coils * 2 * pi, 500);
    x_coordinate_between_offsets = distance_between_offsets * linspace(0, 1, numel(t));
    
    % ratio between x amplitude and y
    ratio_X_div_Y = 0.5;
    
    x = x_coordinate_between_offsets + ratio_X_div_Y * y_amplitude * cos(t);
    y = y_amplitude * sin(t);
    
    coil_positions = [x; y];
    rotation_matrix = [normalvector_AB, null(normalvector_AB')];
    rotated_coil_positions = rotation_matrix * coil_positions;
    h = plot([A(1), offset_A(1) + rotated_coil_positions(1,:), B(1)], ...
         [A(2), offset_A(2) + rotated_coil_positions(2,:), B(2)], 'k');
end

