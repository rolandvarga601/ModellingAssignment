% Simulation of coupled pendulum by Lagrangian mechanics
close all; clear all; clc
% generalized coordinates
syms t dum_
th1 = str2sym('th1(t)');
th2 = str2sym('th2(t)');
q = str2sym('q(t)');

% constants, length, mass, g, geometry
L1 = 5;
L2 = 6;
M = 20;
m_1 = 2;
m_2 = 2;
m_3 = 4;
g = 9.81;
d_0 = 7.5; % rest length spring
% positions and velocities as function of the generalized coordinates

% First mass
x1 = q + L1 * sin(th1);
y1 = -L1 * cos(th1);

% Second mass
x2 = q + L1 * sin(th2); 
y2 = -L1 * cos(th2); 


% kinetic and potential energy
% T = M/2 * (diff(q,t))^2 + m_1/2 * ((diff(x1, t))^2 + (diff(y1, t))^2) +...
%     m_2/2 * ((diff(x2, t))^2 + (diff(y2, t))^2) + ...
%     m_3/2 * ((diff(x3, t))^2 + (diff(y3, t))^2);
T = M/2 * (diff(q,t))^2 + m_1/2 * ((diff(x1, t))^2 + (diff(y1, t))^2) +...
    m_2/2 * ((diff(x2, t))^2 + (diff(y2, t))^2);
k = 2;
% V = m_1 * g * y1 + m_2 * g * y2 +  m_3 * g * y3 +...
%     1/2 * k * (sqrt((x2 - x1)^2 + (y2 - y1)^2) - d_0)^2;
V = m_1 * g * y1 + m_2 * g * y2 +...
    1/2 * k * (sqrt(((x2 - x1)^2 + (y2 - y1)^2))-d_0)^2;

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


% generalized equations of motion
deq_1 = diff(dL_dth1dot, t) - dL_dth1 ;
deq_2 = diff(dL_dth2dot, t) - dL_dth2 ;
deq_3 = diff(dL_dqdot, t) - dL_dq ;


deqs = struct('th1', deq_1, 'th2', deq_2, 'q', deq_3); 
var = [th1;th2;q];

F_Ext=50*cos(2*t);
% F_Ext=0;
leqs = [deqs.th1 == 0; deqs.th2 == 0; deqs.q == F_Ext];

[eqs,vars,m] = reduceDifferentialOrder(leqs,var);

[MassMatrix,ForceMatrix] = massMatrixForm(eqs,vars);
MassMatrix = simplify(MassMatrix);
ForceMatrix = simplify(ForceMatrix);

MM = odeFunction(MassMatrix, vars);
FF = odeFunction(ForceMatrix, vars);
%% Solve non linear ode system without losses
time = linspace(0, 60, 10000);
% initial conditions [th1, th2, q, th1dot, th2dot, qdot]
x_0 = [-pi/3 pi/3 0 0 0 0]; 
res =0.9;
opts=odeset('Mass', MM, 'Stats','on','Events',@myEventsFcn);
x_nl = 0;
counter=0;
te=0;
% 
% [tt, x_nl1,te,ye,ie] = ode45(FF, time, x_0, opts);
% 
% switch ie(end) 
%     
%     case 1
%         ye(4)=-ye(4)*res;
%     case 2
%         ye(5)=-ye(5)*res;
%     case 3
%         ye(4)=-ye(4);
%         ye(5)=-ye(5);
% end
% 
% % Set new initial conditions
% x_0 = ye; 
% timeindex = find(t>te(end),1);
% 
% timecut = time(timeindex:end);
% % time = linspace(te, 60, 10000/(2*te));
% opts=odeset('Mass', MM, 'Stats','on');
% % [tt, x_nl2,te,ye,ie] = ode45(FF, timecut, x_0, opts);
% [tt, x_nl2] = ode45(FF, timecut, x_0, opts);
% x_nl=[x_nl1;x_nl2]; 


timecut=time;
while timecut(1)~=time(end)
    
[tt, x_nl1,te,ye,ie] = ode45(FF, timecut, x_0, opts);
ye = ye(end,:);
switch ie(end) 
    
    case 1
        ye(4)=-ye(4)*res;
    case 2
        ye(5)=-ye(5)*res;
    case 3
        ye(4)=-ye(4);
        ye(5)=-ye(5);
end
% Set new initial conditions
x_0 = ye; 
% timecut = linspace(te, 60, 10000/(2*te));
timeindex = find(time>te(end),1);

timecut = time(timeindex:end);
if counter==0
    x_nl=x_nl1;
    else
    x_nl=[x_nl;x_nl1];    
end
counter=counter+1;
end






% Calculate positions as function of generalized coordinates
X1_nl = x_nl(:, 3) + L1 * sin(x_nl(:, 1));
Y1_nl = -L1 * cos(x_nl(:, 1));
X2_nl = x_nl(:, 3) + L1 * sin(x_nl(:, 2));
Y2_nl = -L1 * cos(x_nl(:, 2));
Q_nl = x_nl(:, 3);
Q_dot_nl = x_nl(:,6);

%% Check Mechanical Power Loss=0 in no losses case 

for i=1:size(X2_nl,1)

X1_dot_nl(i,1)=L1*cos(x_nl(i,1))*x_nl(i,4) + x_nl(i,6);
Y1_dot_nl(i,1)=L1*sin(x_nl(i,1))*x_nl(i,4);

X2_dot_nl(i,1)=L1*cos(x_nl(i,2))*x_nl(i,5) + x_nl(i,6);
Y2_dot_nl(i,1)=L1*sin(x_nl(i,2))*x_nl(i,5);


% kinetic and potential energy

T_nl(i,1) = M/2 * Q_dot_nl(i,1)^2 + m_1/2 * (X1_dot_nl(i,1)^2 + Y1_dot_nl(i,1)^2) +...
    m_2/2 * (X2_dot_nl(i,1)^2 + Y2_dot_nl(i,1)^2) ;

V_nl(i,1) = m_1 * g * Y1_nl(i,1) + m_2 * g * Y2_nl(i,1) + ...
    1/2 * k *(sqrt(((X2_nl(i,1) - X1_nl(i,1))^2 + (Y2_nl(i,1) - Y1_nl(i,1))^2))-d_0)^2;

end

L_nl = T_nl + V_nl;

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
    plot([Q_nl(i), X1_nl(i)], [0, Y1_nl(i)], 'b', 'Linewidth', 2);
    plot(X1_nl(i), Y1_nl(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4 * m_1);
    plot([Q_nl(i), X2_nl(i)], [0, Y2_nl(i)], 'k', 'Linewidth', 2);
    plot(X2_nl(i), Y2_nl(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4 * m_2);
    axis([-20, 20, -15, 15]);
    h = draw_spring_2D([X1_nl(i); Y1_nl(i)], [X2_nl(i); Y2_nl(i)], 10, 0.5);
    drawnow
end
%%
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
function [value,isterminal,direction] = myEventsFcn(t,y)
value = [y(1)+pi/2;y(2)-pi/2;abs(y(1)-y(2))-2*asin(3/(2*5))];
isterminal = [1;1;1];
direction = [0;0;0];

% value = [y(1)+pi/2;y(2)-pi/2];
% isterminal = [1;1];
% direction = [0;0];
end
