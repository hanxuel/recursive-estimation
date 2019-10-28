function [posEst,linVelEst,oriEst,driftEst,...
          posVar,linVelVar,oriVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,driftEst,...
%    posVar,linVelVar,oriVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
%

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % generate random pos
    a = 2*pi*rand();
    r = estConst.StartRadiusBound * sqrt(rand());
    posEst = [r*cos(a),r*sin(a)]; % 1x2 matrix
    linVelEst = [0.0,0.0]; % 1x2 matrix
    oriEst = 2*rand()*estConst.RotationStartBound - estConst.RotationStartBound; % 1x1 matrix
    driftEst = 0.0; % 1x1 matrix
    
    % initial state variance
    %posVar = [0.0,0.0]; % 1x2 matrix
    %linVelVar = [0.0,0.0]; % 1x2 matrix
    %oriVar = 0.0; % 1x1 matrix
    driftVar = 0.0; % 1x1 matrix
    posVar = [estConst.StartRadiusBound^2/4 estConst.StartRadiusBound^2/4]; % 1x2 matrix
    linVelVar = [0 0]; % 1x2 matrix
    oriVar = estConst.RotationStartBound^2/3; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,oriVar,linVelVar,driftVar]);
    % estimator state
    estState.xm = transpose([posEst,oriEst,linVelEst,driftEst]);
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = estState.tm;
estState.tm = tm; % update measurement update time


% get Const from estimator const
Q_d = estConst.DragNoise;
Q_r = estConst.RudderNoise;
Q_b = estConst.GyroDriftNoise;
C_d = estConst.dragCoefficient;
C_r = estConst.rudderCoefficient;
x_a = estConst.pos_radioA(1); y_a = estConst.pos_radioA(2);
x_b = estConst.pos_radioB(1); y_b = estConst.pos_radioB(2);
x_c = estConst.pos_radioC(1); y_c = estConst.pos_radioC(2);
sigma_a2 = estConst.DistNoiseA;
sigma_b2 = estConst.DistNoiseB;
sigma_c2 = estConst.DistNoiseC;
sigma_g2 = estConst.GyroNoise;
sigma_n2 = estConst.CompassNoise;
u_t = actuate(1); 
u_r = actuate(2); 
Q_c = diag([Q_r,Q_d,Q_d,Q_b]);

%combine state and var to get initial variable vector
x = [estState.xm(:); estState.Pm(:)];
% prior update

%solve ode equation for prior update
tspan = [dt tm];
%options = odeset('RelTol',1e-15,'AbsTol',1e-15);
[t,y] = ode45(@(t,y) odefcn(t,y,C_r,C_d,u_r,u_t,Q_c),tspan,x);

%get state and variance at the last step
xp = transpose(y(end,1:6));
pp = reshape(y(end,7:end),[6,6]);

% measurement update
% compute H and M for measurement update
% based on the value of z to determine the shape of H and M
if(isinf(sense(3)))
    H = zeros(4,6);
    H(1,:) = [(xp(1)-x_a)/sqrt((xp(1)-x_a)*(xp(1)-x_a)+(xp(2)-y_a)*(xp(2)-y_a)), (xp(2)-y_a)/sqrt((xp(1)-x_a)*(xp(1)-x_a)+(xp(2)-y_a)*(xp(2)-y_a)),0.0,0.0,0.0,0.0];
    H(2,:) = [(xp(1)-x_b)/sqrt((xp(1)-x_b)*(xp(1)-x_b)+(xp(2)-y_b)*(xp(2)-y_b)), (xp(2)-y_b)/sqrt((xp(1)-x_b)*(xp(1)-x_b)+(xp(2)-y_b)*(xp(2)-y_b)),0.0,0.0,0.0,0.0];
    H(3,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
    H(4,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    
    %compute h
    h = zeros(4,1);
    h(1) = sqrt((xp(1)-x_a)*(xp(1)-x_a)+(xp(2)-y_a)*(xp(2)-y_a));
    h(2) = sqrt((xp(1)-x_b)*(xp(1)-x_b)+(xp(2)-y_b)*(xp(2)-y_b));
    h(3) = xp(3)+xp(6);
    h(4) = xp(3);
    %compute M
    M = eye(4,4);
    R = diag([sigma_a2,sigma_b2,sigma_g2,sigma_n2]);
    z = [sense(1:2),sense(4:5)];
else
    H =zeros(5,6);
    H(1,:) = [(xp(1)-x_a)/sqrt((xp(1)-x_a)*(xp(1)-x_a)+(xp(2)-y_a)*(xp(2)-y_a)), (xp(2)-y_a)/sqrt((xp(1)-x_a)*(xp(1)-x_a)+(xp(2)-y_a)*(xp(2)-y_a)),0.0,0.0,0.0,0.0];
    H(2,:) = [(xp(1)-x_b)/sqrt((xp(1)-x_b)*(xp(1)-x_b)+(xp(2)-y_b)*(xp(2)-y_b)), (xp(2)-y_b)/sqrt((xp(1)-x_b)*(xp(1)-x_b)+(xp(2)-y_b)*(xp(2)-y_b)),0.0,0.0,0.0,0.0];
    H(3,:) = [(xp(1)-x_c)/sqrt((xp(1)-x_c)*(xp(1)-x_c)+(xp(2)-y_c)*(xp(2)-y_c)), (xp(2)-y_c)/sqrt((xp(1)-x_c)*(xp(1)-x_c)+(xp(2)-y_c)*(xp(2)-y_c)),0.0,0.0,0.0,0.0];
    H(4,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
    H(5,:) = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    %compute h
    h = zeros(5,1);
    h(1) = sqrt((xp(1)-x_a)*(xp(1)-x_a)+(xp(2)-y_a)*(xp(2)-y_a));
    h(2) = sqrt((xp(1)-x_b)*(xp(1)-x_b)+(xp(2)-y_b)*(xp(2)-y_b));
    h(3) = sqrt((xp(1)-x_c)*(xp(1)-x_c)+(xp(2)-y_c)*(xp(2)-y_c));
    h(4) = xp(3) + xp(6);
    h(5) = xp(3);
    %compute M
    M = eye(5,5);
    %compute R
    R = diag([sigma_a2,sigma_b2,sigma_c2,sigma_g2,sigma_n2]);
    z = sense;
end

%compute resulting measurement update
K = pp*transpose(H)/(H*pp*transpose(H)+M*R*transpose(M)+eye(size(R))*0.1*rand());
xm = xp + K*(transpose(z)-h);
pm = (eye(6,6)-K*H)*pp;

estState.xm = xm;
estState.Pm = pm;

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(4:5);
oriEst = estState.xm(3);
driftEst = estState.xm(6);

posVar = [estState.Pm(1,1),estState.Pm(2,2)];
linVelVar = [estState.Pm(4,4),estState.Pm(5,5)];
oriVar = estState.Pm(3,3);
driftVar = estState.Pm(6,6);

end

function dydt = odefcn(t,y,Cr,Cd,ur,ut,Qc)

%state vector ode equations
dydt = zeros(42,1);
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = Cr*ur;
dydt(4) = cos(y(3))*(tanh(ut)-Cd*(y(4)*y(4)+y(5)*y(5)));
dydt(5) = sin(y(3))*(tanh(ut)-Cd*(y(4)*y(4)+y(5)*y(5)));
dydt(6) = 0.0;

%compute A = dq/dx
A = zeros(6,6);
A(1,:) = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
A(2,:) = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
A(3,:) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
A(4,:) = [0.0, 0.0, -sin(y(3))*(tanh(ut)-Cd*(y(4)*y(4)+y(5)*y(5))),cos(y(3))*(-Cd*2*(y(4))),cos(y(3))*(-Cd*2*(y(5))), 0.0];
A(5,:) = [0.0, 0.0, cos(y(3))*(tanh(ut)-Cd*(y(4)*y(4)+y(5)*y(5))),sin(y(3))*(-Cd*2*(y(4))), sin(y(3))*(-Cd*2*(y(5))), 0.0];
A(6,:) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
%A = transpose(A);

%compute L = dq/dv
L = zeros(6,4);
L(1,:) = [0.0,0.0,0.0,0.0];
L(2,:) = [0.0,0.0,0.0,0.0];
L(3,:) = [Cr*ur,0.0,0.0,0.0];
L(4,:) = [0.0,-cos(y(3))*Cd*(y(4)*y(4)+y(5)*y(5)),0.0,0.0];%-cos(y(3))*Cd*(y(4)*y(4)+y(5)*y(5)),0.0];
L(5,:) = [0.0,-sin(y(3))*Cd*(y(4)*y(4)+y(5)*y(5)),0.0,0.0];%-sin(y(3))*Cd*(y(4)*y(4)+y(5)*y(5)),0.0];
L(6,:) = [0.0,0.0,0.0,1.0];

%reshape variance matrix
Pp = reshape(y(7:end),[6,6]);
%compute update Pp
dPpdt = A*Pp + Pp*transpose(A) + L*Qc*transpose(L);
dydt(7:end) = dPpdt(:);
end
