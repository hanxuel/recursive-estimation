function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==1, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index, scalar
%                       corresponds to continous time t = k*Ts
%                       If tm==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
%

% Set number of particles:
N_particles = 15000; % obviously, you will need more particles than 10.
K = 0.0050;% roughening noise tuning parameter
d = 3; % dimension of state space

[maxvalue,index_max]=max(estConst.contour);
max_x=maxvalue(1);
max_y=maxvalue(2);
[minvalue,index_min]=min(estConst.contour);
min_x=minvalue(1);
min_y=minvalue(2);
%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    r = sqrt(rand(1,N_particles)).*estConst.d;
    theta = rand(1,N_particles).*2*pi-pi;
    index = rand(1,N_particles);
    postParticles.x_r = zeros(1,N_particles);
    postParticles.y_r = zeros(1,N_particles);
    postParticles.phi = rand(1,N_particles).*estConst.phi_0.*2-estConst.phi_0; % 1xN_particles matrix
    postParticles.x_r(index>0.5) = estConst.pA(1)+ r(index>0.5).*cos(theta(index>0.5));% [1xN_particles matrix, 1xN_particles matrix]
    postParticles.x_r(index<=0.5) = estConst.pB(1)+ r(index<=0.5).*cos(theta(index<=0.5));
    postParticles.y_r(index>0.5) = estConst.pA(2)+ r(index>0.5).*sin(theta(index>0.5));
    postParticles.y_r(index<=0.5) = estConst.pB(2)+ r(index<=0.5).*sin(theta(index<=0.5)); 
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
if (km>0)
% If km > 0, we perform a regular update of the estimator.
% Implement your estimator here!


    % Prior Update:
    v_f=rand(1,N_particles)*estConst.sigma_f-estConst.sigma_f/2;
    priorParticles.x_r=prevPostParticles.x_r+(act(1)+v_f).*cos(prevPostParticles.phi);
    priorParticles.y_r=prevPostParticles.y_r+(act(1)+v_f).*sin(prevPostParticles.phi);
    v_phi=rand(1,N_particles)*estConst.sigma_phi-estConst.sigma_phi/2;
    priorParticles.phi=prevPostParticles.phi+(v_phi+act(2));

    % Posterior Update:
    
    % Obtain each line segment
    [width,~] = size(estConst.contour);
    contour = zeros(size(estConst.contour));
    contour(1:width-1,:) = estConst.contour(2:width,:);
    contour(width,:) = estConst.contour(1,:);
    vertex = zeros(width,4);
    vertex(:,1:2) = estConst.contour;
    vertex(:,3:4) = contour;
    
    % the relationship between measurement and prior position
    lamda = zeros(width,N_particles);
    r_distance = zeros(width,N_particles);
    for i =1:width
        lamda(i,:) = (priorParticles.x_r.*sin(priorParticles.phi)-vertex(i,3).*sin(priorParticles.phi)-priorParticles.y_r.*cos(priorParticles.phi)...
        +vertex(i,4).*cos(priorParticles.phi))./(vertex(i,1).*sin(priorParticles.phi)-vertex(i,3).*sin(priorParticles.phi)-vertex(i,2).*cos(priorParticles.phi)...
        +vertex(i,4).*cos(priorParticles.phi));
        r_distance(i,:) = (lamda(i,:).*(vertex(i,1)-vertex(i,3))+vertex(i,3)-priorParticles.x_r)./cos(priorParticles.phi);
    end
    %distance = zeros(1,N_particles);
    lamda(lamda<=1 & lamda>=0)=1;
    lamda(lamda>1 | lamda<0)=inf;
    r_distance(r_distance<0)=inf;
    distance = min(lamda.*r_distance);
    %for i=1: N_particles
    %    if(~isempty(min(r_distance((lamda(:,i)>=0 & lamda(:,i)<=1)& r_distance(:,i)>=0,i))))
    %        distance(i) = min(r_distance((lamda(:,i)>=0 & lamda(:,i)<=1)& r_distance(:,i)>=0,i));
    %    end
    %end
    % posterior probability
    prob_sense = prob(sens-distance,estConst.epsilon);
    %size(prob_sense)
    if sum(prob_sense) ~=0
        prob_sense = prob_sense/sum(prob_sense);
        distri_prob=cumsum(prob_sense);
        distri_prob=[0,distri_prob];
    % resampling 
        seed = rand(1,N_particles);
        Y=discretize(seed,distri_prob);
        postParticles.x_r = priorParticles.x_r(Y);
        postParticles.y_r = priorParticles.y_r(Y);
        postParticles.phi = priorParticles.phi(Y);
            % roughening
        Ex = range(postParticles.x_r);
        Ey = range(postParticles.y_r);
        Ephi = range(postParticles.phi);
        noise_sigma_x = K*Ex*N_particles^(-1/d);
        noise_sigma_y = K*Ey*N_particles^(-1/d);
        noise_sigma_phi = K*Ephi*N_particles^(-1/d);
        noise_x = noise_sigma_x*randn(1,N_particles);
        noise_y = noise_sigma_y*randn(1,N_particles);
        noise_phi = noise_sigma_phi*randn(1,N_particles);
        postParticles.x_r = postParticles.x_r + 0.5*noise_x;
        postParticles.y_r = postParticles.y_r + 0.5*noise_y;
        postParticles.phi = postParticles.phi + 0.5*noise_phi;
        randomsample = [postParticles.x_r;postParticles.y_r;postParticles.phi];
        randomsample(:,(randomsample(2,:)<min_y|randomsample(2,:)>max_y|randomsample(1,:)<min_x|randomsample(1,:)>max_x))=[];
%         randomsample(:,(randomsample(1,:)<estConst.contour(10,1) & randomsample(2,:)<estConst.contour(10,2)))=[];
%         randomsample(:,randomsample(1,:)<estConst.contour(7,1) & randomsample(2,:)>estConst.contour(7,2))=[];
%         randomsample(:,randomsample(1,:)>estConst.contour(4,1) & randomsample(2,:)>(-randomsample(1,:)+2.5))=[];
        [~,n]=size(randomsample);
        
        doubleramdom=repmat(randomsample,[1,ceil(N_particles/n)]);
        index=randperm(size(doubleramdom,2),N_particles);
        postParticles.x_r = doubleramdom(1,index);
        postParticles.y_r = doubleramdom(2,index);
        postParticles.phi = doubleramdom(3,index);
    else
        disp('no weights');

        postParticles.x_r = priorParticles.x_r;
        postParticles.y_r = priorParticles.y_r;
        postParticles.phi = priorParticles.phi;
        Ex = range(postParticles.x_r);
        Ey = range(postParticles.y_r);
        Ephi = range(postParticles.phi);
        noise_sigma_x = K*Ex*N_particles^(-1/d);
        noise_sigma_y = K*Ey*N_particles^(-1/d);
        noise_sigma_phi = K*Ephi*N_particles^(-1/d);
        noise_x = noise_sigma_x*randn(1,N_particles);
        noise_y = noise_sigma_y*randn(1,N_particles);
        noise_phi = noise_sigma_phi*randn(1,N_particles);
        postParticles.x_r = postParticles.x_r + 1.8*noise_x;
        postParticles.y_r = postParticles.y_r + 1.8*noise_y;
        postParticles.phi = postParticles.phi + 0.8*noise_phi;
        randomsample = [postParticles.x_r;postParticles.y_r;postParticles.phi];
        randomsample(:,(randomsample(2,:)<0|randomsample(2,:)>2|randomsample(1,:)<0|randomsample(1,:)>2))=[];
        randomsample(:,randomsample(1,:)<estConst.contour(10,1) & randomsample(2,:)<estConst.contour(10,2))=[];
        randomsample(:,randomsample(1,:)<estConst.contour(7,1) & randomsample(2,:)>estConst.contour(7,2))=[];
        randomsample(:,randomsample(1,:)>estConst.contour(4,1) & randomsample(2,:)>(-randomsample(1,:)+2.5))=[];
        [~,n]=size(randomsample);
        
        doubleramdom=repmat(randomsample,[1,ceil(N_particles/n)]);
        index=randperm(size(doubleramdom,2),N_particles);
        postParticles.x_r = doubleramdom(1,index);
        postParticles.y_r = doubleramdom(2,index);
        postParticles.phi = doubleramdom(3,index);
    end
    

    
    
end % end estimator
end

function p = prob(w,epsilon)
x = -abs(w);
p = (x>=-3*epsilon & x<-2.5*epsilon).*(1/(2.5*epsilon^2).*x+3/(2.5*epsilon))+(x>=-2.5*epsilon & x<-2*epsilon).*(-1/(2.5*epsilon^2).*x-4/(5*epsilon))...
    +(x>=-2*epsilon & x<=0).*(1/(5*epsilon^2).*x+2/(5*epsilon));
p(isnan(p))=0.0;
end

