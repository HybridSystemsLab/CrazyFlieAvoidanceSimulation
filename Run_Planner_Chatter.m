close all
clear
clear GLOBAL
clc

global HEQOptions;
HEQOptions = odeset('AbsTol', 1e-2,'RelTol',1e-2);

global gamma lambda sigma_max planTime planStep executeTime worldTimeStep obstacleSize FirstRunEject ConstThrust objColStep prevRef ObjPlot55
gamma = 9.8;
lambda = .6;
sigma_max = 0;
planTime = 0.2;
planStep = 0.1;
executeTime = 0.01;
worldTimeStep = 0.1;
obstacleSize = 0.5;
FirstRunEject = 0;
ConstThrust = 0;
objColStep = 1;
prevRef = [];
ObjPlot55 = 0;

global Target TargetRad
Target = [1.5,0,0];
TargetRad = 0.3;

global mobilityRadius usingPointMassVehicle maxAcc Mass InertialTensor MaxMotorThrust d_xy PlanningThrustMult chatterSim
usingPointMassVehicle = 0;
mobilityRadius = planTime;
maxAcc = 23;
Mass = 0.028; %kg
InertialTensor = [16.571710, 0.830806, 0.718277; 0.830806, 16.655602, 1.800197; 0.718277, 1.800197, 29.261652] * 10^-6; % 10^-6kg meters^2
MaxMotorThrust = 0.1597;%newtons % 2.130295*10^-11 *65535^2 + 1.032633*10^-6 *65535 + 5.484560*10^-4
d_xy = (( 2^0.5)/4)*0.092; %meters
PlanningThrustMult = 0.9;
chatterSim = 1;

%Hybrid controller constants
global Qdelta alpha K k_V_0 k_z k_p k_v k k_omega k_b1 k_b2 beta k_q
Qdelta = 0.5; % Qdelta \in (0,1)
alpha = 0.5; % \in (0,1)
K = 1;
k_V_0 = 0.01;
k_z = 0.3;
k_p = 3;
k_v = 6;
k = 3;
k_omega = 40;
k_b1 = 1;
k_b2 = 1;
beta = k_v/4; % beta \in (0,k_v)
k_q = 1;

global spherePts spx spy spz
spherePts = 6;
[spx, spy, spz] = sphere(spherePts);
spx = obstacleSize*reshape(spx(:,1:spherePts),[],1);
spy = obstacleSize*reshape(spy(:,1:spherePts),[],1);
spz = obstacleSize*reshape(spz(:,1:spherePts),[],1);

% p_x_a; p_y_a; p_z_a; v_x_a; v_y_a; v_z_a; R_11-R_31; R_12-R_32; R_13-R_33; omega; z; h; qhat
startq = [0.01 1 0 0]/norm([0.01 1 0 0]);
startR = doubleCover([0.01 1 0 0]/norm([0.01 1 0 0]));
x_a = [0, 0, 0, 0, 0, 0, reshape(startR,[1,9]), zeros(1,3), zeros(1,3), -1, startq];
% p_x_o; v_x_o; p_y_o; v_y_o; p_z_o; v_z_o; sigmaMin; sigmaMax
xi_o_0 = [0.4*obstacleSize/0.5 + 0.01, 0, -0.3*obstacleSize/0.5, 0, -0.0, 0, -sigma_max, sigma_max];
xi_o_1 = [0.4*obstacleSize/0.5 + 0.01, 0, 0.3*obstacleSize/0.5, 0, -0.0, 0, -sigma_max, sigma_max];
% x_o = mat2cell([0,0,xi_o_0], 1, 10);
x_o = [mat2cell([0,0,xi_o_0], 1, 10);mat2cell([0,0,xi_o_1], 1, 10)];
% t,
xi_p = 0;

t_a = 0;
j_a = 0;
r_log = zeros(1,31);

%Run planner and update world
global loopCount
for i=1:300
    %Run planner
    loopCount = i;
    if(i == 1)
        [r,rcost] = setBasedPlanner(x_a(end,:), x_o, @cost);
    else
        x_a(end,:)
        rAsStateTemp = interp1(prevRef(:,1),prevRef(:,2:end), executeTime)
        rAsState = [rAsStateTemp(1:6), rAsStateTemp(16:27)]
        [r,rcost] = setBasedPlanner(rAsState,x_o,@cost);
    end
    r_log = [r_log;(i*executeTime)+r(:,1), r(:,2:end)];
    
    %Update vehicle
    [t_a,j_a,x_a] = updateVehicle(r,t_a,j_a,x_a);
    
    %Update obstacles
    x_o = updateObjs(x_o,x_a);
    
    plotResults(t_a, x_a, x_o, r_log(2:end,:), rcost);
    dist = max(norm(x_a(end,3:5)-Target) - TargetRad, 0);
    if(dist<0)
        break
    end
end

plotResults(t_a, x_a, x_o, r_log(2:end,:), rcost);

function newU = updateObjs(U, x_a)
    global executeTime
    persistent chatter holdcount lastDir
    [nobj,~] = size(U);
    newU = cell(nobj,1);
    if(isempty(chatter) || isempty(holdcount) || isempty(lastDir))
        chatter = 1;
        holdcount = 1;
        lastDir = 0;
    end
    if(lastDir ~= sign(x_a(end,6)))
        chatter = sign(x_a(end,6));
    else
        chatter  = 0;
    end
    for i=1:nobj
        x_o = cell2mat(U(i));
        out = [x_o; x_o(end,1)+executeTime, x_o(end,2), x_o(end,3:6), sign(x_a(end,6))*0.05, x_o(end,8:end)];
        [r,c] = size(out);
        newU(i) = mat2cell(out,r,c);
    end
    lastDir = sign(x_a(end,6));
    
    if(holdcount == 1)
        chatter = 0;
        holdcount = 2;
    elseif(holdcount == 2)
        chatter = -1;
        holdcount = -1;
    elseif(holdcount == 3)
        chatter = 0;
        holdcount = 0;
    elseif(holdcount == 0)
        chatter = 1;
        holdcount = 1;
    else
        chatter = chatter*-1;
    end
end

function [t_out,j_out,x_a_out] = updateVehicle(r,t,j,x_a)
    global executeTime usingPointMassVehicle prevRef
    t_new = 0;
    j_new = 0;
    x_a_new = zeros(1,28);
    if(~isempty(usingPointMassVehicle) & usingPointMassVehicle == 1)
        [t_new,j_new,x_a_new] = pointMassVehicle(r, x_a, executeTime);
    else
        [t_new,j_new,x_a_new] = hybridVehicle(r, x_a(end,:), executeTime);
    end
    
    t_out = [t; t(end)+t_new(2:end)];
    j_out = [j; j(end)+j_new(2:end)];
    x_a_out = [x_a; x_a_new(2:end,:)];
    
    prevRef = r;
end

function [out, decelTerm, hysteresisCost] = cost(traj,plotCost)
    global Target TargetRad maxAcc prevRef executeTime planTime
    if(~exist('plotCost','var'))
        plotCost = 0;
    end
    if plotCost == 1
        exState = traj;
    else
        exState = interp1(traj(:,1),traj(:,3:end),traj(end,1) - planTime + executeTime);
    end
    if(any(isnan(exState)))
        traj(end,1) - planTime + executeTime
        traj
        AngTerm = inf;
        decelTerm = inf;
        out = inf;
        pause
        return
    end
    dist = max(norm(traj(end,3:5)-Target) - TargetRad, 0);
    if(isempty(prevRef))
        hysteresisCost = 0;
    else
        if(plotCost == 1)
            hysteresisCost = 0;
        else
            trajAtEndT = interp1(traj(:,1),traj(:,3:5),planTime - executeTime);
            prevRefTimeMatch = interp1(prevRef(:,1),prevRef(:,2:4),prevRef(1,1) + planTime);
            if(norm(trajAtEndT-prevRefTimeMatch)/norm(trajAtEndT-traj(1,3:5)) > 0.1)
                hysteresisCost = 10*max(norm(prevRef(end,2:4)-Target) - TargetRad, 0);
            else
                hysteresisCost = 0;
            end
        end
    end
    
    R = reshape(exState(end,7:15),[3,3]);
    phi = acos(dot(R*[0;0;1], [0;0;1]));
    phi = min(abs(phi), abs(phi-pi));
    poseTerm = 1;
    if(phi>=pi/2)
        poseTerm = Inf;
    end
    
    decelTerm = 0;
   
    out = dist*poseTerm + decelTerm + hysteresisCost;
%     out = dist*poseTerm + decelTerm;
end

function plotResults(t_a,x_a, x_o,r, rcost)

    disp("Plotting Results")
    global Target TargetRad obstacleSize
    figure(1);
    plot3(x_a(:,1), x_a(:,2), x_a(:,3), 'b');
    hold on;
    scatter3(x_a(:,1), x_a(:,2), x_a(:,3), 'b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    for i = 1:size(x_o)
        x_o_1 = cell2mat(x_o(i));
        plot3(x_o_1(:,3), x_o_1(:,5), x_o_1(:,7), 'r');
        [SX,SY,SZ] = sphere();
        SXT = SX*obstacleSize + x_o_1(1,3);
        SYT = SY*obstacleSize + x_o_1(1,5);
        SZT = SZ*obstacleSize + x_o_1(1,7);
        oSurface(i) = surface(SXT, SYT, SZT);
        set(oSurface(i),'FaceColor',[1 0 0], 'FaceAlpha',0.5,'FaceLighting','gouraud')
    end
    
    [SX,SY,SZ] = sphere();
    SXT = SX*TargetRad + Target(1);
    SYT = SY*TargetRad + Target(2);
    SZT = SZ*TargetRad + Target(3);
    tSurface = surface(SXT, SYT, SZT);
    set(tSurface,'FaceColor',[0 1 0], 'FaceAlpha',0.5,'FaceLighting','gouraud');
    hold off;
    
    figure(10);
    subplot(2,1,1)
    plot3(x_a(:,1), x_a(:,2), x_a(:,3), 'b');
    hold on;
    scatter3(r(:,2), r(:,3), r(:,4), 'b*');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    for i = 1:size(x_o)
        x_o_1 = cell2mat(x_o(i));
        plot3(x_o_1(:,3), x_o_1(:,5), x_o_1(:,7), 'r');
        [SX,SY,SZ] = sphere();
        SXT = SX*obstacleSize + x_o_1(1,3);
        SYT = SY*obstacleSize + x_o_1(1,5);
        SZT = SZ*obstacleSize + x_o_1(1,7);
        oSurface(i) = surface(SXT, SYT, SZT);
        set(oSurface(i),'FaceColor',[1 0 0], 'FaceAlpha',0.5,'FaceLighting','gouraud')
    end
    
    [SX,SY,SZ] = sphere();
    SXT = SX*TargetRad + Target(1);
    SYT = SY*TargetRad + Target(2);
    SZT = SZ*TargetRad + Target(3);
    tSurface = surface(SXT, SYT, SZT);
    set(tSurface,'FaceColor',[0 1 0], 'FaceAlpha',0.5,'FaceLighting','gouraud');
    hold off;
    
    figure(2)
    subplot(3,1,1); plot(t_a, x_a(:,1)); xlabel('t'); ylabel('x'); title('Quad pos')
    subplot(3,1,2); plot(t_a, x_a(:,2)); xlabel('t'); ylabel('y');
    subplot(3,1,3); plot(t_a, x_a(:,3)); xlabel('t'); ylabel('z');
    
    figure(3)
    subplot(3,1,1); plot(t_a, x_a(:,4)); xlabel('t'); ylabel('x');title('Quad vel')
    subplot(3,1,2); plot(t_a, x_a(:,5)); xlabel('t'); ylabel('y');
    subplot(3,1,3); plot(t_a, x_a(:,6)); xlabel('t'); ylabel('z');
    
    figure(4)
    hold on
    persistent numCost
    if(isempty(numCost))
        numCost = 1;
    end
    scatter(numCost, rcost)
    numCost = numCost +1;
    
    figure(5)
    subplot(3,1,1); plot(r(:,1), r(:,2))
    title('Reference Trajectory, includes discarded section')
    subplot(3,1,2); plot(r(:,1), r(:,3))
    subplot(3,1,3); plot(r(:,1), r(:,4))
    
    figure(6)
    subplot(3,1,1); plot(x_o_1(:,1), x_o_1(:,3)); xlabel('t'); ylabel('x'); title('obj pos')
    subplot(3,1,2); plot(x_o_1(:,1), x_o_1(:,5)); xlabel('t'); ylabel('y');
    subplot(3,1,3); plot(x_o_1(:,1), x_o_1(:,7)); xlabel('t'); ylabel('z');
    
    figure(7)
    scatter3(r(:,2), r(:,3), r(:,4),'r');
    title('Reference and Actual Trajectories');
    xlabel('x'); ylabel('y'); zlabel('z');
    
    figure(8)
    plot3(r(:,2), r(:,3), r(:,4),'r');
    title('Reference and Actual Trajectories');
    xlabel('x'); ylabel('y'); zlabel('z');
    hold on;
    plot3(x_a(:,1), x_a(:,2), x_a(:,3), 'b');
    
    figure(9)
    plot3(x_a(:,1),x_a(:,2),x_a(:,3));
    title('Vehicle Pose')
    xlabel('x'); ylabel('y'); zlabel('z');
    hold on;
    numPts = size(x_a);
    for i = 1:min(floor(numPts(1)/50),1):numPts(1)
        pose = reshape(x_a(i,7:15),[3 3])*[0; 0; 0.1];
        plot3([x_a(i,1),pose(1)+x_a(i,1)], [x_a(i,2),pose(2)+x_a(i,2)], [x_a(i,3),pose(3)+x_a(i,3)], 'r')
        pose = reshape(x_a(i,7:15),[3 3])*[0; 0.1; 0];
        plot3([x_a(i,1),pose(1)+x_a(i,1)], [x_a(i,2),pose(2)+x_a(i,2)], [x_a(i,3),pose(3)+x_a(i,3)], 'b')
        pose = reshape(x_a(i,7:15),[3 3])*[0.1; 0; 0];
        plot3([x_a(i,1),pose(1)+x_a(i,1)], [x_a(i,2),pose(2)+x_a(i,2)], [x_a(i,3),pose(3)+x_a(i,3)], 'g')
    end
end

function R = doubleCover(q)
    R = eye(3) + 2*q(1)*S(q(2:4)) + 2*S(q(2:4))^2;
end

function out = S(x)
    out = [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0];
end