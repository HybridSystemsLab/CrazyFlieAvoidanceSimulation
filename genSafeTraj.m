function [phis,rs] = genSafeTraj(x_a, U)
    global usingPointMassVehicle obstacleSize objColStep ObjPlot55 chatterSim
    
    f56 = figure(56);
    clf(f56);
    
    Mset = [];
    spherePts = 6;
    [spx,spy,spz] = sphere(spherePts);
    spx = obstacleSize*reshape(spx(:,1:spherePts),[],1);
    spy = obstacleSize*reshape(spy(:,1:spherePts),[],1);
    spz = obstacleSize*reshape(spz(:,1:spherePts),[],1);
    if(~isempty(usingPointMassVehicle) && usingPointMassVehicle == 1)
        %Generate trajectories
        [Mset,all_rs] = pointMassMobility(x_a);
    else
        [Mset,all_rs] = hybridModelMobility(x_a);
    end
    numSafe = 0;
    numObj = size(U);
    %Remove unsafe sets
    numTraj = size(Mset);
    phis = {};
    rs = {};
    for i = 1:numTraj
        unsafeFlag = 0;
        traj = cell2mat(Mset(i));
        if(chatterSim == 1)
            for j = 1:numObj(1)
                x_o = cell2mat(U(j,1));
                if(any(vecnorm(traj(:,3:5) - [x_o(1,3), x_o(1,5), x_o(1,7)],2,2)<obstacleSize))
                    unsafeFlag = 1;
                end
            end
        else
            for j = 1:numObj(1)
                x_o = cell2mat(U(j,1));
                [~,u_t,~] = unique(x_o(:,1));
                xos = [];
                xos(:,:,1) = x_o(u_t,:);
                % Resample obstacle trajectories to match time stamps of first obstacle trajectory
                for k = 2:4
                    x_o = cell2mat(U(j,k));
                    [~,u_t,~] = unique(x_o(:,1));
                    u_x_o = x_o(u_t,:);
                    if(u_t > 1)
                        xos(:,:,k) = interp1(u_x_o(:,1),u_x_o,xos(:,1,1));
                    else
                        xos(:,:,k) = u_x_o;
                    end
                end

                hullPts = zeros(spherePts*(spherePts-1)*(4*(objColStep+1)),3);
                numXos = size(xos(:,1,1));
                for t = 1:objColStep:numXos(1)-objColStep
                    %Unroll points in sphere arond each point
                    pt_count = 0;
                    % Build convext hull over subset of obstacle samapled points
                    for k = 1:4
                        [sp_num_pts,~] = size(spx);
                        for l = 0:objColStep
                            hullPts(pt_count+1:pt_count+sp_num_pts,:) = [spx,spy,spz]+[xos(t+l,3,k),xos(t+l,5,k),xos(t+l,7,k)];
                            pt_count = pt_count + sp_num_pts;
                        end
                        if(i == 1 && ObjPlot55 == 1)
                            figure(55);
                            hold on;
                            color = ['r','g','b','k'];
                            hullPtsNotOrig = (hullPts(:,1:3) ~= [0 0 0]);
                            scatter3(hullPts(hullPtsNotOrig(:,1),1),hullPts(hullPtsNotOrig(:,1),2),hullPts(hullPtsNotOrig(:,1),3),color(k));
                            figure(56);
                            hold on
                            hullPtsNotOrig = (hullPts(:,1:3) ~= [0 0 0]);
                            scatter3(hullPts(hullPtsNotOrig(:,1),1),hullPts(hullPtsNotOrig(:,1),2),hullPts(hullPtsNotOrig(:,1),3),color(k));
                        end
                    end

                    % Check if any trajectory point is in hull
                    ptsInHull = inhull(traj(:,3:5),hullPts);
                    if(sum(ptsInHull) > 0)
                        unsafeFlag = 1;
                        break
                    end
                end
            end
        end
        
        if unsafeFlag == 0
            % Does not collide with any
            numSafe = numSafe+1;
            phis(numSafe) = Mset(i);
            rs(numSafe) = all_rs(i);
        end
    end
end

% returns path with structure [t,p_x, p_y, p_z, v_x, v_y, v_z, R(3x3), omeag(3x1)]
function [trajs,rs] = pointMassMobility(x_a)
    global planTime planStep maxAcc
%Discritize reachable space
    sampleAngs = 20; %number of angles to test
    sampleAccSteps = 5; %number of distances to test at each angle

    angstep = 2*pi/sampleAngs;
    numTraj = 0;
    trajs = cell(ceil(pi/angstep)*sampleAngs*sampleAccSteps,1);
    rs = cell(ceil(pi/angstep)*sampleAngs*sampleAccSteps,1);
%     figure;
    for i = (-pi/2):angstep:(pi/2)
        for j = angstep:angstep:2*pi
            for k = 1:sampleAccSteps
                r_x = [maxAcc/k * cos(i) * cos(j), maxAcc/k * cos(i) * sin(j), maxAcc/k * sin(i), 0,0,0];
                r = [(planStep:planStep:planTime)', repmat(r_x, [size((planStep:planStep:planTime)'),1])];
                if(any(isnan(r(:))))
                    r
                    r_x
                    pause
                end
                [time,jump,phi] = pointMassVehicle(r, x_a, planTime);
                if(any(isnan(phi(:))))
                    r
                    r_x
                    phi
                    pause
                end
                [row,col] = size([time,jump,phi]);
                numTraj = numTraj+1;
                trajs(numTraj) = mat2cell([time,jump,phi],row,col);
                [row,col] = size(r);
                rs(numTraj) = mat2cell(r,row,col);
            end
        end
    end
end

function [trajs, rs] = hybridModelMobility(x_a)
    global d_xy MaxMotorThrust planTime FirstRunEject ConstThrust PlanningThrustMult
    numMVals = 6; % +1
    
    %Generate reference trajs for vehicle
    persistent TMs i
    if(FirstRunEject == 1)
        disp("first run")
        FirstRunEject = 2;
        TMs = zeros(4,1);
        m1 = 1;
        m2 = 1;
        m3 = 1;
        m4 = 1;
        roll = (-m1-m2+m3+m4)*MaxMotorThrust*d_xy;
        pitch = (+m1-m2-m3+m4)*MaxMotorThrust*d_xy;
        yaw = (m1-m2+m3-m4)*MaxMotorThrust*d_xy;
        thrust = (m1+m2+m3+m4)*MaxMotorThrust;
        i = 1;
        TMs(:,i) = [thrust, roll, pitch, yaw];
        
        pause
    elseif(isempty(TMs)||FirstRunEject == 2)
        disp("genning TM ball")
        FirstRunEject = 0;
        TMs = zeros(4, 2+(numMVals-1)*(numMVals+1)^2);
        i = 2;
        TMs(:,1) = [0; 0; 0; 0];
        for trst=1/(numMVals+2):1/(numMVals+2):(1-1/(numMVals+2))
            diffVals = min(1-trst, trst)/2; % base is over 2
            for m1m3diff=-diffVals:2*diffVals/numMVals:diffVals
                for m2m4diff=-diffVals:2*diffVals/numMVals:diffVals
                    m1 = (trst+m1m3diff)*PlanningThrustMult;
                    m2 = (trst+m2m4diff)*PlanningThrustMult;
                    m3 = (trst-m1m3diff)*PlanningThrustMult;
                    m4 = (trst-m2m4diff)*PlanningThrustMult;
                    if(any([m1 m2 m3 m4] > 1)|| any([m1 m2 m3 m4] < 0))
                        disp('out of motor range')
                        pause
                        continue
                    end
                    roll = (-m1-m2+m3+m4)*MaxMotorThrust*d_xy;
                    pitch = (+m1-m2-m3+m4)*MaxMotorThrust*d_xy;
                    yaw = (m1-m2+m3-m4)*MaxMotorThrust*d_xy;
                    thrust = (m1+m2+m3+m4)*MaxMotorThrust;
                    if(abs(yaw) <= 1e-6)
                        yaw = 0;
                    else
                        disp('bad yaw');
                        pause
                        continue
                    end
                    i = i+1;
                    TMs(:,i) = [thrust; roll; pitch; yaw];
                end
            end
        end
    else
        disp("No gen, using last TM ball")
    end
    disp("finished ref gen");
    numTraj = 0;
    trajs = cell(i,1);
    rs = cell(i,1);
    figure(76)
    hold on;
    for j = 1:i
        % Find reference trajectories
        if(ConstThrust == 1)
            PitchRollTime = planTime/3;
            [t,simJump,r] = simHybridVehicleOnly(x_a(1:18)', TMs(1,j), TMs(2:4,j), [0 PitchRollTime]);
            [t2,simJump2,r2] = simHybridVehicleOnly(r(end,1:end-4)', TMs(1,j), [0; 0; 0], [PitchRollTime planTime]);
            t = [t;t(end)+t2(3:end)];
            simJump = [simJump;simJump(end)+simJump2(3:end)];
            r = [r;r2(3:end,:)];
        else
            [t,simJump,r] = simHybridVehicleOnly(x_a(1:18)', TMs(1,j), TMs(2:4,j), planTime);
        end
        r = r(:,1:18);
        if(t(1) == t(2))
            t = t(2:end);
            simJump = simJump(2:end);
            r = r(2:end,:);
        end
        %Use hybrid sym as traj out put
        phi = r;
        jump = simJump;
        time = t;
        
        dt = gradient(t);
        [~,dv] = gradient(r(:,4:6));
        vdot = dv./dt;
        numrpts = size(t);
        adot = zeros(numrpts(1),3);
        jdot = zeros(numrpts(1),3);
        [~,dw] = gradient(r(:,16:18));
        wdot = dw./dt;
        r = [t, r(:,1:6), vdot, adot, jdot, r(:,7:end), wdot];

        %Save reference and resulting trajectories
        [row,col] = size([time,jump,phi]);
        if(row == 1)
            continue
        end
        numTraj = numTraj+1;
        trajs(numTraj) = mat2cell([time,jump,phi],row,col);
        [row,col] = size(r);
        rs(numTraj) = mat2cell(r,row,col);
    end
end