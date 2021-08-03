function [out, outCost] = setBasedPlanner(x_a, x_o, costFun)
    global planTime sigma_max
    %Generate unsafe sets for obstacles
    xSpin = [sigma_max, sigma_max, -sigma_max, -sigma_max];
    ySpin = [sigma_max, -sigma_max, sigma_max, -sigma_max];
    numObj = size(x_o);
    U = cell(numObj(1),numObj(2)*4);
    for i=1:numObj(1)
        for j = 1:4
            x_o_mat = cell2mat(x_o(i));
            x_o_mat_size = size(x_o_mat);
            x_o_mat(:,9:10) = repmat([xSpin(j), ySpin(j)], [x_o_mat_size(1), 1]);
            [bb_t, bb_j, bb_x] = bouncingBallModel(x_o_mat(end,end-7:end), planTime);
            bballOut = [bb_t, bb_j, bb_x];
            [row, col] = size(bballOut);
            U(i,j) = mat2cell(bballOut,row,col);
        end
    end
    
    %Generate safe reachable map
    [phi,r] = genSafeTraj(x_a,U);
    f77 = figure(77);
    clf(f77)
    numRef = size(r);
    hold on;
    for i = 1:numRef(2)
        ref = cell2mat(r(i));
        plot3(ref(:,2),ref(:,3),ref(:,4));
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
    daspect([1 1 1])
    
    %Find "optimal" trajectory
    [~,numTraj] = size(phi);
    if(numTraj == 0)
        error("No safe trajectories, phi:",phi)
    end
    path = cell2mat(phi(1));
    closestTraj = [1, costFun(path)];
    for i=2:numTraj
        path = cell2mat(phi(i));
        [phir, ~] = size(path);
        if(phir == 1)
            error('bad phi')
        end
        cost = costFun(path);
        if(closestTraj(2) > cost)
            closestTraj = [i,cost];
        end
    end
    
    out = cell2mat(r(closestTraj(1)));
    outCost = closestTraj(2);
end