% x_a = [p_x_a; p_y_a; p_z_a; v_x_a; v_y_a; v_z_a; R_11-R_31; R_12-R_32; R_13-R_33; omega]
function [t,j,out] =  pointMassVehicle(r, x_a, simTime)
    global HEQOptions QUAD_REFERENCE
    QUAD_REFERENCE = [0,r(1,2:end);r];
    TSPAN = [0, simTime];
    JSPAN = [0, 1];
    [t,j,out_new] =  HyEQsolver(@f,@g,@C,@D,[x_a(end,1:18),zeros(1,7)],TSPAN,JSPAN,1,HEQOptions);
    if(any(isnan(out_new(:))))
        disp("Fail in pointMassVehicle, found NAN")
        r
        x_a
        [t,j,out_new]
        pause
    end
    out = out_new(:,1:end-7);
end

function out = f(x)
    global QUAD_REFERENCE
    if(x(end) < QUAD_REFERENCE(end,1))
        r = interp1(QUAD_REFERENCE(:,1),QUAD_REFERENCE(:,2:end),x(end))';
    else
        r = QUAD_REFERENCE(end,2:end)';
    end
    out = [x(4:6); r(1:3) - [0;0;9.8]; zeros(12,1); x(end-3:end-1); r(1:3); 1];
end

function out = g(x)
    out = x;
end

function out = C(x)
    out = 1;
end

function out = D(x)
    out = 0;
end