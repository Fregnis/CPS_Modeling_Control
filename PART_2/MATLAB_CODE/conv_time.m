function [tc, tc_sim] = conv_time(x_tilde,threshold,sim_time)
%CONV_TIME Calculate convergence time as the time at which the estimation 
% error goes below a given threshold

% if the algorithm doesn't provide a solution it means the system doesn't
% converge within the simulation time
% Convergence time is set to sim_time to represent an "infinite" time
tc = sim_time(size(sim_time,1));
tc_sim = size(sim_time,1);

for t = 1:size(x_tilde,3)
    done = 1;
    for j = 1:size(x_tilde,2)
        if abs(x_tilde(1,j,t)) > threshold 
            % if there is still at least a node k with relative x_tilde(1,k,t) 
            % not below threshold the convergence is not reached
            done = 0;
        end
    end
    if (done == 1)
        tc_sim = t;
        tc = sim_time(t);
        return;
    end
end

return;

