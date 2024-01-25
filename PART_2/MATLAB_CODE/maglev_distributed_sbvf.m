clear all
close all
clc

% Agents parameters 

A = [0 1;880.87 0];

B = [0; -9.9453];

C = [708.27 0];

D = zeros(1,2);

%% -- CHOOSE NETWORK TOPOLOGY -- %%
topology = 4;
[N,Ad,Aug,D1,L,G,Gdiag] = network_topology(topology);

%% Randomize initial conditions of follower nodes
% in range [-5, 5]

x0 = rand([2 N]).* 10 - 5;

%% Fixed initial conditions of follower nodes used for simulations in the reports
x0 = [ 2.0936   -2.2397    1.5510   -3.8100    4.5974    0.8527
    2.5469    1.7970   -3.3739   -0.0164   -1.5961   -2.7619 ];

%% -- CHOOSE LEADER NODE STEADY STATE BEHAVIOUR - %%
% run only desired section

%% LEADER NODE: Steady state constant

R0 = 2;

% State feedback gain
K0 = place(A,B,[0 -20]);
A0 = A-B*K0;

x0_leader = [R0 0]';

%% LEADER NODE: Steady state ramp

R0 = 1; % slope

% State feedback gain
K0 = acker(A,B,[0 0]);
A0 = A-B*K0;

x0_leader = [0 R0]';

%% LEADER NODE: Steady state sinusoid

w0 = 1;
R0 = 1; % amplitude, needs to have same name R0 for the last part of the script to work properly

% State feedback gain
K0 = place(A,B,[w0*1i -w0*1i]);
A0 = A-B*K0;

x0_leader = [R0 0]';


%% SBVF and OBSERVER parameters SIMULATION
% This section is used to calculate the values of the SBVF algorithm and
% observers, both distributed and local

%% --- CHOOSE ONE OF THE FOLLOWING ---
% Choose only one of the three following sections depending on what you
% want to do:
% - iterate for different values of ratio R/Q and then create plots
% - iterate for different values of the coupling gain c and then create
%   plots
% - just generate values (1 iteration) in order to run last section of this
%   script

%% ITERATE for different values of ratio R/Q
% uncomment line 152 and comment line 154
weights = [1/10 1/8 1/6 1/4 1/2 1:2:10];
c_gain = 2*ones(size(weights));
runs = 1;

tc_distributed_weights = zeros(size(weights));
tc_local_weights = zeros(size(weights));
peak_command_dist = zeros(size(weights));
peak_command_local = zeros(size(weights));

%% ITERATE for different values of the coupling gain c
% IT IS NECESSARY TO STRENGTHEN THE LOCAL OBSERVERS EIGENVALUES
% uncomment line 154 and comment line 152

c_gain = 1:8;
weights = 5*ones(size(c_gain));
runs = 1;

tc_distributed_weights = zeros(size(c_gain));
tc_local_weights = zeros(size(c_gain));
peak_command_dist = zeros(size(c_gain));
peak_command_local = zeros(size(c_gain));

%% NO ITERATION
% uncomment line 152 and comment line 154

tc_distributed_weights = 0;
tc_local_weights = 0;
peak_command_dist = 0;
peak_command_local = 0;
ratio = 5;
weights = ratio;
c_gain = 2;
runs = 1;

%% SBVF and OBSERVER parameters CALCULATION
% Always run this

for k = 1:runs
    for z = 1:size(weights,2)
        
        % --- SBVF PARAMETERS ---
        
        lambda = eig(L+G);

        cmin = 1/(2*min(real(lambda)));
        c = c_gain(z)*cmin;

        % Q R weights
        q = 1;
        r = q/weights(z);
        Q = q*eye(2); % input weight
        R = r;        % state weight

        % Distributed controller riccati equation
        P = are(A0,B*R^-1*B',Q);
        K = R^-1 * B' * P;

        % Verify conditions
        % Tracking error convergence
        Ac = kron(eye(N),A0)-kron(c*(L+G),B*K);
        eig_Ac = eig(Ac)

        
        % --- OBSERVERS PARAMETERS ---

        % DISTRIBUTED observer gain (riccati equation)
        P_dist = are(A0',C'*R^-1*C,Q);
        F_dist = P_dist * C' * R^-1;

        % LOCAL observer gain (riccati equation)
        % P = are(A0',C'*R^-1*C,Q);
        % F_local = P * C' * R^-1;
        
        % VALUES FOR SIMULATIONS 
        F_local = place(A0',C',[-4 -2])';
        % VALUES FOR COUPLING GAIN SIMULATIONS
        % F_local = place(A0',C',[-6 -3])';

        % LEADER observer as standard Luenberger observer 
        % xdot = Ax + Bu + c*L_leader*y_tilde
        L_leader = place(A0',C',[-5 -4])';

        % IDEAL observer to obtain a perfect extreme performance in order to
        % produce the ideal steady state behaviour to be used in computations
        % L_ideal = place(A0',C',

        Tsim = 30;

        k = 1;
        std_dev(k) = 1;
        threshold = 5/100 * R0;
        
        % Distributed
        out_dist = sim('sbvf_model_2019b_distributed.slx');
        tc1 = conv_time(out_dist.x1_tilde.Data,threshold,out_dist.tout);
        tc2 = conv_time(out_dist.x1_tilde.Data,threshold,out_dist.tout);
        tc_distributed_weights(z) = tc_distributed_weights(z) + max(tc1,tc2);
        % Peak command u value for second state (convergence to constant)
        peak_command_dist(z) = peak_command_dist(z) + max(out_dist.command.Data(2,1,:));
        
        % Local
        out_local = sim('sbvf_model_2019b_local.slx');
        tc1 = conv_time(out_local.x1_tilde_local.Data,threshold,out_local.tout);
        tc2 = conv_time(out_local.x1_tilde_local.Data,threshold,out_local.tout);
        tc_local_weights(z) = tc_local_weights(z) + max(tc1,tc2);
        % Peak command u value for second state (convergence to constant)
        peak_command_local(z) = peak_command_local(z) + max(out_local.command.Data(2,1,:));

    end
end

tc_distributed_weights = tc_distributed_weights / runs;
tc_local_weights = tc_local_weights / runs;
peak_command_local = peak_command_local / runs;
peak_command_dist = peak_command_dist / runs;

[ weights', c_gain', tc_distributed_weights', tc_local_weights']
fprintf("R/Q ratio     c/cmin     distributed     local\n");

%% PLOTS for R/Q variations
% run this ONLY IF you previously chose "ITERATE for different values of
% R/Q ratio

% PLOT Q/R
figure(1);
hold on;
grid on;
title("Convergence time for variation of the weights R and Q");
xlabel("R/Q ratio");
ylabel("Convergence time");
p1 = plot(weights, tc_distributed_weights);
p2 = plot(weights, tc_local_weights);
ymin = min([min(tc_local_weights) min(tc_distributed_weights)]);
ymax = max([max(tc_local_weights) max(tc_distributed_weights)]);
legend([p1 p2],["Distributed" "Local"]);
axis([weights(1), weights(size(weights,2)), ymin-1 ymax+1]);

% PLOT PEAK COMMAND
figure(2)
hold on;
grid on;
title("Peak command value for different R/Q ratio");
xlabel("R/Q ratio");
ylabel("Peak command value");
p1 = semilogx(weights, peak_command_dist);
p2 = semilogx(weights, peak_command_local);
ymin = min([min(peak_command_local) min(peak_command_dist)]);
ymax = max([max(peak_command_local) max(peak_command_dist)]);
legend([p1 p2],["Distributed" "Local"]);
axis([weights(1), weights(size(weights,2)), ymin-1 ymax+1]);
set(gca, 'XScale', 'log');

%% PLOTS for coupling gain variations (c/cmin)
% run this ONLY IF you previously chose "ITERATE for different values of
% coupling gain c

% PLOT c_gain
figure(1);
hold on;
grid on;
title("Convergence time for variation of c/cmin");
xlabel("c/cmin");
ylabel("Convergence time");
p1 = plot(c_gain, tc_distributed_weights);
p2 = plot(c_gain, tc_local_weights);
ymin = min([min(tc_local_weights) min(tc_distributed_weights)]);
ymax = max([max(tc_local_weights) max(tc_distributed_weights)]);
legend([p1 p2],["Distributed" "Local"]);
axis([c_gain(1), c_gain(size(c_gain,2)), ymin-1, ymax+1]);

% PLOT PEAK COMMAND
figure(2)
hold on;
grid on;
title("Peak command value for different c/cmin ratio");
xlabel("c/cmin ratio");
ylabel("Peak command value");
p1 = plot(c_gain, peak_command_dist);
p2 = plot(c_gain, peak_command_local);
ymin = min([min(peak_command_local) min(peak_command_dist)]);
ymax = 300; % it diverges so would go to 10^300 and the graphic is not readable anymore
legend([p1 p2],["Distributed" "Local"]);
axis([c_gain(1), c_gain(size(c_gain,2)), ymin-1, ymax]);

%% DISTRIBUTED OBSERVER SIMULINK PROJECT
% State estimation error convergence
% Open simulink

open ('sbvf_model_2019b_distributed.slx')

%% LOCAL OBSERVER SIMULINK PROJECT
% State estimation error convergence
% Open simulink

open ('sbvf_model_2019b_local.slx')

%% SIMPLE SIMULATION 
% Use this if you just want to run a simulation with small realistic errors
% just to analyze and compare the behaviours of the two kind of systems
% Can run this ONLY IF you previously chose "NO ITERATION"

clear k;
clear std_dev;
clear tc_distributed;
clear tc_local;

threshold = 5/100 * R0;
k = 1;
std_dev(k) = 2;

% Distributed
out_dist = sim('sbvf_model_2019b_distributed.slx');
[tc1, tc1_sim] = conv_time(out_dist.x1_tilde.Data,threshold,out_dist.tout);
[tc2, tc2_sim] = conv_time(out_dist.x1_tilde.Data,threshold,out_dist.tout);
tc_distributed = max(tc1,tc2);
tc_dist_sim = max(tc1_sim,tc2_sim);

% Local
out_local = sim('sbvf_model_2019b_local.slx');
[tc1, tc1_sim] = conv_time(out_local.x1_tilde_local.Data,threshold,out_local.tout);
[tc2, tc2_sim] = conv_time(out_local.x1_tilde_local.Data,threshold,out_local.tout);
tc_local = max(tc1,tc2);
tc_local_sim = max(tc1_sim,tc2_sim);

[ std_dev', tc_distributed, tc_local]
fprintf("std_dev     distributed     local\n");
[ mean(mean(abs(out_dist.x2_diff.Data(:,:,tc_dist_sim:end)))), mean(mean(abs(out_local.x2_diff.Data(:,:,tc_local_sim:end)))) ]
fprintf("diff_ideal_dist      diff_ideal_local\n");

%% CONVERGENCE TIME and MEAN DIFFERENCE with respect to ideal behaviour
% comparisons while varying noise's standard deviation on output measurements
% Can run this ONLY IF you previously chose "NO ITERATION"


% Convergence time, threshold calculated as 5% * amplitude of steady state
% behaviour of the signal
threshold = 5/100 * R0;

std_dev = 0:1:10;
samples = size(std_dev,2);
tc_distributed = zeros(samples,1);
tc_local = zeros(samples,1);
diff_mean_dist_x1 = zeros(samples,1);
diff_mean_dist_x2 = zeros(samples,1);
diff_mean_local_x1 = zeros(samples,1);
diff_mean_local_x2 = zeros(samples,1);

runs = 1;
for z = 1:runs
    for k = 1:samples
        % Calculate errors

        % Distributed
        out_dist = sim('sbvf_model_2019b_distributed.slx');
        [tc1, tc_dist_sim1] = conv_time(out_dist.x1_tilde.Data,threshold,out_dist.tout);
        [tc2, tc_dist_sim2] = conv_time(out_dist.x1_tilde.Data,threshold,out_dist.tout);
        tc_distributed(k) = tc_distributed(k) + max(tc1,tc2);
        tc_dist_sim = max(tc_dist_sim1,tc_dist_sim2);
        diff_mean_dist_x1(k) = diff_mean_dist_x1(k) + mean(mean(abs(out_dist.x1_diff.Data(:,:,tc_dist_sim:end))));
        diff_mean_dist_x2(k) = diff_mean_dist_x2(k) + mean(mean(abs(out_dist.x2_diff.Data(:,:,tc_dist_sim:end))));

        % Local
        out_local = sim('sbvf_model_2019b_local.slx');
        [tc1, tc_local_sim1] = conv_time(out_local.x1_tilde_local.Data,threshold,out_local.tout);
        [tc2, tc_local_sim2] = conv_time(out_local.x1_tilde_local.Data,threshold,out_local.tout);
        tc_local(k) = tc_local(k) + max(tc1,tc2);
        tc_local_sim = max(tc_local_sim1,tc_local_sim2);
        diff_mean_local_x1(k) = diff_mean_local_x1(k) + mean(mean(abs(out_local.x1_diff.Data(:,:,tc_local_sim:end))));
        diff_mean_local_x2(k) = diff_mean_local_x2(k) + mean(mean(abs(out_local.x2_diff.Data(:,:,tc_local_sim:end))));

    end
end

% calculate average above all runs, to derive the general behaviour and set
% aside occasional outliers due to specific values of noise generated by
% casualty

tc_distributed = tc_distributed / runs;
tc_local = tc_local / runs;
diff_mean_local_x1 = diff_mean_local_x1 / runs;
diff_mean_local_x2 = diff_mean_local_x2 / runs;
diff_mean_dist_x1 = diff_mean_dist_x1 / runs;
diff_mean_dist_x2 = diff_mean_dist_x2 / runs;

[ std_dev', tc_distributed, tc_local, diff_mean_dist_x1, diff_mean_dist_x2, diff_mean_local_x1, diff_mean_local_x2]
fprintf("std_dev     distributed     local    mean_diff_dist_x1     mean_diff_dist_x2     mean_diff_local_x1     mean_diff_local_x2\n");

%% PLOTS
% create plots, can run this ONLY IF you chose NO ITERATION before and then
% run previous section

% PLOT Convergence time with respect to std deviation on same plot
f1 = figure(1)
title("Convergence time");
xlabel("Standard deviation");
ylabel("Convergence time");
hold on;
grid on;
p1 = plot(std_dev, tc_distributed);
p2 = plot(std_dev, tc_local);
ymin = min([min(tc_local) min(tc_distributed)]);
ymax = max([max(tc_local) max(tc_distributed)]);
legend([p1 p2],["Distributed" "Local"]);
axis([std_dev(1), std_dev(size(std_dev,2)), ymin-0.1, ymax+0.1]);

% PLOT Difference with respect to ideal behaviour
f2 = figure(2)
title("Steady state mean difference to ideal behaviour");
xlabel("Standard deviation");
ylabel("Difference to ideal");
hold on;
grid on;
p1 = plot(std_dev, diff_mean_local_x1);
p2 = plot(std_dev, diff_mean_local_x2);
p3 = plot(std_dev, diff_mean_dist_x1);
p4 = plot(std_dev, diff_mean_dist_x2);
legend([p1 p2 p3 p4],["Local x2" "Local x1" "Distributed x2" "Distributed x1"]);
