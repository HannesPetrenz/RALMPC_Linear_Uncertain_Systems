%% ========================================================================
%   Robust Adaptive Learning Model Predictive Control (RALMPC) 
%   for Linear Systems with Parametric Uncertainties
%
%   This script initializes system parameters, constraints, cost functions, 
%   and disturbance models for implementing Robust Adaptive Learning MPC.
%
%   Author: Hannes Petrenz
%   Date: 02.04.2025 
%   Paper: Robust MPC for uncertain linear systems - Combining model adaptation and iterative learning
%   Description:  
%     - Defines system dynamics with parametric uncertainties  
%     - Configures disturbance sets and cost functions  
%     - Initializes state constraints and input constraints  
%     - Prepares hypercube parameter set for adaptive estimation 
%   Steps:
%     1. Perform offline computation for RALMPC to obtain the required parameters.
%     2. Compute the initial feasible solution.
%     3. Solve the robust optimal solution based on true system parameters.
%     4. Initialize and solve the Robust Adaptive MPC (RAMPC).
%     5. Solve RALMPC and RLMPC for various prediction horizons (N) and 
%        save the results.
%
% ========================================================================
%   Copyright Notice
%   ----------------
%MIT License
%
%Copyright (c) 2025 Hannes Petrenz and University of California, Berkeley
%
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.
%% ====================== Clear Workspace & Initialize =======================
clear; clc; close all;

%% ====================== General Settings =======================
% Number of iterations for adaptive learning
numberitertions = 20;  

% Array of prediction horizons for RALMPC
N_array_RALMPC = [6, 8, 12, 18]; 

% Prediction horizon for RAMPC
N_RAMPC = 25; 

% Prediction horizon for RAMPC
N_OS = 100;   
%% ====================== System Definition =======================
% Sampling time
T_s = 0.1;  

% Uncertain physical parameters (mass-spring-damper system)
m = 1;       % Mass
k = 1;       % Spring constant (uncertain parameter p1)
c = 0.2;     % Damping coefficient (uncertain parameter p2)

% True system parameters (actual but unknown)
k_star = 0.5;  
c_star = 0.3;

% Correct system matrices (for reference)
A_star = [1 T_s; -T_s*k_star/m 1 - T_s*c_star/m];
B_star = [0; T_s/m];

% System dimensions
n = size(A_star, 1);  % Number of states
m = size(B_star, 2);  % Number of inputs
p = 2;                % Number of uncertain parameters

%% ====================== Parametric Uncertainty Model =======================
% Scaling factor for uncertainty (bounded within [-1, 1])
parametric_unc = 0.5; 

% Uncertainty bounds
p1_unc = parametric_unc * (k/m);
p2_unc = parametric_unc * (c/m);

% Affine uncertainty description (A matrices)
A_0 = [1 T_s; -T_s*k/m 1 - T_s*c/m];  
A_1 = p1_unc * [0 0; T_s 0];  
A_2 = p2_unc * [0 0; 0 T_s];  

% Input matrices
B_0 = [0; T_s/m];  
B_1 = zeros(n, m);  
B_2 = zeros(n, m);  

%% ====================== Disturbance Model =======================
% Maximum disturbance
d_max = 0.2;  

% Rescaled disturbance for state-space representation
w_max = d_max/m * T_s;  

% Polyhedral disturbance set
W_V = [[0, w_max]; [0, -w_max]];
W = Polyhedron(W_V);

% Polyhedral constraints (H_w * w <= h_w)
H_w = [0, 1; 0, -1; 1, 0; -1, 0];
h_w = [w_max; w_max; 0; 0];

%% ====================== Cost Function Definition =======================
% State cost matrix
Q = [1 0; 0 1e-2];  

% Control effort cost matrix
R = 1e-1;  

%% ====================== Hypercube Parameter Set =======================
% Hypercube for parameter set B_p
HB_p = [1 0; -1 0; 0 1; 0 -1];
hB_p = 0.5 * [1; 1; 1; 1];  

% Create polyhedron representation
B_p = Polyhedron(HB_p, hB_p);

%% ====================== Constraints Definition =======================
% State and input constraints (F * x + G * u <= l)
F = [1/4.1 0; 1/-0.2 0; 0 1/5; 0 1/-5; zeros(2,2)];  
G = [zeros(4,1); 1/15; 1/-15];  
L = [F, G];  
l = ones(size(L, 1), 1);  

%% ====================== Initial Conditions =======================
% Initial states 
x0 = [4; 0];  
s0 = 0;  

%% ---------------- Offline Computation: Robust Adaptive Learning MPC ----------------
% This section precomputes necessary parameters for the Robust Adaptive Learning MPC.
% These parameters are used in the optimization problem.

fprintf("\n Step 1: Performing Offline Computation for RALMPC...\n");
[Theta_HC0, theta_bar0, eta_0, K, P, mu, H, rho_theta0, L_B, d_bar, c, c_max, ...
    L_cost, l_maxsteady, s_steady] = offlineComputation(Q, R, F, G, B_p, W, ...
    A_0, A_1, A_2, B_0, B_1, B_2, m, n, p, L, l);
fprintf(" Offline Computation Completed.\n\n");

%% ---------------- Compute Initial Solution ----------------
% The initial solution provides an initial sample set and worst-case cost to go.

fprintf(" Step 2: Computing Initial Feasible Solution...\n");
[SS_0, J_wc_0, X_bar_inital, S_inital, J_0] = get_Initalsolution( ...
    x0, s0, theta_bar0, B_p, eta_0, rho_theta0, L_B, d_bar, c, c_max, ...
    H, A_0, A_1, A_2, B_0, B_1, B_2, K, Q, R, P, F, G, m, n, p, ...
    L_cost, l_maxsteady, s_steady);
fprintf(" Initial Solution Computed.\n\n");

%% ---------------- Solve the Robust Optimal Solution ----------------
% This step solves the **robust optimal control problem** 
% with access to the true parameters.

fprintf("Step 3: Solving Robust Optimal Solution...\n");
[X_bar_OS, V_OS, S_OS, J_OS, X_OS, U_OS, time_OS] = solve_OS( ...
    x0, A_star, B_star, Q, R, P, K, F, G, d_bar, c_max, c, H, n, m, N_OS, W_V, T_s);
fprintf("Robust Optimal Solution Computed.\n\n");

%% ---------------- Solve the Robust Adaptive MPC (RAMPC) ----------------
% This section implements the **Robust Adaptive MPC (RAMPC)**, which adapts to 
% unknown parameters iteratively using set membership estimation.
% Reference: "Linear Robust Adaptive Model Predictive Control: 
% Computational Complexity and Conservatism."

fprintf("Step 4: Initializing Robust Adaptive MPC (RAMPC)...\n");
% RAMPC Settings
adpation = true; % Enable adaptation to unknown parameters
M = 10;         % Number of iterations for moving window hypercube update

% Initialize parameter uncertainty set (Delta) for adaptation
for i = 1:M  
    Delta{i} = Theta_HC0; % Initial uncertainty set from offline computation
end
fprintf("RAMPC Initialization Complete.\n\n");

% Solve the RAMPC for multiple iterations
fprintf("Step 5: Solving RAMPC Iterations...\n");
for i = 1:numberitertions
    fprintf("   Iteration %2d/%d... ", i, numberitertions);
    tStart = cputime; % Start CPU time measurement
    
    % Solve RAMPC optimization problem for the current iteration
    [X_hat_OL_RAMPC{i}, X_bar_OL_RAMPC{i}, V_OL_RAMPC{i}, S_OL_RAMPC{i}, J_RAMPC{i}, ...
        X_RAMPC{i}, U_RAMPC{i}, time_RAMPC{i}, Theta_HC_RAMPC{i}, t_solve_RAMPC{i}] = ...
        solve_RAMPC(x0, A_0, A_1, A_2, B_0, B_1, B_2, A_star, B_star, H_w, h_w, W_V, ...
        Q, R, P, K, F, G, d_bar, L_B, c_max, c, H, B_p, n, m, N_RAMPC, p, ...
        Theta_HC0, theta_bar0, eta_0, rho_theta0, Delta, mu, T_s, adpation);
    
    % Store computation time for performance evaluation
    t_cpu_RAMPC{i} = cputime - tStart;
    
    fprintf(" Completed in %.2f seconds.\n", t_cpu_RAMPC{i});
end
%% ---------------- Solve RALMPC and RLMPC for Different N ----------------
% This section iterates over different prediction horizons (N) and solves 
% the Robust Learning MPC (RLMPC) and Robust Adaptive Learning MPC (RALMPC).

fprintf("\nStep 6: Solving RALMPC and RLMPC for Different Prediction Horizons (N)...\n");

for i = 1:length(N_array_RALMPC)
    N = N_array_RALMPC(i);
    fprintf("\n➡️  Running Simulations for N = %d...\n", N);
    
    %---------------- Solve the RLMPC ----------------
    % RLMPC (Robust Learning MPC) is solved first without adaptation.
    
    fprintf("   Solving RLMPC (No Adaptation)...\n");
    
    % RLMPC Settings
    M = 0;           % No adaptation
    adpation = false;
    
    % Initialize parameter uncertainty set (Delta)
    Delta{1} = 0;
    
    % Set up the initial sample set
    SS = SS_0;
    Q_func = J_wc_0;
    
    % Solve the optimization problem
    [X_RLMPC, U_RLMPC, S_RLMPC, J_RLMPC, X_bar_OL_RLMPC, V_OL_RLMPC, ...
        S_OL_RLMPC, J_wc_RLMPC, time_RLMPC, Theta_HC_RLMPC, t_cpu_RLMPC] = ...
        solve_RALMPC(x0, SS, Q_func, A_0, A_1, A_2, B_0, B_1, B_2, ...
        A_star, B_star, H_w, h_w, W_V, Q, R, K, F, G, d_bar, L_B, c, H, ...
        B_p, n, m, N, p, Theta_HC0, theta_bar0, eta_0, rho_theta0, Delta, ...
        L_cost, l_maxsteady, T_s, numberitertions, adpation);
    
    fprintf("   RLMPC Solved Successfully!\n\n");

    %---------------- Solve the RALMPC ----------------
    % RALMPC (Robust Adaptive Learning MPC) includes adaptation to unknown parameters.
    
    fprintf("   Solving RALMPC (With Adaptation)...\n");
    
    % RALMPC Settings
    M = 10;           % Number of adaptation iterations
    adpation = true;
    
    % Initialize parameter uncertainty set (Delta) for adaptation
    for j = 1:M  
        Delta{j} = Theta_HC0;
    end
    
    % Set up the initial sample set
    SS = SS_0;
    Q_func = J_wc_0;
    
    % Solve the convex robust safe set implementation
    [X_RALMPC, U_RALMPC, S_RALMPC, J_RALMPC, X_bar_OL_RALMPC, V_OL_RALMPC, ...
        S_OL_RALMPC, J_wc_RALMPC, time_RALMPC, Theta_HC_RALMPC, t_cpu_RALMPC, ...
        t_solve_RALMPC] = ...
        solve_RALMPC(x0, SS, Q_func, A_0, A_1, A_2, B_0, B_1, B_2, ...
        A_star, B_star, H_w, h_w, W_V, Q, R, K, F, G, d_bar, L_B, c, H, ...
        B_p, n, m, N, p, Theta_HC0, theta_bar0, eta_0, rho_theta0, Delta, ...
        L_cost, l_maxsteady, T_s, numberitertions, adpation);
    
    fprintf("   RALMPC Solved Successfully!\n");

    %---------------- Save Results to File ----------------
    % Save all computed results in a .mat file for further analysis.
    saveResults(numberitertions, N, H, X_bar_inital, J_wc_0, J_0, S_inital, ...
                     X_RALMPC, X_RAMPC, X_OS, Theta_HC0, Theta_HC_RAMPC, ...
                     Theta_HC_RALMPC, J_RAMPC, J_RALMPC, J_OS, t_cpu_RAMPC, ...
                     t_cpu_RALMPC, t_solve_RALMPC, t_solve_RAMPC, T_s)
end

fprintf("\nAll RALMPC and RLMPC Simulations Completed Successfully!\n");


%% Help function
function saveResults(numberitertions, N, H, X_bar_inital, J_wc_0, J_0, S_inital, ...
                     X_RALMPC, X_RAMPC, X_OS, Theta_HC0, Theta_HC_RAMPC, ...
                     Theta_HC_RALMPC, J_RAMPC, J_RALMPC, J_OS, t_cpu_RAMPC, ...
                     t_cpu_RALMPC, t_solve_RALMPC, t_solve_RAMPC, T_s)
    % Get the current directory (assumes script is in RALMPC_LTI_Paper)
    baseDir = pwd;  

    % Define a random path before RALMPC_LTI_Paper (you can modify this logic if needed)
    randomFolder = strcat(baseDir, filesep, 'RandomFolder', num2str(randi([1000, 9999])), filesep);

    % Check if the random folder exists, if not, create it
    if ~exist(randomFolder, 'dir')
        mkdir(randomFolder);
        fprintf('Created directory: %s\n', randomFolder);
    else
        fprintf('Directory already exists: %s\n', randomFolder);
    end
    
    % Now create or check the 'Data' directory inside 'RALMPC_LTI_Paper'
    dataDir = fullfile(baseDir, 'Data');
    
    if ~exist(dataDir, 'dir')
        mkdir(dataDir);
        fprintf('Created directory: %s\n', dataDir);
    else
        fprintf('Directory already exists: %s\n', dataDir);
    end

    % Create the file name with iteration number and N value
    filename = fullfile(randomFolder, sprintf("data_Iter%d_N%d.mat", numberitertions, N));

    % Save the variables to the .mat file
    save(filename, "H", "X_bar_inital", "J_wc_0", "J_0", "S_inital", ...
        "X_RALMPC", "X_RAMPC", "X_OS", "Theta_HC0", "Theta_HC_RAMPC", ...
        "Theta_HC_RALMPC", "J_RAMPC", "J_RALMPC", "J_OS", "N", ...
        "numberitertions", "t_cpu_RAMPC", "t_cpu_RALMPC", "t_solve_RALMPC", "t_solve_RAMPC", "T_s");

    % Display the confirmation message
    fprintf("   Results saved to: %s\n", filename);
end