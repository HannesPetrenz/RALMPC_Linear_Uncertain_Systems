%% ========================================================================
%   Data Visualization Script - generatePlot.m
%
%   This script generates plots based on given data or mathematical 
%   functions. Modify the script to customize the appearance and 
%   behavior of the generated plots.
%
%   Author: Hannes Petrenz
%   Date: 03.04.2025 
%   Description:  
%     - Loads or defines data for visualization  
%     - Configures plot settings (labels, titles, grid, etc.)  
%     - Generates and displays the plot  
%     - Saves the plot if specified  
%   Steps:
%     1. Define or import the data.  
%     2. Set up figure properties and plotting parameters.  
%     3. Generate the plot using MATLAB's built-in functions.  
%     4. Customize axes, labels, legends, and other elements.  
%     5. Save the plot if needed.  
%% Clear Workspace & Initialize
clc;
close all;

% Set the directory path and search for all .mat files

[outputDir, mat] = getOrCreatePlotsDir();
%% Initialize variables for storing results
Numberofiteration = []; % Store the number of iterations for each dataset
X_RALMPC_N = {}; % Store X_RALMPC data for each dataset
J_RALMPC_N_all = {}; % Store all J_RALMPC values for each dataset
J_RLMPC_N = {}; % Store the last J_RLMPC value for each dataset
J_RALMPC_N = {}; % Store the last J_RALMPC value for each dataset
T_cpu_RALMPC = {}; % Store the average CPU time for RALMPC
T_cpu_RLMPC = {}; % Store the average CPU time for RLMPC

%% Loop through all data sets to load and process
for q = 1:length(mat)
    % Load the current dataset
    load(fullfile(mat(q).folder, mat(q).name));
    
    % Store relevant variables from each dataset
    X_RALMPC_N{q} = X_RALMPC;
    J_RALMPC_N_all{q} = J_RALMPC;
    J_RLMPC_N{q} = J_RLMPC{end};  % Last value of J_RLMPC
    J_RALMPC_N{q} = J_RALMPC{end}; % Last value of J_RALMPC
    
    % Store the number of iterations and compute average CPU times
    Numberofiteration(q) = N;
    T_cpu_RALMPC{q} = mean(vertcat(t_cpu_RALMPC{:})); % Average CPU time for RALMPC
    T_cpu_RLMPC{q} = mean(vertcat(t_cpu_RLMPC{:}));   % Average CPU time for RLMPC
end

% Compute the overall average CPU time for RAMPC
T_cpu_RAMPC = mean(vertcat(t_cpu_RAMPC{:}));

% Sort the datasets based on the number of iterations (in descending order)
[~, id] = sort(Numberofiteration, 'descend');

% Reorganize the data based on the sorted order
Numberofiteration = Numberofiteration(id); % Reorder the iteration counts
k = 1;
for i = id
    % Store reordered data
    X_RALMPC_N_{k} = X_RALMPC_N{i};
    J_RALMPC_N_{k} = J_RALMPC_N{i};
    J_RLMPC_N_{k} = J_RLMPC_N{i};
    T_cpu_RALMPC_{k} = T_cpu_RALMPC{i};
    T_cpu_RLMPC_{k} = T_cpu_RLMPC{i};
    J_RALMPC_N_all_{k} = J_RALMPC_N_all{i};
    
    k = k + 1; % Increment the index
end

% Final assignment to variables after reordering
X_RALMPC_N = X_RALMPC_N_;
J_RALMPC_N = J_RALMPC_N_;
J_RLMPC_N = J_RLMPC_N_;
T_cpu_RALMPC = T_cpu_RALMPC_;
T_cpu_RLMPC = T_cpu_RLMPC_;
J_RALMPC_N_all = J_RALMPC_N_all_;
%% Plot 1: Comparison Computation time over N
% MATLAB script for a bar plot comparing computation times of RAMPC and RALMPC
% The plot includes colorblind-friendly colors, a legend, and exports the figure to PDF and EPS.

% Ensure valid indices for the given values of N
[isMember, indices] = ismember([18, 12, 8, 6], Numberofiteration);

% Check if all desired indices were found
if any(~isMember)
    fprintf('Warning: Some specified iteration numbers (18, 12, 8, 6) were not found in Numberofiteration.\n');
end

% Data preparation: 
% Combine reference (RAMPC) and comparison (RALMPC) values for plotting
values = [T_cpu_RAMPC, cell2mat(T_cpu_RALMPC(indices))]; 

% Generate corresponding category labels for the x-axis (N values)
categories = cellstr(['N=25', compose("N=%d", Numberofiteration(indices))]);

% Define colorblind-friendly colors for the bars
refColor = [55, 126, 184] / 255;  % Dark Blue (#377EB8) for RAMPC
compColor = [230, 159, 0] / 255;  % Orange (#E69F00) for RALMPC

% Create a new figure for the bar plot
figure;
b = bar(values, 'FaceColor', 'flat'); % Create the bar plot without categorical alignment

% Assign colors to the bars: RAMPC gets a different color, RALMPC shares a common color
b.CData(1, :) = refColor;  % Color for RAMPC (first bar)
b.CData(2:end, :) = repmat(compColor, length(values)-1, 1);  % Color for RALMPC (remaining bars)

% Set the x-axis ticks and labels
xticks(1:length(categories));  % Set the ticks to match the number of categories
xticklabels(categories);  % Assign the dynamic category labels to the x-ticks

% Label the y-axis and adjust the font size for better readability
ylabel('Computation Time [s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on;  % Display grid lines for better visibility

% Set axis properties
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Add a manual legend with color patches
hold on;
h1 = patch(nan, nan, refColor);  % Create a patch for RAMPC
h2 = patch(nan, nan, compColor);  % Create a patch for RALMPC
legend([h1, h2], {'RAMPC', 'RALMPC'}, 'Location', 'best', 'Interpreter', 'latex');  % Create the legend

% Export the plot to PDF and EPS formats for use in LaTeX documents
pdfFile = fullfile(outputDir, 'PaperPlot_ComputationTime.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_ComputationTime.eps');

% Save the figure in both formats
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Command-line feedback: Notify the user about the saved files
fprintf('Plot saved as: %s\n', pdfFile);
fprintf('Plot saved as: %s\n', epsFile);

%% Plot 2: Comparison Cost over N
% MATLAB script for a bar plot comparing costs of RAMPC and RALMPC
% The plot includes colorblind-friendly colors, a horizontal line for optimal cost,
% a legend, and exports the figure to PDF and EPS.

% Ensure valid indices for the given values of N
[isMember, indices] = ismember([18, 12, 8, 6], Numberofiteration);

% Check if all desired indices were found
if any(~isMember)
    fprintf('Warning: Some specified iteration numbers (18, 12, 8, 6) were not found in Numberofiteration.\n');
end

% Data preparation: 
% Combine reference (RAMPC) and comparison (RALMPC) cost values for plotting
values = [cell2mat(J_RAMPC(end)), cell2mat(J_RALMPC_N(indices))]; 

% Generate corresponding category labels for the x-axis (N values)
categories = cellstr(['N=25', compose("N=%d", Numberofiteration(indices))]);

% Define colorblind-friendly colors for the bars
refColor = [55, 126, 184] / 255;  % Dark Blue (#377EB8) for RAMPC
compColor = [230, 159, 0] / 255;  % Orange (#E69F00) for RALMPC

% Create a new figure for the bar plot
figure;
b = bar(values, 'FaceColor', 'flat'); % Create the bar plot without categorical alignment

% Assign colors to the bars: RAMPC gets a different color, RALMPC shares a common color
b.CData(1, :) = refColor;  % Color for RAMPC (first bar)
b.CData(2:end, :) = repmat(compColor, length(values)-1, 1);  % Color for RALMPC (remaining bars)

% Set the x-axis ticks and labels
xticks(1:length(categories));  % Set the ticks to match the number of categories
xticklabels(categories);  % Assign the dynamic category labels to the x-ticks

% Add a horizontal line representing the optimal cost
yline(J_OS, '--', 'Optimal Cost', 'Color', '#0072B2', ...
    'Interpreter', 'latex', 'FontSize', 12, 'LineWidth', 2);

% Label the y-axis and adjust the font size for better readability
ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 14);

% Set the y-axis limits for better visualization of the data
ylim([100, 165]);

% Display grid lines for better plot visibility
grid on;

% Set axis properties for improved readability
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Add a manual legend with color patches
hold on;
h1 = patch(nan, nan, refColor);  % Create a patch for RAMPC
h2 = patch(nan, nan, compColor);  % Create a patch for RALMPC
legend([h1, h2], {'RAMPC', 'RALMPC'}, 'Location', 'northeast', 'Interpreter', 'latex');  % Create the legend

% Export the plot to PDF and EPS formats for use in LaTeX documents

pdfFile = fullfile(outputDir, 'PaperPlot_Cost.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_Cost.eps');

% Save the figure in both formats
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Command-line feedback: Notify the user about the saved files
fprintf('Plot saved as: %s\n', pdfFile);
fprintf('Plot saved as: %s\n', epsFile);
%% Plot 3: State \( x_1 \) Over Time
% This script generates a time series plot comparing x1 trajectories for RAMPC and RALMPC.
% The plot uses colorblind-friendly colors and is formatted for LaTeX.

% Set global font size for consistency across all plots
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontSize', 16);

% Find the index corresponding to iteration N=12
index = find(Numberofiteration == 12, 1);

% Check if N=12 exists in Numberofiteration
if isempty(index)
    fprintf('Warning: N=12 not found in Numberofiteration. Plotting skipped.\n');
    return;
end

% Extract state trajectories for RAMPC and RALMPC at the identified index
X_RALMPC_iteration = X_RALMPC_N{index};
state_RAMPC_iteration = X_RAMPC{end}; % Final state of RAMPC iteration

% Generate the time vector
time = 0:T_s:T_s*(length(state_RAMPC_iteration)-1);

% Define colorblind-friendly colors (Color Universal Design - CUD)
X_Init_Color = [0.2, 0.2, 0.8];    % Dark Blue for Initial Trajectory
RAMPC_Color = [0.0, 0.45, 0.7];     % Blue for RAMPC
OS_Color = [0.8, 0.4, 0.0];         % Orange for OS
RALMPC_First_Color = [0.6, 0.8, 0.2]; % Light Green for RALMPC (First Iteration)
RALMPC_Last_Color = [0.8, 0.17, 0.55]; % Magenta for RALMPC (Last Iteration)

% Create the figure
figure(3);
hold on;

% Ensure X_bar_initial exists before plotting
if exist('X_bar_inital', 'var') && ~isempty(X_bar_inital)
    plot(time, X_bar_inital(1, 1:length(X_OS)), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, ...
        'Color', X_Init_Color, 'DisplayName', 'Initial Trajectory');
else
    fprintf('Warning: X_bar_inital is missing or empty. Skipping Initial Trajectory plot.\n');
end

% Plot RAMPC (solid blue line)
plot(time, state_RAMPC_iteration(1,:), '--', 'LineWidth', 2, 'DisplayName', 'RAMPC', 'Color', RAMPC_Color);

% Plot OS (solid line with markers)
plot(time, X_OS(1,:), '-*', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'OS', 'Color', OS_Color);

% Plot RALMPC first iteration (light green)
state_RALMPC_first = X_RALMPC_iteration{1};
plot(time, state_RALMPC_first(1,:), '-', 'LineWidth', 2, 'DisplayName', 'RALMPC (First It.)', 'Color', RALMPC_First_Color);

% Plot RALMPC last iteration (magenta)
state_RALMPC_last = X_RALMPC_iteration{end};
plot(time, state_RALMPC_last(1,:), '-', 'LineWidth', 2, 'DisplayName', 'RALMPC (Last It.)', 'Color', RALMPC_Last_Color);

% Add a horizontal reference line at y = 0 (excluded from legend)
yline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Add the new horizontal line at y = -0.2 (included in legend)
yline(-0.2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Lower constraint $x_1$', 'HandleVisibility', 'on');

% Set axis labels with LaTeX formatting
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 18);
grid on;
xlim([0, 5.9]);

% Set axis properties for better readability
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex');

% Add a legend with an appropriate location
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 16);

% Stop adding more plots to the figure
hold off;

% Define output directory and file names
pdfFile = fullfile(outputDir, 'PaperPlot_timex1.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_timex1.eps');

% Export the plot to both formats
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Command-line feedback to confirm file saves
fprintf('Plot successfully saved:\n - PDF: %s\n - EPS: %s\n', pdfFile, epsFile);

%% Plot 4: Cost per Iteration
% This script generates a bar plot for cost reduction over iterations.
% It highlights the first few iterations and the final iteration with a gap in between.

% Create figure
figure(4);

% Ensure valid index for N=12 in Numberofiteration
index = find(Numberofiteration == 12, 1);

if isempty(index)
    fprintf('Warning: N=12 not found in Numberofiteration. Skipping plot.\n');
    return;
end

% Extract cost values for RALMPC over iterations
J_RALMPC_iteration = J_RALMPC_N_all{index}; 

% Convert cell array to numerical array
J_values = cell2mat(J_RALMPC_iteration);

% Ensure there are enough iterations for meaningful selection
numBars = length(J_values);
if numBars < 5
    fprintf('Warning: Not enough iterations to plot (found only %d). Adjusting selection.\n', numBars);
    selectedIndices = 1:numBars; % Take all available
else
    selectedIndices = [1:4, numBars]; % Select first 4 and last
end

% Extract the selected cost values
selectedValues = J_values(selectedIndices);

% Define the initial and optimal costs
if exist('J_0', 'var')
    initialCost = J_0;
else
    fprintf('Warning: J_0 (Initial Cost) is not defined. Using first iteration value instead.\n');
    initialCost = selectedValues(1);
end

if exist('J_OS', 'var')
    optimalCost = J_OS;
else
    fprintf('Warning: J_OS (Optimal Cost) is not defined. Skipping horizontal line.\n');
    optimalCost = NaN; % Placeholder
end

% Define x-axis positions and labels
xPositions = [0, 1:4, 5.5, 7]; % Position 5.5 represents the skipped section
xLabels = ["Initial", string(1:4), "...", string(numBars)]; % Labels with "..." for missing middle values
barValues = [initialCost, selectedValues]; % Include initial cost in data

% Define color scheme
compColor = [230, 159, 0] / 255; % Colorblind-friendly orange (#E69F00)

% Create the bar plot (excluding dots)
b = bar(xPositions([1, 2:5, 7]), barValues, 'FaceColor', 'flat'); % Enable color assignment
b.CData = compColor; % Assign colors to bars

% Adjust x-axis ticks and labels
xticks(xPositions);
xticklabels(xLabels);

% Set axis limits and grid
ylim([130, 170]);
grid on;

% Add a horizontal line for the optimal cost (if defined)
if ~isnan(optimalCost)
    yline(optimalCost, '--', 'Optimal Cost', 'Color', '#0072B2', ...
        'Interpreter', 'latex', 'FontSize', 12, 'LineWidth', 2);
end

% Add labels and title
xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 14);

% Improve appearance
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Define output directory and filenames

pdfFile = fullfile(outputDir, 'PaperPlot_Costiteration.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_Costiteration.eps');

% Export plot to PDF and EPS
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Display confirmation message
fprintf('Plot successfully saved:\n - PDF: %s\n - EPS: %s\n', pdfFile, epsFile);


%% Figure 5: 2D State Trajectory Plot for CDC
% This script plots the state trajectories for different control approaches.

% Create figure
figure(5); clf; hold on;

% Find the index of iteration N=12
index = find(Numberofiteration == 12, 1);
if isempty(index)
    fprintf('Warning: N=12 not found in Numberofiteration. Skipping plot.\n');
    return;
end

% Extract states for RALMPC and RAMPC iterations
if exist('X_RALMPC_N', 'var') && length(X_RALMPC_N) >= index
    X_RALMPC_iteration = X_RALMPC_N{index};
else
    fprintf('Warning: X_RALMPC_N is missing or too short.\n');
    return;
end

if exist('X_RAMPC', 'var') && ~isempty(X_RAMPC)
    state_RAMPC_iteration = X_RAMPC{end};
else
    fprintf('Warning: X_RAMPC is missing or empty.\n');
    return;
end

% Define colorblind-friendly colors
Polyhedron_Color = [0.8, 0.8, 1];   % Light blue for polyhedron
X_Init_Color = [0.2, 0.2, 0.8];     % Dark blue for initial trajectory
RAMPC_Color = [0.0, 0.45, 0.7];     % Blue (RAMPC)
OS_Color = [0.8, 0.4, 0.0];         % Orange (OS)
RALMPC_First_Color = [0.6, 0.8, 0.2]; % Light Green (RALMPC First Iteration)
RALMPC_Last_Color = [0.8, 0.17, 0.55]; % Magenta (RALMPC Last Iteration)

% Plot Polyhedron constraints (if variables exist)
if exist('H', 'var') && exist('S_inital', 'var') && exist('X_bar_inital', 'var')
    for k = 1:length(X_bar_inital)
        X_poly = Polyhedron(H, S_inital(k) * ones(size(H,1),1) + H * X_bar_inital(:,k));
        h = plot(X_poly, 'Color', Polyhedron_Color, 'LineStyle', '--', 'LineWidth', 1.5, 'FaceAlpha', 0.3);
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
    end
else
    fprintf('Warning: H, S_inital, or X_bar_inital is missing. Skipping Polyhedron plotting.\n');
end

% Plot initial trajectory
if exist('X_bar_inital', 'var')
    plot(X_bar_inital(1,:), X_bar_inital(2,:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, ...
        'Color', X_Init_Color, 'DisplayName', 'Initial Trajectory');
else
    fprintf('Warning: X_bar_inital is missing. Skipping initial trajectory plot.\n');
end

% Plot RAMPC trajectory
plot(state_RAMPC_iteration(1,:), state_RAMPC_iteration(2,:), '--', 'LineWidth', 2, ...
    'DisplayName', 'RAMPC', 'Color', RAMPC_Color);

% Plot OS trajectory
if exist('X_OS', 'var')
    plot(X_OS(1,:), X_OS(2,:), '-*', 'LineWidth', 2, 'MarkerSize', 6, ...
        'DisplayName', 'OS', 'Color', OS_Color);
else
    fprintf('Warning: X_OS is missing. Skipping OS plot.\n');
end

% Plot RALMPC first iteration (excluded from legend)
h1 = plot(X_RALMPC_iteration{1}(1,:), X_RALMPC_iteration{1}(2,:), '-', 'LineWidth', 2, ...
    'Color', RALMPC_First_Color);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

% Plot RALMPC last iteration (included in legend)
plot(X_RALMPC_iteration{end}(1,:), X_RALMPC_iteration{end}(2,:), '-', 'LineWidth', 2, ...
    'DisplayName', 'RALMPC (Last It.)', 'Color', RALMPC_Last_Color);

% Define rectangle limits
x1_min = -0.2; x1_max = 4.1;
x2_min = -5; x2_max = 5;
width = x1_max - x1_min;
height = x2_max - x2_min;

% Draw the rectangle (excluded from legend)
rectangle('Position', [x1_min, x2_min, width, height], 'EdgeColor', 'k', ...
    'LineWidth', 2.5, 'HandleVisibility', 'off');

% Set axis limits
xlim([-0.4, 4.3]); 
ylim([-5.3, 5.3]);

% Formatting
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
title('State Trajectory Comparison', 'Interpreter', 'latex', 'FontSize', 16);
grid on; 

% Improve appearance
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Add a legend (Polyhedron and first RALMPC iteration excluded)
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

hold off;

% Define output directory and filenames
pdfFile = fullfile(outputDir, 'PaperPlot_x1-x2.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_x1-x2.eps');

% Export plot to PDF and EPS
set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Display confirmation message
fprintf('Plot successfully saved:\n - PDF: %s\n - EPS: %s\n', pdfFile, epsFile);


%% Figure 6: Combined Bar Plot: Computation Time & Cost (Side-by-Side, Black Y-Labels & Cost Limit)
figure; hold on;

% Ensure valid indices exist
targetIterations = [18, 12, 8, 6];
[isMember, indices] = ismember(targetIterations, Numberofiteration);

if ~all(isMember)
    fprintf('Warning: Some target iterations not found in Numberofiteration.\n');
end

% Filter only found indices
validIndices = indices(isMember);
validIterations = Numberofiteration(validIndices);

% Ensure T_cpu_RALMPC exists
if exist('T_cpu_RALMPC', 'var') && length(T_cpu_RALMPC) >= max(validIndices)
    compTimeValues = [T_cpu_RAMPC, cell2mat(T_cpu_RALMPC(validIndices))];
else
    fprintf('Warning: T_cpu_RALMPC is missing or too short.\n');
    compTimeValues = NaN(1, length(validIndices) + 1);
end

% Ensure cost values exist
if exist('J_RAMPC', 'var') && exist('J_RALMPC_N', 'var') && length(J_RALMPC_N) >= max(validIndices)
    costValues = [cell2mat(J_RAMPC(end)), cell2mat(J_RALMPC_N(validIndices))];
else
    fprintf('Warning: Cost data is missing.\n');
    costValues = NaN(1, length(validIndices) + 1);
end

% Define x-axis category labels
categories = [{'N=25'}, compose("N=%d", validIterations)];

% Colorblind-friendly colors
refColor = [55, 126, 184] / 255;  % Dark Blue (#377EB8) for RAMPC
compColor = [230, 159, 0] / 255;  % Orange (#E69F00) for RALMPC

% Number of categories
numGroups = length(categories);
xPos = 1:numGroups; % X-axis positions

% ---- Plot Computation Time on Left Y-Axis ----
yyaxis left;
b1 = bar(xPos - 0.2, compTimeValues, 0.4, 'FaceColor', 'flat', 'EdgeColor', 'none'); % Left-shift bars
b1.CData(1, :) = refColor; % RAMPC (Computation Time) in blue
b1.CData(2:end, :) = repmat(compColor, length(compTimeValues)-1, 1); % RALMPC (Computation Time) in orange

ylabel('Computation Time [s]', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k'); % Y-label in black
set(gca, 'YColor', 'k'); % Left Y-axis color set to black

% ---- Plot Cost on Right Y-Axis ----
yyaxis right;
b2 = bar(xPos + 0.2, costValues, 0.4, 'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 2.5); % Right-shift bars
b2.CData(1, :) = refColor; % RAMPC (Cost) in blue
b2.CData(2:end, :) = repmat(compColor, length(costValues)-1, 1); % RALMPC (Cost) in orange

ylabel('Cost', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'k'); % Y-label in black
set(gca, 'YColor', 'k'); % Right Y-axis color set to black
ylim([100, 165]); % Limit the right Y-axis

% ---- Formatting ----
xticks(xPos);
xticklabels(categories);
grid on;
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% ---- Add a horizontal line for the optimal cost ----
if exist('J_OS', 'var')
    yline(J_OS, '--', 'Optimal Cost', 'Color', '#0072B2', ...
        'Interpreter', 'latex', 'FontSize', 12, 'LineWidth', 2);
else
    fprintf('Warning: J_OS is missing. Skipping optimal cost line.\n');
end

% ---- Create Corrected Legend ----
hold on;
h1 = bar(nan, nan, 'FaceColor', refColor, 'EdgeColor', 'none'); % RAMPC (Time)
h2 = bar(nan, nan, 'FaceColor', compColor, 'EdgeColor', 'none'); % RALMPC (Time)
h3 = bar(nan, nan, 'FaceColor', refColor, 'EdgeColor', 'k', 'LineWidth', 2.5); % RAMPC (Cost)
h4 = bar(nan, nan, 'FaceColor', compColor, 'EdgeColor', 'k', 'LineWidth', 2.5); % RALMPC (Cost)
h5 = plot(nan, nan, '--', 'Color', '#0072B2', 'LineWidth', 2); % Optimal Cost Line

legend([h1, h2, h3, h4, h5], ...
    {'RAMPC (Time)', 'RALMPC (Time)', 'RAMPC (Cost)', 'RALMPC (Cost)', 'Optimal Cost'}, ...
    'Location', 'southwest', 'Interpreter', 'latex');

% ---- Export as PDF and EPS ----
pdfFile = fullfile(outputDir, 'PaperPlot_Combined.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_Combined.eps');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Confirmation message
fprintf('Plot successfully saved:\n - PDF: %s\n - EPS: %s\n', pdfFile, epsFile);
%% Plot 7: Theta 

figure; hold on;

% Define colorblind-friendly colors from the CUD palette
color_HC0 = [86, 180, 233] / 255;   % Blue (Theta_0^{HC})
color_RALMPC = [230, 159, 0] / 255; % Orange (Theta_{RALMPC})
color_RAMPC = [0, 158, 115] / 255;  % Green (Theta_{RAMPC})

% ---- Check Data Existence Before Plotting ----
plotHandles = []; % Store handles for legend

% Plot Theta_0^{HC} - Blue (Ensuring it exists)
if exist('Theta_HC0', 'var') && ~isempty(Theta_HC0)
    h1 = plot(Theta_HC0, 'LineWidth', 2, 'Color', color_HC0, 'LineStyle', '-', 'DisplayName', 'Theta_0^{HC}');
    plotHandles = [plotHandles, h1]; % Add to legend
else
    fprintf('Warning: Theta_HC0 is missing. Skipping plot.\n');
end

% Plot Theta_{RALMPC} - Orange (Ensuring it exists)
if exist('Theta_HC_RALMPC', 'var') && ~isempty(Theta_HC_RALMPC) && ~isempty(Theta_HC_RALMPC{end})
    h2 = plot(Theta_HC_RALMPC{end}{end}, 'LineWidth', 2, 'Color', color_RALMPC, 'LineStyle', '-', 'DisplayName', 'Theta_{RALMPC}');
    plotHandles = [plotHandles, h2]; % Add to legend
else
    fprintf('Warning: Theta_HC_RALMPC is missing. Skipping plot.\n');
end

% Plot Theta_{RAMPC} - Green (Ensuring it exists)
if exist('Theta_HC_RAMPC', 'var') && ~isempty(Theta_HC_RAMPC) && ~isempty(Theta_HC_RAMPC{end})
    h3 = plot(Theta_HC_RAMPC{end}{end}, 'LineWidth', 2, 'Color', color_RAMPC, 'LineStyle', '-', 'DisplayName', 'Theta_{RAMPC}');
    plotHandles = [plotHandles, h3]; % Add to legend
else
    fprintf('Warning: Theta_HC_RAMPC is missing. Skipping plot.\n');
end

% ---- Axis Labels ----
xlabel('$\theta_1$', 'Interpreter', 'latex', 'FontSize', 14);  
ylabel('$\theta_2$', 'Interpreter', 'latex', 'FontSize', 14);  

% ---- Grid and Axis Formatting ----
grid on;
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% ---- Create Legend Only If Data Exists ----
if ~isempty(plotHandles)
    legend(plotHandles, ...
        {'$\Theta_0^{HC}$', '$\Theta_{RALMPC}^{HC}$', '$\Theta_{RAMPC}^{HC}$'}, ...
        'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);
else
    fprintf('Warning: No data was plotted. Legend omitted.\n');
end

% ---- Export as PDF and EPS ----
pdfFile = fullfile(outputDir, 'PaperPlot_Theta.pdf');
epsFile = fullfile(outputDir, 'PaperPlot_Theta.eps');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpdf', pdfFile);
print('-depsc', epsFile);

% Confirmation message
fprintf('Plot successfully saved:\n - PDF: %s\n - EPS: %s\n', pdfFile, epsFile);

hold off;


%% Help functions 
function [outputDir, matFiles] = getOrCreatePlotsDir()
    % Get the current directory (assumes script is in RALMPC_LTI_Paper)
    baseDir = pwd;  

    % Define the 'Plots' directory path
    outputDir = fullfile(baseDir, 'Plots');

    % Check if 'Plots' directory exists, if not, create it
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
        fprintf('Created directory: %s\n', outputDir);
    else
        fprintf('Directory already exists: %s\n', outputDir);
    end

    % Define the 'Data' directory path
    dataDir = fullfile(baseDir, 'Data');

    % Check if 'Data' directory exists
    if exist(dataDir, 'dir')
        % List all .mat files in the 'Data' directory
        matFiles = dir(fullfile(dataDir, '*.mat'));
        
        if isempty(matFiles)
            fprintf('No .mat files found in: %s\n', dataDir);
        else
            fprintf('Found %d .mat file(s) in: %s\n', length(matFiles), dataDir);
        end
    else
        fprintf('Data directory does not exist: %s\n', dataDir);
        matFiles = []; % Return empty if 'Data' folder is missing
    end
end