clear all; clc;

%% Positions in 2D
positionAPs = [1,    -9,    14;
               -9,    10,   -2];
positionSTA = [4;
               3];

% Estimated distances in 2D
distances_est = vecnorm(positionAPs - positionSTA, 2, 1);

% Estimate the user's position using trilateration
if randi([0, 1])
    % either use least_squares_solution
    positionSTA_est = least_squares_solution(positionAPs, distances_est);
else
    % or use trilateration2D_solver
    positionSTA_est = trilateration2D_solver(positionAPs, distances_est);
end

% Plot the positions of the true/estimated STA and APs in the xy-plane
figure
plot_Nodes(positionAPs,positionSTA,positionSTA_est,distances_est)
xlim([-40 35]), ylim([-40 35]);
grid on
title('2D Trilateration')

%% 3D
% Positions in 3D
positionAPs_3D = [ 1, -9,  14,  5;
                  -9, 10,  -2,  5;
                   0,  0,  0,   3];  % 3 axes and 4 APs
positionSTA_3D = [4;
                  3;
                  2];

% Estimated distances in 3D
distances_est_3D = vecnorm(positionAPs_3D - positionSTA_3D, 2, 1);

% Estimate the user's position using trilateration
positionSTA_est_3D = trilateration3D_solver(positionAPs_3D, distances_est_3D);

% Plot the positions of the true/estimated STA and APs in the xyz-space 
figure
plot_Nodes_3D(positionAPs_3D, positionSTA_3D, positionSTA_est_3D, distances_est_3D)
xlim([-12 12]), ylim([-12 12]), zlim([-5 10])
grid on
title('3D Trilateration')

%% Local functions for 2D
% 2D: Find the estimated position based on the least-squares solution of Qz=b
function sol = least_squares_solution(positionAPs, distances_est)
dim = 2; % dim = 2 means 2D, while dim = 3 means 3D
num_APs = numel(distances_est); % number of APs
Q = zeros(num_APs-1, dim);
b = zeros(num_APs-1, 1);
for i=1:num_APs-1
    Q(i,:) = [2*(positionAPs(1,1) - positionAPs(1,i+1)), 2*(positionAPs(2,1) - positionAPs(2,i+1))];
    b(i) = positionAPs(1,1)^2 - positionAPs(1,i+1)^2 + positionAPs(2,1)^2 - positionAPs(2,i+1)^2 + distances_est(i+1)^2 - distances_est(1)^2;
end
sol = ((Q.'*Q)^(-1)) * (Q.' * b);
end

% 2D: Find the estimated position
function estimated_position = trilateration2D_solver(AP_positions, distances)
    N = size(AP_positions, 2);
    if N < 3
        error('At least 3 APs are required for 2D trilateration.');
    end

    % Use AP1 as reference
    x1 = AP_positions(1,1);
    y1 = AP_positions(2,1);
    d1 = distances(1);

    A = zeros(N-1, 2);
    b = zeros(N-1, 1);

    for i = 2:N
        xi = AP_positions(1,i);
        yi = AP_positions(2,i);
        di = distances(i);

        A(i-1, :) = 2 * [xi - x1, yi - y1];
        b(i-1) = d1^2 - di^2 + (xi^2 - x1^2) + (yi^2 - y1^2);
    end

    % Solve A x = b
    estimated_position = A \ b;
end


% 2D: Plot nodes
function plot_Nodes(positionAPs,positionSTA,positionSTA_est,distances_est)
plot(positionAPs(1,:),positionAPs(2,:),'kv','LineWidth',2) % Plot APs
hold on
plot(positionSTA(1),positionSTA(2),'kp','LineWidth',2) % Plot the user
hold on
plot(positionSTA_est(1),positionSTA_est(2),'bs','LineWidth',1.5,'MarkerSize',12); % Plot the estimated position of the user
hold on;
num_APs = size(positionAPs,2);
angles = 0:2*pi/360:2*pi;
for i = 1:num_APs
    x = distances_est(i) * cos(angles) + positionAPs(1,i);
    y = distances_est(i) * sin(angles) + positionAPs(2,i);
    plot(x,y,'k--','LineWidth',0.2); % plot trilateration circles
    hold on;
end
legend({'APs','User','Estimated user','Trilateration circles'}, ...
        'Location','northwest','FontSize',11)
xlabel('x-axis (m)'),ylabel('y-axis (m)');
axis equal
end

%% Local functions for 3D 
% 3D: Find the estimated position
function estimated_position = trilateration3D_solver(AP_positions, distances)
    N = size(AP_positions, 2);
    if N < 4
        error('At least 4 APs are required for 3D trilateration.');
    end

    % Use AP1 as reference
    x1 = AP_positions(1,1); y1 = AP_positions(2,1); z1 = AP_positions(3,1);
    d1 = distances(1);

    A = zeros(N-1, 3);
    b = zeros(N-1, 1);

    for i = 2:N
        xi = AP_positions(1,i);
        yi = AP_positions(2,i);
        zi = AP_positions(3,i);
        di = distances(i);

        A(i-1, :) = 2 * [xi - x1, yi - y1, zi - z1];
        b(i-1) = d1^2 - di^2 + (xi^2 - x1^2) + (yi^2 - y1^2) + (zi^2 - z1^2);
    end

    % Solve A x = b
    estimated_position = A \ b;
end

% 3D: Plot nodes
function plot_Nodes_3D(positionAPs, true_pos, est_pos, distances)
    % Plot APs
    plot3(positionAPs(1,:), positionAPs(2,:), positionAPs(3,:), 'kv', 'LineWidth', 2);
    hold on

    % Plot the user
    plot3(true_pos(1), true_pos(2), true_pos(3), 'kp', 'LineWidth', 2);

    % Plot the estimated position of the user
    plot3(est_pos(1), est_pos(2), est_pos(3), 'bs', 'MarkerSize', 12, 'LineWidth', 1.5);

    % Draw spheres around each AP
    [X, Y, Z] = sphere(30); % 30 = smoother
    for i = 1:size(positionAPs,2)
        r = distances(i);
        x = r * X + positionAPs(1,i);
        y = r * Y + positionAPs(2,i);
        z = r * Z + positionAPs(3,i);
        surf(x, y, z, ...
            'FaceAlpha', 0.1, ...
            'EdgeColor', 'none', ...
            'FaceColor', [1, 0.6, 0.2]);  % transparent light orange
        hold on
    end

    % Format plot
    legend({'APs', 'User', 'Estimated user', 'Trilateration spheres'}, ...
            'Location', 'northeast', 'FontSize', 10);
    xlabel('x-axis (m)'), ylabel('y-axis (m)'), zlabel('z-axis (m)');
    axis equal
    view(3)
    rotate3d on
end
