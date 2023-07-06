%% Init
clc; clear all; close all;
data_dir = 'E:\KAIST_VDC\sunmyung\calibration2\data\230619\';
% frame_id = '1687161014968505';
frame_id = '1687160853967675';

lidar_target = load(fullfile(data_dir, strcat(frame_id, '_lidar.mat')));
point_cloud = lidar_target.point_cloud;
pointCloudData = squeeze(point_cloud);


% Extract x, y, and z coordinates from point cloud
x = pointCloudData(:,1);
y = pointCloudData(:,2);
z = pointCloudData(:,3);

% Fit a plane to the point cloud using RANSAC
maxDistance = 0.01; % Maximum distance from point to plane to be considered an inlier
[model, inliers] = pcfitplane(pointCloud([x,y,z]), maxDistance);
plane = select(pointCloud([x,y,z]), inliers);
% pcshow(plane)
hold on
% Extract inlier points
pointsInPlane = [x(inliers), y(inliers), z(inliers)];

% rng(5,'twister');
% X = mvnrnd([0 0 0], [1 .2 .7; .2 1 0; .7 0 1],50);
plot3(x(inliers),y(inliers),z(inliers),'b.');
grid on;

[min_x, min_x_index] = min(x(inliers))
[max_x, max_x_index] = max(x(inliers))
[min_z, min_z_index] = min(z(inliers))
[max_z, max_z_index] = max(z(inliers))

axis square
view(-9,12);

[coeff,score,roots] = pca(pointsInPlane);
basis = coeff(:,1:2)

normal = coeff(:,3)

pctExplained = roots' ./ sum(roots)

[n,p] = size(pointsInPlane);
meanX = mean(pointsInPlane,1);
Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';

min_x_point = Xfit(min_x_index)
max_x_point = Xfit(max_x_index)
min_z_point = Xfit(min_z_index)
max_z_point = Xfit(max_z_index)

min_x_index_z = z(min_x_index);
max_x_index_z = z(max_x_index);

Xfit(min_x_index, 1)
Xfit(min_x_index, 2)
residuals = pointsInPlane - Xfit;

error = abs((pointsInPlane - repmat(meanX,n,1))*normal);
sse = sum(error.^2)

[xgrid,ygrid] = meshgrid(linspace(min(pointsInPlane(:,1)),max(pointsInPlane(:,1)),5), ...
                         linspace(min(pointsInPlane(:,2)),max(pointsInPlane(:,2)),5));
zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);

hold on

above = (pointsInPlane-repmat(meanX,n,1))*normal < 0;
below = ~above;
nabove = sum(above);
X1 = [pointsInPlane(above,1) Xfit(above,1) nan*ones(nabove,1)];
X2 = [pointsInPlane(above,2) Xfit(above,2) nan*ones(nabove,1)];
X3 = [pointsInPlane(above,3) Xfit(above,3) nan*ones(nabove,1)];
plot3(X1',X2',X3','-', pointsInPlane(above,1),pointsInPlane(above,2),pointsInPlane(above,3),'o', 'Color','b');
nbelow = sum(below);
X1 = [pointsInPlane(below,1) Xfit(below,1) nan*ones(nbelow,1)];
X2 = [pointsInPlane(below,2) Xfit(below,2) nan*ones(nbelow,1)];
X3 = [pointsInPlane(below,3) Xfit(below,3) nan*ones(nbelow,1)];
plot3(X1',X2',X3','-', pointsInPlane(below,1),pointsInPlane(below,2),pointsInPlane(below,3),'o', 'Color','r');

plot3(Xfit(min_x_index, 1), Xfit(min_x_index, 2), Xfit(min_x_index, 3), 'g.', 'MarkerSize', 15)
plot3(Xfit(max_x_index, 1), Xfit(max_x_index, 2), Xfit(max_x_index, 3), 'g.', 'MarkerSize', 15)
plot3(Xfit(max_z_index, 1), Xfit(max_z_index, 2), Xfit(max_z_index, 3), 'g.', 'MarkerSize', 15)

hold on
set(gca,'color','w')

% Define the given points
Point1 = [Xfit(min_x_index, 1), Xfit(min_x_index, 2), Xfit(min_x_index, 3)];
Point2 = [Xfit(max_x_index, 1), Xfit(max_x_index, 2), Xfit(max_x_index, 3)];
% Point3 = [Xfit(max_z_index, 1), Xfit(max_z_index, 2), Xfit(max_z_index, 3)];
% Point4 = [Xfit(max_z_index, 1), Xfit(max_z_index, 2), Xfit(max_z_index, 3)];
% planeNormal = [normal(1), normal(2), normal(3)];  % Normal vector of Plane 1
% planePoint = [Xfit(100, 1), Xfit(100, 2), Xfit(100, 3)];  % A point on Plane 1

mid_h = (min_x_index_z + max_x_index_z) / 2;

% find out missing part of checkerboard (upper or lower)
if ((max(z(inliers)) - mid_h) >= (mid_h - min(z(inliers))))
    fprintf("missing lower part of checkerboard")
    Point3 = [Xfit(max_z_index, 1), Xfit(max_z_index, 2), Xfit(max_z_index, 3)];
else
    fprintf("missing upper part of checkerboard")
    Point3 = [Xfit(min_z_index, 1), Xfit(min_z_index, 2), Xfit(min_z_index, 3)];
end

Point4 = findIntersection(Point1, Point3, Point2);
Point4

% Display the coordinates of Point4
plot3(Point4(1), Point4(2), Point4(3), 'g.', 'MarkerSize', 12)
hold on

xlim([-1.5 0.5])% xlim([min(x(inliers)) max(x(inliers))])
ylim([4.5 6.5])% ylim([min(y(inliers)) max(y(inliers))])
zlim([-1 1])

%% 3d plane to 2d
% Input
N = [normal(1), normal(2), normal(3)]; % Normal vector of the plane in 3D
P0 = [Xfit(20, 1), Xfit(20, 2), Xfit(20, 3)] % Point on the plane

% Translation: Move P0 to the origin
T = eye(4); % Identity matrix of size 4
T(1:3, 4) = -P0; % Translate the origin

% Rotation: Align normal vector N with the z-axis
theta = acos(N(3) / norm(N)); % Angle between normal vector and z-axis
axis2 = cross(N, [0, 0, 1]); % Axis of rotation
axis2 = axis2 / norm(axis2); % Normalized rotation axis
K = [0 -axis2(3) axis2(2); axis2(3) 0 -axis2(1); -axis2(2) axis2(1) 0]; % Cross product matrix
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2; % Rodrigues' rotation formula
R = [R, [0; 0; 0]; 0 0 0 1]; % Convert to homogeneous transformation

% Homogeneous Transformation
H = R * T; % The transformation matrix

% To apply the transformation to a point P in the plane, do:
% P_transformed = H * [P0'; 1] % Homogeneous coordinates
% P_transformed = P_transformed(1:3) % Discard the fourth dimension

% Project onto the 2D plane by ignoring the z-coordinate
% P_2D = P_transformed(1:2)

% Assuming your matrix of 3D points is P (size 801x3)
P = Xfit; % Replace with your actual points

% Convert to homogeneous coordinates
P_homogeneous = [P, ones(size(P, 1), 1)]; % Add a column of ones

% Apply the transformation to all points
P_transformed = (H * P_homogeneous')'; % Transpose needed for matrix multiplication

% Discard the fourth dimension and get back to 3D
P_transformed = P_transformed(:, 1:3);

% Project onto the 2D plane by ignoring the z-coordinate
P_2D = P_transformed(:, 1:2)

% Assuming P_2D is your matrix of 2D points
x = P_2D(:, 1); % Extract x-coordinates
y = P_2D(:, 2); % Extract y-coordinates

figure; % Create a new figure
[sorted_x, indices] = sort(x, 'descend');

% Select the top 20 points
x_max_top_20_x = sorted_x(1:60)
x_max_top_20_y = y(indices(1:60))

% Plot the points
figure; % Create a new figure
plot(x, y, '.'); % Plot all points
hold on;
% plot(x_max_top_20_x, x_max_top_20_y, 'ro', 'MarkerSize', 6); % Plot the top 20 points in red
xlabel('x'); % Label x-axis
ylabel('y'); % Label y-axis
title('Projected 2D points'); % Title for the plot
hold on
linearCoef = polyfit(x_max_top_20_x,x_max_top_20_y,1);
linearFit = polyval(linearCoef,x_max_top_20_x);
% plot(x_max_top_20_x, linearFit,'r-')
% xlabel('Weight'); ylabel('Proportion');

P_combined = [y, x];

% Sort rows first by y in descending order, then by x in ascending order
P_sorted = sortrows(P_combined, [-1, 2]);

% Extract the sorted x and y coordinates
y_sorted = P_sorted(:, 2);
x_sorted = P_sorted(:, 1);

P_sorted

[min_y_point2, min_y_index2] = min(P_2D(:, 2))
[max_x_point2, max_x_index2] = max(P_2D(:, 1))
[min_x_point2, min_x_index2] = min(P_2D(:, 1))

plot(P_2D(min_y_index2, 1), P_2D(min_y_index2, 2),'go')
plot(P_2D(max_x_index2, 1), P_2D(max_x_index2, 2),'go')
plot(P_2D(min_x_index2, 1), P_2D(min_x_index2, 2), 'go')
hold on

% line1
m1 = (P_2D(min_y_index2, 2) - P_2D(max_x_index2, 2)) / (P_2D(min_y_index2, 1) - P_2D(max_x_index2, 1));
b1 = P_2D(min_y_index2, 2) - m1 * P_2D(min_y_index2, 1);

% line2
m2 = (P_2D(min_y_index2, 2) - P_2D(min_x_index2, 2)) / (P_2D(min_y_index2, 1) - P_2D(min_x_index2, 1));
b2 = P_2D(min_y_index2, 2) - m2 * P_2D(min_y_index2, 1);


line_x1 = [P_2D(max_x_index2, 1) P_2D(min_y_index2, 1)]
line_y1 = [P_2D(max_x_index2, 2) P_2D(min_y_index2, 2)]

line_x2 = [P_2D(min_x_index2, 1) P_2D(min_y_index2, 1)]
line_y2 = [P_2D(min_x_index2, 2) P_2D(min_y_index2, 2)]

plot(line_x1, line_y1, 'r-')
hold on
plot(line_x2, line_y2, 'r-')
hold on

distance1 = [];
distance2 = [];

P_2D_1 = [];
P_2D_2 = [];

[numRows,numCols] = size(P_2D)
for i = 1:numRows
    temp_distance1 = abs(m1*P_2D(i, 1) - P_2D(i, 2) + b1) / sqrt(m1^2 + 1);
    temp_distance2 = abs(m2*P_2D(i, 1) - P_2D(i, 2) + b2) / sqrt(m2^2 + 1);
    distance1 = cat(1, distance1, [temp_distance1, i]);
    distance2 = cat(1, distance2, [temp_distance2, i]);
end

sorted_distance1 = sortrows(distance1, 1)
sorted_distance2 = sortrows(distance2, 1)

for i = 1:40
    plot(P_2D(sorted_distance1(i, 2), 1), P_2D(sorted_distance1(i, 2), 2), 'ro')
    sorted_distance1(i, 2)
    P_2D(sorted_distance1(i, 2), :)
    P_2D_1 = cat(1, P_2D_1, P_2D(sorted_distance1(i, 2), :));
    hold on
    plot(P_2D(sorted_distance2(i, 2), 1), P_2D(sorted_distance2(i, 2), 2), 'bo')
    sorted_distance2(i, 2)
    P_2D(sorted_distance2(i, 2), :)
    P_2D_2 = cat(1, P_2D_2, P_2D(sorted_distance2(i, 2), :));
    hold on
end

P_2D_1
P_2D_2

% Fit a line to the first cluster using polyfit
polyfit1 = polyfit(P_2D_1(:, 1), P_2D_1(:, 2), 1);

% Calculate the perpendicular to this line
m_perpendicular = -1/polyfit1(1);

% Define the error function to be minimized for the second cluster
% The line has the form y = m_perpendicular*x + b
% We want to find the value of b that minimizes the sum of squared offsets
errorFun = @(b) sum((P_2D_2(:, 2) - (m_perpendicular*P_2D_2(:, 1) + b)).^2);

% Use fminunc to find the value of b that minimizes the error function
b0 = 0;  % initial guess for b
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');
b_opt = fminunc(errorFun, b0, options);

% Now, the equations of the two lines are:
perpendicular1 = polyfit1(1)*P_2D_1(:, 1) + polyfit1(2)
perpendicular2 = m_perpendicular*P_2D_2(:, 1) + b_opt

plot(P_2D_1(:, 1), perpendicular1, 'g-')
plot(P_2D_2(:, 1), perpendicular2, 'g-')

hold on

multiplication_a = polyfit1(1)*m_perpendicular


[x_intersect, y_intersect] = line_intersection(polyfit1(1), polyfit1(2), m_perpendicular, b_opt);
fprintf('The intersection point is (%f, %f)\n', x_intersect, y_intersect);

plot(x_intersect, y_intersect, 'g.', 'MarkerSize', 20)


[np_line1_x, np_line1_y] = point_on_line_at_distance(polyfit1(1), polyfit1(2), x_intersect, y_intersect, 0.59)
fprintf('The new point on line1 is (%f, %f)\n', np_line1_x, np_line1_y);

[np_line2_x, np_line2_y] = point_on_line_at_distance(m_perpendicular, b_opt, x_intersect, y_intersect, -0.59)
fprintf('The new point on line2 is (%f, %f)\n', np_line2_x, np_line2_y);


plot(np_line1_x, np_line1_y, 'g.', 'MarkerSize', 20)
plot(np_line2_x, np_line2_y, 'g.', 'MarkerSize', 20)
%% find intersection point in 3d
function intersection_point = findIntersection(point1, point2, point3)
    % Extract coordinates of point1
    x1 = point1(1); y1 = point1(2); z1 = point1(3);
    
    % Extract coordinates of point2
    x2 = point2(1); y2 = point2(2); z2 = point2(3);
    
    % Extract coordinates of point3
    x3 = point3(1); y3 = point3(2); z3 = point3(3);
    
    % Calculate the direction vectors of line1 and line2
    direction_line1 = [x2 - x1, y2 - y1, z2 - z1];
    direction_line2 = [x3 - x2, y3 - y2, z3 - z2];
    
    % Scale the direction vector of line2 to pass through point1
    scale_factor = norm(direction_line1) / norm(direction_line2);
    direction_line4 = direction_line2 * scale_factor;
    
    % Calculate the intersection point
    intersection_point = point1 + direction_line4;
end

%% find intersection point in 2d
function [x_intersect, y_intersect] = line_intersection(a1, b1, a2, b2)
    % Check if the lines are parallel
    if a1 == a2
        error('The lines are parallel and do not intersect.');
    end

    % Calculate the intersection point
    x_intersect = (b2 - b1) / (a1 - a2);
    y_intersect = a1*x_intersect + b1;
end

%% find a new point far from intersect point by 0.59
function [x2, y2] = point_on_line_at_distance(a, b, x1, y1, d)
    % Calculate direction vector of the line
    dirVec = [1, a];

    % Normalize the direction vector
    dirVec = dirVec / norm(dirVec);

    % Scale the direction vector by d
    dirVec = d * dirVec;

    % Calculate the new point
    x2 = x1 + dirVec(1);
    y2 = y1 + dirVec(2);

    % Check if the new point is on the line
    if abs(y2 - (a*x2 + b)) > 1e-6
        % If not, we moved in the wrong direction, so subtract the direction vector instead
        x2 = x1 - dirVec(1);
        y2 = y1 - dirVec(2);
    end
end