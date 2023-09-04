function [transformed_data,transform_matrix, Rotation_GT, Translation_GT] = transformation(data,x,y,z,translation)
%% Data must be always 3xN
s = size(data);
if s(1) ~= 3
    data = data';
end

% Convert the angle in radians
% Trig functions expect their arguments to be in radians, not in degrees
x = deg2rad(x);
y = deg2rad(y);
z = deg2rad(z);

Rx = [1      0      0;
    0 cos(x) -sin(x);
    0 sin(x) cos(x)];
Ry = [cos(y)  0 sin(y);
    0       1      0;
    -sin(y) 0 cos(y)];
Rz = [cos(z) -sin(z) 0;
    sin(z) cos(z)  0;
    0      0       1];

% Rotation matrix
transform_matrix = Rz*Ry*Rx;      

%  Composition --> Rotation matrix + Translation
transform_matrix = [transform_matrix(1,1),transform_matrix(1,2),transform_matrix(1,3),translation(1);
                    transform_matrix(2,1),transform_matrix(2,2),transform_matrix(2,3),translation(2);
                    transform_matrix(3,1),transform_matrix(3,2),transform_matrix(3,3),translation(3);
                      0                       0                       0                  1]; 

% Ground Truth Rotation / Translation
Rotation_GT = transform_matrix(1:3,1:3); % Rotation Ground Truth Value
Translation_GT = transform_matrix(1:3,4); % Translation Ground Truth Value

cols = size(data,2); % How many columns
rows_one = ones(1,cols); % Create a row of ones of length equal to the 
% number of columns 

data = [data;rows_one];    % Add a last row of ones to make it of
% dimension 4 x .., so it can be multiplied with 4x4 transform matrix

transformed_data = transform_matrix * data;
transformed_data = transformed_data(1:3,:);  % Get rid of ones (last row)
% and return the 3D coordinates

end