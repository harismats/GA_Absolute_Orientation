function [New_Rotor, New_Rotor_Reverse, Rnew] = ...
    CM_FGmatrices(F_matrix,G_matrix)

% F and G are 3x3 matrices
% MUST BE DOUBLE OTHERWISE IT IS NOT WORKING, i.e. Left and right operands
% must have components of the same class
F_matrix = double(F_matrix);
G_matrix = double(G_matrix);

%% TRICKY - Sometimes the Clifford Toolbox instead of really small values,
% e.g. 10e-12 gives 0. When a 0 exists and we form the matrix F, the
% Clifford Toolbox calculates 1 element (when we have two 0s) or 2 elemetns
% (when we have one 0) in vector1_Fclif, or vector2_Fclif, or
% vector2_Fclif. So instead on X1e1 + X2e2 + X3e3 (X1,X2,X3 numbers), it
% gives, e.g. X2e2 + X3e3 or X1e1 + X3e3, or X2e2,..etc. Thus, for
% computational purposes we set it after the calculation as:
if any(F_matrix==0)
    F_matrix(F_matrix==0) = 10e-30;
end

%% Vectors - F matrix, G matrix
% Convert them to Clifford

% Matrix F
vector1_Fclif = F_matrix(1,1) * e1 + F_matrix(2,1) * e2 + F_matrix(3,1)*e3;
vector2_Fclif = F_matrix(1,2) * e1 + F_matrix(2,2) * e2 + F_matrix(3,2)*e3;
vector3_Fclif = F_matrix(1,3) * e1 + F_matrix(2,3) * e2 + F_matrix(3,3)*e3;

% Matrix G
vector1_Gclif = G_matrix(1,1) * e1 + G_matrix(2,1) * e2 + G_matrix(3,1)*e3;
vector2_Gclif = G_matrix(1,2) * e1 + G_matrix(2,2) * e2 + G_matrix(3,2)*e3;
vector3_Gclif = G_matrix(1,3) * e1 + G_matrix(2,3) * e2 + G_matrix(3,3)*e3;

%% Calculate Reciprocals
vector1_Fclif_up = - (1 / abs(wedge(vector1_Fclif, vector2_Fclif, vector3_Fclif))) * wedge(vector2_Fclif, vector3_Fclif) * wedge(e1,e2,e3);
vector2_Fclif_up =   (1 / abs(wedge(vector1_Fclif, vector2_Fclif, vector3_Fclif))) * wedge(vector1_Fclif, vector3_Fclif) * wedge(e1,e2,e3);
vector3_Fclif_up = - (1 / abs(wedge(vector1_Fclif, vector2_Fclif, vector3_Fclif))) * wedge(vector1_Fclif, vector2_Fclif) * wedge(e1,e2,e3);

%% TRICKY - It happens very rarely the Clifford Multivector Toolbox to give
% values +/- Inf for the Reciprocals. In this case replace the +/- Inf with
% very large numbers --> SOLVED: BECAUSE OF THE RANK (=1)
% vector1_FC_up
coef_vector1_FC_up = cell2mat(coefficients(vector1_Fclif_up));

if any(coef_vector1_FC_up == Inf) || any(coef_vector1_FC_up == -Inf)
    coef_vector1_FC_up(coef_vector1_FC_up == Inf) = 10e30;
    coef_vector1_FC_up(coef_vector1_FC_up == -Inf) = -10e30;
    coef_vector1_FC_up = [coef_vector1_FC_up, 10e30];
    vector1_Fclif_up = coef_vector1_FC_up(1) *e1 + coef_vector1_FC_up(2) *e2 +...
        coef_vector1_FC_up(3) *e3;
end


%% Calculate Invariants
Inv1 = vector1_Fclif_up * vector1_Gclif + vector2_Fclif_up * vector2_Gclif + vector3_Fclif_up * vector3_Gclif;
Inv2 = wedge(vector2_Fclif_up,vector1_Fclif_up) * wedge(vector1_Gclif, vector2_Gclif) + wedge(vector3_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector3_Gclif) + wedge(vector3_Fclif_up, vector2_Fclif_up) * wedge(vector2_Gclif, vector3_Gclif);
Inv3 = wedge(vector3_Fclif_up, vector2_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector2_Gclif, vector3_Gclif);
Inv_all = Inv1 + Inv2 + Inv3;

%% Rotor
Rotor_reverse = 1 + Inv_all; % the rotor reverse
Rotor_reverse_normalized = unit(Rotor_reverse); % the rotor reverse normalized

New_Rotor = reverse(Rotor_reverse_normalized); % the rotor (normalized)
New_Rotor_Reverse = Rotor_reverse_normalized;

R_Rtilde_3D = New_Rotor * New_Rotor_Reverse; % check that is = 1

%% Forming the columns of the matrix
f1_F = New_Rotor * e1 * New_Rotor_Reverse;
f2_F = New_Rotor * e2 * New_Rotor_Reverse;
f3_F = New_Rotor * e3 * New_Rotor_Reverse;

%% Extract coefficients
coefficients_f1_F = coefficients(f1_F);
coefficients_f1_F = cell2mat(coefficients_f1_F);

coefficients_f2_F = coefficients(f2_F);
coefficients_f2_F = cell2mat(coefficients_f2_F);

coefficients_f3_F = coefficients(f3_F);
coefficients_f3_F = cell2mat(coefficients_f3_F);

% In the last iteration it might happen that the coefficients have only one
% value, i.e. '1', INSTEAD OF 3.
if length(coefficients_f1_F) && length(coefficients_f2_F) && length(coefficients_f3_F)  == 1
    try
        % Keep the coefficients of interest
        F1 = coefficients_f1_F(1:3);
        F2 = coefficients_f2_F(1:3);
        F3 = coefficients_f3_F(1:3);
    catch
        F1 = coefficients_f1_F(1);
        F2 = coefficients_f2_F(1);
        F3 = coefficients_f3_F(1);
        % Form matrix F (IDENTITY MATRIX AS THE TRANSFORMATION HAS BEEN
        % RECOVERED
        F = [F1 0 0; 0 F2 0; 0 0 F3];
    end
    % F1--> 1 element, F2 --> 2 elements, F3 --> 2 elements
elseif length(coefficients_f1_F) == 1
    F1 = coefficients_f1_F(1);
    F2 = coefficients_f2_F(1:2);
    F3 = coefficients_f3_F(1:2);
    F = [F1 0 0; 0 F2; 0 F3];
    % F1--> 2 element, F2 --> 3 elements, F3 --> 3 elements
elseif length(coefficients_f1_F) == 2
    F1 = coefficients_f1_F(1:2);
    F2 = coefficients_f2_F(1:3);
    F3 = coefficients_f3_F(1:3);
    F_matrix = [F1 0; F2; F3];
    % F1--> 3 elements, F1 --> 2 elements, F3 --> 2 elements
elseif length(coefficients_f2_F) == 1
    F2 = coefficients_f2_F(1);
    F1 = coefficients_f1_F(1:2);
    F3 = coefficients_f3_F(1:2);
    F = [F1 0; 0 F2 0; 0 F3];
else
    F1 = coefficients_f1_F(1:3);
    F2 = coefficients_f2_F(1:3);
    F3 = coefficients_f3_F(1:3);

    F = [F1;F2;F3]; % by row
end

%% Trasnpose it
Rnew = F;

end