%% Both points must be 3xN

function F = matrix_F_optimized(u, v)
%% Checking size of data
[rows_u, ~] = size(u);
[rows_v, ~] = size(v);

if rows_u ~=3
    u=u';
end

if rows_v ~=3
    v=v';
end

%% Covnvert to double
u = double(u);
v = double(v);

%% Subtract the mean from each point set as this is what GA method does
u = u - mean(u,2);
v = v - mean(v,2);

%% Vector containing e1,e2,e3
vec = [e1,e2,e3];

%% Multiplication of u and v (COLUMN) - Convert it to Clifford
clif_u = vec * u;
clif_v = vec * v;

%% Forming Matrix F --> F_{jk} = sum_{i=1}^{i=N} (u_i . e_j)(v_i . e_k)
%% 1st Column
F11a = bsxfun(@scalar_product, clif_u(1,:), e1);
F11a_coef = cell2mat(coefficients(F11a));
F11b = bsxfun(@scalar_product, clif_v(1,:),e1);
F11b_coef = cell2mat(coefficients(F11b));
if isempty(F11a_coef) || isempty(F11b_coef)
    F11 = 0;
else 
    F11 = sum(F11a_coef.*F11b_coef);
end

F21a = bsxfun(@scalar_product, clif_u(1,:), e2);
F21a_coef = cell2mat(coefficients(F21a));
F21b = bsxfun(@scalar_product, clif_v(1,:), e1);
F21b_coef = cell2mat(coefficients(F21b));
if isempty(F21a_coef) || isempty(F21b_coef)
    F21 = 0;
else 
    F21 = sum(F21a_coef.*F21b_coef);
end

F31a = bsxfun(@scalar_product, clif_u(1,:),e3);
F31a_coef = cell2mat(coefficients(F31a));
F31b = bsxfun(@scalar_product, clif_v(1,:),e1);
F31b_coef = cell2mat(coefficients(F31b));
if isempty(F31a_coef) || isempty(F31b_coef)
    F31 = 0;
else 
    F31 = sum(F31a_coef.*F31b_coef);
end

%% 2nd Column
F12a = bsxfun(@scalar_product, clif_u(1,:),e1);
F12a_coef = cell2mat(coefficients(F12a));
F12b = bsxfun(@scalar_product, clif_v(1,:),e2);
F12b_coef = cell2mat(coefficients(F12b));
if isempty(F12a_coef) || isempty(F12b_coef)
    F12 = 0;
else 
    F12 = sum(F12a_coef.*F12b_coef);
end

F22a = bsxfun(@scalar_product, clif_u(1,:),e2);
F22a_coef = cell2mat(coefficients(F22a));
F22b = bsxfun(@scalar_product, clif_v(1,:),e2);
F22b_coef = cell2mat(coefficients(F22b));
if isempty(F22a_coef) || isempty(F22b_coef)
    F22 = 0;
else 
    F22 = sum(F22a_coef.*F22b_coef);
end

F32a = bsxfun(@scalar_product, clif_u(1,:),e3);
F32a_coef = cell2mat(coefficients(F32a));
F32b = bsxfun(@scalar_product, clif_v(1,:),e2);
F32b_coef = cell2mat(coefficients(F32b));
if isempty(F32a_coef) || isempty(F32b_coef)
    F32 = 0;
else 
    F32 = sum(F32a_coef.*F32b_coef);
end

%% 3rd Column
F13a = bsxfun(@scalar_product, clif_u(1,:),e1);
F13a_coef = cell2mat(coefficients(F13a));
F13b = bsxfun(@scalar_product, clif_v(1,:),e3);
F13b_coef = cell2mat(coefficients(F13b));
if isempty(F13a_coef) || isempty(F13b_coef)
    F13 = 0;
else 
    F13 = sum(F13a_coef.*F13b_coef);
end

F23a = bsxfun(@scalar_product, clif_u(1,:),e2);
F23a_coef = cell2mat(coefficients(F23a));
F23b = bsxfun(@scalar_product, clif_v(1,:),e3);
F23b_coef = cell2mat(coefficients(F23b));
if isempty(F23a_coef) || isempty(F23b_coef)
    F23 = 0;
else 
    F23 = sum(F23a_coef.*F23b_coef);
end

F33a = bsxfun(@scalar_product, clif_u(1,:),e3);
F33a_coef = cell2mat(coefficients(F33a));
F33b = bsxfun(@scalar_product, clif_v(1,:),e3);
F33b_coef = cell2mat(coefficients(F33b));
if isempty(F33a_coef) || isempty(F33b_coef)
    F33 = 0;
else 
    F33 = sum(F33a_coef.*F33b_coef);
end

% Matrix F - by column
F = [F11, F12, F13; F21, F22, F23; F31, F32, F33];

end