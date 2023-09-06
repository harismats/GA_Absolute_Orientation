function G = matrix_G_optimized(v)
%% Checking size of data
[rows_v, ~] = size(v);

if rows_v ~=3
    v=v';
end
%% Covnvert to double
v = double(v);

%% Subtract the mean from point set v as this is what GA method does
v = v - mean(v,2);

%% Vector containing e1,e2,e3
vec = [e1,e2,e3];

%% Multiplication of v (COLUMN) - Convert it to Clifford
clif_v = vec * v;

%% Forming Matrix G --> G_{jk} = sum_{i=1}^{i=N} (v_i . e_j)(v_i . e_k)
%% 1st Column
G11a = bsxfun(@scalar_product, clif_v(1,:),e1);
G11b = bsxfun(@scalar_product, clif_v(1,:),e1);
G11 = sum(G11a.*G11b);

G21a = bsxfun(@scalar_product, clif_v(1,:),e2);
G21b = bsxfun(@scalar_product, clif_v(1,:),e1);
G21 = sum(G21a.*G21b);

G31a = bsxfun(@scalar_product, clif_v(1,:),e3);
G31b = bsxfun(@scalar_product, clif_v(1,:),e1);
G31 = sum(G31a.*G31b);

%% 2nd Column
G12a = bsxfun(@scalar_product, clif_v(1,:),e1);
G12b = bsxfun(@scalar_product, clif_v(1,:),e2);
G12 = sum(G12a.*G12b);

G22a = bsxfun(@scalar_product, clif_v(1,:),e2);
G22b = bsxfun(@scalar_product, clif_v(1,:),e2);
G22 = sum(G22a.*G22b);

G32a = bsxfun(@scalar_product, clif_v(1,:),e3);
G32b = bsxfun(@scalar_product, clif_v(1,:),e2);
G32 = sum(G32a.*G32b);

%% 3rd Column

G13a = bsxfun(@scalar_product, clif_v(1,:),e1);
G13b = bsxfun(@scalar_product, clif_v(1,:),e3);
G13 = sum(G13a.*G13b);

G23a = bsxfun(@scalar_product, clif_v(1,:),e2);
G23b = bsxfun(@scalar_product, clif_v(1,:),e3);
G23 = sum(G23a.*G23b);

G33a = bsxfun(@scalar_product, clif_v(1,:),e3);
G33b = bsxfun(@scalar_product, clif_v(1,:),e3);
G33 = sum(G33a.*G33b);

% Matrix G - by column
G = coefficients([G11, G12, G13; G21, G22, G23; G31, G32, G33]);
G = cell2mat(G);

end