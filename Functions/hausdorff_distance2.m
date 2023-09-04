% Gives the SAME RESULTS AS hausdorff_distance.m
function distance = hausdorff_distance2(A, B)

D = pdist2(A,B);
dist_ab = max(min(D,[],2)); % Directed from point a to point b
dist_ba = max(min(D)); % Directed from point b to point a
distance = max([dist_ab,dist_ba]);
end