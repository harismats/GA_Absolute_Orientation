function result = mean_absolute_error(X,Y)
result = sum(sum(abs(X-Y))) / length(X(:));
end

%% or 2nd WAY
% Calculate First the absolute error
% score = abs(X(:)-Y(:)); % the formula is score = abs(actual(:)-prediction(:))
% AND then
% result = mean(score);

%% or 3rd WAY using MATLAB's in-built function
% result = mae(X-Y); % the formula is mae(Actual - Predicted)