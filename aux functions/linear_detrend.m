function [X_yfit] = linear_detrend(X, t)
% applies linear detrend on the first dimension of X
% X_yfit: detrended signal

if size(X,1) ~= size(t,1)
    disp(['Warning: X first dimension and t should have the same size.'])
    return
end


if ndims(X) == 2;
    for j = 1:size(X,2) % channel
        p = polyfit(t,X(:,j),1);
        yfit = polyval(p,t);
        X_yfit(:,j) = X(:,j)-yfit;
    end
elseif ndims(X) == 3;
    for i = 1:size(X,2)
        for j = 1:size(X,3) % channel
            p = polyfit(t,X(:,i,j),1);
            yfit = polyval(p,t);
            X_yfit(:,i,j) = X(:,i,j)-yfit;
        end
    end
end

