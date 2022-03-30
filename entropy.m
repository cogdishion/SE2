%called by discrete_nmi.m

function h = entropy(varargin)
    if length(varargin) == 1
        h = entropySolo(varargin{1});
        return
    end

    c1 = varargin{1};
    type = varargin{2};
    c2 = varargin{3};

    switch type
      case 'given'
        h = entropyGiven(c1, c2);
      case 'and'
        h = entropyJoint(c1, c2);
    end
end

function h = entropySolo(x)
    p = probabilityDistribution(x);
    h = entropyBase(p);
end

function p = probabilityDistribution(x)
    nk = sum(x, 1);
    n = size(x, 1);
    p = nk / n;
    p = [p; (1 - p)];
end

function h = entropyBase(p)
    h = p .* log2(p);
    h(isnan(h)) = 0;
    h = - sum(h, 1);
end

function h = entropyGiven(x, y)
    h = zeros(1, size(x, 2));
    hy = entropy(y);
    for k = 1:size(x, 2)
        hxk = abs(entropy(x(:, k), 'and', y) - hy);
        h(k) = min(hxk);
    end
    h = mean(h ./ entropy(x));
end

function h = entropyJoint(x, y)
    h = entropyBase(jointDistribution(x, y));
end

function p = jointDistribution(x, y)
    k = size(y, 2);

    xy = x & y;
    nkX = sum(x, 1);
    nkY = sum(y, 1);
    nkXY = sum(xy, 1);
    n = size(x, 1);

    p = zeros(4, k);
    p(1, :) = nkXY;
    p(2, :) = nkX - nkXY;
    p(3, :) = nkY - nkXY;
    p(4, :) = n - sum(x | y, 1);
    p = p / n;
end
