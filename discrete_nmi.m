function I = discrete_nmi(x, y)
%NMI calculate normalized mutual information between two partitions
%   I = NMI(x, y) computes the normalized mutual information of
%   partitions x and y. Partitions should be a vector of community ids
%   and both partitions should have the same number of nodes.

    [x, y] = membership2cover(x, y);
    hx = entropy(x);
    hy = entropy(y);
    I = (hx + hy - entropy(x, y));
    I = I / (sqrt(hx * hy));

    if isnan(I)
        % If x and y have no communities in common, I will be nan.
        I = 0;
    end
end

function [c1, c2] = membership2cover(x, y)
% Convert a list of community ids to a cover. Sets the first community
% to 1.
    assert(min(size(x)) == 1, ...
           "x most be a vector of community ids. Received a matrix instead.")
    assert(min(size(y)) == 1, ...
           "y most be a vector of community ids. Received a matrix instead.")
    assert(length(x) == length(y), ...
           "Partitions x and y most have the same number of nodes." + ...
           "Recieved %d nodes for x and %d nodes for y.", ...
           length(x), length(y))

    if (size(x, 2) ~= 1)
        x = x';
    end

    if (size(y, 2) ~= 1)
        y = y';
    end

    communityNumberingOffset = min(min(x), min(y)) - 1;
    x = x - communityNumberingOffset;
    y = y - communityNumberingOffset;
    nCommunities = max(max(x), max(y));

    c1 = sparse(1:length(x), x, ones(length(x), 1), length(x), nCommunities);
    c2 = sparse(1:length(y), y, ones(length(y), 1), length(y), nCommunities);
end

function h = entropy(varargin)
    if length(varargin) == 1
        h = entropySolo(varargin{1});
        return
    else
        x = varargin{1};
        y = varargin{2};
        h = entropyJoint(x, y);
    end
end

function h = entropySolo(x)
    p = probabilityDistribution(x);
    h = entropyBase(p);
end

function p = probabilityDistribution(x)
    p = mean(x, 1);
end

function h = entropyBase(p)
    p = nonzeros(p);
    h = -p' * log2(p);
end

function h = entropyJoint(x, y)
    h = entropyBase(jointDistribution(x, y));
end

function p = jointDistribution(x, y)
    p = (x' * y) / size(x, 1);
end