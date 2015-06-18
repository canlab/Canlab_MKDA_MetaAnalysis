function d = dist_to_surface(coords, vertices)
% Euclidean distance of each of a set of [x y z] coords from a surface defined by a list of [x y z] vertices
%
% e.g.,
% coords = [-22 0 -12];
% vertices = vertices = get(p(1), 'Vertices');
%
% coords = [-22 0 -12; 30 30 30; 42 12 3];
% d = dist_to_surface(coords, vertices)
%
% Tor Wager, 3/2013

for i = 1:size(coords, 1)
    
    m = bsxfun(@minus, coords(i, :), vertices);
    
    m = sum(m .^ 2, 2);                            % sum of squared distances
    
    [d(i, 1), wh_vertex] = min(m);                 % min squared distance
    
end

d = d .^ .5;  % min distance for each point

end % function

