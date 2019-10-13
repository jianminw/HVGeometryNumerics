

function Vq = interpOnS1andTime(T, X, V, Tq, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    %disp(X)
    %disp(size(newX))
    %disp(size(newV))
    Vq = interp2( T, newX, newV, Tq, Xq, 'spline');
end