% A small script to test out the efficiency of using circshift for finite
% differences compared to matrix multiplication. 

function testing()
%{
    k = 100;
    results = zeros(3, k);
    for i = 1:k
        n = 10 * i;
        results(1, i) = n;
        f1 = @() useMultiply(n);
        results(2, i) = timeit(f1);
        f2 = @() useCircshift(n);
        results(3, i) = timeit(f2);
    end
    hold on
    plot( 1:k, results(2, :), 'r' )
    plot( 1:k, results(3, :), 'b' )
%}
end

function useMultiply(n)
    dx = circshift(eye(n), 1) - eye(n);
    for i = 1:100
        f = rand(n);
        fdiff = dx * f;
    end
end

function useCircshift(n)
    for i = 1:100
        f = rand(n);
        fdiff = circshift(f, 1) - f;
    end
end

function Vq = interpOnS1andTime(T, X, V, Tq, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    %disp(X)
    Vq = interp2( T, newX, newV, Tq, Xq );
end