
function Action = ComputeAction(Path)

% compute vx, and take L2 norms. 

% Use balanced derivatives. 
% From some testing, circshift seems to have better performance than 
% matrix multiplication. 

vdiff = circshift(Path.v, -1) - circshift(Path.v, 1);

vx = vdiff * (size(Path.v, 1) - 1) / 2;

Action = L2(Path.v) + L2(Path.z) + L2(vx);
%disp(L2(Path.v))
%disp(L2(Path.z))
%disp(L2(vx))

end

function y = L2(f)
y = mean( f .* f, 'all');
end