function result = even_RL(int)
% Check if a number is even
%
% Input:
%   int - Integer or array of integers to test
%
% Output:
%   result - Logical: true if even, false if odd

    result = not(rem(int,2));
end