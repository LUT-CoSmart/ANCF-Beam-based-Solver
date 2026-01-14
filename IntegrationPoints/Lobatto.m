function [x,w] = Lobatto(n)

% http://www.mathworks.com/matlabcentral/newsreader/view_thread/11507
% MKM: It is found from matlabcentral and it is not
% verified carefully...!!!!!!!!! Do it later!


% Weights w and nodes x of the n-point Lobatto quadrature formula.

if n < 2
   error('Number of nodes must be greater than one.')
end
s = 2/(n*(n-1));
lP = LegendP(n-1);
coeff = derp(lP);
x = sort(roots(coeff));
rt = polyval(lP,x);
w = 1./rt.^2;
w = s*[1 w(:)' 1];
x = [-1 x(:)' 1];

w=w(:);
x=x(:);


function P = LegendP(n)

% Coefficients P of the nth Legengre polynomial.
% They are stored in the decreasing order of powers.

p0 = 1;
p1 = [1 0];
if n == 0
   P = p0;
elseif n == 1
   P = p1;
else
   for k=2:n
      P = [(2*k-1)/k*p1 0] - [0 0 (k-1)/k*p0];
      p0 = p1;
      p1 = P;
   end
end


function dp = derp(p)

% Derivative dp of an algebraic polynomial that is
% represented by its coefficients p. They must be stored
% in the descending order of powers.

n = length(p) - 1;
p = p(:)';
dp = p(1:n).*(n:-1:1);
k = find(dp ~= 0);
if ~isempty(k)
   dp = dp(k(1):end);
else
   dp = 0;
end