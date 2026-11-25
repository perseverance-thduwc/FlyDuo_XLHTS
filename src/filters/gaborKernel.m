function g = gaborKernel(theta, f, p)
% GABORKERNEL  Kernel Gabor 2D cho van tay.
%   g = gaborKernel(theta, f, p)
%   theta : huong doc van (rad). f : tan so (chu ky/pixel).
%   p     : struct tu defaultGaborParams (can sigmaScale, thetaOffset).
%
%   g(x,y) = exp(-(x^2+y^2)/(2 sigma^2)) * cos(2 pi f * x_n)
%   x_n    : toa do theo phap tuyen cua van (song mang chay ngang van).
%   sigma  = sigmaScale * 5/(12 f);  kich thuoc s = 1 + 2*ceil(3 sigma)
%   Chuan hoa: trung binh 0, chuan L2 bang 1.

if nargin < 3 || isempty(p), p = defaultGaborParams(); end
assert(isscalar(f) && f > 0, 'gaborKernel: f phai la so duong (chu ky/pixel).');

sigma = p.sigmaScale * 5 / (12 * f);
half  = ceil(3 * sigma);                 % s = 1 + 2*half
[x, y] = meshgrid(-half:half, -half:half);

th = theta + p.thetaOffset;
% u = (cos th, sin th) doc van ; n = (-sin th, cos th) vuong goc van
xn = -x * sin(th) + y * cos(th);

g = exp(-(x.^2 + y.^2) / (2 * sigma^2)) .* cos(2 * pi * f * xn);
g = g - mean(g(:));
g = g / norm(g(:));
end

