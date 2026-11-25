function g = gaborKernel(theta, f, p)
% GABORKERNEL  Kernel Gabor 2D cho van tay.
%   g = gaborKernel(theta, f, p)
%   theta : goc phuong vuong goc van (rad), theo docs/CONVENTIONS.md muc 2.
%   f     : tan so (chu ky/pixel), so duong.
%   p     : struct (cfg.filters.gabor); can sigmaK, thetaOffset. Thieu -> mac dinh.
%   g     : kernel vuong, canh le, double.
%
%   xp = x*cos(th) + y*sin(th)    vuong goc van : song cos nam o day
%   yp = -x*sin(th) + y*cos(th)   doc theo van  : chi con bao Gauss
%   g  = exp(-(x^2+y^2)/(2 sigma^2)) .* cos(2 pi f xp)
%   sigma = sigmaK / f ;  kich thuoc s = 1 + 2*ceil(3 sigma)
%   Chuan hoa: trung binh 0, chuan L2 bang 1.

if nargin < 2
    error('gaborKernel:args', 'Cach dung: g = gaborKernel(theta, f [, p]).  Vi du: g = gaborKernel(0, 0.1);');
end
if nargin < 3 || isempty(p), p = struct(); end
p = fillDefaults(p);
assert(isscalar(f) && f > 0, 'gaborKernel: f phai la so duong (chu ky/pixel).');

sigma = p.sigmaK / f;
half  = ceil(3 * sigma);                 % s = 1 + 2*half
[x, y] = meshgrid(-half:half, -half:half);

th = theta + p.thetaOffset;
xp = x * cos(th) + y * sin(th);

g = exp(-(x.^2 + y.^2) / (2 * sigma^2)) .* cos(2 * pi * f * xp);
g = g - mean(g(:));
g = g / norm(g(:));
end

function p = fillDefaults(params)
p = defaultGaborParams();
f = fieldnames(params);
for i = 1:numel(f)
    p.(f{i}) = params.(f{i});
end
end
