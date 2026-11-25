function bank = gaborBank(p)
% GABORBANK  Dung ngan hang nTheta x nFreq kernel Gabor.
%   bank = gaborBank(params)      params = cfg.filters.gabor (co the thieu truong)
%   bank.kernels{iT,iF} : kernel
%   bank.thetas : 1 x nTheta (rad, [0,pi))   bank.freqs : 1 x nFreq (chu ky/pixel)
%   bank.gain(iT,iF)    : do loi kernel voi song sin cung (theta,f), bien do 1.
%                         (kernel chuan L2 nen do loi ~ sigma ~ 1/f; filterGabor chia
%                          cho gain de dap ung khong phu thuoc tan so.)

if nargin < 1 || isempty(p), p = struct(); end
d = defaultGaborParams();
fn = fieldnames(p);
for i = 1:numel(fn), d.(fn{i}) = p.(fn{i}); end
p = d;

nT = p.nTheta;  nF = p.nFreq;
bank.thetas = (0:nT-1) * pi / nT;                       % [0, pi)
if nF == 1
    bank.freqs = mean(p.freqRange);
else
    bank.freqs = linspace(p.freqRange(1), p.freqRange(2), nF);
end

bank.kernels = cell(nT, nF);
bank.gain    = zeros(nT, nF);
for iF = 1:nF
    f = bank.freqs(iF);
    for iT = 1:nT
        g = gaborKernel(bank.thetas(iT), f, p);
        bank.kernels{iT, iF} = g;

        h = (size(g, 1) - 1) / 2;
        [x, y] = meshgrid(-h:h, -h:h);
        th = bank.thetas(iT) + p.thetaOffset;
        xp = x * cos(th) + y * sin(th);
        bank.gain(iT, iF) = sum(g(:) .* cos(2 * pi * f * xp(:)));
    end
end
end
