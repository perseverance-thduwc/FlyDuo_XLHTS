function I_enh = filterGabor(I_noisy, mask, theta, freq, params)
% FILTERGABOR  Loc ngu canh Gabor (ngan hang nTheta x nFreq).
%   I_enh = filterGabor(I_noisy, mask, theta, freq, params)
%
%   I_noisy : anh double [0,1]            mask : logical, cung kich thuoc
%   theta   : ban do huong (rad, [0,pi)), phuong VUONG GOC van (CONVENTIONS muc 2)
%   freq    : ban do tan so (chu ky/pixel)
%   params  : cfg.filters.gabor (co the thieu truong -> lay mac dinh); co the dua
%             san params.bank (tu gaborBank) de khong phai dung lai moi lan goi.
%   I_enh   : double, cung kich thuoc I_noisy, khong NaN.
%
%   Moi pixel trong mask dung kernel gan nhat voi (theta, freq) cua no.
%   Cai dat: loc ca anh voi tung kernel dang duoc dung, chi giu pixel tuong ung.
%   theta/freq khong hop le (NaN/Inf) duoc thay bang gia tri dai dien trong mask.

if nargin < 4
    error('filterGabor:args', ['Cach dung: I_enh = filterGabor(I_noisy, mask, theta, freq [, params]).' newline ...
          'Vi du: I_enh = filterGabor(rand(64), true(64), zeros(64), 0.1*ones(64));']);
end
if nargin < 5 || isempty(params), params = struct(); end
p = fillDefaults(params);

I    = double(I_noisy);
mask = logical(mask);
theta = double(theta);
if isfield(params, 'bank'), bank = params.bank; else, bank = gaborBank(p); end
if ~isfield(bank, 'gain'), bank.gain = ones(size(bank.kernels)); end
nT = p.nTheta;  nF = p.nFreq;

% ---- theta khong hop le -> trung binh vong (tren 2*theta) trong mask
tv = isfinite(theta) & mask;
if any(tv(:))
    th0 = 0.5 * atan2(mean(sin(2*theta(tv))), mean(cos(2*theta(tv))));
else
    th0 = 0;
end
theta(~isfinite(theta)) = th0;

% ---- luong tu hoa (theta, f) -> chi so kernel
iT = mod(round(mod(theta, pi) / (pi / nT)), nT) + 1;      % 1..nT, quay vong

fq = double(freq);
valid = isfinite(fq) & mask;
if any(valid(:)), fq(~isfinite(fq)) = median(fq(valid)); else, fq(~isfinite(fq)) = mean(p.freqRange); end
fq = min(max(fq, bank.freqs(1)), bank.freqs(end));
if nF == 1
    iF = ones(size(fq));
else
    iF = round((fq - bank.freqs(1)) / (bank.freqs(end) - bank.freqs(1)) * (nF - 1)) + 1;
end
idx = iT + (iF - 1) * nT;                                  % chi so tuyen tinh vao bank.kernels

% ---- loc
out = zeros(size(I));
for k = unique(idx(mask))'
    R = imfilter(I, bank.kernels{k}, p.padding, 'corr');
    if p.gainComp, R = R / bank.gain(k); end
    sel = mask & (idx == k);
    out(sel) = R(sel);
end

% ---- dau ra
switch lower(p.outputMode)
    case 'raw'
        I_enh = out;
        I_enh(~mask) = 0;
    otherwise   % 'rescale'
        v = sort(out(mask));
        if isempty(v)
            I_enh = p.bgValue * ones(size(I));
            return
        end
        n  = numel(v);
        lo = v(max(1, round(p.clipPct(1) / 100 * n)));
        hi = v(min(n, max(1, round(p.clipPct(2) / 100 * n))));
        if hi <= lo, hi = lo + eps; end
        I_enh = (out - lo) / (hi - lo);
        I_enh = min(max(I_enh, 0), 1);
        I_enh(~mask) = p.bgValue;
end
end

function p = fillDefaults(params)
p = defaultGaborParams();
f = fieldnames(params);
for i = 1:numel(f)
    if ~strcmp(f{i}, 'bank'), p.(f{i}) = params.(f{i}); end
end
end
