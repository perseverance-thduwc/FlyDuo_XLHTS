%% test_gabor.m - kiem thu nhanh Gabor
% (A) kiem thu don vi bang assert   (B) anh van vong tron + oracle
% (C) kiem tra do loi khi f thay doi (D) ve ngan hang 16 x 9 -> docs/figures
% Quy uoc theta: phuong vuong goc van (docs/CONVENTIONS.md muc 2).
% Yeu cau: Image Processing Toolbox (imfilter), MATLAB R2020a+ (exportgraphics).
clear; clc; close all; rng(1);

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
% Them thu muc chua ham vao path: neu co src/filters (cau truc repo) thi dung,
% neu khong thi dung chinh thu muc chua file nay (khi 5 file nam chung mot cho).
srcDir = fullfile(here, '..', 'src', 'filters');
if exist(srcDir, 'dir')
    addpath(srcDir);
    figDir = fullfile(here, '..', 'docs', 'figures');
else
    addpath(here);
    figDir = fullfile(here, 'figures');
end
if ~exist(figDir, 'dir'), mkdir(figDir); end

p    = defaultGaborParams();
bank = gaborBank(p);

%% (A) Kiem thu don vi ---------------------------------------------------
assert(isequal(size(bank.kernels), [16 9]),  'bank phai la 16 x 9');
assert(numel(bank.thetas) == 16 && numel(bank.freqs) == 9);
assert(all(bank.thetas >= 0 & bank.thetas < pi), 'theta phai nam trong [0, pi)');

for iF = 1:p.nFreq
    for iT = 1:p.nTheta
        g  = bank.kernels{iT, iF};
        gr = rot90(g, 2);
        assert(abs(sum(g(:))) < 1e-10,        'tong he so phai bang 0');
        assert(abs(norm(g(:)) - 1) < 1e-10,   'chuan L2 phai bang 1');
        assert(size(g,1) == size(g,2) && mod(size(g,1), 2) == 1, 'kernel vuong, canh le');
        assert(max(abs(g(:) - gr(:))) < 1e-12, 'kernel phai doi xung tam (pha khong)');
    end
end

% kich thuoc theo f (bang trong bao cao)
assert(size(gaborKernel(0, 0.2,  p), 1) == 15);
assert(size(gaborKernel(0, 0.1,  p), 1) == 27);
assert(size(gaborKernel(0, 1/14, p), 1) == 37);

% huong: theta = 0 -> van thang dung, song bien thien theo x
f = 0.1;  g0 = gaborKernel(0, f, p);  g90 = gaborKernel(pi/2, f, p);
h = (size(g0,1)-1)/2;  [X, ~] = meshgrid(-h:h, -h:h);
P = cos(2*pi*f*X);                                % van thang dung
assert(sum(g0(:).*P(:)) > 10 * abs(sum(g90(:).*P(:))), 'sai huong kernel (quy uoc goc)');

% do chon loc huong (in ra de dua vao bao cao)
fprintf('Do chon loc huong (f=%.2f):\n', f);
r0 = sum(g0(:).*P(:));
for a = [0 30 60 90]
    ga = gaborKernel(a*pi/180, f, p);
    fprintf('  lech %2d deg: %.3f\n', a, sum(ga(:).*P(:)) / r0);
end

% quay vong: theta ~ pi phai cung kernel voi theta = 0
Ir = rand(64);  M = true(64);  pb = struct('bank', bank);
o1 = filterGabor(Ir, M, zeros(64),          0.1*ones(64), pb);
o2 = filterGabor(Ir, M, (pi-1e-6)*ones(64), 0.1*ones(64), pb);
assert(isequal(o1, o2), 'theta ~ pi phai quay ve kernel 1');

% theta / f NaN khong duoc lam hong ham
thN = zeros(64);  thN(5:9, 5:9) = NaN;
fN  = 0.1*ones(64);  fN(20:25, 20:25) = NaN;
oN  = filterGabor(Ir, M, thN, fN, pb);
assert(all(isfinite(oN(:))), 'NaN trong theta/f lam hong dau ra');
assert(all(oN(:) >= 0 & oN(:) <= 1), 'dau ra phai trong [0,1]');
disp('(A) Kiem thu don vi: DAT');

%% (B) Anh van vong tron (oracle) ----------------------------------------
N = 256; f0 = 0.10;
[X, Y] = meshgrid(1:N, 1:N);
dx = X - N/2;  dy = Y - N/2;
r  = hypot(dx, dy);

I_clean = 0.5 + 0.5 * cos(2*pi*f0*r);
theta   = mod(atan2(dy, dx), pi);            % phap tuyen van = phuong ban kinh (dx, dy)
freq    = f0 * ones(N);
mask    = r < 110;

sigmaN  = 0.20;
I_noisy = min(max(I_clean + sigmaN * randn(N), 0), 1);

tic; I_enh = filterGabor(I_noisy, mask, theta, freq, struct('bank', bank)); t1 = toc;
I_bad = filterGabor(I_noisy, mask, theta, freq, struct('thetaOffset', pi/2));   % sai quy uoc

psnrM = @(a, b) 10*log10(1 / mean((a(mask) - b(mask)).^2));
psnrN = psnrM(I_noisy, I_clean);  psnrG = psnrM(I_enh, I_clean);  psnrB = psnrM(I_bad, I_clean);
fprintf('PSNR anh nhieu        : %.2f dB\n', psnrN);
fprintf('PSNR Gabor (dung goc) : %.2f dB   (%.2f s)\n', psnrG, t1);
fprintf('PSNR Gabor (lech pi/2): %.2f dB\n', psnrB);
assert(psnrG > psnrN + 8,  'Gabor phai cai thien >= 8 dB tren anh vong tron');
assert(psnrB < psnrN,      'lech quy uoc phai te hon anh nhieu');
disp('(B) Anh vong tron: DAT');

fig1 = figure('Position', [50 50 1200 330]);
subplot(1,4,1); imshow(I_clean);  title('Sach');
subplot(1,4,2); imshow(I_noisy);  title(sprintf('Nhieu AWGN %.2f  (%.1f dB)', sigmaN, psnrN));
subplot(1,4,3); imshow(I_enh);    title(sprintf('Gabor, dung quy uoc (%.1f dB)', psnrG));
subplot(1,4,4); imshow(I_bad);    title(sprintf('Gabor, lech pi/2 (%.1f dB)', psnrB));
set(findall(fig1, 'Type', 'axes'), 'Toolbar', []);   % tat thanh cong cu truoc khi xuat
drawnow;
exportgraphics(fig1, fullfile(figDir, 'gabor_ring_test.png'), 'Resolution', 200);

%% (C) Do loi khong phu thuoc f (gainComp) --------------------------------
% Hai nua anh: van ngang f=1/14 va f=0.2, cung bien do 0.5 quanh 0.5.
M2 = 128;  [~, Yb] = meshgrid(1:M2, 1:M2);
Ia = 0.5 + 0.5*cos(2*pi*(1/14)*Yb);   Ib = 0.5 + 0.5*cos(2*pi*0.2*Yb);
Iab = [Ia, Ib];  Mab = true(M2, 2*M2);
fab = [ (1/14)*ones(M2), 0.2*ones(M2) ];  thab = (pi/2) * ones(M2, 2*M2);
rawOn  = filterGabor(Iab, Mab, thab, fab, struct('outputMode','raw', 'gainComp', true));
rawOff = filterGabor(Iab, Mab, thab, fab, struct('outputMode','raw', 'gainComp', false));
c = 30:M2-30;                                   % bo vung sat bien
ampOn  = [max(abs(rawOn (c, c)),[],'all'),  max(abs(rawOn (c, M2+c)),[],'all')];
ampOff = [max(abs(rawOff(c, c)),[],'all'),  max(abs(rawOff(c, M2+c)),[],'all')];
fprintf('Bien do dap ung (f=1/14 | f=0.2): gainComp tat = %.2f | %.2f ; bat = %.2f | %.2f\n', ampOff, ampOn);
assert(abs(ampOn(1)/ampOn(2) - 1) < 0.1,  'gainComp bat: bien do phai xap xi nhau');
assert(ampOff(1)/ampOff(2) > 1.5,         'gainComp tat: bien do lech theo f (xem bao cao)');
disp('(C) Do loi theo f: DAT');

%% (D) Ve ngan hang 16 huong x 9 tan so ------------------------------------
fig2 = figure('Position', [20 20 1700 900], 'Color', 'w');
tl = tiledlayout(p.nFreq, p.nTheta, 'TileSpacing', 'none', 'Padding', 'compact');
for iF = 1:p.nFreq
    for iT = 1:p.nTheta
        nexttile; g = bank.kernels{iT, iF};
        imagesc(g, [-1 1] * max(abs(g(:)))); axis image; colormap(gray);
        set(gca, 'XTick', [], 'YTick', []);
        if iF == 1, title(sprintf('%.1f\\circ', bank.thetas(iT)*180/pi), 'FontSize', 7); end
        if iT == 1, ylabel(sprintf('%.3f', bank.freqs(iF)), 'FontSize', 7); end
    end
end
title(tl, 'Ngan hang Gabor 16 huong x 9 tan so (cot: huong, hang: f chu ky/pixel)');
set(findall(fig2, 'Type', 'axes'), 'Toolbar', []);
drawnow;
exportgraphics(fig2, fullfile(figDir, 'gabor_bank_16x9.png'), 'Resolution', 200);
fprintf('Da luu hinh vao %s\n', figDir);

fprintf('PASS test_gabor\n');
