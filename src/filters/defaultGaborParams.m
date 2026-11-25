function p = defaultGaborParams()
% DEFAULTGABORPARAMS  Tham so mac dinh cho nhanh Gabor.
%   p = defaultGaborParams()
%   Moi truong co the bi ghi de bang struct 'params' (cfg.filters.gabor)
%   truyen vao filterGabor; truong thieu se lay gia tri o day.
%
% QUY UOC GOC (theo docs/CONVENTIONS.md muc 2):
%   theta = goc cua phuong VUONG GOC voi duong van (huong gradient), do tu
%   truc +x quay ve +y, radian, lay mod pi.
%   He toa do anh: x = cot (sang phai), y = hang (xuong duoi).
%   Toa do xoay: xp = x*cos(theta) + y*sin(theta)  (song cos nam o day).
%   theta = 0: van thang dung ; theta = pi/2: van nam ngang.
%   Neu mot module khac dung quy uoc lech di mot goc a thi dat thetaOffset = a.

p.nTheta      = 16;          % so huong, cach deu trong [0, pi)
p.nFreq       = 9;           % so muc tan so, cach deu
p.freqRange   = [1/15, 1/5]; % [fmin fmax] chu ky/pixel (= cfg.freq.range)
p.sigmaK      = 5/12;        % sigma = sigmaK / f  (= cfg.filters.gabor.sigmaK)
p.thetaOffset = 0;           % cong them vao theta khi dung kernel
p.gainComp    = true;        % chia dap ung cho do loi cua kernel -> bien do khong phu thuoc f
p.outputMode  = 'rescale';   % 'rescale' -> [0,1] ; 'raw' -> dap ung tho cua bo loc
p.clipPct     = [1 99];      % phan vi cat khi rescale (phan tram)
p.bgValue     = 1;           % gia tri ngoai mat na (1 = nen trang, van toi)
p.padding     = 'symmetric'; % kieu dem bien cho imfilter
end
