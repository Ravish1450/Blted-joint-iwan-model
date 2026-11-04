clear; clc; close all;

%% Fig. 7: natural frequency and damping vs modal displacement amplitude
%  (QSMA, NP-Iwan, 4-parameter Iwan and full FE decay)
%  INPUTS : qr_curve_2614.txt, iwan_decays_ds.csv, fe_modal_decay.csv
%  OUTPUT : fig7_freevib_matlab.png
%  Frequency and damping from each decay: Hilbert transform + straight-line fits

%% SETTINGS
QR_FILE     = 'qr_curve_2614.txt';
MACRO_IDX   = 533;
IWAN_FILE   = 'iwan_decays_ds.csv';
FE_FILE     = 'fe_modal_decay.csv';
FE_DT       = 1e-5;
DLN         = 0.1;        % FREEVIB settings
MIN_CYCLES  = 4;
TRIM_CYCLES = 3;
A_MIN_IWAN  = 2.5e-5;
Q_MIN_PLOT  = 2e-6;
OUT         = 'fig7';

%% QSMA curve
[qq, fq, zq, f0] = qsma_curves(QR_FILE, MACRO_IDX);
fprintf('QSMA f0 = %.3f Hz\n', f0);

%% Decay signals
iw  = dlmread(IWAN_FILE, ',');
sig = {};
if exist(FE_FILE, 'file')
    fe = dlmread(FE_FILE, ',', 1, 0);
    [tf, iu] = unique(fe(:,1));
    tu = (tf(1):FE_DT:tf(end)).';
    sig(end+1,:) = {'Dynamic Full FE', tu, interp1(tf, fe(iu,2), tu), ':o', [0 0.7 0]};
end
sig(end+1,:) = {'NP-Iwan',          iw(:,1), iw(:,2), '--', [0 0 1]};
sig(end+1,:) = {'4-parameter Iwan', iw(:,1), iw(:,3), '-.', [1 0 0]};

ns  = size(sig,1);
res = cell(ns,3);                                  
for i = 1:ns
    [A, f, z] = hilbert_freevib(sig{i,2}, sig{i,3}, f0, DLN, MIN_CYCLES, TRIM_CYCLES);
    if ~strcmp(sig{i,1}, 'Dynamic Full FE')
        keep = A >= A_MIN_IWAN; A = A(keep); f = f(keep); z = z(keep);
    end
    res(i,:) = {A, f, z};
end

%% Plot: (a) frequency, (b) damping
m    = qq > Q_MIN_PLOT;
qsma = {fq, zq};
ylab = {'Natural Frequency (Hz)', 'Damping'};
lloc = {'southwest', 'northwest'};
figure('Position', [100 100 1150 450]);
for p = 1:2
    subplot(1,2,p);
    semilogx(qq(m), qsma{p}(m), 'k-', 'LineWidth', 2, 'DisplayName', 'QSMA'); hold on
    for i = 1:ns
        semilogx(res{i,1}, res{i,p+1}, sig{i,4}, 'Color', sig{i,5}, 'LineWidth', 1.5, ...
                 'MarkerSize', 5, 'DisplayName', sig{i,1});
    end
    xlabel('Modal Displacement Amplitude [kg^{1/2}m]', 'FontWeight', 'bold');
    ylabel(ylab{p}, 'FontWeight', 'bold');
    xlim([1e-6 2.5e-3]); grid on; box on; legend('Location', lloc{p});
    set(gca, 'FontWeight', 'bold', 'FontSize', 11);
    title(sprintf('(%c)', 'a'+p-1), 'FontWeight', 'normal');
end
print(gcf, [OUT '.png'], '-dpng', '-r200');

%% Values at chosen amplitudes
fprintf('\nValues at chosen amplitudes (f in Hz, zeta)\n');
for A0 = [5e-5 1e-4 2e-4 5e-4 1e-3 1.5e-3]
    fprintf('A = %.1e  QSMA: f %.2f z %.3e', A0, value_at(A0,qq,fq), value_at(A0,qq,zq));
    for i = 1:ns
        fprintf(' | %s: f %.2f z %.3e', sig{i,1}(1:4), value_at(A0,res{i,1},res{i,2}), ...
                value_at(A0,res{i,1},res{i,3}));
    end
    fprintf('\n');
end

%% FUNCTIONS
function z = hilb(x)
x = x(:); n = numel(x); X = fft(x); h = zeros(n,1);
if mod(n,2) == 0
    h([1 n/2+1]) = 1; h(2:n/2) = 2;
else
    h(1) = 1; h(2:(n+1)/2) = 2;
end
z = ifft(X.*h);
end

function [A, f, zeta] = hilbert_freevib(t, x, f_guess, dln, min_cycles, trim_cycles)
t = t(:); x = x(:); dt = t(2) - t(1);
z   = hilb(x);
lnA = log(abs(z));
th  = unwrap(angle(z));
Ns  = round(1/f_guess/dt);
ntrim = trim_cycles*Ns;
idx = (ntrim+1):(numel(t)-ntrim);
tt = t(idx); lnA = lnA(idx); th = th(idx); n = numel(tt);
lnm = cummin(lnA);
levels = lnm(1):-dln/2:lnm(end);
A = []; f = []; zeta = [];
for k = 1:numel(levels)-2
    s = find(lnm <= levels(k), 1);
    e = find(lnm <= levels(k+2), 1);
    if isempty(s) || isempty(e), break; end
    c = round((s+e)/2);
    ncyc = max(min_cycles, round((e-s)/Ns));
    h = floor(ncyc*Ns/2);
    if c-h < 1 || c+h > n, continue; end
    sl = (c-h):(c+h-1);
    tc = mean(tt(sl));
    pa = polyfit(tt(sl)-tc, lnA(sl), 1);
    pw = polyfit(tt(sl)-tc, th(sl), 1);
    sg = -pa(1); wd = pw(1); wn = sqrt(wd^2 + sg^2);
    A(end+1,1) = exp(pa(2)); f(end+1,1) = wn/(2*pi); zeta(end+1,1) = sg/wn;
end
end

function [qm, f, z, f0] = qsma_curves(qr_file, macro_idx)
d  = dlmread(qr_file, ',', 1, 0);
a  = d(:,3); q = abs(d(:,4));
Kt = diff(a)./diff(q);
tail = Kt(macro_idx:end);
Kest = median(Kt(end-4:end));
Kinf = mean(tail(abs(tail-Kest)/Kest < 0.01));
qm = q(1:macro_idx); am = a(1:macro_idx);
fnl = am - Kinf*qm;
D = 8*cumtrapz(qm, fnl) - 4*qm.*fnl;
z = D./(2*pi*qm.*am);
f = sqrt(am./qm)/(2*pi);
ok = qm > 0;
qm = qm(ok); f = f(ok); z = z(ok);
f0 = sqrt(Kt(1))/(2*pi);
end

function y0 = value_at(A0, A, y)
if isempty(A) || A0 < min(A) || A0 > max(A)
    y0 = NaN; return
end
[As, o] = sort(A(:)); ys = y(o);
[As, iu] = unique(As); ys = ys(iu);
y0 = interp1(As, ys, A0);
end