% Fig. 8: natural frequency and damping vs modal displacement amplitude of the
% non-parametric Iwan model for different slider-spacing ratios alpha,
% with QSMA and the full FE decay.
% INPUTS : qr_curve_2614.txt, iwan_alpha_sweep.csv, fe_modal_decay.csv
% OUTPUTS: fig8.png, fig8_data.csv

clear; clc; close all;

%% SETTINGS
QR_FILE     = 'qr_curve_2614.txt';
MACRO_IDX   = 533;
SWEEP_FILE  = 'iwan_alpha_sweep.csv';
FE_FILE     = 'fe_modal_decay.csv';
FE_DT       = 1e-5;
ALPHAS      = [0.5 0.75 1 1.25 1.5];
DLN         = 0.10;
MIN_CYCLES  = 4;
TRIM_CYCLES = 3;
A_MIN_IWAN  = 2.5e-5;
Q_MIN_PLOT  = 2e-6;
OUT         = 'fig8_freevib';

%% QSMA curve and NP-Iwan decays for each alpha
[qq, fq, zq, f0] = qsma_curves(QR_FILE, MACRO_IDX);
sw = dlmread(SWEEP_FILE, ',');
t  = sw(:,1);
if size(sw,2)-1 ~= numel(ALPHAS)
    error('%s has %d decay columns, ALPHAS has %d', SWEEP_FILE, size(sw,2)-1, numel(ALPHAS));
end

na  = numel(ALPHAS);
res = cell(na,3);                                   
for i = 1:na
    [A, f, z] = hilbert_freevib(t, sw(:,i+1), f0, DLN, MIN_CYCLES, TRIM_CYCLES);
    keep = A >= A_MIN_IWAN;
    res(i,:) = {A(keep), f(keep), z(keep)};
end

%% Full FE decay
haveFE = exist(FE_FILE, 'file');
if haveFE
    fe = dlmread(FE_FILE, ',', 1, 0);
    [tf, iu] = unique(fe(:,1));
    tu = (tf(1):FE_DT:tf(end)).';
    [Afe, ffe, zfe] = hilbert_freevib(tu, interp1(tf, fe(iu,2), tu), f0, DLN, MIN_CYCLES, TRIM_CYCLES);
    fev = {ffe, zfe};
end

%% FIGURE: (a) frequency, (b) damping
styles = {'--','-.',':','--','-.'};
cols   = [0.49 0.18 0.56; 0 0.75 0.75; 0 0 1; 0.93 0.69 0.13; 1 0 0];
m      = qq > Q_MIN_PLOT;
qsma   = {fq, zq};
xpos   = [0.10 0.585];
ylab   = {'Natural Frequency (Hz)', 'Damping'};
lloc   = {'southwest', 'northwest'};

f8 = figure('Units', 'centimeters', 'Position', [2 2 17 7.5], 'Color', 'w');
set(f8, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 17 7.5], 'PaperSize', [17 7.5]);

for p = 1:2
    ax = subplot(1,2,p);
    ax.Position = [xpos(p) 0.25 0.38 0.67];
    semilogx(qq(m), qsma{p}(m), 'k-', 'LineWidth', 3, 'DisplayName', 'QSMA'); hold on
    if haveFE
        semilogx(Afe, fev{p}, ':', 'Color', [0 0.8 0], 'LineWidth', 4, 'DisplayName', 'Dynamic Full FE');
    end
    for i = 1:na
        semilogx(res{i,1}, res{i,p+1}, styles{i}, 'Color', cols(i,:), 'LineWidth', 1.8, ...
                 'DisplayName', sprintf('NP-Iwan, \\alpha = %g', ALPHAS(i)));
    end
    xlabel('Modal Displacement Amplitude [kg^{1/2}m]', 'FontWeight', 'bold');
    ylabel(ylab{p}, 'FontWeight', 'bold');
    lg = legend('Location', lloc{p}, 'FontSize', 6.5);
    lg.ItemTokenSize = [14 8];
    text(ax, 0.5, -0.27, sprintf('(%c)', 'a'+p-1), 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12);
    xlim([1e-6 2.5e-3]); grid on; box on;
    set(ax, 'FontWeight', 'bold', 'FontSize', 8);
    ax.XLabel.FontSize = 8;  ax.YLabel.FontSize = 8;
    if p == 2, ylim([0 5.2e-3]); end
end

%% SAVE
print(f8, [OUT '.png'], '-dpng', '-r600');

fid = fopen([OUT '_data.csv'], 'w');
fprintf(fid, 'alpha,amplitude,frequency_Hz,zeta\n');
for i = 1:na
    fprintf(fid, '%g,%.6e,%.6f,%.6e\n', [ALPHAS(i)*ones(numel(res{i,1}),1) res{i,1} res{i,2} res{i,3}].');
end
fclose(fid);

%% Damping at chosen amplitudes 
fprintf('damping at chosen amplitudes\n        A       QSMA ');
fprintf('  a=%-6g', ALPHAS); fprintf('\n');
for A0 = [5e-5 1e-4 2e-4 5e-4 1e-3 1.5e-3]
    fprintf('%9.1e %.3e', A0, value_at(A0,qq,zq));
    for i = 1:na
        fprintf(' %.3e', value_at(A0,res{i,1},res{i,3}));
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
Ns  = round(1/f_guess/dt);                 % samples per cycle
ntrim = trim_cycles*Ns;                    % Hilbert end effects
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