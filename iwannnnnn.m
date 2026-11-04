% This script is for QSMA, non-parametric and 4-paramtric iwan model
clear; clc; close all;
%%
set(0, 'DefaultAxesFontWeight', 'bold', 'DefaultTextFontWeight', 'bold', ...
       'DefaultLegendFontWeight', 'bold');   

INFILE        = 'qr_curve_2614.txt';
N_SLIDERS     = 70;                            
ALPHA_RATIO   = 1.0;                          
KNOT_INTERVAL = 26;                            % selected data points
ALPHA_SWEEP   = [0.5, 0.75, 1.0, 1.25, 1.5];   


%% 1. LOAD
d = readmatrix(INFILE, 'Delimiter', ',', 'NumHeaderLines', 1);  

beta = d(:,3);  q = d(:,4);  tip = d(:,5);
m = q > 0;
beta = beta(m);  q = q(m);  tip = tip(m);
fprintf('QSMA points: %d\n', numel(q));

%% 2. MACROSLIP, omega_0^2, omega_inf^2, phi_max
Kt = diff(beta) ./ diff(q);         % Slope of every segment

% macroslip: estimate K_inf from the tail of the sweep, then find the first
% segment whose slope has settled onto K_inf (within 1%). The segment before
% it contains the stiffness drop, so phi_max is the point before that segment.

kinf_est = median(Kt(end-4:end));                    % The last 5 slopes are on the macroslip plateau, so their median is a good estimate of K∞. 
on_plat  = abs(Kt - kinf_est) / kinf_est < 0.01;     % first segment that has settled onto that slope, 
if ~any(on_plat)
    error('macroslip plateau not found - extend ALPHA_MAX');
end
i_plat = find(on_plat, 1, 'first');       % 1-based index of first plateau segment, A segment counts as “on the plateau” if its slope is within 1% of kinf_est. i_plat is the first such segment
j      = i_plat - 2;                      % so that phi_max = q(j+1) = q(i_plat-1)
w02    = Kt(1);                           % slope of the first segment
winf2  = mean(Kt(j+2:end));
phimax = q(j+1);                       

fprintf('%5s %14s %14s %18s %12s\n', 'point', 'q_r', 'beta', 'slope', 'drop_%');
for i = 1:length(Kt)
    fprintf('%5d %14.6e %14.4f %18.2f %12.6f\n', ...
        i, q(i), beta(i), Kt(i), (Kt(1)-Kt(i))/Kt(1)*100);
end

fprintf('\n');
fprintf('omega_0^2   = %.6e   f_0   = %.4f Hz\n', w02,   sqrt(w02)  / (2*pi));
fprintf('omega_inf^2 = %.6e   f_inf = %.4f Hz\n', winf2, sqrt(winf2)/ (2*pi));
fprintf('phi_max     = %.6e\n', phimax);
fprintf('tip at macroslip = %.4f mm\n', tip(j+1)*1000);

% Segment i_plat runs from point i_plat to i_plat+1. 
% The segment before it, i_plat-1, still contains the last of the stiffness drop. So phimax is taken as q(i_plat-1), 
% the point at the start of that last transition segment, not the start of the plateau itself.

%% 3. f_nl (Eq. 15) and K_T,f
f_nl = beta - winf2 * q;
K_Tf = (f_nl(2) - f_nl(1)) / (q(2) - q(1));
fprintf('\n');
fprintf('K_T,f       = %.6e\n', K_Tf);


%% 4. MICROSLIP BRANCH ONLY
sel      = q <= q(j+1);                  
q_m = q(sel);  f_m = f_nl(sel);
fprintf('microslip points: %d of %d\n', numel(q_m), numel(q));


%% 5. KNOT INTERVAL SCAN  (paper's acceptance criterion)
fprintf('\n');
fprintf('knot interval scan (alpha = %g):\n', ALPHA_RATIO);
fprintf('  ki knots  neg          S     S/K_T,f      chi  verdict\n');
chi_ok = [];
for ki = 2:50
    if floor(numel(q_m) / ki) < 3
        continue;
    end
    [pn, dn, r, nk] = distribution(ALPHA_RATIO, ki, q_m, f_m, N_SLIDERS, phimax);
    S  = K_Tf - sum(r .* dn);
    nn = sum(r(1:end-1) < 0);                        
    sc = sum(diff(sign(diff(r(1:end-1)))) ~= 0);
    ok = (nn == 0) && (0.0 < S) && (S < K_Tf) && (sc == 1);

    cs = NaN;
    [~, imin] = min(r(1:end-1));
    dd = (1:imin)';
    if nn == 0 && all(r(dd) > 0) && numel(dd) > 1
        pf = polyfit(log(pn(dd)), log(r(dd)), 1);
        cs = pf(1);
    end
    if ok && isfinite(cs)
        chi_ok(end+1) = cs; %#ok<SAGROW>
    end
    if ok, verdict = 'accept'; else, verdict = 'reject'; end
    fprintf('  %2d %5d %4d  %.4e  %7.4f  %7.4f  %s\n', ...
        ki, nk, nn, S, S/K_Tf, cs, verdict);
end
if ~isempty(chi_ok)
    fprintf('  chi over accepted intervals: %.4f to %.4f  (mean %.4f)\n', ...
        min(chi_ok), max(chi_ok), mean(chi_ok));
end

%% 6. CHOSEN FIT
[phi, dphi, rho, nknot] = distribution(ALPHA_RATIO, KNOT_INTERVAL, q_m, f_m, N_SLIDERS, phimax);
sum_rho = sum(rho .* dphi);
S       = K_Tf - sum_rho;

fprintf('\n');
fprintf('chosen: knot interval %d -> %d knots\n', KNOT_INTERVAL, nknot);
fprintf('  sum rho_n dphi_n = %.6e\n', sum_rho);
fprintf('  S                = %.6e\n', S);
fprintf('  K_T,rho          = %.6e   (Eq. 23, = K_T,f)\n', sum_rho + S);
fprintf('  S / K_T,f        = %.4f\n', S / K_Tf);
fprintf('  negative rho     = %d of %d\n', sum(rho(1:end-1) < 0), numel(rho) - 1);


%% 7. FOUR-PARAMETER IWAN FIT  (Eqs. 25-27)
chi = NaN; R = NaN; F_S = NaN; bta = NaN;
[~, imin] = min(rho(1:end-1));                           % bottom of the U 
dec = (1:imin)';                                         % falling part only
if all(rho(dec) > 0)                                     % log needs positive values
    A   = [log(phi(dec)), ones(numel(dec), 1)];          % Shetty’s Eq.27 the straight-line fit.
    sol = pinv(A) * log(rho(dec));                       
    chi = sol(1);                                        % slope = chi
    R = exp(sol(2));                                     % R = e^(intercept)
    KT_d = R * phimax^(chi+1) / (chi+1);
    S4   = K_Tf - KT_d;                                  % the leftover stiffness

    % converting [R, S, χ, φmax] into [F_S, K_T, χ, β].
    % χ stays the same nothing to do. Three new numbers are calculated with Segalman’s formulas:

    F_S  = R * phimax^(chi+2) / (chi+2) + S4 * phimax;   % Eq. 24 segalsman,force when everything slips
    bta  = S4 / KT_d;                                    % size of the lump compared with the smooth part
    fprintf('\n');
    fprintf('4-parameter Iwan\n');
    fprintf('chi = %.12f   beta = %.12f\n', chi, bta);
    fprintf('  R    = %.4g \n', R);             
    fprintf('  beta = %8.4f \n', bta);  
    fprintf('  F_S  = %8.2f  \n', F_S);
    fprintf('  K_T  = %.4e    \n', KT_d + S4);
    fprintf('  S4    = %.4e    \n', S4);
else
    fprintf('\n');
    fprintf('4-parameter fit skipped: rho not positive on decreasing branch\n');
end

if isfinite(chi)
    rho4 = R * phi .^ chi;
else
    rho4 = NaN(size(phi));
end
   

%% 8. WRITE
fh = fopen('rho_of_phi.txt', 'w');
fprintf(fh, 'phi,dphi,rho_nonparam,rho_4param\n');
fprintf(fh, '%.18e,%.18e,%.18e,%.18e\n', [phi, dphi, rho, rho4]');
fclose(fh);

fh = fopen('backbone.txt', 'w');
fprintf(fh, 'q_r,beta,f_nl\n');
fprintf(fh, '%.18e,%.18e,%.18e\n', [q, beta, f_nl]');
fclose(fh);


%% 9. FIGURES
f3 = figure('Position', [100 100 1100 450]);

a31 = subplot(1, 2, 1);
a31.Position = [0.09 0.22 0.38 0.72];
plot(q, beta, 'k-', 'LineWidth', 1.4);
xl = xlim; xlim([0 xl(2)]);              
yl = ylim; ylim([0 yl(2)]);
xlabel('Modal Displacement q_r  [kg^{1/2}m]');
ylabel('Modal Force \beta  [kg^{3/2}ms^{-2}]');
panel_label(a31, '(a)');
sci_x(a31);
box on; grid on;                         

a32 = subplot(1, 2, 2);
a32.Position = [0.57 0.22 0.38 0.72];
plot(q, f_nl, 'b-', 'LineWidth', 1.4); hold on;
xl = xlim; xlim([0 xl(2)]);
yl = ylim; ylim([0 yl(2)]);
[~, jm] = min(abs(q - phimax));
plot(q(jm), f_nl(jm), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r', ...
     'MarkerSize', 11, 'LineWidth', 2);
xlabel('Modal Displacement q_r  [kg^{1/2}m]');
ylabel('Nonlinear Modal Force f_{nl}  [kg^{3/2}ms^{-2}]');
panel_label(a32, '(b)');
sci_x(a32);
box on; grid on;                        


drawnow;
pos = a32.Position;  xl = xlim(a32);  yl = ylim(a32);
x2f = @(x) pos(1) + pos(3) * (x - xl(1)) / (xl(2) - xl(1));
y2f = @(y) pos(2) + pos(4) * (y - yl(1)) / (yl(2) - yl(1));
ARROW_LEN = 0.075;          
GAP       = 0.016;          
xa   = x2f(q(jm));                     
yh   = y2f(f_nl(jm)) - GAP;             
yt0  = yh - ARROW_LEN;                  
annotation('arrow', [xa xa], [yt0 yh], 'Color', 'k', 'LineWidth', 1.5, ...
           'HeadStyle', 'vback2', 'HeadLength', 7, 'HeadWidth', 7);
annotation('textbox', [xa-0.06, yt0-0.045, 0.12, 0.04], 'String', 'macroslip', ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 10);
bold_all(f3);
print(f3, 'fig3.png', '-dpng', '-r300');
save_pdf(f3, 'fig3.pdf');

%% Figure 4 : distribution function
f4 = figure('Position', [100 100 1100 450]);

ZOOM_X   = [-0.5e-5 7e-5];                            
HEADROOM = 6;                
                             
sty_col = [166 77 255; 0 208 224; 0 0 255; 230 168 0; 230 0 0] / 255;
sty_ls  = {'--', '-', '-', '-.', ':'};
sty_lw  = [1.5, 1.5, 2.2, 1.5, 2.2];

a41 = subplot(1, 2, 1);
a41.Position = [0.09 0.22 0.39 0.72];
leg = cell(1, numel(ALPHA_SWEEP));
Psw = cell(1, numel(ALPHA_SWEEP));
Rsw = cell(1, numel(ALPHA_SWEEP));
for ia = 1:numel(ALPHA_SWEEP)
    al = ALPHA_SWEEP(ia);
    [p, dn, r, ~] = distribution(al, KNOT_INTERVAL, q_m, f_m, N_SLIDERS, phimax);
    Psw{ia} = p;  Rsw{ia} = abs(r);
    semilogy(p, abs(r), 'Color', sty_col(ia,:), 'LineStyle', sty_ls{ia}, ...
             'LineWidth', sty_lw(ia)); hold on;
    leg{ia} = sprintf('\\alpha = %g', al);
end
xlabel('Slider Strength \phi  [kg^{1/2}m]');
ylabel('Slider Distribution Function \rho(\phi)');
panel_label(a41, '(a)');
lg41 = legend(a41, leg, 'FontSize', 9);
box on;
a41.FontWeight = 'bold';
xlim(a41, [-0.015*phimax, 1.01*phimax]);  
style_shetty(a41);                         
a41.XTick = (0:2:16) * 1e-4;              
a41.XAxis.Exponent = -4;                   
allr = []; zz = [];
for ia = 1:numel(ALPHA_SWEEP)
    pp = Psw{ia}(:);  rr = Rsw{ia}(:);
    allr = [allr; rr(rr > 0)];                                    
    in = pp >= ZOOM_X(1) & pp <= ZOOM_X(2) & rr > 0;
    zz = [zz; rr(in)];                                            
end
if isempty(zz)                      
    zz = allr;
end
ymax = max(allr);  ylo = min(allr) / 1.2;  yhi = ymax * HEADROOM;
ylim(a41, [ylo yhi]);
ZOOM_Y = [min(zz) / 1.08, max(zz) * 1.08];    
set(a41, 'XLim', a41.XLim, 'YLim', a41.YLim);     

logticks_tex(a41, 1);          
plot(a41, [ZOOM_X(1) ZOOM_X(2) ZOOM_X(2) ZOOM_X(1) ZOOM_X(1)], ...
          [ZOOM_Y(1) ZOOM_Y(1) ZOOM_Y(2) ZOOM_Y(2) ZOOM_Y(1)], ...
          'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
drawnow;
pos  = a41.Position;
frac = @(y) (log10(y) - log10(ylo)) / (log10(yhi) - log10(ylo)); 
yb   = frac(max(ymax, ZOOM_Y(2)) * 1.25);        
top  = 0.97;                                    
lg41.Units = 'normalized';
lw = lg41.Position(3);  lh = lg41.Position(4);
if lh > (top - yb) * pos(4)
    warning('increase HEADROOM');
end
lg41.Position = [pos(1) + 0.97*pos(3) - lw, pos(2) + top*pos(4) - lh, lw, lh];
inset_w = min(0.30*pos(3), lg41.Position(1) - (pos(1) + 0.10*pos(3)) - 0.02*pos(3));
ain = axes('Position', [pos(1) + 0.10*pos(3), pos(2) + (yb + 0.02)*pos(4), ...
                        inset_w,              (top - yb - 0.02)*pos(4)]);
hold(ain, 'on');
for ia = 1:numel(ALPHA_SWEEP)
    semilogy(ain, Psw{ia}, Rsw{ia}, 'Color', sty_col(ia,:), ...
             'LineStyle', sty_ls{ia}, 'LineWidth', sty_lw(ia));
end
set(ain, 'YScale', 'log');
xlim(ain, ZOOM_X);  ylim(ain, ZOOM_Y);
box(ain, 'on');  ain.LineWidth = 2;  ain.FontSize = 8;  ain.FontWeight = 'bold';
style_shetty(ain);
ain.XTick = []; ain.YTick = [];           
ain.YAxis.MinorTick = 'off';

a42 = subplot(1, 2, 2);
a42.Position = [0.58 0.22 0.39 0.72];
semilogy(phi, rho,  'b-',  'LineWidth', 2); hold on;
semilogy(phi, rho4, 'r--', 'LineWidth', 2);
xlabel('Slider Strength \phi  [kg^{1/2}m]');
ylabel('Slider Distribution Function \rho(\phi)');
panel_label(a42, '(b)'); legend({'non-parametric Iwan model', '4-parameter Iwan model'}, ...
                     'FontSize', 9, 'Location', 'north');
box on;
a42.FontWeight = 'bold';
allv = [rho(rho > 0); rho4(isfinite(rho4) & rho4 > 0)];
ylim(a42, [min(allv)/1.3, max(allv)*1.6]);          
style_shetty(a42);
a42.XTick = (0:2:18) * 1e-4;
a42.XAxis.Exponent = -4;
logticks_tex(a42, [1 2 3 4 6]);                     
bold_all(f4);
print(f4, 'fig4.png', '-dpng', '-r300');
fprintf('\nwritten: fig3.png, fig4.png\n');


%% LOCAL FUNCTIONS

function [pn, dn] = slider_grid(al, N, phimax)
    if al == 1                                             % If α = 1, all blocks are equal, so d1 = phimax/70.
        d1 = phimax / N;                               
    else                                                                                                         
        d1 = phimax * (al - 1) / (al^N - 1);               % Eq. 17,  %If α ≠ 1, the blocks grow geometrically, and d1 comes from Eq. 17.
    end
    ju = [0; cumsum(d1 * al .^ (0:N-1)')];                 % Eq. 18, N+1 edges
    pn = [ju(1:end-1) + diff(ju) / 2.0; phimax];           % N centres + phi_max
    dn = [diff(ju); 0.0];                                  % N widths + 0
end

function [pn, dn, rho, nk] = distribution(al, ki, q_m, f_m, N, phimax)
% Eqs. 20-21 on the microslip branch.
    kk = (1:ki:numel(q_m))';
    kk(end) = numel(q_m);                    
    xk = [0.0; q_m(kk)];                    
    fk = [0.0; f_m(kk)];
    pp = spline(xk, fk);                     
    [pn, dn] = slider_grid(al, N, phimax);

    % second derivative of the spline at pn
    [breaks, coefs] = unmkpp(pp);
    b   = breaks(:);
    idx = arrayfun(@(x) find(b <= x, 1, 'last'), pn);
    idx = min(idx, numel(b) - 1);
    t   = pn - b(idx);
    rho = -(6 * coefs(idx, 1) .* t + 2 * coefs(idx, 2));
    nk  = numel(xk);
end

function sci_x(ax)
    mx = max(abs(ax.XLim));
    if mx > 0
        ax.XAxis.Exponent = floor(log10(mx));
    end
end

function logticks_tex(ax, cands)
    yl  = ax.YLim;
    dec = floor(log10(yl(1))) : ceil(log10(yl(2)));
    tk  = [];  lab = {};
    for dd = dec
        for c = cands
            v = c * 10^dd;
            if v >= yl(1) && v <= yl(2)
                tk(end+1) = v; 
                if c == 1
                    lab{end+1} = sprintf('10^{%d}', dd); 
                else
                    lab{end+1} = sprintf('%g\\times10^{%d}', c, dd); 
                end
            end
        end
    end
    if numel(tk) >= 2
        ax.YTick = tk;
        ax.YTickLabel = lab;
        ax.TickLabelInterpreter = 'tex';
    end
    try ax.YAxis.Exponent = 0; catch, end
end

function bold_all(fig)
    for ax = findall(fig, 'Type', 'axes')'
        ax.FontWeight = 'bold';
        ax.XLabel.FontWeight = 'bold';
        ax.YLabel.FontWeight = 'bold';
        ax.Title.FontWeight  = 'bold';
    end
    set(findall(fig, 'Type', 'legend'), 'FontWeight', 'bold');
    set(findall(fig, 'Type', 'text'),   'FontWeight', 'bold');
end

function style_shetty(ax)
    ax.TickDir = 'in';
    ax.LineWidth = 0.8;
    ax.Layer = 'top';
    grid(ax, 'on');
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.15;
    ax.XMinorGrid = 'off';  ax.YMinorGrid = 'off';
    ax.MinorGridLineStyle = 'none';
    ax.XAxis.MinorTick = 'off';
    ax.YAxis.MinorTick = 'off';
end

function save_pdf(fig, name)
    drawnow;
    fig.Units = 'inches';
    pos = fig.Position;
    fig.PaperUnits = 'inches';
    fig.PaperSize = pos(3:4);
    fig.PaperPosition = [0 0 pos(3:4)];
    print(fig, name, '-dpdf', '-vector');
end

function panel_label(ax, s)
    text(ax, 0.5, -0.19, s, 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
         'FontSize', 12, 'FontWeight', 'bold');
end


%% Normalised non-parametric density and power-law fit
phin = phi(1:end-1);  rhon = rho(1:end-1);
[~, i_min] = min(rhon);
FS_y = f_nl(j+1);                                  
fprintf('chi = %.4f, R = %.4g\n', chi, R);
fprintf('F_S/(K_T*phi_max) = %.4f\n', FS_y/(K_Tf*phimax));
fprintf('S4/K_T            = %.4f\n', S4/K_Tf);
figure; hold on; grid on;
plot(phin/phimax, rhon*phimax/K_Tf, 'b-o', 'DisplayName', 'non-parametric');
plot(phin(1:i_min)/phimax, R*phin(1:i_min).^chi*phimax/K_Tf, 'b--', 'DisplayName', 'power-law fit');
set(gca, 'YScale', 'log'); legend show;
xlabel('\phi / \phi_{max}'); ylabel('\rho \phi_{max} / K_T');