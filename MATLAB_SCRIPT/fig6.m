clear; clc; close all;
%% fig6 modal displacement vs time,
% INPUTS:
%   - iwan_decays_ds.csv 

% WHAT IT DOES:
%   1. joins the FE pieces by TOTAL time 
%   2. checks the overlap between consecutive pieces (same time computed twice)
%   3. converts tip U2 to modal displacement: q = (U2 - U2_preload) / phi_tip
%   4. plots FE and Iwan in ONE panel with a zoom inset

% OUTPUTS: fig6.png, fe_modal_decay.csv

%% 
FE_PIECES  = {'u2_job57.rpt', 'u2_job57r5.rpt', 'u2_job57r6.rpt', ...
              'u2_job57r9.rpt', 'u2_job57r10.rpt', 'u2_job57r11.rpt', 'u2_job57r12.rpt', ...
              'u2_job57r13.rpt', 'u2_job57r14.rpt', 'u2_job57r15.rpt'};

IWAN_FILE  = 'iwan_decays_ds.csv';   
T_RELEASE  = 2.0;                    
U2_PRELOAD = 5.32251e-7;             
PHI_TIP    = 0.388016;               
ZOOM       = [0.13 0.16];           
FIG_W      = 9;                     
FIG_H      = 7.5;                    
FS         = 8;                      
LW         = 1.5;                   
OUT        = 'fig6_decay';

set(0, 'DefaultAxesFontWeight', 'bold', 'DefaultTextFontWeight', 'bold');

%% 1) read and join FE pieces

np_ = numel(FE_PIECES);
T = cell(np_,1);  U = cell(np_,1);
for k = 1:np_
    [T{k}, U{k}] = read_rpt(FE_PIECES{k});
end
t_fe = [];  u_fe = [];
for k = 1:np_
    if k < np_, t_end = T{k+1}(1); else, t_end = inf; end
    keep = T{k} < t_end;
    t_fe = [t_fe; T{k}(keep)];   u_fe = [u_fe; U{k}(keep)];          
    fprintf('%20s: total time %.6f -> %.6f  (used up to %.6f), %d points\n', ...
            FE_PIECES{k}, T{k}(1), T{k}(end), min(t_end, T{k}(end)), nnz(keep));
end

%% 2) overlap check between consecutive pieces
for k = 1:np_-1
    t1 = T{k}; u1 = U{k}; t2 = T{k+1}; u2 = U{k+1};
    lo = max(t1(1), t2(1));  hi = min(t1(end), t2(end));
    if hi > lo
        m   = t2 >= lo & t2 <= hi;
        u1i = interp_u(t1, u1, t2(m));
        amp = max(abs(u1i));
        err = max(abs(u1i - u2(m)));
        fprintf('overlap %s vs %s: %.5f-%.5f, max |diff| = %.3e m = %.3f%% of max |U2|\n', ...
                FE_PIECES{k}, FE_PIECES{k+1}, lo, hi, err, 100*err/amp);
        fprintf('   value at restart time %.6f: first piece %.6e, second piece %.6e\n', ...
                t2(1), interp_u(t1, u1, t2(1)), u2(1));
    end
end

%% 3) convert to modal displacement
tau_fe = t_fe - T_RELEASE;    
q_fe   = (u_fe - U2_PRELOAD) / PHI_TIP;
fid = fopen('fe_modal_decay.csv', 'w');         
fprintf(fid, 'time_from_release_s,q_modal\n');  
fprintf(fid, '%.18e,%.18e\n', [tau_fe q_fe].'); 
fclose(fid);
fprintf('FE decay: %.5f -> %.5f s, q0 = %.6e\n', tau_fe(1), tau_fe(end), q_fe(1));

%% 4) plot 
iw = dlmread(IWAN_FILE, ',');
[t_iw, iu] = unique(iw(:,1));             
q_np = iw(iu, 2);
TMAX = 0.5;                                % Keep 0 to 0.5 s for both curves.
m  = t_iw <= TMAX;
mf = tau_fe <= TMAX;

f6 = figure('Units', 'centimeters', 'Position', [2 2 FIG_W FIG_H], 'Color', 'w');
set(f6, 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 FIG_W FIG_H], 'PaperSize', [FIG_W FIG_H]);
ax = axes(f6); hold(ax, 'on'); box(ax, 'on');
h1 = plot(ax, t_iw(m),   q_np(m),   '-',  'Color', 'b', 'LineWidth', LW);
h2 = plot(ax, tau_fe(mf), q_fe(mf), '--', 'Color', 'r', 'LineWidth', LW);
qmax = max(abs([q_np(m); q_fe(mf)]));
YL   = 1.25 * qmax * [-1 1];              
xlim(ax, [0 TMAX]);  ylim(ax, YL);
set(ax, 'XTick', 0:0.1:TMAX, 'FontSize', FS, 'FontWeight', 'bold', 'LineWidth', 0.75);
grid(ax, 'on');  try set(ax, 'GridAlpha', 0.15); catch, end    
xlabel(ax, 'Time [s]', 'FontSize', FS, 'FontWeight', 'bold');
ylabel(ax, 'Modal Displacement [kg^{1/2}m]', 'FontSize', FS, 'FontWeight', 'bold');
sci_y(ax);

zi = t_iw >= ZOOM(1) & t_iw <= ZOOM(2);
zr = 1.15 * max(abs(q_np(zi)));
rectangle(ax, 'Position', [ZOOM(1), -zr, diff(ZOOM), 2*zr], 'EdgeColor', 'k', 'LineWidth', 1.2);

lg = legend(ax, [h1 h2], {'non-parametric Iwan model', 'Dynamic analysis of full FE model'}, ...
            'Location', 'southeast', 'FontSize', FS - 1, 'FontWeight', 'bold');


drawnow;
pos  = get(ax, 'Position');
x0   = 0.40;                                         
env  = max(abs(q_np(t_iw >= x0*TMAX & t_iw <= TMAX)));
yb   = (env - YL(1)) / diff(YL) + 0.04;             
ain  = axes(f6, 'Position', [pos(1) + x0*pos(3), pos(2) + yb*pos(4), ...
                             (0.98 - x0)*pos(3),   (0.97 - yb)*pos(4)]);
hold(ain, 'on'); box(ain, 'on');
zf = tau_fe >= ZOOM(1) & tau_fe <= ZOOM(2);
plot(ain, t_iw(zi),   q_np(zi),   '-',  'Color', 'b', 'LineWidth', LW);
plot(ain, tau_fe(zf), q_fe(zf),   '--', 'Color', 'r', 'LineWidth', LW);
xlim(ain, ZOOM);  ylim(ain, [-zr zr]);
set(ain, 'XTick', [], 'YTick', [], 'LineWidth', 1.2);

save_fig(f6, [OUT '']);
fprintf('inset window: %.3f to %.3f s\n', ZOOM(1), ZOOM(2));
fprintf('saved %s.png/ and fe_modal_decay.csv\n', OUT);


%% FUNCTIONS 
function [t, y] = read_rpt(fname)
    fid = fopen(fname, 'r');
    if fid < 0, error('cannot open %s', fname); end
    t = [];  y = [];
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        v = sscanf(ln, '%f');
        if numel(v) == 2 && numel(strsplit(strtrim(ln))) == 2
            t(end+1,1) = v(1);  y(end+1,1) = v(2);              
        end
    end
    fclose(fid);
    [t, o] = sort(t);  y = y(o);
end

function v = interp_u(t, u, tq)
    [tt, iu] = unique(t);
    v = interp1(tt, u(iu), tq, 'linear', 'extrap');
end

function sci_y(ax)
    try
        yl = get(ax, 'YLim');  mx = max(abs(yl));
        if mx > 0, ax.YAxis.Exponent = floor(log10(mx)); end
    catch
    end
end

function save_fig(fig, name)
    try
        exportgraphics(fig, [name '.png'], 'Resolution', 600);
        exportgraphics(fig, [name '.pdf'], 'ContentType', 'vector');
    catch
        print(fig, [name '.png'], '-dpng', '-r600');
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, [name '.pdf'], '-dpdf');
    end
end
