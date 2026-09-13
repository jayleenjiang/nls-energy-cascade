% log-log fit + heatmaps, batch over sites
% configurable (q, p) file pair, loops over sites, adds a residual panel.
% Source of all x / R^2 / T_0 numbers in the report.
%
%--- n=25 main (NESS vs Eq6) ---
% sites  = [6 12 18];
% fn_q   = @(j) sprintf('simd_n25_j%d.hist', j);      % NESS (q)
% fn_p   = @(j) sprintf('simd_n25_eq_j%d.hist', j);   % Eq (p_6)
% tag    = 'n25';

% --- n=50 (NESS vs Eq6) ---
% sites  = [12 24 36];
% fn_q   = @(j) sprintf('simd_n50_j%d.hist', j);
% fn_p   = @(j) sprintf('simd_n50_eq_j%d.hist', j);
% tag    = 'n50';

% --- n=100 (NESS vs Eq6) ---
% sites  = [24 48 72];
% fn_q   = @(j) sprintf('n100_j%d.hist', j);
% fn_p   = @(j) sprintf('n100_eq_j%d.hist', j);
% tag    = 'n100';

% --- (T_L,T_R) = (8,4) vs Eq6 baseline ---
% sites  = [6 12 18];
% fn_q   = @(j) sprintf('simd_T8_4_j%d.hist', j);
% fn_p   = @(j) sprintf('simd_n25_eq_j%d.hist', j);
% tag    = 'T84';

% --- (T_L,T_R) = (9,3) vs Eq6 baseline ---
% sites  = [6 12 18];
% fn_q   = @(j) sprintf('simd_T9_3_j%d.hist', j);
% fn_p   = @(j) sprintf('simd_n25_eq_j%d.hist', j);
% tag    = 'T93';

% % --- control 1: NESS vs p_7.12 ---
sites = [6 12 18];
fn_q  = @(j) sprintf('simd_n25_j%d.hist', j);
fn_p  = @(j) sprintf('simd_eq_T712_j%d.hist', j);
tag   = 'C712';

% --- control 2 (null): p_7.12 vs p_6 ---
% sites = [6 12 18];
% fn_q  = @(j) sprintf('simd_eq_T712_j%d.hist', j);
% fn_p  = @(j) sprintf('simd_n25_eq_j%d.hist', j);
% tag   = 'fam';


NB = 80; NTOT = NB^3;
T_eq = 6;
theta = 20;          % theta-bin for heatmap slice
idx  = (0:NTOT-1)';                              % bin centers (for box fit)
Ia_c = (floor(idx/6400)       + 0.5) * 4/NB;
Ib_c = (mod(floor(idx/80),NB) + 0.5) * 4/NB;

for s = 1:numel(sites)
    j = sites(s);
    A = readmatrix(fn_q(j), 'FileType','text', 'CommentStyle','#');
    B = readmatrix(fn_p(j), 'FileType','text', 'CommentStyle','#');
    A = A(~any(isnan(A),2), :);
    B = B(~any(isnan(B),2), :);
    assert(size(A,1) == NTOT && size(B,1) == NTOT, 'row count mismatch');
    Q = A(:,4);                          % q counts
    P = B(:,4);                          % p counts

    %% ---- WLS fit + 误差棒 (weights = counts = inverse variance of log-density) ----
    % 全支撑 (count>50)
    m  = (Q > 50) & (P > 50);
    qn = Q / sum(Q);  pn = P / sum(P);
    u  = log(pn(m));  v = log(qn(m));  w = Q(m);
    [b, se] = lscov([ones(size(u)) u], v, w);     % b = [c; x]
    x_fit = b(2);  se_x = se(2);
    T0    = T_eq / x_fit;          se_T0  = T0 * se_x / abs(x_fit);
    cf    = [b(2) b(1)];                            % [x c] 给下面 heatmap 的 polyval 用
    vhat  = [ones(size(u)) u]*b;  vbar = sum(w.*v)/sum(w);
    R2w   = 1 - sum(w.*(v-vhat).^2) / sum(w.*(v-vbar).^2);   % weighted R^2

    % 核心 box (I_a,I_b<2.5)
    bx = (Ia_c<2.5) & (Ib_c<2.5) & (Q>0) & (P>0);
    ub = log(pn(bx)); vb = log(qn(bx)); wb = Q(bx);
    [bb, seb] = lscov([ones(size(ub)) ub], vb, wb);
    xb   = bb(2);  se_xb = seb(2);
    T0box= T_eq / xb;             se_T0b = T0box * se_xb / abs(xb);
    vhb  = [ones(size(ub)) ub]*bb; vbb = sum(wb.*vb)/sum(wb);
    R2box= 1 - sum(wb.*(vb-vhb).^2) / sum(wb.*(vb-vbb).^2);

    fprintf(['%s a=%2d | full: x=%.4f±%.4f  T0=%.3f±%.3f  R2w=%.5f', ...
        ' | box: x=%.4f±%.4f  T0=%.3f±%.3f  R2w=%.5f\n'], ...
        tag, j, x_fit, se_x, T0, se_T0, R2w, xb, se_xb, T0box, se_T0b, R2box);


    %% heatmaps: data | rescaled baseline | residual
    Ppred = exp(polyval(cf, log(pn)));
    Ppred(P == 0) = 0;
    nlQ = -log(qn);     nlQ(Q == 0) = nan;
    nlP = -log(Ppred);  nlP(P == 0) = nan;
    sQ = reshape(nlQ(theta:NB:NTOT), NB, NB);
    sP = reshape(nlP(theta:NB:NTOT), NB, NB);
    figure('Position', [50 50 1500 420]);
    subplot(1,3,1);
    imagesc(sQ); colorbar; axis square;
    title(sprintf('a=%d: -log q (NESS)', j));
    subplot(1,3,2);
    imagesc(sP); colorbar; axis square;
    title(sprintf('-log p_6^x (Eq rescaled, x=%.3f)', x_fit));
    subplot(1,3,3);
    imagesc(sP - sQ); colorbar; axis square;
    title('residual (data - pred)');
    sgtitle(sprintf('(%s), a=%d: T_0 = %.2f, R^2 = %.4f', tag, j, T0, R2w));
end