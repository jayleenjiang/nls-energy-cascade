NB = 80; NTOT = NB^3; T_eq = 6; theta = 20;
idx  = (0:NTOT-1)';
Ia_c = (floor(idx/6400)       + 0.5) * 4/NB;
Ib_c = (mod(floor(idx/80),NB) + 0.5) * 4/NB;

configs = {
    4,  @(j)sprintf('n15_j%d.hist',j),      @(j)sprintf('n15_eq_j%d.hist',j),      'n15';
    6,  @(j)sprintf('simd_n25_j%d.hist',j), @(j)sprintf('simd_n25_eq_j%d.hist',j), 'n25';
    12, @(j)sprintf('simd_n50_j%d.hist',j), @(j)sprintf('simd_n50_eq_j%d.hist',j), 'n50';
    };

figure('Position',[60 60 1600 480]);
for s = 1:3
    j      = configs{s,1};
    fn_q   = configs{s,2};
    fn_p   = configs{s,3};
    tag    = configs{s,4};

    A = readmatrix(fn_q(j),'FileType','text','CommentStyle','#');
    B = readmatrix(fn_p(j),'FileType','text','CommentStyle','#');
    A = A(~any(isnan(A),2),:);  B = B(~any(isnan(B),2),:);
    Q = A(:,4);  P = B(:,4);
    qn = Q/sum(Q);  pn = P/sum(P);

    mask = (Q > 50) & (P > 50);
    cf = polyfit(log(pn(mask)), log(qn(mask)), 1);
    Ppred = polyval(cf, log(pn));
    sQ = -log(qn);  sP = -Ppred;
    resid = sP - sQ;
    resid(Q < 50 | P == 0) = NaN;
    R = reshape(resid(theta:NB:NTOT), NB, NB);
    ax = (0:NB-1)*4/NB;

    subplot(1,3,s);
    mesh(ax, ax, R);
    xlabel('I_b'); ylabel('I_a'); zlabel('log q - fit');
    title(sprintf('%s, a=%d: residual (x=%.3f)', tag, j, cf(1)));
end