
NB = 80; NTOT = NB^3; T_eq = 6; theta = 20;   
idx  = (0:NTOT-1)';
Ia_c = (floor(idx/6400)       + 0.5) * 4/NB;
Ib_c = (mod(floor(idx/80),NB) + 0.5) * 4/NB;

% sites = [4 7 11]; fn_q=@(j)sprintf('n15_j%d.hist',j);
% fn_p=@(j)sprintf('n15_eq_j%d.hist',j); tag='n15';

sites = [6 12 18]; fn_q=@(j)sprintf('simd_n25_j%d.hist',j);
fn_p=@(j)sprintf('simd_n25_eq_j%d.hist',j); tag='n25';

% sites = [12 24 36]; fn_q=@(j)sprintf('simd_n50_j%d.hist',j);
% fn_p=@(j)sprintf('simd_n50_eq_j%d.hist',j); tag='n50';

% sites = [24 48 72]; fn_q=@(j)sprintf('n100_j%d.hist',j);
%  fn_p=@(j)sprintf('n100_eq_j%d.hist',j);   tag='n100';

for s = 1:numel(sites)
    j = sites(s);
    A = readmatrix(fn_q(j),'FileType','text','CommentStyle','#');
    B = readmatrix(fn_p(j),'FileType','text','CommentStyle','#');
    A = A(~any(isnan(A),2),:);  B = B(~any(isnan(B),2),:);
    Q = A(:,4);  P = B(:,4);
    qn = Q/sum(Q);  pn = P/sum(P);

    % % WLS box fit 
    % bx = (Ia_c<2.5)&(Ib_c<2.5)&(Q>0)&(P>0);
    % bb = lscov([ones(nnz(bx),1) log(pn(bx))], log(qn(bx)), Q(bx));
    % cf = [bb(2) bb(1)]; % [x c]

    % resid = log q - (x log p + c) = sP - sQ
    mask = (Q > 50) & (P > 50);
    cf = polyfit(log(pn(mask)), log(qn(mask)), 1);
    Ppred = polyval(cf, log(pn));
    sQ = -log(qn);  sP = -Ppred;
    resid = sP - sQ;                              % data - pred

    resid(Q < 50 | P == 0) = NaN;

    R = reshape(resid(theta:NB:NTOT), NB, NB);
    ax = (0:NB-1)*4/NB;                         

    figure('Position',[60 60 620 480]);
    mesh(ax, ax, R);
    xlabel('I_b'); ylabel('I_a'); zlabel('log q - fit');
    title(sprintf('%s, a=%d: residual (x=%.3f)', tag, j, cf(1)));
end