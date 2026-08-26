filename = 'simd_long_j12.hist';
NB = 80;
fid = fopen(filename, 'r');

A = [];
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || startsWith(line, '#'), continue; end
    nums = sscanf(line, '%f');
    if numel(nums) == 4
        A = nums(:).';
        rest = textscan(fid, '%f %f %f %f');
        A = [A; [rest{1} rest{2} rest{3} rest{4}]];
        break;
    end
end

fclose(fid);
hist3d = zeros(NB, NB, NB);
for r = 1:size(A,1)
    hist3d(A(r,1)+1, A(r,2)+1, A(r,3)+1) = A(r,4);
end
dI  = 4/NB;  dth = 2*pi/NB;
I_c  = dI  * ((1:NB) - 0.5);
th_c = -pi + dth * ((1:NB) - 0.5);
P3 = hist3d / (sum(hist3d(:)) * dI * dI * dth);

% cos = -1, 0, 1, 0
target_thetas = [-pi, -pi/2, 0, pi/2];   
ns = numel(target_thetas);
I_max_plot = 2.5;
keep = (I_c < I_max_plot);
[Iam, Ibm] = ndgrid(I_c(keep), I_c(keep));
min_count = 20;

figure('Position', [50 50 1800 900]);
Hmax = 0;
for s = 1:ns
    [~, it] = min(abs(th_c - target_thetas(s)));
    cv = cos(th_c(it));
    cnt = hist3d(keep, keep, it);
    Ps  = P3(keep, keep, it);
    nlP = nan(size(cnt));
    nlP(cnt >= min_count) = -log(Ps(cnt >= min_count));

    % H_2mode = (1/2)(Ia+Ib)^2 - (1/4)(Ia^2+Ib^2) + Ia Ib cos(theta)
    H = 0.5*(Iam+Ibm).^2 - 0.25*(Iam.^2+Ibm.^2) + Iam.*Ibm*cv;
    Hmax = max(Hmax, max(H(:)));

    subplot(2, ns, s);
    surf(Iam, Ibm, nlP, 'EdgeColor', 'none');
    xlabel('I_a'); ylabel('I_b'); zlabel('-log P');
    title(sprintf('data: cos\\theta = %.2f', cv));
    view(45, 30); grid on;

    subplot(2, ns, s + ns);
    surf(Iam, Ibm, H, 'EdgeColor', 'none');
    xlabel('I_a'); ylabel('I_b'); zlabel('H_{2mode}');
    title(sprintf('H_{2mode}: cos\\theta = %.2f', cv));
    view(45, 30); grid on;
end

for s = 1:ns
    subplot(2, ns, s + ns); zlim([0, Hmax]);
end
