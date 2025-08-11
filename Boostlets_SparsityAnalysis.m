%BOOSTLETS_SPARSITYANALYSIS 
% -----------------------------------------------------------------------------------
%       This code can be used to reproduce the sparsity analyses in the paper:
%   Zea, Laudato, Andén, "A boostlet transform for wave-based acoustic signal 
%   processing in space-time", arXiv:2403.11362v2, 2025. 
% -----------------------------------------------------------------------------------
% EXAMPLE OF USAGE: 
%   To reproduce the plots in Figure 7, run the command: 
%
%       > Boostlets_SparsityAnalysis;
%
%   The output produces a figure (Figure 1) corresponding to Figure 7 in the paper. 
%   A table including the l1-norm of the 10.0000 largest coefficients, corresponding 
%   to Table I in the manuscript, is output in the Command Window. 
% -----------------------------------------------------------------------------------
% DEPENDENCIES
%       Folders:    'Datasets'          measured acoustic fields in three rooms [1]
%                   'Toolboxes'         curvelets, wave atoms, and shearlets 
% 
%   The curvelet toolbox (CurveLab-2.1.3) is taken from [2], the wave atoms toolbox 
%   (WaveAtom-1.1.1) from [3], and the shearlet toolbox (FFST) from [4]. 
% -----------------------------------------------------------------------------------
% REFERENCES: 
% [1] E. Zea, 'Compressed sensing of impulse responses in rooms of unknown properties 
%     and contents', J. Sound Vib 459, 114871 (2019). 
% [2] E. Candès, L. Demanet, D. Donoho, L. Ying, 'Fast Discrete Curvelet Transforms'
%     Multiscale Modeling & Simulation 5(3), 861–899 (2006). 
% [3] L. Demanet, L. Ying, 'Wave atoms and sparsity of oscillatory patterns', Appl. 
%     Comput. Harmon. Anal. 23(3), 368–387 (2007)
% [4] S. Häuser, G. Steidl, 'Fast finite shearlet transform: a tutorial', ArXiv 
%     1202.1773, 1-41 (2012)
% -----------------------------------------------------------------------------------
% Code: E. Zea
% Code history: 
% - Version 001 [March 15, 2024]
% - Version 002 [August 11, 2025]: Renamed repo and included sparse
%                                  reconstruction errors
% -----------------------------------------------------------------------------------
% CONTACT: Elias Zea (zea@kth.se)
%          Marcus Wallenberg Laboratory for Sound and Vibration Research
%          KTH Royal Institute of Technology
%          Teknikringen 8
%          10044 Stockholm, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%               BEGIN CODE...               %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; 

% add dependencies
addpath(genpath('dependencies'));
fprintf('Dependencies added to MATLAB path successfully.\n');

% choose rooms
room = {'Balder','Freja','Munin','Munin'};

% starting time samples
T_start = [0,0,0,2500];

% no. decomposition scales (wavelets, shearlets, boostlets)
L = 3; % integer between [2,4] — otherwise, some systems do not work!

% sampling parameters and 2D space-time coordinates
N = 100; % hardwired from experimental data
N_max = 1e4; % limit to X largest coefficients
eval_array = unique(round(logspace(0,log10(N_max),49)));
T = N;
M = N;
dx = 0.03;
fs = 11250;
x = linspace(-M*dx/2,M*dx/2,M);

% Preallocate arrays to store l1- and l2-norm results
l1_norms = zeros(length(T_start), 6); % 6 decomposition methods
l2_norms = zeros(length(T_start), 6);

% Preallocate boostlet coefficients
BT = [];
plot_count = 1;

for splt = 1:length(T_start)
    %% Load acoustic field
    load(['dependencies/Datasets/' room{splt} 'RIR.mat'])
    y = out.image(1+T_start(splt):T+T_start(splt),:);

    % define space-time coordinates
    t = 1e3*linspace(0+T_start(splt)/fs,(T+T_start(splt))/fs,T); % in ms
    [XX,TT] = meshgrid(x,t);

    % pad each row to next power‐of‐two for wave atoms
    T_ext = 2^nextpow2(N);
    y_ext = zeros(T_ext,T_ext);
    dN = (T_ext - N)/2;    % here (128 − 100)/2 = 14
    for ii = 1:N
        y_ext(ii,:) = lpbp1D(y(ii,:), T_ext, 16);
    end

    %% COMPUTE SPARSITY (L1-norm) & RECONSTRUCTION (L2-norm)
    %% DAUBECHIES45 WAVELETS
    [WT, L_db] = wavedec2(y, L, 'db45');
    WT_flat    = WT(:);
    [~, id_db] = sort(abs(WT_flat), 'descend');
    nEval      = numel(eval_array);
    l1_norm_db = zeros(nEval,1);
    RMSRE_db   = zeros(nEval,1);
    for k = 1:nEval
        m = eval_array(k);
        % top-m magnitudes
        sel_db = WT_flat(id_db(1:m));       
        sel_db = sel_db ./ norm(sel_db);
        l1_norm_db(k) = norm(sel_db,1);
        % reconstruct from top-m only
        rec_flat      = zeros(size(WT_flat));
        rec_flat(id_db(1:m)) = WT_flat(id_db(1:m));
        WT_rec        = reshape(rec_flat, size(WT));
        y_rec_db      = waverec2(WT_rec, L_db, 'db45');
        % compute RMSRE
        RMSRE_db(k)   = computeRMSRE(y, y_rec_db);
    end

    %% MEYER WAVELETS
    [WT2, L_Meyer] = wavedec2(y, L, 'dmey');
    WT2_flat       = WT2(:);
    [~, id_meyer]  = sort(abs(WT2_flat), 'descend');
    nEval          = numel(eval_array);
    l1_norm_meyer  = zeros(nEval,1);
    RMSRE_meyer    = zeros(nEval,1);
    for k = 1:nEval
        m = eval_array(k);
        % top-m magnitudes
        sel_mey = WT2_flat(id_meyer(1:m));
        sel_mey = sel_mey ./ norm(sel_mey);
        l1_norm_meyer(k) = norm(sel_mey, 1);
        % reconstruct from top-m only
        rec_flat = zeros(size(WT2_flat));
        rec_flat(id_meyer(1:m)) = WT2_flat(id_meyer(1:m));
        rec_WT2 = reshape(rec_flat, size(WT2));
        % — inverse transform & crop back to N×N
        y_rec_meyer = waverec2(rec_WT2, L_Meyer, 'dmey');
        % — compute RMSRE
        RMSRE_meyer(k) = computeRMSRE(y, y_rec_meyer);
    end

    %% CURVELETS
    CT = fdct_wrapping(y,0,2,L);
    [flatC, bandList] = flattenCurveletCoeffs(CT);
    [~,id_clet]  = sort(abs(flatC),'descend');
    l1_norm_clet = zeros(numel(eval_array),1);
    RMSRE_clet   = zeros(numel(eval_array),1);
    for k = 1:numel(eval_array)
        m = eval_array(k);
        % top-m magnitudes
        sel_clet = flatC(id_clet(1:m));
        sel_clet = sel_clet./norm(sel_clet);
        l1_norm_clet(k) = norm(sel_clet,1);
        % reconstruct from top-m only
        recFlat = zeros(size(flatC));
        recFlat(id_clet(1:m)) = flatC(id_clet(1:m));
        CT_rec = reconstructCurveletCoeffs(recFlat,bandList,CT);
        y_rec_clet = ifdct_wrapping(CT_rec);
        % compute RMSRE
        RMSRE_clet(k) = computeRMSRE(y,y_rec_clet);
    end

    %% SHEARLETS
    [ST, shearlet_bases] = shearletTransformSpect(y, L);
    ST_flat = ST(:);
    [~, id_shear] = sort(abs(ST_flat), 'descend');
    nEval         = numel(eval_array);
    l1_norm_shear = zeros(nEval,1);
    RMSRE_shear   = zeros(nEval,1);
    for k = 1:nEval
        m = eval_array(k);
        % top-m magnitudes
        sel_shear = ST_flat(id_shear(1:m));
        sel_shear = sel_shear ./ norm(sel_shear);
        l1_norm_shear(k) = norm(sel_shear, 1);
        % reconstruct from top-m only
        rec_flat = zeros(size(ST_flat));
        rec_flat(id_shear(1:m)) = ST_flat(id_shear(1:m));
        ST_rec = reshape(rec_flat, size(ST));
        y_rec_shear = inverseShearletTransformSpect(ST_rec, shearlet_bases);
        % compute RMSRE
        RMSRE_shear(k) = computeRMSRE(y, y_rec_shear);
    end

    %% WAVE ATOMS
    pat = 'p'; tp = 'complex';
    c = fwa2(y_ext,pat,tp);
    % flatten all coeffs into one vector
    watCoeffs = [];
    for is = 1:size(c,1)
        for dir = 1:size(c,2)
            for ii = 1:size(c{is,dir},1)
                for jj = 1:size(c{is,dir},2)
                    watCoeffs = [watCoeffs; c{is,dir}{ii,jj}(:)];
                end
            end
        end
    end
    [~,id_wat] = sort(abs(watCoeffs),'descend');
    l1_norm_wat   = zeros(numel(eval_array),1);
    RMSRE_wat     = zeros(numel(eval_array),1);
    for k = 1:numel(eval_array)
        m = eval_array(k);
        % top-m magnitudes
        sel_wat = watCoeffs(id_wat(1:m));
        sel_wat = sel_wat./norm(sel_wat);
        l1_norm_wat(k) = norm(sel_wat,1);
        % reconstruct from top-m only
        recFlat = zeros(size(watCoeffs));
        recFlat(id_wat(1:m)) = watCoeffs(id_wat(1:m));
        % rebuild c‐structure
        cc = 1;
        c_rec = c;
        for is = 1:size(c,1)
            for dir = 1:size(c,2)
                for ii = 1:size(c{is,dir},1)
                    for jj = 1:size(c{is,dir},2)
                        sz = numel(c{is,dir}{ii,jj});
                        block = recFlat(cc:cc+sz-1);
                        c_rec{is,dir}{ii,jj} = reshape(block, size(c{is,dir}{ii,jj}));
                        cc = cc + sz;
                    end
                end
            end
        end
        y_rec_wat = iwa2(c_rec,pat,tp);
        y_rec_wat = y_rec_wat(1 : N,  dN+1 : dN+N);
        % compute RMSRE
        RMSRE_wat(k) = computeRMSRE(y,y_rec_wat);
    end

    %% compute boostlet coefficients and sort in descending amplitude
    [BT, phi] = ffbt(y, L);               
    BT_flat = BT(:);
    [~, id_boost] = sort(abs(BT_flat), 'descend');
    nEval         = numel(eval_array);
    l1_norm_boost = zeros(nEval,1);
    RMSRE_boost   = zeros(nEval,1);
    for k = 1:nEval
        m = eval_array(k);
        % top-m magnitudes
        sel_boost = BT_flat(id_boost(1:m));
        sel_boost = sel_boost ./ norm(sel_boost);              
        l1_norm_boost(k) = norm(sel_boost,1);
        % reconstruct from top-m only
        rec_flat = zeros(size(BT_flat));
        rec_flat(id_boost(1:m)) = BT_flat(id_boost(1:m));
        BT_rec = reshape(rec_flat, size(BT));
        y_rec_boost = iffbt(BT_rec, phi);
        % compute RMSRE
        RMSRE_boost(k) = computeRMSRE(y, y_rec_boost);
    end

    %% ——— Gather l1 norms ———
    l1_norms(splt,:) = [ ...
        l1_norm_db(end), ...
        l1_norm_meyer(end), ...
        l1_norm_clet(end), ...
        l1_norm_shear(end), ...
        l1_norm_wat(end), ...
        l1_norm_boost(end) ];

    %% ——— Gather relative reconstruction errors ———
    [~, idx_n] = min(abs(eval_array - 1e3));
    l2_norms(splt,:) = [ ...
        RMSRE_db(idx_n), ...
        RMSRE_meyer(idx_n), ...
        RMSRE_clet(idx_n), ...
        RMSRE_shear(idx_n), ...
        RMSRE_wat(idx_n), ...
        RMSRE_boost(idx_n) ];

    %% ——— Build selected coefficient‐magnitude vectors ———
    selected_wavelet_coeffs_db    = abs(sel_db);
    selected_wavelet_coeffs_meyer = abs(sel_mey);
    selected_waveats_coeffs       = abs(sel_wat);
    selected_curvelet_coeffs      = abs(sel_clet);
    selected_shearlet_coeffs      = abs(sel_shear);
    selected_boostlet_coeffs      = abs(sel_boost);


    %% plot fields and approximation errors
    figure(1);
    % subplot: acoustic field in space-time
    subplot(3,length(T_start),plot_count);
    surf(XX,TT,y,'edgecolor','none');
    view(0,90);
    if splt == 3
        xlabel('$x$ (m)','Interpreter','latex');
    end
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$t$ (ms)','Interpreter','latex');
    colormap gray; axis tight;
    cb = colorbar('NorthOutside');
    ttl = title(cb,'[Pa]','Interpreter', 'latex', 'Position', [70 50 0]);
    cb.TickLabelInterpreter = 'latex';
    set(gca,'fontsize',25,'TickLabelInterpreter', 'latex');
    text(-0.3, 2, ['(' char(96 + splt) ')'], 'Color', 'black', ...
        'FontSize', 30, 'Units', 'normalized');

    % subplot: coefficient decays
    figure(1);
    subplot(3,length(T_start),length(T_start)+plot_count);
    fig = loglog(...
        1:N_max, selected_wavelet_coeffs_db,    'x:', ...
        1:N_max, selected_wavelet_coeffs_meyer, ':', ...
        1:N_max, selected_curvelet_coeffs,      '*--', ...
        1:N_max, selected_shearlet_coeffs,      '-.', ...
        1:N_max, selected_waveats_coeffs,       '+--', ...
        1:N_max, selected_boostlet_coeffs,      '-'  ...
        );
    % style each line
    fig(1).LineWidth = 2.5; fig(1).MarkerSize = 2; fig(1).Color = [0.15,0.82,0.79];
    fig(2).LineWidth = 2.5; fig(2).MarkerSize = 3; fig(2).Color = [0,0,1];
    fig(3).LineWidth = 2.5; fig(3).MarkerSize = 3; fig(3).Color = [0.53,0.48,1];
    fig(4).LineWidth = 2.5; fig(4).MarkerSize = 2; fig(4).Color = [0,0.45,0.74];
    fig(5).LineWidth = 2.5; fig(5).MarkerSize = 3; fig(5).Color = [0.59,0.12,0.75];
    fig(6).LineWidth = 2.5;                        fig(6).Color = [0.89,0.07,0.52];
    grid on;
    xlabel('$n$','Interpreter','latex','FontSize',14);
    set(gca,'fontsize',25,'TickLabelInterpreter', 'latex');
    axis tight;
    xticks([1,10,100,1000,10000]);
    text(-0.3, 1.3, ['(' char(100 + splt) ')'], 'Color', 'black', ...
        'FontSize', 30, 'Units', 'normalized');
    if splt == 3
        xlabel('$n$','Interpreter','latex');
    end
    xlabel('$n$','Interpreter','latex');
    if splt == 1
        ylabel('Coeff. magnitude','Interpreter','latex');
    end

    % subplot: approximation errors
    figure(1);
    subplot(3,length(T_start),2*length(T_start)+plot_count);
    fig = semilogx(...
        eval_array, RMSRE_db,    'x:', ...
        eval_array, RMSRE_meyer, ':', ...
        eval_array, RMSRE_clet,      '*--', ...
        eval_array, RMSRE_shear,      '-.', ...
        eval_array, RMSRE_wat,       '+--', ...
        eval_array, RMSRE_boost,      '-'  ...
        );
    % style each line
    fig(1).LineWidth = 2.5; fig(1).MarkerSize = 2; fig(1).Color = [0.15,0.82,0.79];
    fig(2).LineWidth = 2.5; fig(2).MarkerSize = 3; fig(2).Color = [0,0,1];
    fig(3).LineWidth = 2.5; fig(3).MarkerSize = 3; fig(3).Color = [0.53,0.48,1];
    fig(4).LineWidth = 2.5; fig(4).MarkerSize = 2; fig(4).Color = [0,0.45,0.74];
    fig(5).LineWidth = 2.5; fig(5).MarkerSize = 3; fig(5).Color = [0.59,0.12,0.75];
    fig(6).LineWidth = 2.5;                        fig(6).Color = [0.89,0.07,0.52];
    grid on;
    xlabel('$n$','Interpreter','latex','FontSize',14);
    set(gca,'fontsize',25,'TickLabelInterpreter', 'latex');
    axis tight;
    xticks([1,10,100,1000,10000]);
    text(-0.3, 1.33, ['(' char(104 + splt) ')'], 'Color', 'black', ...
        'FontSize', 30, 'Units', 'normalized');
    if splt == 3
        xlabel('$n$','Interpreter','latex');
    end
    xlabel('$n$','Interpreter','latex');
    if splt == 1
        ylabel('$\| f - f_n \|_2^2 / \| f \|_2^2$','Interpreter','latex');
    end
    % horizontal legend below the third row
    hL = legend( ...
        {'Daubechies45','Meyer','Curvelets','Shearlets','Wave atoms','Boostlets'}, ...
        'Interpreter','latex', ...
        'FontSize',25, ...
        'Orientation','horizontal', ...
        'Units','normalized', ...
        'Position',[0.17, 0.01, 0.80, 0.05] ...
        );
    legend boxoff;

    plot_count = plot_count + 1;
end

% store Figure7
print(gcf, 'Figure7.eps', '-depsc2', '-r600');

% Method and scenario names
method_names   = {'Daubechies45', 'Meyer', 'Curvelets', 'Shearlets', 'Wave atoms', 'Boostlets'};
scenario_names = {'Field (a)', 'Field (b)', 'Field (c)', 'Field (d)'};

% Create tables
l1_table = array2table(l1_norms, 'RowNames', scenario_names, 'VariableNames', method_names);
l2_table = array2table(100*l2_norms, 'RowNames', scenario_names, 'VariableNames', method_names);

% Display
disp('L1-norms (top 10 000 coeffs):');
disp(l1_table);
disp('RMSRE in % (top 1 000 coeffs):');
disp(l2_table);

% remove dependencies...
rmpath(genpath('dependencies'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                END CODE...                %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALLBACK FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RMSRE = computeRMSRE(y,y_rec)
RMSRE = ( norm(y(:) - y_rec(:), 2) )^2  / norm(y(:),2)^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y_rec_boost = iffbt(BT,phi)
    y_rec_boost = real(sum(ifft2(ifftshift(BT.*conj(phi))),3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BT,phi] = ffbt(y,L)
    Nk = size(y,1);
    Nw = size(y,2);
    y_hat = fftshift(fft2(y));
    cc = 1; % redundancy counter (cc = 1 means scaling function coeffs)
    % set dilation and boost levels
    a_grid = 2.^(linspace(0,L-1,L));
    theta_grid = linspace(-pi/2,pi/2,7); % hard-wired, should change if L changes!
    % first, scaling function
    far_or_near = 0;
    boost_type = 1;
    a_j = L-1;
    phi(:,:,cc) = genBoostlet(Nk,Nw,a_j,0,far_or_near,boost_type);
    BT(:,:,cc) = reshape(y_hat.*phi(:,:,cc),Nk,Nw,1);
    % then, decompose through boostlets
    far_or_near = [0,1];
    boost_type = 2;
    for fnfn = far_or_near
        for aa = 1:length(a_grid)
            for thth = 1:length(theta_grid)
                % generate boostlet functions
                a_j = a_grid(aa); % dilation level
                theta_j = theta_grid(thth); % boost level
                % generate boostlet function, compute (cc+1)-th boostlet coefficient and store
                cc = cc + 1;
                phi(:,:,cc) = genBoostlet(Nk,Nw,a_j,theta_j,fnfn,boost_type);
                BT(:,:,cc) = reshape(y_hat.*phi(:,:,cc),Nk,Nw,1);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [phi,KX,OM] = genBoostlet(Nk,Nw,a_j,theta_j,far_or_near,boost_type)
% create Cartesian wavenumber-frequency space
kx = linspace(-1,1,Nk);
om = linspace(-1,1,Nw);
[KX,OM] = meshgrid(kx,om);

% preallocate boosted/dilated points
KX_atheta = zeros(size(KX));
OM_atheta = zeros(size(OM));

% define boost/dilation matrix
M_a_theta = [ a_j*cosh(theta_j) -a_j*sinh(theta_j);
             -a_j*sinh(theta_j)  a_j*cosh(theta_j)];

% first boost and dilate:
for ii = 1:Nk
    for jj = 1:Nw
        boosted_points = M_a_theta*[KX(ii,jj); OM(ii,jj)];
        KX_atheta(ii,jj) = boosted_points(1);
        OM_atheta(ii,jj) = boosted_points(2);
    end
end

% apply diffeo to boosted/dilated points
[Ad,Th] = computeDiffeo(OM_atheta,KX_atheta,far_or_near);

switch boost_type
    case 1 % scaling function
        phi = MeyerScalingFun(Ad);
    case 2 % boostlet functions
        PHI_1 = MeyerWaveletFun(Ad);
        PHI_2 = MeyerScalingFun(Th);
        phi = PHI_1.*PHI_2;
end
phi = reshape(phi,Nk,Nw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ad,Th] = computeDiffeo(OM,KX,far_or_near)
if ~far_or_near
    % move Cartesian grid to diffeomorphic (wavelet) space
    Ad = sqrt(OM.^2-KX.^2); Ad = Ad(:);
    Th = atanh(KX./OM); Th = Th(:);
else
    % move Cartesian grid to diffeomorphic (wavelet) space
    Ad = sqrt(KX.^2-OM.^2); Ad = Ad(:);
    Th = atanh(OM./KX); Th = Th(:);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi_1 = MeyerWaveletFun(x)
xB1 = 1/3;
xB2 = 2*xB1;
xB3 = 4*xB1;

int1 = find((abs(x) >= xB1) & (abs(x) < xB2)); 
int2 = find((abs(x) >= xB2) & (abs(x) < xB3));

psi_1 = zeros(numel(x),1);
psi_1(int1) = sin( pi/2*meyeraux( abs(x(int1))/xB1-1 ) );
psi_1(int2) = cos( pi/2*meyeraux( abs(x(int2))/xB2-1 ) );
psi_1 = reshape(psi_1,size(x,1),size(x,2)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi_2 = MeyerScalingFun(a)
aB1 = 1/6;
aB2 = 2*aB1;

int1 = find(abs(a) < aB1); 
int2 = find((abs(a) >= aB1) & (abs(a) < aB2));

psi_2 = zeros(numel(a),1);
psi_2(int1) = ones(size(int1));
psi_2(int2) = cos( pi/2*meyeraux( abs(a(int2))/aB1-1 ) );
psi_2 = reshape(psi_2,size(a,1),size(a,2)); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signalOut = lpbp1D( signalIn, Ntot, order )
%LPBP1D - 1D Linear Predictive Border Padding
%   This function extrapolates linear microphone array data by means of
%   designing and applying filters with AR prediction coefficients. 
%
%   INPUTS: 
%       signalIn  - 1D signal in a vector
%       Ntot      - no. of total samples (after extrapolation)
%       order     - order of the filter ("no. of Fourier peaks")
%
%   OBS: The difference in samples between Ntot and the length of 
%   signalIn MUST be an even number!
%
%   OUTPUTS:
%       signalOut - extrapolated signal
%
%   EXAMPLE: 
%   Apply a 4-th order LPBP extrapolation to the next power-of-two samples
%   >> p_ext = lpbp1D( p, 2^(nextpow(length(p)), 4 );
%
% Based on the reference: 
%   [1] R. Scholte et al. Truncated aperture extrapolation for Fourier-based 
%       near-field acoustic holography by means of border-padding, JASA 125(6) 
%       pp. 3844-3854 (2009)
%
% Code: E. Zea
% Version: 001
% Date: 2018-10-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L  = length(signalIn); 
dN = 0.5*( Ntot - L ); % difference in samples between in and out

if mod(dN,2) ~= 0
    error('The extrapolated number of samples must be even!');
end

signalOut = zeros( Ntot, 1 ); % preallocate output signal
signalOut(dN+1:Ntot-dN,:) = signalIn; % allocate input signal in the middle

signal_extrapolated_flip = flipud(signalOut); % needs flip first!
    
% backwards extrapolation
a_back = arburg(signal_extrapolated_flip(dN+1:Ntot-dN,1),order);
Z_back = filtic(1,a_back,signal_extrapolated_flip(Ntot-dN-(0:(order-1)),1));
signalOut(1:dN,1) = fliplr(filter(1,a_back,zeros(1,dN),Z_back));

% forward extrapolation
a_forw = arburg(signalOut(dN+1:Ntot-dN,1),order);
Z_forw = filtic(1,a_forw,signalOut(Ntot-dN-(0:(order-1)),1));
signalOut(Ntot-dN+1:Ntot,1) = filter(1,a_forw,zeros(1,dN),Z_forw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flatCoeffs,bandList] = flattenCurveletCoeffs(CT)
% Flattens a CurveLab fdct_wrapping output into a vector + a map
flatCoeffs = [];
bandList = {};
for i = 1:length(CT)
    for j = 1:length(CT{i})
        bandList{end+1} = [i,j];
        block = CT{i}{j}(:);
        flatCoeffs = [flatCoeffs; block];
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CTrec = reconstructCurveletCoeffs(flatCoeffsRec,bandList,CTtemplate)
% Rebuilds a CurveLab cell‐array from a flattened vector
CTrec = CTtemplate;
pos = 1;
for b = 1:numel(bandList)
    i = bandList{b}(1);
    j = bandList{b}(2);
    sz = size(CTtemplate{i}{j});
    n  = prod(sz);
    CTrec{i}{j} = reshape(flatCoeffsRec(pos:pos+n-1), sz);
    pos = pos + n;
end
end
