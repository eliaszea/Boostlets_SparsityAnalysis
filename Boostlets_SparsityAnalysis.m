%BOOSTLETS_SPARSITYANALYSIS Code to reproduce sparsity analysis
% -----------------------------------------------------------------------------------
%       This code can be used to reproduce the sparsity analysis in the paper:
%   E. Zea, M. Laudato, J. Andén, "A boostlet transform for wave-based acoustic signal 
%   processing in space--time", arXiv:2403.11362, DOI: 10.48550/arXiv.2403.11362, 2025.
% -----------------------------------------------------------------------------------
% EXAMPLE OF USAGE: 
%   To reproduce the sparsity analysis in Figure 6 of the paper, run the command: 
%
%       > Boostlets_SparsityAnalysis;
%
%   The output produces a figure (Figure 1) corresponding to Figure 6 in the paper. 
%   A table including the l1-norm of the 20.0000 largest coefficients, corresponding 
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
% - Version 002 [July 20, 2025]: Renamed repo
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
S = 2; % integer between [2,4] — otherwise, some systems do not work!

% sampling parameters and 2D space-time coordinates
N = 100;
N_max = 2e4; % limit to 20.000 largest coefficients
T = N;    % no. time samples
M = N;    % no. mic positions
dx = 0.03;   % mic spacing (m)
fs = 11250;   % sampling frequency (Hz)
x = linspace(-M*dx/2,M*dx/2,M);    % x-axis coordinates

% Preallocate arrays to store l1-norm results
l1_norms = zeros(length(T_start), 6); % 6 decomposition methods

% Preallocate boostlet coefficients
BT = [];

plot_count = 1;
for splt = 1:length(T_start)
    % load acoustic field
    load(['dependencies/Datasets/' room{splt} 'RIR.mat'])
    y = out.image(1+T_start(splt):T+T_start(splt),:);
    
    % take 2D fft
    Y_hat = fftshift(fftshift(fft2(y),1),2)/(M*T);

    %% COMPUTE DECOMPOSITIONS AND SORT INTO DECAYING AMPLITUDE

    % Debauchies45 wavelet coefficients and sort in descending amplitude
    [WT, ~] = wavedec2(y, S, 'db45');   
    wavelet_coeffs_db = sort(abs(WT(:)),'descend');
    selected_wavelet_coeffs_db = wavelet_coeffs_db(1:N_max);
    selected_wavelet_coeffs_db = selected_wavelet_coeffs_db/sqrt(sum(selected_wavelet_coeffs_db(:).^2));

    % Meyer wavelet coefficients and sort in descending amplitude
    [WT2, ~] = wavedec2(y, S, 'dmey');  
    wavelet_coeffs_meyer = sort(abs(WT2(:)),'descend');
    selected_wavelet_coeffs_meyer = wavelet_coeffs_meyer(1:N_max);
    selected_wavelet_coeffs_meyer = selected_wavelet_coeffs_meyer/sqrt(sum(selected_wavelet_coeffs_meyer(:).^2));

    % wave atom coefficients and sort in descending amplitude
    % (OBS! needs LPBP extrapolation to 2^n)
    pat = 'p';
    tp = 'directional';
    T_ext = 128;
    y_ext = zeros(T_ext,T_ext);
    for ii = 1:T
        y_ext(ii,:) = lpbp1D(y(ii,:),T_ext,16);
    end
    c = fwa2(y_ext,pat,tp);
    ns = size(c,1);
    cc = 1;
    for is=1:ns
        for dir = [1,2]
            nw = size(c{is,dir},1);
            nb = size(c{is,dir}{nw,nw},1);
            for ia = 1:ceil(nw/2)
                for ib = 1:ceil(nw/4)
                    c{is,2}{ia,ib}(nb/2+1,nb/2+1) = 1;
                    waveatom = iwa2(c,pat,tp);
                    wave_hat = fftshift(fft2(waveatom));
                    WAT(:,:,cc) = reshape(fftshift(fft2(y_ext)).*wave_hat,T_ext,T_ext,1);
                    cc = cc + 1;
                end
            end
        end
    end
    waveats_coeffs = sort(abs(WAT(:)),'descend');
    selected_waveats_coeffs = waveats_coeffs(1:N_max);
    selected_waveats_coeffs = selected_waveats_coeffs/sqrt(sum(selected_waveats_coeffs(:).^2));

    % boostlet coefficients and sort in descending amplitude
    cc = 1; % redundancy counter (cc = 1 means scaling function coeffs)
    % set dilation and boost levels
    a_grid = 2.^(linspace(0,S-1,S));
    theta_grid = linspace(-pi/2,pi/2,7); % hard-wired, should change if S changes!
    % first, scaling function
    [phi,KX,OM] = genBoostlet(N,S-1,0,0,1);
    BT(:,:,cc) = reshape(Y_hat.*phi,N,N,1);
    % then, decompose through boostlets
    far_or_near = [0,1];
    boost_type = 2;
    for fnfn = far_or_near
        for aa = 1:length(a_grid)
            for thth = 1:length(theta_grid)
                % generate boostlet function
                [phi,KX,OM] = genBoostlet(N,a_grid(aa),theta_grid(thth),fnfn,boost_type);
                % compute (cc+1)-th boostlet coefficient and store
                cc = cc + 1;
                BT(:,:,cc) = reshape(Y_hat.*phi,N,N,1);
            end
        end
    end
    boostlet_coeffs = sort(abs(BT(:)),'descend');
    selected_boostlet_coeffs = boostlet_coeffs(1:N_max);
    selected_boostlet_coeffs = selected_boostlet_coeffs/sqrt(sum(selected_boostlet_coeffs(:).^2));

    % compute curvelet coefficients and sort in descending amplitude
    CT = fdct_wrapping(y, 0, 2, S+1);
    coeff_vector = [];
    % Iterate over each element of CT
    for i = 1:length(CT)
        % Iterate over each element of the nested cell array
        for j = 1:length(CT{i})
            % Convert the nested cell array into a vector and concatenate it with coeff_vector
            coeff_vector = [coeff_vector CT{i}{j}(:)'];
        end
    end
    curvelet_coeffs = sort(abs(coeff_vector(:)),'descend');
    selected_curvelet_coeffs = curvelet_coeffs(1:N_max);
    selected_curvelet_coeffs = selected_curvelet_coeffs/sqrt(sum(selected_curvelet_coeffs(:).^2));

    % compute shearlet coefficients and sort in descending amplitude
    ST = shearletTransformSpect(y,S);
    shearlet_coeffs = sort(abs(ST(:)),'descend');
    selected_shearlet_coeffs = shearlet_coeffs(1:N_max);
    selected_shearlet_coeffs = selected_shearlet_coeffs/sqrt(sum(selected_shearlet_coeffs(:).^2));

    % compute l1 norms
    l1_db = norm(selected_wavelet_coeffs_db(:),1);
    l1_meyer = norm(selected_wavelet_coeffs_meyer(:),1);
    l1_wavat = norm(selected_waveats_coeffs(:),1);
    l1_clets = norm(selected_curvelet_coeffs(:),1);
    l1_slets = norm(selected_shearlet_coeffs(:),1);
    l1_blets = norm(selected_boostlet_coeffs(:),1);
    % Store results in the array
    l1_norms(splt, :) = [l1_db, l1_meyer, l1_wavat, l1_clets, l1_slets, l1_blets];

    % plot fields and coefficient decays
    t = 1e3*linspace(0+T_start(splt)/fs,(T+T_start(splt))/fs,T); % in ms
    [XX,TT] = meshgrid(x,t);
    figure(1);
    % subplot: acoustic field in space-time
    subplot(2,length(T_start),plot_count);
    surf(XX,TT,y,'edgecolor','none');
    view(0,90);
    if splt == 3
        xlabel('$x$ (m)','Interpreter','latex');
    end
    xlabel('$x$ (m)','Interpreter','latex');
    ylabel('$t$ (ms)','Interpreter','latex');
    colormap gray; axis tight;
    cb = colorbar('NorthOutside'); 
    ttl = title(cb,'[Pa]','Interpreter', 'latex', 'Position', [-30 25 0]);
    cb.TickLabelInterpreter = 'latex';
    set(gca,'fontsize',25,'TickLabelInterpreter', 'latex');
    text(0.02, 0.95, ['(' char(96 + splt) ')'], 'Color', 'white', ...
        'FontSize', 30, 'FontWeight', 'bold', 'Units', 'normalized');
    % subplot: coefficient decay
    subplot(2,length(T_start),length(T_start)+plot_count);
    fig = loglog(1:N_max,selected_wavelet_coeffs_db,'x:', ...
                 1:N_max,selected_wavelet_coeffs_meyer,':', ...
                 1:N_max,selected_waveats_coeffs,'+--', ...
                 1:N_max,selected_shearlet_coeffs,'-.', ...
                 1:N_max,selected_curvelet_coeffs,'*--', ...
                 1:N_max,selected_boostlet_coeffs,'-');
    fig(1).LineWidth = 2; fig(1).MarkerSize = 10; fig(1).Color = [0.15,0.82,0.79];
    fig(2).LineWidth = 4; fig(2).Color = [0,0.45,0.74];
    fig(3).LineWidth = 1.5; fig(3).MarkerSize = 8; fig(3).Color = [0,0,1]; 
    fig(4).LineWidth = 3; fig(4).Color = [0.53,0.48,1]; 
    fig(5).LineWidth = 3; fig(5).MarkerSize = 2; fig(5).Color = [0.59,0.12,0.75]; 
    fig(6).LineWidth = 5; fig(6).Color = [0.89,0.07,0.52];
    set(gca,'fontsize',25,'TickLabelInterpreter', 'latex');
    axis tight;
    text(0.02, 0.95, ['(' char(100 + splt) ')'], 'Color', 'black', ...
        'FontSize', 30, 'FontWeight', 'bold', 'Units', 'normalized');
    ylim([5e-5,5e-1]);
    xticks([1,10,100,1000,10000])
    yticks([1e-4,1e-3,1e-2,1e-1])
    if splt == 3 
        xlabel('No. coefficients','Interpreter','latex');
    end
    xlabel('No. coefficients','Interpreter','latex');
    if splt == 1
        ylabel('Coeff. magnitude','Interpreter','latex');
        legend('Daubechies45','Meyer','Wave atoms','Shearlets','Curvelets','Boostlets', ...
               'Orientation','vertical','Location','best','NumColumns',1,'Interpreter', 'latex','fontsize',20);
        legend boxoff;
    end
    plot_count = plot_count + 1;
end

% Create table with l1 norms
method_names = {'Daubechies45', 'Meyer', 'Wave atoms', 'Curvelets', 'Shearlets', 'Boostlets'};
scenario_names = {'Field (a)', 'Field (b)', 'Field (c)', 'Field (d)'};
results_table = array2table(l1_norms, 'RowNames', scenario_names, 'VariableNames', method_names);

% Display the table
disp('L1 Norms:');
disp(results_table);

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
function [phi,KX,OM] = genBoostlet(N,a_j,theta_j,far_or_near,boost_type)
% create Cartesian wavenumber-frequency space
om = linspace(-1,1,N);
kx = linspace(-1,1,N);
[KX,OM] = meshgrid(kx,om);

% preallocate boosted/dilated points
KX_atheta = zeros(size(KX));
OM_atheta = zeros(size(OM));

% define boost/dilation matrix
M_a_theta = [ a_j*cosh(theta_j) -a_j*sinh(theta_j);
             -a_j*sinh(theta_j)  a_j*cosh(theta_j)];

% first boost and dilate:
for ii = 1:N
    for jj = 1:N
        boosted_points = M_a_theta*[KX(ii,jj); OM(ii,jj)];
        KX_atheta(ii,jj) = boosted_points(1);
        OM_atheta(ii,jj) = boosted_points(2);
    end
end

% apply diffeo to boosted/dilated points
[Ad,Th] = computeDiffeo(OM_atheta,KX_atheta,far_or_near);

switch boost_type
    case 1 % scaling function
        phi = MeyerScalingFun(N,Ad);
    case 2 % boostlet functions
        PHI_1 = MeyerWaveletFun(N,Ad);
        PHI_2 = MeyerScalingFun(N,Th);
        phi = PHI_1.*PHI_2;
end
phi_max = max(abs(phi(:)));
phi = phi/phi_max;
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
function psi_1 = MeyerWaveletFun(N,a)
aB1 = 1/3;
aB2 = 2*aB1;
aB3 = 4*aB1;

int1 = find((abs(a) >= aB1) & (abs(a) < aB2)); 
int2 = find((abs(a) >= aB2) & (abs(a) < aB3));

psi_1 = zeros(numel(a),1);
psi_1(int1) = 1/sqrt(2*pi)*exp( 1i*a(int1)/2 ).*...
                         sin( pi/2*meyeraux( abs(a(int1))/aB1-1 ) );
psi_1(int2) = 1/sqrt(2*pi)*exp( 1i*a(int2)/2 ).*...
                         cos( pi/2*meyeraux( abs(a(int2))/aB2-1 ) );
psi_1 = reshape(psi_1,N,N); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi_2 = MeyerScalingFun(N,theta)
thetaB1 = 1/6;
thetaB2 = 2*thetaB1;

int1 = find(abs(theta) < thetaB1); 
int2 = find((abs(theta) >= thetaB1) & (abs(theta) < thetaB2));

psi_2 = zeros(numel(theta),1);
psi_2(int1) = 1/sqrt(2*pi)*ones(size(int1));
psi_2(int2) = 1/sqrt(2*pi)*cos( pi/2*meyeraux( abs(theta(int2))/thetaB1-1 ) );
psi_2 = reshape(psi_2,N,N); 
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
