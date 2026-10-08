%% Beamspace_PARAMING_M3_MSE_vs_SNR.m
% Extended beamspace-PARAMING demonstration for M = 3 targets.
% This is part of the work done in iSEE-6G. 
% UNIVERSITY OF PIRAEUS RESEARCH CENTER
% CORRESDPONDING AUTHOR: Harris K. Armeniakos, harmen@unipi.gr



% Uses all 64 x 64 beam pairs.
% Beamspace channel model:
%
%   H_B[n] = sum_{m=1}^{M} alpha_m exp(-j2*pi*n*Delta_f*tau_m)
%            f_R(Omega_R,m) f_T(Omega_T,m)^H + noise
%
% The algorithm:
%   1) Generate beamspace CSI for M targets.
%   2) Stack all beam pairs over subcarriers.
%   3) Build a frequency-Hankel matrix.
%   4) Estimate delays using the PARAMING/ESPRIT eigenvalue step.
%   5) Remove delays by LS.
%   6) Estimate 3D AoA/AoD by known-pattern matching.
%   7) Estimate complex gains.
%   8) Repeat over SNR and Monte Carlo trials.
%   9) Plot MSE vs SNR curves.
%
% Important:
%   This is not antenna-domain PARAMING.
%   The spatial phase-regression step is replaced by beam-pattern matching.

clear; 

%% ===============================================================
%  1) System parameters
% ===============================================================

Br = 64;                 % RX beams
Bt = 64;                 % TX beams
Nb = Br * Bt;            % all beam pairs = 4096

NP = 64;                 % number of OFDM subcarriers
Delta_f = 960e3;         % subcarrier spacing [Hz]
Tmax = 1 / Delta_f;      % maximum unambiguous delay

M = 3;                   % number of targets/paths

% Frequency-Hankel parameters
MP = 8;                  % frequency window length
KN = NP - MP + 1;

if KN <= M
    error('Increase NP or reduce MP. Need KN > M.');
end

%% ===============================================================
%  2) 64-beam 3D grids
% ===============================================================

% 64 beams = 16 azimuth beams x 4 elevation beams.
% You may change this layout according to your real codebook.

Nphi_beams   = 16;
Ntheta_beams = 4;

phi_min_deg   = -60;
phi_max_deg   =  60;
theta_min_deg = -30;
theta_max_deg =  30;

rxBeamDirs = create_beam_grid(Nphi_beams, Ntheta_beams, ...
    phi_min_deg, phi_max_deg, theta_min_deg, theta_max_deg);

txBeamDirs = create_beam_grid(Nphi_beams, Ntheta_beams, ...
    phi_min_deg, phi_max_deg, theta_min_deg, theta_max_deg);

if size(rxBeamDirs,1) ~= Br || size(txBeamDirs,1) ~= Bt
    error('Beam-grid size does not match Br and Bt.');
end

% Synthetic beam-pattern widths.
% These define the smoothness of the known beam pattern.
phi_spacing_deg   = (phi_max_deg - phi_min_deg)/(Nphi_beams-1);
theta_spacing_deg = (theta_max_deg - theta_min_deg)/(Ntheta_beams-1);

sigma_phi_deg   = 0.75 * phi_spacing_deg;
sigma_theta_deg = 0.75 * theta_spacing_deg;

%% ===============================================================
%  3) True target parameters
% ===============================================================

% Delays must lie in [0, 1/Delta_f) to avoid ambiguity.
tau_true = [1.03; 1.022; 1.04] * 1e-6;

if any(tau_true >= Tmax)
    error('Some true delays exceed the unambiguous delay 1/Delta_f.');
end

% RX AoA angles
phi_R_true_deg   = [-35.38;  12.21;  43.89];
theta_R_true_deg = [-18.32;   5.56;  21.22];

% TX AoD angles
phi_T_true_deg   = [-28.55;  24.87;  50.65];
theta_T_true_deg = [ 19.11; -22.32;  24.51];

% Complex gains
alpha_true = [ ...
    1.00 * exp(1j*0.30); ...
    0.80 * exp(1j*1.10); ...
    0.65 * exp(1j*2.00) ];

trueParams.tau    = tau_true;
trueParams.phiR   = phi_R_true_deg;
trueParams.thetaR = theta_R_true_deg;
trueParams.phiT   = phi_T_true_deg;
trueParams.thetaT = theta_T_true_deg;
trueParams.alpha  = alpha_true;

%% ===============================================================
%  4) Angular search dictionaries
% ===============================================================

% Coarser grid gives faster execution but larger angular MSE floor.
% Use 1 degree first. Then reduce to 0.5 degree if needed.

phi_search_deg   = phi_min_deg:1:phi_max_deg;
theta_search_deg = theta_min_deg:1:theta_max_deg;

[FR_dict, angleGrid_R] = build_pattern_dictionary(rxBeamDirs, ...
    phi_search_deg, theta_search_deg, sigma_phi_deg, sigma_theta_deg);

[FT_dict, angleGrid_T] = build_pattern_dictionary(txBeamDirs, ...
    phi_search_deg, theta_search_deg, sigma_phi_deg, sigma_theta_deg);

%% ===============================================================
%  5) Generate clean beamspace channel for the M targets
% ===============================================================

[Hbar_clean, trueAtoms] = generate_clean_beamspace_channel( ...
    trueParams, rxBeamDirs, txBeamDirs, ...
    sigma_phi_deg, sigma_theta_deg, NP, Delta_f, Br, Bt);

% Hbar_clean size:
%   Br*Bt x NP = 4096 x 64

%% ===============================================================
%  6) Single-SNR sanity check
% ===============================================================

SNR_check_dB = 0;

Hbar_noisy = add_awgn_to_channel(Hbar_clean, SNR_check_dB);

estCheck = estimate_beamspace_paraming( ...
    Hbar_noisy, Br, Bt, MP, M, Delta_f, ...
    FR_dict, angleGrid_R, FT_dict, angleGrid_T);

idxBest = best_permutation(estCheck, trueParams);

estCheck = reorder_estimates(estCheck, idxBest);

fprintf('\n================ SINGLE-SNR SANITY CHECK ================\n');
fprintf('SNR = %.1f dB\n', SNR_check_dB);

for m = 1:M
    fprintf('\n--- Target %d ---\n', m);

    fprintf('True tau      = %.4e s\n', trueParams.tau(m));
    fprintf('Estimated tau = %.4e s\n', estCheck.tau(m));
    fprintf('Delay error   = %.4e s\n', estCheck.tau(m)-trueParams.tau(m));

    fprintf('True RX phi      = %.2f deg\n', trueParams.phiR(m));
    fprintf('Estimated RX phi = %.2f deg\n', estCheck.phiR(m));
    fprintf('RX phi error     = %.2f deg\n', ...
        angle_diff_deg(estCheck.phiR(m), trueParams.phiR(m)));

    fprintf('True RX theta      = %.2f deg\n', trueParams.thetaR(m));
    fprintf('Estimated RX theta = %.2f deg\n', estCheck.thetaR(m));
    fprintf('RX theta error     = %.2f deg\n', ...
        estCheck.thetaR(m)-trueParams.thetaR(m));

    fprintf('True TX phi      = %.2f deg\n', trueParams.phiT(m));
    fprintf('Estimated TX phi = %.2f deg\n', estCheck.phiT(m));
    fprintf('TX phi error     = %.2f deg\n', ...
        angle_diff_deg(estCheck.phiT(m), trueParams.phiT(m)));

    fprintf('True TX theta      = %.2f deg\n', trueParams.thetaT(m));
    fprintf('Estimated TX theta = %.2f deg\n', estCheck.thetaT(m));
    fprintf('TX theta error     = %.2f deg\n', ...
        estCheck.thetaT(m)-trueParams.thetaT(m));

    fprintf('True |alpha|      = %.4f\n', abs(trueParams.alpha(m)));
    fprintf('Estimated |alpha| = %.4f\n', abs(estCheck.alpha(m)));
end

Hbar_est_check = reconstruct_beamspace_channel(estCheck, NP, Delta_f, Nb);
NMSE_check = norm(Hbar_clean - Hbar_est_check, 'fro')^2 / ...
             norm(Hbar_clean, 'fro')^2;

fprintf('\nNMSE at %.1f dB = %.4e\n', SNR_check_dB, NMSE_check);
fprintf('=========================================================\n');

%% ===============================================================
%  7) Monte Carlo MSE vs SNR
% ===============================================================

% SNR_dB_vec = -20:5:30;
% 
% Nmc = 50;       % Increase to 100 or 200 for smoother paper-quality curves.
% 
% delayMSE = zeros(size(SNR_dB_vec));
% angleMSE = zeros(size(SNR_dB_vec));
% gainMSE  = zeros(size(SNR_dB_vec));
% chanNMSE = zeros(size(SNR_dB_vec));
% 
% fprintf('\n================ MONTE CARLO SIMULATION ================\n');
% fprintf('M = %d targets, Br = Bt = %d, all beam pairs = %d\n', M, Br, Nb);
% fprintf('Monte Carlo trials per SNR = %d\n\n', Nmc);
% 
% for iSNR = 1:numel(SNR_dB_vec)
% 
%     SNR_dB = SNR_dB_vec(iSNR);
% 
%     delayErrAccum = 0;
%     angleErrAccum = 0;
%     gainErrAccum  = 0;
%     nmseAccum     = 0;
% 
%     fprintf('Running SNR = %4.1f dB ... ', SNR_dB);
% 
%     for imc = 1:Nmc
% 
%         Hbar_noisy = add_awgn_to_channel(Hbar_clean, SNR_dB);
% 
%         est = estimate_beamspace_paraming( ...
%             Hbar_noisy, Br, Bt, MP, M, Delta_f, ...
%             FR_dict, angleGrid_R, FT_dict, angleGrid_T);
% 
%         % Resolve permutation ambiguity.
%         idxBest = best_permutation(est, trueParams);
%         est = reorder_estimates(est, idxBest);
% 
%         % Delay MSE
%         delayErr = est.tau - trueParams.tau;
%         delayErrAccum = delayErrAccum + mean(abs(delayErr).^2);
% 
%         % Angular MSE in deg^2
%         dPhiR   = angle_diff_deg(est.phiR, trueParams.phiR);
%         dThetaR = est.thetaR - trueParams.thetaR;
% 
%         dPhiT   = angle_diff_deg(est.phiT, trueParams.phiT);
%         dThetaT = est.thetaT - trueParams.thetaT;
% 
%         angleErrors = [dPhiR(:); dThetaR(:); dPhiT(:); dThetaT(:)];
%         angleErrAccum = angleErrAccum + mean(angleErrors.^2);
% 
%         % Gain MSE
%         gainErrAccum = gainErrAccum + mean(abs(est.alpha - trueParams.alpha).^2);
% 
%         % Channel reconstruction NMSE
%         Hbar_est = reconstruct_beamspace_channel(est, NP, Delta_f, Nb);
% 
%         nmseAccum = nmseAccum + ...
%             norm(Hbar_clean - Hbar_est, 'fro')^2 / ...
%             norm(Hbar_clean, 'fro')^2;
%     end
% 
%     delayMSE(iSNR) = delayErrAccum / Nmc;
%     angleMSE(iSNR) = angleErrAccum / Nmc;
%     gainMSE(iSNR)  = gainErrAccum  / Nmc;
%     chanNMSE(iSNR) = nmseAccum     / Nmc;
% 
%     fprintf('done.\n');
% end
% 
% fprintf('========================================================\n');
% 
% %% ===============================================================
% %  8) Plot MSE vs SNR curves
% % ===============================================================
% 
% figure;
% semilogy(SNR_dB_vec, delayMSE, '-o', 'LineWidth', 1.6);
% grid on;
% xlabel('SNR [dB]');
% ylabel('Delay MSE [s^2]');
% title('Beamspace-PARAMING: Delay MSE vs SNR');
% 
% figure;
% semilogy(SNR_dB_vec, angleMSE, '-s', 'LineWidth', 1.6);
% grid on;
% xlabel('SNR [dB]');
% ylabel('Angular MSE [deg^2]');
% title('Beamspace-PARAMING: Angular MSE vs SNR');
% 
% figure;
% semilogy(SNR_dB_vec, gainMSE, '-^', 'LineWidth', 1.6);
% grid on;
% xlabel('SNR [dB]');
% ylabel('Gain MSE');
% title('Beamspace-PARAMING: Complex Gain MSE vs SNR');
% 
% figure;
% semilogy(SNR_dB_vec, chanNMSE, '-d', 'LineWidth', 1.6);
% grid on;
% xlabel('SNR [dB]');
% ylabel('Channel NMSE');
% title('Beamspace-PARAMING: Channel NMSE vs SNR');

% %% ===============================================================
% %  9) Plot beam-pair energy for visual inspection
% % ===============================================================
% 
% Pbeam = sum(abs(Hbar_clean).^2, 2);
% PbeamMat = reshape(Pbeam, Br, Bt);
% 
% figure;
% imagesc(10*log10(PbeamMat / max(PbeamMat(:))));
% colorbar;
% xlabel('TX beam index');
% ylabel('RX beam index');
% title('Clean normalized beam-pair energy [dB]');
% axis xy;

%% =================================================================
%  Local functions
% =================================================================

function beamDirs = create_beam_grid(Nphi, Ntheta, ...
                                     phiMin, phiMax, thetaMin, thetaMax)

    phiGrid   = linspace(phiMin, phiMax, Nphi);
    thetaGrid = linspace(thetaMin, thetaMax, Ntheta);

    beamDirs = zeros(Nphi*Ntheta, 2);

    idx = 0;
    for itheta = 1:Ntheta
        for iphi = 1:Nphi
            idx = idx + 1;
            beamDirs(idx,1) = phiGrid(iphi);
            beamDirs(idx,2) = thetaGrid(itheta);
        end
    end
end


function f = beam_pattern_response(beamDirs, phiDeg, thetaDeg, ...
                                   sigmaPhiDeg, sigmaThetaDeg)
% Synthetic complex beam pattern F.
%
% Replace this function with your actual measured/simulated complex pattern:
%
%   f(q) = F_q(phiDeg, thetaDeg), q = 1,...,B.
%
% For example, if your measured pattern table is F_meas(q, iphi, itheta),
% this function should interpolate that table at (phiDeg, thetaDeg).

    phiCenters   = beamDirs(:,1);
    thetaCenters = beamDirs(:,2);

    dphi   = angle_diff_deg(phiDeg, phiCenters);
    dtheta = thetaDeg - thetaCenters;

    % Smooth main-lobe-like magnitude response
    amp = exp(-0.5*((dphi./sigmaPhiDeg).^2 + ...
                    (dtheta./sigmaThetaDeg).^2));

    % Deterministic complex phase pattern
    phase = exp(1j*2*pi*( ...
        0.30*sind(phiDeg).*cosd(phiCenters) + ...
        0.20*sind(thetaDeg).*cosd(thetaCenters)));

    f = amp .* phase;

    % Normalize to remove arbitrary pattern gain scaling
    if norm(f) > 0
        f = f / norm(f);
    end
end


function [Fdict, angleGrid] = build_pattern_dictionary(beamDirs, ...
                                                       phiSearch, thetaSearch, ...
                                                       sigmaPhiDeg, sigmaThetaDeg)

    Nphi   = numel(phiSearch);
    Ntheta = numel(thetaSearch);

    B = size(beamDirs,1);
    Nang = Nphi*Ntheta;

    Fdict = zeros(B, Nang);
    angleGrid = zeros(Nang, 2);

    idx = 0;
    for itheta = 1:Ntheta
        for iphi = 1:Nphi
            idx = idx + 1;

            phiDeg   = phiSearch(iphi);
            thetaDeg = thetaSearch(itheta);

            angleGrid(idx,:) = [phiDeg, thetaDeg];

            Fdict(:,idx) = beam_pattern_response(beamDirs, ...
                phiDeg, thetaDeg, sigmaPhiDeg, sigmaThetaDeg);
        end
    end

    colNorms = vecnorm(Fdict,2,1);
    colNorms(colNorms == 0) = 1;
    Fdict = Fdict ./ colNorms;
end


function [Hbar_clean, atoms] = generate_clean_beamspace_channel( ...
    trueParams, rxBeamDirs, txBeamDirs, ...
    sigmaPhiDeg, sigmaThetaDeg, NP, Delta_f, Br, Bt)

    M = numel(trueParams.tau);
    Nb = Br * Bt;

    Hbar_clean = zeros(Nb, NP);

    atoms.fR = zeros(Br, M);
    atoms.fT = zeros(Bt, M);
    atoms.b  = zeros(Nb, M);

    for m = 1:M

        fR = beam_pattern_response(rxBeamDirs, ...
            trueParams.phiR(m), trueParams.thetaR(m), ...
            sigmaPhiDeg, sigmaThetaDeg);

        fT = beam_pattern_response(txBeamDirs, ...
            trueParams.phiT(m), trueParams.thetaT(m), ...
            sigmaPhiDeg, sigmaThetaDeg);

        b = kron(conj(fT), fR);

        atoms.fR(:,m) = fR;
        atoms.fT(:,m) = fT;
        atoms.b(:,m)  = b;

        for nIdx = 1:NP
            n = nIdx - 1;

            delayPhase = exp(-1j*2*pi*n*Delta_f*trueParams.tau(m));

            Hbar_clean(:,nIdx) = Hbar_clean(:,nIdx) + ...
                trueParams.alpha(m) * delayPhase * b;
        end
    end
end


function Hbar_noisy = add_awgn_to_channel(Hbar_clean, SNR_dB)

    signalPower = mean(abs(Hbar_clean(:)).^2);
    noisePower  = signalPower / (10^(SNR_dB/10));

    
    noise = sqrt(noisePower/2) * ...
        (randn(size(Hbar_clean)) + 1j*randn(size(Hbar_clean)));

    Hbar_noisy = Hbar_clean + noise;
end


function est = estimate_beamspace_paraming( ...
    Hbar, Br, Bt, MP, M_est, Delta_f, ...
    FR_dict, angleGrid_R, FT_dict, angleGrid_T)

    [Nb, NP] = size(Hbar);

    if Nb ~= Br*Bt
        error('Hbar row dimension must be Br*Bt.');
    end

    KN = NP - MP + 1;

    %% Build frequency-Hankel matrix
    Hhankel = build_frequency_hankel(Hbar, MP);

    H1 = Hhankel(:,1:KN-1);
    H2 = Hhankel(:,2:KN);

    %% Rank-M SVD
    [U,S,V] = svd(H1, 'econ');

    U_M = U(:,1:M_est);
    S_M = S(1:M_est,1:M_est);
    V_M = V(:,1:M_est);

    %% Eigenvalue step
    Tmat = S_M \ (U_M' * H2 * V_M);

    gamma_hat = eig(Tmat);

    tau_hat = -angle(gamma_hat) / (2*pi*Delta_f);
    tau_hat = mod(real(tau_hat), 1/Delta_f);
    tau_hat = tau_hat(:);

    %% Delay removal
    Ctau_hat = build_delay_matrix(NP, Delta_f, tau_hat);

    Yhat = Hbar * conj(Ctau_hat) / (Ctau_hat.' * conj(Ctau_hat));

    %% Angle and gain estimation
    phiR_hat   = zeros(M_est,1);
    thetaR_hat = zeros(M_est,1);

    phiT_hat   = zeros(M_est,1);
    thetaT_hat = zeros(M_est,1);

    alpha_hat = zeros(M_est,1);
    b_hat_all = zeros(Nb, M_est);

    for m = 1:M_est

        y_m = Yhat(:,m);
        Ymat_m = reshape(y_m, Br, Bt);

        % Since Ymat_m ≈ alpha_m fR_m fT_m^H,
        % its dominant left/right singular vectors approximate fR_m and fT_m.
        [Uy,~,Vy] = svd(Ymat_m, 'econ');

        u1 = Uy(:,1);
        v1 = Vy(:,1);

        % RX pattern matching
        scoreR = abs(FR_dict' * u1).^2;
        [~,idxR] = max(scoreR);

        phiR_hat(m)   = angleGrid_R(idxR,1);
        thetaR_hat(m) = angleGrid_R(idxR,2);

        % TX pattern matching
        scoreT = abs(FT_dict' * v1).^2;
        [~,idxT] = max(scoreT);

        phiT_hat(m)   = angleGrid_T(idxT,1);
        thetaT_hat(m) = angleGrid_T(idxT,2);

        fR_hat = FR_dict(:,idxR);
        fT_hat = FT_dict(:,idxT);

        b_hat = kron(conj(fT_hat), fR_hat);

        alpha_hat(m) = (b_hat' * y_m) / (b_hat' * b_hat);

        b_hat_all(:,m) = b_hat;
    end

    est.tau    = tau_hat;
    est.phiR   = phiR_hat;
    est.thetaR = thetaR_hat;
    est.phiT   = phiT_hat;
    est.thetaT = thetaT_hat;
    est.alpha  = alpha_hat;
    est.bHat   = b_hat_all;
end


function Hhankel = build_frequency_hankel(Hbar, MP)

    [Nb, NP] = size(Hbar);

    KN = NP - MP + 1;

    if KN <= 1
        error('MP must satisfy 1 < MP < NP.');
    end

    Hhankel = zeros(MP*Nb, KN);

    for q = 1:MP
        rowRange = (q-1)*Nb + (1:Nb);
        colRange = q:(q+KN-1);

        Hhankel(rowRange,:) = Hbar(:,colRange);
    end
end


function Ctau = build_delay_matrix(NP, Delta_f, tauVec)

    tauVec = tauVec(:);
    M = numel(tauVec);

    n = (0:NP-1).';

    Ctau = zeros(NP, M);

    for m = 1:M
        Ctau(:,m) = exp(-1j*2*pi*n*Delta_f*tauVec(m));
    end
end


function idxBest = best_permutation(est, trueParams)
% Finds the permutation of estimated targets that best matches the true ones.
% For M=3 this brute-force search is simple and robust.

    M = numel(trueParams.tau);

    P = perms(1:M);

    bestCost = inf;
    idxBest = 1:M;

    tauScale = 1e-6;
    angScale = 10;

    for ip = 1:size(P,1)

        idx = P(ip,:);

        dtau = (est.tau(idx) - trueParams.tau) / tauScale;

        dphiR = angle_diff_deg(est.phiR(idx), trueParams.phiR) / angScale;
        dtheR = (est.thetaR(idx) - trueParams.thetaR) / angScale;

        dphiT = angle_diff_deg(est.phiT(idx), trueParams.phiT) / angScale;
        dtheT = (est.thetaT(idx) - trueParams.thetaT) / angScale;

        cost = sum(abs(dtau).^2) + ...
               sum(abs(dphiR).^2) + sum(abs(dtheR).^2) + ...
               sum(abs(dphiT).^2) + sum(abs(dtheT).^2);

        if cost < bestCost
            bestCost = cost;
            idxBest = idx;
        end
    end
end


function estOut = reorder_estimates(estIn, idx)

    estOut = estIn;

    estOut.tau    = estIn.tau(idx);
    estOut.phiR   = estIn.phiR(idx);
    estOut.thetaR = estIn.thetaR(idx);
    estOut.phiT   = estIn.phiT(idx);
    estOut.thetaT = estIn.thetaT(idx);
    estOut.alpha  = estIn.alpha(idx);
    estOut.bHat   = estIn.bHat(:,idx);
end


function Hbar_est = reconstruct_beamspace_channel(est, NP, Delta_f, Nb)

    M = numel(est.tau);

    Hbar_est = zeros(Nb, NP);

    for m = 1:M

        for nIdx = 1:NP
            n = nIdx - 1;

            delayPhase = exp(-1j*2*pi*n*Delta_f*est.tau(m));

            Hbar_est(:,nIdx) = Hbar_est(:,nIdx) + ...
                est.alpha(m) * delayPhase * est.bHat(:,m);
        end
    end
end


function d = angle_diff_deg(a, b)
% Wrapped difference a-b in degrees in [-180,180).

    d = mod((a - b) + 180, 360) - 180;
end


