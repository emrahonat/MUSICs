%% High-Resolution Direction Finding on Uniform Linear Arrays: 
%% A Comparative Evaluation of MUSIC Algorithm Derivatives 
%
% 19.08.2025
% Dr. Emrah Onat
% 
% You can use my MATLAB codes as long as you cite my article below.
%
% Onat, E. (2025). High-Resolution Direction Finding on Uniform Linear Arrays: A Comparative Evaluation of MUSIC Algorithm Derivatives. Gazi Üniversitesi Fen Bilimleri Dergisi Part C: Tasarım ve Teknoloji, 13(3), 1122-1136. https://doi.org/10.29109/gujsc.1696029
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clc
close all
clear all

%angles = [30]/180*pi; % DOA
w = [pi/4]'; % % Frequency of incident wavefronts
M = 8; % Number of array elements
D = length(w); % Number of incident wavefronts
lambda = 100; % Wavelength of incident wavefronts
%spacing = lambda/2; % Spacing between array elements of ULA
% spacing = 10; % Spacing between array elements of ULA
theta = -90:1:90; % Range of thetas to simulate for
% theta = 0:1:180; % Range of thetas to simulate for
num_trials = 10000; % Number of trials for each snapshot value

M_values = 8:2:32;
spacing_values = [10 50];%, 100/9, 100/8, 100/7, 100/6, 20, 25, 100/3, 50];
snapshot_values = [100]; % Different snapshot values
SNR_values = -15; %-20:5:30;

execution_times = zeros(1, num_trials);
colors = {'r', 'g', 'b'}; % Colors
sure_tut=1;

grafik = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
grafik_ucaroot = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
grafik_root = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
grafik_fbss = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
grafik_improved = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));

ind = 0;
time_MUSIC = [];
time_ROOT = [];
time_FBSS = [];
time_IMP = [];
time_UCA = [];

RMSE_total_music = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_total_UCA = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_total_ROOT = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_total_fbss = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_total_improved = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));

RMSE_avg_music = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_avg_UCA = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_avg_ROOT = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_avg_fbss = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));
RMSE_avg_improved = zeros(length(M_values), length(spacing_values), length(snapshot_values), length(SNR_values));

for ii = 1:length(M_values)

    M = M_values(ii);

    for spc = 1:length(spacing_values)

        spacing = spacing_values(spc);

        for snap_idx = 1:length(snapshot_values)
           
            snapshot = snapshot_values(snap_idx); % Number of instants of time
            color_index = mod(snap_idx - 1, length(colors)) + 1; % Color index

        
            bos_MUSIC = zeros(length(-20:5:30),1);
            bos_FBSS = zeros(length(-20:5:30),1);
            bos_Root = zeros(length(-20:5:30),1);
            bos_Imp = zeros(length(-20:5:30),1);
            bos_UCA = zeros(length(-20:5:30),1);
        
            for snr_idx = 1:length(SNR_values) % Different SNR values
%                 SNR = -20 + (snr_idx - 1) * 5;
                SNR = SNR_values(snr_idx);
                        
        
                for trial = 1:num_trials
                    ind = ind +1;
                    %% Signal Generation %%
        
                    sinyal = (rand(1) * 2 - 1);
        
                    % Angle Scope
                    min_angle_deg = -80;
                    max_angle_deg = 80;
        
                    % Random angle generation (btw 10 - 80 degree)
                    first_random_angle_deg = (max_angle_deg - min_angle_deg) * rand(1) + min_angle_deg;
%                     disp(['degree:', num2str(first_random_angle_deg)]);

                    % Degree to radian
                    first_random_angle_rad = deg2rad(first_random_angle_deg);
                                
                    %angles = first_random_angle;
                    angles=first_random_angle_rad;
                    w = 3e8/lambda; 
        
                    A = zeros(D,M);
                    for k=1:D
                        A(k,:) = exp(-1i*2*pi*spacing*sin(angles(k))/lambda*(0:M-1)); % Steering/Mode vectors
                    end
                    
                    A = A';
                    F = 2*exp(1j*(w*(1:snapshot))); % Incident signals
                    X = A*F;
                    X = X+awgn(X,SNR); % Observed signal
                    S = X*X'; % Covariance
        
                    %% Traditional MUSIC
                    t1=tic;
                    Pmusic = MuSiC(S,M,D,lambda,spacing,theta); % The function P_{MU}(theta)
                    [pks,locs] = findpeaks(Pmusic,theta,'SortStr','descend','Annotate','extents');
                    t2=toc(t1);
                    time_MUSIC(ind)=t2;
        
                    % locs control
                    if length(locs) >= D
                        estimated_angles = sort(locs(1:D));
                    else
                        estimated_angles = sort(locs);
                    end
                    
                    if isempty(estimated_angles)
                        disp('Empty');
                        bos_MUSIC(snr_idx) = bos_MUSIC(snr_idx) + 1;
                    else
                        angle_derece = angles*180/pi;
                        RMSE = (estimated_angles - angle_derece).^2;
                        RMSE_total_music(ii,spc,snap_idx,snr_idx) = RMSE_total_music(ii,spc,snap_idx,snr_idx) + RMSE;
                    end
                    estimated_angles_MUSIC = estimated_angles;
%                     disp(['MUSIC:', num2str(estimated_angles_MUSIC)]);
             
                    %%  ROOT MUSIC %%
                    t1=tic;
                    doa_estimate = -1*root_music_doa(X, D, spacing/lambda);
                    t2=toc(t1);
                    time_ROOT(ind)=t2;
        
                    if isreal(doa_estimate)
                        angle_derece = angles*180/pi;
                        RMSE = (doa_estimate - angle_derece).^2;
                        RMSE_total_ROOT(ii,spc,snap_idx,snr_idx) = RMSE_total_ROOT(ii,spc,snap_idx,snr_idx) + RMSE;
                    else
                        disp('Subspace imaginary');
                        bos_Root(snr_idx) = bos_Root(snr_idx) + 1;
                    end
                    estimated_angles_ROOT = doa_estimate;
        
                    %% FBSS MUSIC %%
                    t1=tic;            
                    M0 = 5; % Size
                    L_fb = floor(M / M0); % Number
                    Pmusic = FBSS_MuSiC(S, M, L_fb, lambda, theta, D,spacing); % P_{MU}(theta) 
                    [pks,locs] = findpeaks(Pmusic,theta,'SortStr','descend','Annotate','extents');
                    t2=toc(t1);
                    time_FBSS(ind)=t2;
        
                    % locs control
                    if length(locs) >= D
                        estimated_angles = sort(locs(1:D));
                    else
                        estimated_angles = sort(locs);
                    end
                    
                    if isempty(estimated_angles)
                        disp('Empty');
                        bos_FBSS(snr_idx) = bos_FBSS(snr_idx) + 1;
                    else
                        angle_derece = angles*180/pi;
                        RMSE = (estimated_angles - angle_derece).^2;
                        RMSE_total_fbss(ii,spc,snap_idx,snr_idx) = RMSE_total_fbss(ii,spc,snap_idx,snr_idx) + RMSE;
                    end
                    estimated_angles_FBSS = estimated_angles;
%                     disp(['MUSIC FBSS:', num2str(estimated_angles_FBSS)]);
        
                    %% Improved MUSIC %%
                    t1=tic;
                    J=fliplr(eye(M));
                    S=S+J*conj(S)*J;
                    Pmusic = MuSiC(S,M,D,lambda,spacing,theta); % The function P_{MU}(theta)
                    [pks,locs] = findpeaks(Pmusic,theta,'SortStr','descend','Annotate','extents');
                    t2=toc(t1);
                    time_IMP(ind)=t2;
                   
                    % size of locs control
                    if length(locs) >= D
                        estimated_angles = sort(locs(1:D));
                    else
                        estimated_angles = sort(locs);
                    end
                    
        
                    if isempty(estimated_angles)
                        disp('Empty');
                        bos_Imp(snr_idx) = bos_Imp(snr_idx) + 1;
                    else
                        angle_derece = angles*180/pi;
                        RMSE = (estimated_angles - angle_derece).^2;
                        RMSE_total_improved(ii,spc,snap_idx,snr_idx) = RMSE_total_improved(ii,spc,snap_idx,snr_idx) + RMSE;
                    end
                    estimated_angles_Imp = estimated_angles;
%                     disp(['MUSIC Imp:', num2str(estimated_angles_Imp)]);
        
                    %% UCA ROOT MUSIC  %%%%
                   
                    %input data
                    %theta1 = rand(1)*180 * (pi / 180);
                    theta1 = first_random_angle_rad + (pi/2);
        
                    %angle control
                    
                    aci_1 = theta1 / (pi / 180); 
                    %aci_1 = first_random_angle_deg;
            
                    data1 = sign(2 * rand(1, snapshot) - 1) + 1i * sign(2 * rand(1, snapshot) - 1);
                   
                    % Signal to Noise ratio
                    % transmitted signals.
                    s1 = sqrt(10^(SNR / 10)) * data1;
                    
                    t1=tic;
                    % Direction of Arrival (DOA) for three uncorrelated sources
                    % array response vector
                    m = M / 2;
                    i = 0:M;
                    A1 = exp(-1i * (i - m) * (theta1));  
                    
                   
                    % (DFT) Submatrix
                    Fv = zeros(M + 1);
                    for h = 0:M
                        for p = 1:M+1
                            Fv(h + 1, p) = exp((1i * 2 * pi * (h - M / 2) * (p - 1)) / M);
                        end
                    end
                    F1 = (1 / sqrt(M)) * Fv;
                    % Diagonal matrix
                    Jo = zeros(M + 1, 1);
                    for i = 0:M
                        Jo(i + 1) = 1 / (sqrt(M) * 1i^(i - M / 2) * besseli(i - M / 2, 1.6 * 2 * pi));
                    end
                    J = diag(Jo);
                    % the observation vectors from the two subarrays
                    Z = eye(8);
                    u1 = A1' * s1 ; 
                    % add noise to the observation vector
                    n1 = sqrt(0.5) * randn(size(u1)) + 1i * sqrt(0.5) * randn(size(u1));
                    Zr1 = u1 + n1; % the observation data
                    h = (M / 2) + 1;
                    K = M + 1;
                    % Constructing the Toeplitz matrix
                    Ztop = zeros(h, h);
                    for s = 0:(h - 1)
                        Ztop(s + 1, 1:h) = Zr1(h - s:K - s, 1);
                    end
                    % Compute eigendecomposition and order by descending eigenvalues
        
                    %D => L 
                    [Q0, L] = eig(Ztop * Ztop');
                    [lambda1, index] = sort(abs(diag(L)));
                    %lambda -> lambda_uca
                    lambda_uca = lambda1(h:-1:1);
                    Q = Q0(:, index(h:-1:1));
                    % Compute pseudo-spectrum
                    AA = zeros(2 * h - 1, 1);
                    for n = 1:(h - 2)
                        AA = AA + conv(Q(:, h - (n - 1)), conj(Q(h:-1:1, h - (n - 1))));
                    end
                    r_A = roots(AA);
                    [i_min] = find(abs(r_A) < 1);
                    r_A_min = r_A(i_min);
                    freq = abs(angle(r_A_min) / (pi / 180));
            
                    % angle 1
                    aci1_diff = abs(freq - aci_1);
                    [~, min_index] = min(aci1_diff);
                    en_kucuk_aci1 = freq(min_index);
            
                   
                    theta_median_1 = en_kucuk_aci1;
        
                    t2=toc(t1);
                    time_UCA(ind)=t2;
                    
                    if ~isnan(theta_median_1)
                        RMSE_UCA = (theta_median_1 - aci_1)^2;
                        RMSE_total_UCA(ii,spc,snap_idx,snr_idx) =  RMSE_UCA + RMSE_total_UCA(ii,spc,snap_idx,snr_idx);
                    end
                    estimated_angles_UCA = theta_median_1;
                end

                %%Traditional MUSIC%%
                RMSE_total_music(ii,spc,snap_idx,snr_idx)=RMSE_total_music(ii,spc,snap_idx,snr_idx)/(num_trials-bos_MUSIC(snr_idx));
                RMSE_avg_music(ii,spc,snap_idx,snr_idx)=sqrt(RMSE_total_music(ii,spc,snap_idx,snr_idx));
                grafik(ii,spc,snap_idx,snr_idx) = mean(RMSE_avg_music(ii,spc,snap_idx,snr_idx)); 
            
                %%UCA ROOT MUSIC%%
                RMSE_total_UCA(ii,spc,snap_idx,snr_idx) = RMSE_total_UCA(ii,spc,snap_idx,snr_idx) / num_trials;
                RMSE_avg_UCA(ii,spc,snap_idx,snr_idx) = sqrt(RMSE_total_UCA(ii,spc,snap_idx,snr_idx));
                grafik_ucaroot(ii,spc,snap_idx,snr_idx) = mean(RMSE_avg_UCA(ii,spc,snap_idx,snr_idx)); 
            
                %%ROOT MUSIC%%
                RMSE_total_ROOT(ii,spc,snap_idx,snr_idx)=RMSE_total_ROOT(ii,spc,snap_idx,snr_idx)/(num_trials-bos_Root(snr_idx));
                RMSE_avg_ROOT(ii,spc,snap_idx,snr_idx)=sqrt(RMSE_total_ROOT(ii,spc,snap_idx,snr_idx));
                grafik_root(ii,spc,snap_idx,snr_idx) = mean(RMSE_avg_ROOT(ii,spc,snap_idx,snr_idx)); 
                
                %%FBSS MUSIC%%
                RMSE_total_fbss(ii,spc,snap_idx,snr_idx)=RMSE_total_fbss(ii,spc,snap_idx,snr_idx)/(num_trials-bos_FBSS(snr_idx));
                RMSE_avg_fbss(ii,spc,snap_idx,snr_idx)=sqrt(RMSE_total_fbss(ii,spc,snap_idx,snr_idx));
                grafik_fbss(ii,spc,snap_idx,snr_idx) = mean(RMSE_avg_fbss(ii,spc,snap_idx,snr_idx)); 
            
                %%Improved MUSIC%%
                RMSE_total_improved(ii,spc,snap_idx,snr_idx)=RMSE_total_improved(ii,spc,snap_idx,snr_idx)/(num_trials-bos_Imp(snr_idx));
                RMSE_avg_improved(ii,spc,snap_idx,snr_idx)=sqrt(RMSE_total_improved(ii,spc,snap_idx,snr_idx));
                grafik_improved(ii,spc,snap_idx,snr_idx) = mean(RMSE_avg_improved(ii,spc,snap_idx,snr_idx)); 
                        
            end
        end
    end


end

ele_graph = 1;
if ele_graph == 1
    
    rmse_figure = figure;
    
    semilogy(M_values, grafik(:,1,1,1), '-o', 'LineWidth', 1.5, 'Color', 'g');
    hold on;
    semilogy(M_values, grafik_improved(:,1,1,1), '->', 'LineWidth', 1.5, 'Color', 'k');
    hold on;
    semilogy(M_values, grafik_fbss(:,1,1,1), '-+', 'LineWidth', 1.5, 'Color', 'r');
    hold on;
    semilogy(M_values, grafik_root(:,1,1,1), '-x', 'LineWidth', 1.5, 'Color', 'b');
    hold on;
    

    semilogy(M_values, grafik(:,2,1,1), '--o', 'LineWidth', 1.5, 'Color', 'g');
    hold on;
    semilogy(M_values, grafik_improved(:,2,1,1), '-->', 'LineWidth', 1.5, 'Color', 'k');
    hold on;
    semilogy(M_values, grafik_fbss(:,2,1,1), '--+', 'LineWidth', 1.5, 'Color', 'r');
    hold on;
    semilogy(M_values, grafik_root(:,2,1,1), '--x', 'LineWidth', 1.5, 'Color', 'b');
    grid on;
        
    xlabel('Number of Antennas (M)');
    ylabel('RMSE (deg)');
    title(['2A - RMSE Values when d = λ/2 & λ/10; SNR = ' num2str(SNR_values) 'dB ; N = ' num2str(snapshot_values)]);
    legend('Traditional MUSIC (d = λ/10)', 'Improved MUSIC (d = λ/10)', 'FBSS MUSIC (d = λ/10)', 'ROOT MUSIC (d = λ/10)', 'Traditional MUSIC (d = λ/2)', 'Improved MUSIC (d = λ/2)', 'FBSS MUSIC (d = λ/2)', 'ROOT MUSIC (d = λ/2)');
        
end

name = string(date);
saveas(rmse_figure,name+"_rmse_M_d.eps")
saveas(rmse_figure,name+"_rmse_M_d.fig")

%% Figure 11

avg_duration_MUSIC = mean(time_MUSIC);
avg_duration_ROOT = mean(time_ROOT);
avg_duration_FBSS = mean(time_FBSS);
avg_duration_IMP = mean(time_IMP);
avg_duration_UCA = mean(time_UCA);

colormap(1,:) = [1 0 0]; 
colormap(2,:) = [0 1 0];
colormap(3,:) = [0 0 1];
colormap(4,:) = [0 1 1];
colormap(5,:) = [1 0 1];
colormap(6,:) = [1 1 0];
colormap(7,:) = [0 0 0];
colormap(8,:) = [1 1 1];
colormap(9,:) = [0 0.4470 0.7410];
colormap(10,:) = [0.8500 0.3250 0.0980];
colormap(11,:) = [0.9290 0.6940 0.1250];
colormap(12,:) = [0.4940 0.1840 0.5560];
colormap(13,:) = [0.4660 0.6740 0.1880];
colormap(14,:) = [0.3010 0.7450 0.9330];
colormap(15,:) = [0.6350 0.0780 0.1840];

Durations = [round(avg_duration_MUSIC,4);round(avg_duration_ROOT,4);round(avg_duration_FBSS,4);round(avg_duration_IMP,4)];
centers = 1:4;
Fig_time_figure = figure;
colors = [colormap(2,:); colormap(7,:); colormap(1,:); colormap(3,:)];

hBar = bar(centers,Durations,'hist');
set(hBar,'FaceVertexCData',colors);
labels = [string('Traditional MUSIC');string('Root MUSIC');string('FBSS MUSIC');string('Improved MUSIC')];
axis([0 5 0 4e-3])
grid on

x = get(hBar,'XData');
y = get(hBar,'YData');
ygap = 0.0001;  %// Specify vertical gap between the bar and label
ylimits = get(gca,'YLim');

for i = 1:1:length(x) %// Loop over each bar
    xpos = x(2,i)+0.5;        %// Set x position for the text label
    ypos = y(2,i) + ygap; %// Set y position, including gap
    htext = text(xpos,ypos,labels{i});          %// Add text label
    set(htext,'VerticalAlignment','bottom', 'HorizontalAlignment','center')
end

ylabel('Duration (second)');
title('Durations of the Algorithms');

name = string(date);
saveas(Fig_time_figure,name+"_Fig11_durations.eps")
saveas(Fig_time_figure,name+"_Fig11_durations.fig")


