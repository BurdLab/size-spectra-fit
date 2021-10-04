% 
% Matlab script to provide estimates for the slopes of particle data
%
% NOTES:
% This will also be set up as a Matlab Livescript version that will include
% documentation etc.
%
% HISTORY:
%   11-06-20: Refactor from proof of concept code (ABB)
%   25-01-21: Added Monte Carlo error calculations (ABB)
%
% AUTHORS:
%   Adrian Burd, University of Georgia, Athens, GA, USA.

% Clear workspace to avoid any hangovers from previous runs
close all
clear variables

% Define power law values for POC as a function of particle diameter, based
% on Alldredge. 
poc_param(1) = 5.4106e-5;  % POC vs d prefactor
poc_param(2) = 1.404;      % POC vs d exponent

% Define size class boundaries. The only free parameter here is the upper
% size limit of the large particles. Vary this to see if it makes a
% difference to the calculated slopes.
d_s_lo = 1;       % Lower diameter limit of small particles [microns]
d_s_hi = 51;      % Upper diameter limit of small particles [microns]
d_l_lo = 51;      % Lower diameter limit of large particles [microns]
d_l_hi = 10000;   % Upper diameter limit of large particles [microns]

size_ranges = [d_s_lo d_s_hi d_l_lo d_l_hi];

% Read in the data file
data_file_name = '../GA03_GP16_GN01_LSF_SSF_SPM_Lam_std_updated_20210720.csv';
data           = readmatrix(data_file_name);
n_data         = length(data);

% Separate out the small and large particle data
large_ptcle_data = data(:,7);
large_ptcle_err  = data(:,8);
small_ptcle_data = data(:,9);
small_ptcle_err  = data(:,10);

ptcle_data = [large_ptcle_data large_ptcle_err small_ptcle_data small_ptcle_err];
indx_nan   = isnan(ptcle_data);
nan_sum    = sum(indx_nan,2);

% Set up storage: column 1 = prefactor, column 2 = prefactor err, 3 =
% slope, 4 = slope err
poc_power_law_params = zeros(n_data, 4);
n_power_law_params   = zeros(n_data, 4);
 
% Loop through the data checking for NaNs. If there are none, then fit
% those data to a single power law and save the slopes
my_options = optimoptions('fsolve', 'Display', 'final-detailed', 'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, 'FunctionTolerance', 1.0e-15, 'OptimalityTolerance', 1.0e-12, ...
    'StepTolerance', 1.0e-10);
%my_options = optimoptions('fsolve', 'Display', 'final-detailed', ...
%    'SpecifyObjectiveGradient', true);

x0         = [1.0; 1.2];    % Initial guesses
 
exitflag_store = zeros(n_data,1);

n_trials = 1000;
 
for i_data = 1 : n_data
 
    disp(['Data point ' num2str(i_data) ' out of ' num2str(n_data)])
    
    if nan_sum(i_data) < 1
        
        poc_prefac = zeros(n_trials,1);
        poc_slopes = zeros(n_trials,1);
        n_prefac   = zeros(n_trials,1);
        n_slopes   = zeros(n_trials,1);
        
        mc_obs = zeros(n_trials,2);
        
        for i_trial = 1 : n_trials
           
            obsdata = [small_ptcle_data(i_data) + small_ptcle_err(i_data)*randn(1), ...
                        large_ptcle_data(i_data) + large_ptcle_err(i_data)*randn(1)];
             
            x = fsolve(@(x)TestFunc1(x, obsdata, size_ranges), x0, my_options);
            
            poc_prefac(i_trial) = x(1);
            poc_slopes(i_trial) = x(2);

            [x, fval, exitflag] = fsolve(@(x)TestFunc2(x, obsdata, poc_param, size_ranges), x0, my_options);
%            exitflag
        
            if exitflag >= 0
                n_prefac(i_trial) = x(1);
                n_slopes(i_trial) = x(2);
            else
                n_prefac(i_trial) = NaN;
                n_slopes(i_trial) = NaN;
            end

            mc_obs(i_trial,:) = obsdata;
            
        end
        
        if (i_data == 280)
            small_particles_obs = [small_ptcle_data(280) small_ptcle_err(280)];
            large_particles_obs = [large_ptcle_data(280) large_ptcle_err(280)];
            save('data_280.mat', 'poc_prefac', 'poc_slopes', 'n_prefac', 'n_slopes', ...
                'small_particles_obs', 'large_particles_obs', 'mc_obs');
        end
        if (i_data == 120)
            small_particles_obs = [small_ptcle_data(120) small_ptcle_err(120)];
            large_particles_obs = [large_ptcle_data(120) large_ptcle_err(120)];
            save('data_120.mat', 'poc_prefac', 'poc_slopes', 'n_prefac', 'n_slopes', ...
                'small_particles_obs', 'large_particles_obs', 'mc_obs');
        end
        if (i_data == 850)
            small_particles_obs = [small_ptcle_data(850) small_ptcle_err(850)];
            large_particles_obs = [large_ptcle_data(850) large_ptcle_err(850)];
            save('data_850.mat', 'poc_prefac', 'poc_slopes', 'n_prefac', 'n_slopes', ...
                'small_particles_obs', 'large_particles_obs', 'mc_obs');
        end
        if  (i_data == 561)
            small_particles_obs = [small_ptcle_data(561) small_ptcle_err(561)];
            large_particles_obs = [large_ptcle_data(561) large_ptcle_err(561)];
            save('data_561.mat', 'poc_prefac', 'poc_slopes', 'n_prefac', 'n_slopes', ...
                'small_particles_obs', 'large_particles_obs', 'mc_obs');
        end

        if  (i_data == 600)
            small_particles_obs = [small_ptcle_data(600) small_ptcle_err(600)];
            large_particles_obs = [large_ptcle_data(600) large_ptcle_err(600)];
            save('data_600.mat', 'poc_prefac', 'poc_slopes', 'n_prefac', 'n_slopes', ...
                'small_particles_obs', 'large_particles_obs', 'mc_obs');
        end

        if  (i_data == 620)
            small_particles_obs = [small_ptcle_data(620) small_ptcle_err(620)];
            large_particles_obs = [large_ptcle_data(620) large_ptcle_err(620)];
            save('data_620.mat', 'poc_prefac', 'poc_slopes', 'n_prefac', 'n_slopes', ...
                'small_particles_obs', 'large_particles_obs', 'mc_obs');
        end

        [poc_prefac, poc_slopes, n_prefac, n_slopes] = ...
            CleanModelResults(poc_prefac, poc_slopes, n_prefac, n_slopes);
       
        poc_power_law_params(i_data,1) = mean(poc_prefac,'omitnan');
        poc_power_law_params(i_data,2) = std(poc_prefac,'omitnan');
        poc_power_law_params(i_data,3) = mean(poc_slopes,'omitnan');
        poc_power_law_params(i_data,4) = std(poc_slopes,'omitnan');
        n_power_law_params(i_data,1) = mean(n_prefac,'omitnan');
        n_power_law_params(i_data,2) = std(n_prefac,'omitnan');
        n_power_law_params(i_data,3) = mean(n_slopes,'omitnan');
        n_power_law_params(i_data,4) = std(n_slopes,'omitnan');
    
    
    end
end

%
% Plotting
%

figure(1)
subplot(2,1,1)
errorbar(large_ptcle_data, large_ptcle_err, 'o')
xlabel('Data Number')
ylabel('Large SPM [\mu{}g L^{-1}]')
title('Large (> 51 \mu{}m) Particle Concentration Data')

subplot(2,1,2)
errorbar(small_ptcle_data, small_ptcle_err, 'o')
xlabel('Data Number')
ylabel('Small SPM [\mu{}g L^{-1}]')
title('Small (> 51 \mu{}m) Particle Concentration Data')

orient tall
print('Large_Particle_Data.pdf', '-dpdf', '-fillpage')

figure(2)
subplot(2,1,1)
errorbar(poc_power_law_params(:,1), poc_power_law_params(:,2), 'o')
xlabel('Data Number')
ylabel('POC Power Law Prefactor')
title('POC Power Law Prefactor')

subplot(2,1,2)
errorbar(poc_power_law_params(:,3), poc_power_law_params(:,4), 'o')
xlabel('Data Number')
ylabel('POC Power Law Slope')
title('POC Power Law Slope')

orient tall
print('POC_Fit_Parameters.pdf', '-dpdf', '-fillpage')

figure(3)
subplot(2,1,1)
errorbar(n_power_law_params(:,1), n_power_law_params(:,2), 'o')
xlabel('Data Number')
ylabel('N Power Law Prefactor')
title('N Power Law Prefactor')

subplot(2,1,2)
errorbar(n_power_law_params(:,3), n_power_law_params(:,4), 'o')
xlabel('Data Number')
ylabel('N Power Law Slope')
title('N Power Law Slope')

orient tall
print('N_Fit_Parameters.pdf', '-dpdf', '-fillpage')


% Save results
output_data    = [data poc_power_law_params n_power_law_params];
output_headers = ["Cruise" "Station", "Lat", "Long", "Depth", "Gtnum_all", "POC_Large", ...
    "POC_Large_err", "POC_Small", "POC_Small_err",...
    "p", "p_err", "q", "q_err", "a", "a_err", "b", "b_err"];
filename = 'SpectraResults.csv';

[fid, msg] = fopen(filename, 'wt');
if fid < 0
    error('Could not open file "%s" because "%s%"', fid, msg)
end

fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', output_headers');

tmp = output_data';

fprintf(fid, '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', tmp(:));
fclose(fid)

save("run5000")


% 
% % Sanity Checks
% eps_poc_1 = zeros(n_data, 2);
% eps_poc_2 = zeros(n_data, 2);
% 
% for i_data = 1 : n_data
%     
%     if (~isnan(large_ptcle_data(i_data)) && ~isnan(small_ptcle_data(i_data)))
%         
%         obsdata = [small_ptcle_data(i_data), large_ptcle_data(i_data)];
%     
%         tmp = TestFunc1(poc_power_law_params(i_data,:), obsdata, size_ranges);
%         eps_poc_1(i_data, :) = tmp'./obsdata;
%         
%         if exitflag_store(i_data) > 0
%             tmp = TestFunc2(n_power_law_params(i_data,:), obsdata, poc_param, size_ranges);
%             eps_poc_2(i_data, :) = tmp'./obsdata;
%         else
%             eps_poc_2(i_data, :) = NaN;
%         end
%         
%     end
%     
% 
% 
% end
% eps_poc_1(isnan(large_ptcle_data),:) = NaN;
% eps_poc_1(isnan(small_ptcle_data),:) = NaN;
% 
% % Plot results of sanity check
% figure(1)
% subplot(2,1,1)
% plot([1 : n_data]', eps_poc_1(:, 1), 'r', [1 : n_data]', eps_poc_1(:,2), 'b')
% xlabel('Data number')
% ylabel('Relative Error in POC')
% title('Relative Error for POC spectrum')
% subplot(2,1,2)
% plot([1 : n_data]', eps_poc_2(:, 1), 'r', [1 : n_data]', eps_poc_2(:,2), 'b')
% xlabel('Data number')
% ylabel('Relative Error in POC')
% title('Relative Error for number spectrum')
% 
% % Save results
% output_data    = [data poc_power_law_params n_power_law_params];
% output_headers = ["Cruise" "Station", "Lat", "Long", "Depth", "POC Large", "POC Small", ...
%     "p", "q", "a", "b"];
% filename = 'SpectraResults.csv';
% 
% [fid, msg] = fopen(filename, 'wt');
% if fid < 0
%     error('Could not open file "%s" because "%s%"', fid, msg)
% end
% 
% fprintf(fid, '%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', output_headers');
% 
% tmp = output_data';
% 
% fprintf(fid, '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n', tmp(:));
% fclose(fid)
