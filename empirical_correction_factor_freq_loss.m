% empirical_correction_factor_freq_loss.m
% PL 17.05.2006
%
% See Aubinet et al. "The Euroflux methodology" for details.
%
% Calculate the transfer function for the tube + IRGA response empirically from
% observed data. The tube + IRGA act as a low pass filter, reducing the HF
% components of the observed signal.
% This should be calculated over a sunny period of several hours.
%
tic;
% Example with data from N2O TGA
%
% DAT_FILE = 'H:\TGA_Data_2002_2005\May_2004\05311114.DAT';   % sunny day, June 1 2004
% DAT_FILE = 'H:\TGA_Data_2002_2005\May_2004\05281625.DAT';   % sunny day, June 1 2004
% DAT_FILE = 'H:\TGA_Data_2002_2005\June_2004\06131350.DAT';  % sunny day, June 14t 2004
% DAT_FILE = 'F:\data\N2O_closed\10Hz_n2o\TGA_Data\April_2003\04231103.DAT'; % good day, identified by algorithm. 23/4/2003
DAT_FILE = 'F:\data\N2O_closed\10Hz_n2o\TGA_Data\July 2004\07281203.DAT'; % good day, identified by algorithm. 29/7/2004
HDR_FILE = strrep(DAT_FILE,'.DAT','.HDR');
MET_DATA_LOC = 'F:\data\Donoughmore\30min_co2\raw_data\';
[tga_period, fid, fields] = simple_read_tga_data(DAT_FILE);
% 11:00 to 14:30 on 29/5/2004 is sunny continuously

%% Check that stationarity and fetch criteria satisfied.

dn_start=datenum([2004 7 29 11 30 00]); % actual start of good data
%dn_end=datenum([2003 4 23 22 59 00]);


[t1vec  t_end] = n2o_get_start_end_times(DAT_FILE); % read data start time from hdr file
t1=datenum(t1vec);
%t1=datenum([2004 6 13 13 50 44]);
timebase=t1+(0:(size(tga_period,1)-1))./24./60./60./10;

% restrict range to 2^17 values = approx 3.5 hours

range = find(timebase > dn_start);
range = range(1:2^17);
timebase=timebase(range);       % truncate teh timebase
dn_end=timebase(end);
tga_period=tga_period(range,:); % truncate the tga data series

%% get the corresponding met values
met_dn=datevec(timebase(1));
met_year=met_dn(1);
met_2digit_yr_str=num2str(mod(met_year,1000),'%02.0f'); % e.g. generate '04' from 2004
met_filename=[MET_DATA_LOC,num2str(met_year),'\T0_1_',num2str(is_leap_year(met_year)+365),'_',met_2digit_yr_str,'.dat'];
[f m s]=read_cork_data_new(met_filename);
metrange=intersect(find(m.datenum>dn_start), find(m.datenum<=dn_end));

%% need to do all the corrections, rotations etc.
% this section taken fro m n2o_process_raw.m
% [4] co-ordinate rotation
% 4.1 rotate into prevailing wind direction (make v=0)
n_periods=1;
n_cols=34;
i=1;
results=ones(n_cols,n_periods).*NaN;
[u_rot,v_rot,angle_yaw] = rotate_(tga_period(:,5),tga_period(:,6)); % slava rotation
% 4.2 v should now be 0. rotate pitch:
[u_rot,w_rot,angle_pitch] = rotate_(u_rot, tga_period(:,7));
% check that v is close to zero:
meanv = nanmean(v_rot);
if (abs(meanv) >  10^(-8))
    warning(['<v> after yaw rotation was ', num2str(meanv)]);
end
% despike the n2o data. remove all points 3 stdev's outside mean
[cn_despiked n_spike cn_bar_orig cn_std_orig]=despike_series(tga_period(:,1),3);
tga_period(:,1)=cn_despiked;
%% despike the sonic temp
[Ts_despiked n_Ts_spike Ts_bar_orig Ts_std_orig]=despike_series(tga_period(:,8),3);
tga_period(:,8)=Ts_despiked;
results(33,i)=n_Ts_spike;

% record N2O conc (before detrending):
results(27,i)=nanmean(tga_period(:,1));

% remove mean N2O conc
tga_period(:,1)=tga_period(:,1)-results(27,i);

% form the covariance input matrix: cn | u | v | w | Ts
cov_inputs = [tga_period(:,1).' ; u_rot.' ; v_rot.' ; w_rot.' ; Ts_despiked.'].';

% need to mask the NaNs before cov: cn/u_rot/v_rot/w_rot/T
non_nan=find(~isnan(sum(cov_inputs,2)));

% only perform lag correl & cov() if there are sufficient non-NaN elements:
% [5] lag correlation of N2O conc and w:
% & [6] covariance calculations: and application of lag correction
n_subints=64; % 64 subints = 2^11 / 2^7 (approx 5 minute subints)
[cov_results_temp peak_lag peak_strength cn_w_corrs subint_covs]=...
    n2o_calc_lag_covariance(cov_inputs(non_nan,:), 1, 1, 35, 9,n_subints );
% calc_lag: lag calculation *ON*
% min_lag: 1
% max_lag: 35
% default_lag 9
% n_subints =1 ( subintervals )

%% apply the calculated lag correction to the RAW cn data
tga_period(:,1)=circshift(tga_period(:,1), -1.*peak_lag); % minus => sshift N2O FORWARD in time
tga_period(end-peak_lag+1:end,1)=NaN; % remove "wrapped around" values
results(2,i)=cov_results_temp(1,2); %cn/u cov
results(3,i)=cov_results_temp(1,3); %cn/v cov
results(4,i)=cov_results_temp(1,4); %cn/w cov
results(5,i)=cov_results_temp(2,3); %u/v cov
results(6,i)=cov_results_temp(2,4); %u/w cov
results(7,i)=cov_results_temp(3,4); %v/w cov
results(29,i)=cov_results_temp(5,1); % Ts/cn cov
results(30,i)=cov_results_temp(5,2); % Ts/u cov
results(31,i)=cov_results_temp(5,3); % Ts/v cov
results(32,i)=cov_results_temp(5,4); % Ts/w cov

% std devs are the square roots of the elements of the cov matrix on the
% diagonal: [cn u v w Ts]
cov_diag = (diag(cov_results_temp)).^(0.5);
results(8,i)=cov_diag(1);
results(9,i)=cov_diag(2);
results(10,i)=cov_diag(3);
results(11,i)=cov_diag(4);
results(28,i)=cov_diag(5); % Ts std
% compute u*:
%results(17,i)=sqrt(-1.*results(6,i)); % Formula replaced. uw may
%                                           be positive sometimes.
results(17,i)=( (results(6,i)).^2 + (results(7,i)).^2 )^(0.25); % (uw^2 + vw^2)^(1./4)
results(18,i)=cov_results_temp(4,5); % w/T cov
% convert <wc> to flux:
% c is in [ppm_V] ; w in [m.s-1]
results(20,i)=1e3.*1820.*results(4,i);
% 1820 g.m-3 = density of N2O at std temp, press. from Todd. add T dependence?
% 1e3 = product of ppm (1e-6) and g->ng (1e9) conversion.
% results are therefore in ng N2O m-2 s-1.
%
results(23,i)=peak_lag;
% convert peak_strength to a correlation:
results(26,i)=peak_strength./cov_diag(1)./cov_diag(4);
% <<>> PL: this block moved here from outside else 07.12.2005 :
results(12,i)=nanmean(tga_period(:,1));
results(13,i)=nanmean(u_rot);
results(14,i)=nanmean(v_rot);
results(15,i)=nanmean(w_rot);
results(16,i)=nanmean(tga_period(:,8));
% obukhov length: (without water vapour term)
results(19,i)= -1.*((results(17,i)).^3.*(results(16,i)+273.15))./ ...
    (0.4.*9.81.*results(18,i));

%  add misc. processing info to the results:
results(21,i)=angle_yaw;
results(22,i)=angle_pitch;

results(24,i)=numel(non_nan);
results(25,i)=n_spike; % record number of skpikes

% end of section taken from n2o_process_raw.m
%%% -------------------------------------------

%% simple normalisation factor calculation: ratio of covs
non_nan2=intersect(non_nan,find(~isnan(tga_period(:,1))));
cc_wn=cov(tga_period(non_nan2,1),cov_inputs(non_nan2,4)); % use shifted cn series
cc_wn=cc_wn(2,1);
cc_wT = cov_results_temp(4,5); %wT
NT_NS_1=cc_wT./cc_wn;

%% simple normalisation factor calcultion: ratio of sdevs
NT_NS_2 = cov_diag(5)./cov_diag(1);

%% calculate the normalisation factors
% Aubinet equation [34]
%M_max=2.^3; % maximum freq scale
tic; % cn | u | v | w | Ts
[D_w_T wrbar Tbar wr Tr] = multires_decomp(cov_inputs(:,4),cov_inputs(:,5)); 
toc;
tic;
[D_w_cn cnbar Tbar cnr Tr] = multires_decomp(cov_inputs(:,4),tga_period(:,1));
toc;

%%
for j=1:17
    nt_ns_3(j)=sum(D_w_T(1:j))./sum(D_w_cn(1:j));
end
figure(3);
set(gcf, 'name', 'ratio of cospectra (unnormalised)');
plot(nt_ns_3, 'k-s');
xlabel({'m (mode) \rightarrow', 'f'});
grid on;

%% integration:
%delta_f=2.^(17:-1:1);
delta_f=10.*(2.^(1:1:17))./2.^17; 
% ... Progressively higher pass filter, from m=0 (just leaves all zeros)
% to m=M-1 (single mean removed, preserves 10 Hz).
figure(4);
set(gcf,'name','normalisation factor calculation');
semilogx(delta_f,cumsum(nt_ns_3.*delta_f),'m-^');
grid on;
xlabel('f (Hz)');
ylabel('ratio of integrated cospectra <wT>/<ws>');

%% calculate transfer function.
NT_NS_assumed=5000; % seems reasonable from above estimates
TF_exp=nt_ns_3./NT_NS_assumed;
TF_exp(find(TF_exp<0))=0;
figure(5);
clf;
set(gcf,'name','transfer function fit');
loglog(delta_f,TF_exp(1:17),'k-s');
beta_fit=nlinfit(delta_f,TF_exp(1:17),'exp_lp_filter_tf',0.1);
TF_fit=exp_lp_filter_tf(beta_fit,delta_f);
hold on;
plot(delta_f,TF_fit,'r:o');
legend({'experimental TF',['Fitted TF fco=',num2str(beta_fit)]});
grid on;

%% calculate correction factor
corr_factor=sum((D_w_cn(1:17)).'.*delta_f)./sum(TF_fit.*delta_f);

save corr_factor_02102007.mat D_w_T wrbar Tbar wr Tr D_w_cn cc_wT cc_wn NT_NS_1 NT_NS_2 nt_ns_3 DAT_FILE peak_lag delta_f beta_fit corr_factor;


limits = 10:5000:numel(range);
corr_factor = ones(size(limits)).*NaN;
for i_limit=1:numel(limits)
    N_t = cov(tga(range(1:limits(i_limit)), 7), tga(range(1:limits(i_limit)) ,8)); % <w'T'>
    N_s = cov(tga(range(1:limits(i_limit)), 7), tga(range(1:limits(i_limit)), 1)); % <w'c'>
    corr_factor(i_limit) = (N_t(2,1))./(N_s(2,1));
end

%% plotting results
figure(1);
clf;

% plot w
subplot(2,3,1);
plot(timebase(range),tga(range,7));
ylabel('w [ms-1]');
grid on;
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,15));
xlabel(['time, ',datestr(xt(1),24)]);

% plot T
subplot(2,3,2);
plot(timebase(range),tga(range,8));
ylabel('T [deg C]');
grid on;
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,15));
xlabel(['time, ',datestr(xt(1),24)]);

% plot cN2O
subplot(2,3,3);
plot(timebase(range),tga(range,1));
ylabel('N2O conc. [ppm_V]');
grid on;
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,15));
xlabel(['time, ',datestr(xt(1),24)]);


% plot met (Rswin or PAR)
subplot(2,3,4);
plot(m.datenum(metrange),m.Rswin.data(metrange),'b-o');
ylabel(['R_{SW,in} [',m.Rswin.unit,']']);
grid on;
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,15));
xlabel(['time, ',datestr(xt(1),24)]);


% plot corr factor
subplot(2,3,5);
plot(limits,corr_factor,'k-o');
xlabel('integration limit [# samples]');
ylabel('correction factor [unitless]');
grid on;

% plot TF
subplot(2,3,6);
plot(limits,corr_factor(end).*N_s./N_t,'k-o');
xlabel('integration limit [# samples]');
ylabel('transfer function [unitless]');
grid on;

%% finish
mark;
toc;
