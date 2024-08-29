%% identify_n2o_correction_factor_period.m
% PL 31.05.2007
% Identify a suitable period for calculating the TGA transfer function and
% correction factors (see Aubinet et al. "The Euroflux methodoology").
%
% Period must be ~3 hours and :
% (a) continuously sunny
% (b) stationarity constraints satisfied
% (c) fetch constraints satisfied
%

%% general setup
FLUX_PATH='F:\data\N2O_closed\analysed\5min_mats\';
%FLUX_PATH='H:\TGA_Data_2002_2005\analysed\5min_mats\';
MET_DATA_LOC='F:\data\Donoughmore\30min_co2\raw_data\';

n_x_subplots=2;
n_y_subplots=4;
subplot_handles=cell(n_x_subplots.*n_y_subplots,1);

met_year=2004; % year to examine
rswin_threshold=650; % [W m-2], threshold for "sunny" conditions

%% read the met data
met_2digit_yr_str=num2str(mod(met_year,1000),'%02.0f'); % e.g. generate '04' from 2004
met_filename=[MET_DATA_LOC,num2str(met_year),'\T0_1_',num2str(is_leap_year(met_year)+365),'_',met_2digit_yr_str,'.dat'];
[f m s]=read_cork_data_new(met_filename);

%% plot Rsinw data
figure(1);
clf;
subplot_handles{1}=subplot(n_y_subplots,n_x_subplots,1);
plot(m.datenum,m.Rswin.data,'b-+');
hold on;
grid on;
ylim([-100 1.1.*max(m.Rswin.data)]);

%% scan 3 hr moving window across data to identify periods that satisfy the
% Rswin threshold.
windowsize=6; % 6 x 30 min values
good_intervals=[]; % initialise results array
for i_win=1:(numel(m.Rswin.data)-windowsize+1)
    windowrange=i_win:(i_win+windowsize-1);
    if (sum(m.Rswin.data(windowrange)>=rswin_threshold) == windowsize)
        good_intervals=[good_intervals i_win];
    end
end

%% indicate found intervals on pllot
figure(1);
plot(m.datenum((good_intervals)),m.Rswin.data((good_intervals)),'ro');
xlim([m.datenum(good_intervals(1)) m.datenum(good_intervals(end))]);
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,19));
xlabel(['date, ',num2str(met_year)]);
ylabel({'Rswin','[Wm-2]'});

%% for stationarity, need to look at the 5 minute covariances.
%
pc_difference=zeros(size(good_intervals))+NaN; % // initialise
wind_max=zeros(size(good_intervals));
wind_min=zeros(size(good_intervals));

for i_int=1:numel(good_intervals)
    %for i_int=3:3 % TEST ONLY!
    figure(1);

    curr_datenum=datenum(m.datenum(good_intervals(i_int)));
    curr_datevec=datevec(curr_datenum);
    curr_month_str=lower(datestr(curr_datevec,'mmm'));
    curr_n2o_filename=[FLUX_PATH,'n2o_proc_',curr_month_str,'_',num2str(met_year),'.mat'];
    n2o=load(curr_n2o_filename);
    n2o_range=intersect(find(n2o.all_results(1,:) >= m.datenum(good_intervals(i_int))),...
        find(n2o.all_results(1,:) < m.datenum(windowsize+good_intervals(i_int))) );
    if (numel(n2o_range)>0)

        %% set title of current event
        subplot(n_y_subplots,n_x_subplots,1);
        [year day_of_year]=datenum_2_jd(curr_datenum);
        day_of_year=floor(day_of_year);
        title(['Current time :',datestr(n2o.all_results(1,n2o_range(1))),...
            ' = jd ',num2str(day_of_year)]);


        %% plotting of w'N2O' and w'T' 5 min covs
        subplot_handles{2}=subplot(n_y_subplots,n_x_subplots,2);
        plot(n2o.all_results(1,n2o_range),n2o.all_results(4,n2o_range),'ro');
        hold on;
        plot(n2o.all_results(1,n2o_range),n2o.all_results(18,n2o_range),'bo'); % col 18 = wT
        %% calculate the 3 hour cov by combinign 5min cov's
        [combined_cov combined_mean_w combined_mean_x]=...
            combine_covs( n2o.all_results(4,n2o_range),...
            n2o.all_results(15,n2o_range),...   % mean vertical wind speed
            n2o.all_results(12,n2o_range),...   % mean n2o conc
            n2o.all_results(24,n2o_range));      % number non-nan points
        [combined_wT combined_mean_w combined_mean_T]=...
            combine_covs( n2o.all_results(18,n2o_range),...
            n2o.all_results(15,n2o_range),...   % mean vertical wind speed
            n2o.all_results(16,n2o_range),...   % mean Tsonic
            n2o.all_results(24,n2o_range));      % number non-nan points
        %% plot the combined cov with a line
        line([n2o.all_results(1,n2o_range([1 end]))], [combined_cov combined_cov],...
            'color','r');
        line([n2o.all_results(1,n2o_range([1 end]))], [combined_wT combined_wT],...
            'color','b');
        ylabel('cov []');
        legend({'wN2O','wT','wN2O comb.','wT comb.'},'orientation',...
            'horizontal','location','southwest');
        title('raw covs and combined covs');


        mean_cov=mean(n2o.all_results(4,n2o_range)); % col 4: wN2O
        mean_wT=mean(n2o.all_results(18,n2o_range)); % col 18: wT
        %% plot all the subinterval (5 minute) covs
        %plot(n2o.all_results(1,n2o_range),n2o.all_results(34:39,n2o_range),'k.');

        % percent difference is the absolute difference divided by the
        % smaller of the two values , times 100
        pc_difference(i_int)    = 100.*abs((mean_cov-combined_cov)./abs(min([mean_cov combined_cov])));
        pc_difference_wT(i_int) = 100.*abs((mean_wT-combined_wT)./abs(min([mean_wT combined_wT])));

        grid on;
        %% plot stationarity

        subplot_handles{3}=subplot(n_y_subplots,n_x_subplots,3);
        plot(n2o.all_results(1,n2o_range(18)),pc_difference(i_int),'ro');
        hold on;
        plot(n2o.all_results(1,n2o_range(18)),pc_difference_wT(i_int),'b+');
        ylabel({'stationarity','difference [%]'});
        ylim([0 65]);
        set(gca,'ytick',[0 30 60]);
        grid on;
        legend({'wN2O','wT'});




        %% plot wind dir (by yaw rot'n angle)
        subplot_handles{4}=subplot(n_y_subplots,n_x_subplots,4);
        wind_dir=180./pi.*asin(sin(n2o.all_results(21,n2o_range)));
        plot(n2o.all_results(1,n2o_range),wind_dir,'k-+');
        hold on;
        wind_max(i_int)=max(wind_dir);
        wind_min(i_int)=min(wind_dir);
        ylabel('yaw rotation [deg]');
        ylim([-180 180]);
        grid on;

        %% plot u*
        subplot_handles{5}=subplot(n_y_subplots,n_x_subplots,5);
        plot(n2o.all_results(1,n2o_range),n2o.all_results(17,n2o_range),'k-+');
        hold on;
        ylabel('u* [m s-1]');
        ylim([0 1.5]);
        grid on;

        %% plot stability (z/L)
        subplot_handles{6}=subplot(n_y_subplots,n_x_subplots,6);
        z=6; % measurement height
        plot(n2o.all_results(1,n2o_range),z./n2o.all_results(19,n2o_range),'k-+');
        hold on;
        grid on;
        ylabel('stability z/L []');
        %% plot ITC measures
        %         subplot(n_y_subplots,n_x_subplots,7);
        %         itc_meas=n2o.all_results(9,n2o_range)./ ...
        %             n2o.all_results(17,n2o_range);

        %% apply same properties (xlim)to all subplots except first (Rswin
        %% plot)
        for i_sp=1:numel(subplot_handles)
            if (~isempty(subplot_handles{i_sp})) % check subplot is used
                % xlimits
                set(subplot_handles{i_sp},'xlim',...
                    [n2o.all_results(1,n2o_range(1)) n2o.all_results(1,n2o_range(end))]);
                % xlabel
                set(get(subplot_handles{i_sp},'xlabel'),'string',num2str(met_year));
                % xticks
                xt=get(subplot_handles{i_sp},'xtick');
                set(subplot_handles{i_sp},'xticklabel',datestr(xt,'HHh dd/mm'));
                set(get(subplot_handles{i_sp},'xlabel'),'string',['date, ',num2str(met_year)]);

            end
        end
    end
    set(gcf,'name',[num2str(i_int)]);
    disp(['(',num2str(i_int),'/',num2str(numel(good_intervals)),')']);
    mb=msgbox(['(',num2str(i_int),'/',num2str(numel(good_intervals)),')'],'continue');
    set(mb,'position', [321.5 100.667 125 51.5]);
    waitfor(mb);
end
%% relabel cov() plot
subplot(n_y_subplots,n_x_subplots,2);
xlim([m.datenum(good_intervals(1)) m.datenum(good_intervals(end))]);
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,19));
xlabel(['date, ',num2str(met_year)]);



%% relabel wind dir plot
subplot(n_y_subplots,n_x_subplots,4);
xlim([m.datenum(good_intervals(1)) m.datenum(good_intervals(end))]);
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,19));
hold on;
grid on;
%% overplot with vertical lines min/max
for i_int=1:numel(good_intervals);
    line([m.datenum(good_intervals(i_int)) m.datenum(good_intervals(i_int)) ],...
        [wind_min(i_int) wind_max(i_int)],'color','g');
end

%% relabel u* plot
subplot(n_y_subplots,n_x_subplots,5);

xlim([m.datenum(good_intervals(1)) m.datenum(good_intervals(end))]);
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,19));
xlabel(['date, ',num2str(met_year)]);
hold on;
grid on;


%% relabel stability plot
subplot(n_y_subplots,n_x_subplots,6);

xlim([m.datenum(good_intervals(1)) m.datenum(good_intervals(end))]);
xt=get(gca,'xtick');
set(gca,'xticklabel',datestr(xt,19));
xlabel(['date, ',num2str(met_year)]);

hold on;
grid on;

% qc codes explained in Foken et al. (1997)
% 0% < stat < 30%   :   good,       qc=0
% 30% < stat < 60%   :  acceptable, qc=1
% 60% < stat < 100%   :   bad,        qc=2
stationarity_qc=floor(pc_difference./30);
stationarity_qc_wT=floor(pc_difference_wT./30);


%% finsish
mark;
