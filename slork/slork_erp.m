shsc = filterStruct(events,'src_correct == 1 & ismember(type,varargin{1}) & ismember(testType,varargin{2})',{'TARG_PRES'},{'side'});
thsc = filterStruct(events,'src_correct == 1 & ismember(type,varargin{1}) & ismember(testType,varargin{2})',{'TARG_PRES'},{'task'});

shsi = filterStruct(events,'src_correct == 0 & ismember(type,varargin{1}) & ismember(testType,varargin{2})',{'TARG_PRES'},{'side'});
thsi = filterStruct(events,'src_correct == 0 & ismember(type,varargin{1}) & ismember(testType,varargin{2})',{'TARG_PRES'},{'task'});

cr = filterStruct(events,'src_correct == 1 & ismember(type,varargin{1})',{'LURE_PRES'});

durationMS = 2000;      % signal time length in samples
offsetMS = -200;           % offset at which to start in samples
bufferMS = 1000;        % buffer (needed for filtering or resampling)
                        %   default is 0
filtfreq = [58 62];     % Filter freq (depends on type, see buttfilt)
                        %   default is []
filttype = 'stop';      % Filter type (see buttfilt)
filtorder = 1;          % Filter order (see buttfilt)
resampledRate = 200;    % Sample rate of the returned data
RelativeMS = [-200 0];  % Range for use with the relative subtraction

channel = 52;

eeg_shsc = gete_ms(channel,shsc,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate,RelativeMS);
eeg_thsc = gete_ms(channel,thsc,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate,RelativeMS);
eeg_shsi = gete_ms(channel,shsi,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate,RelativeMS);
eeg_thsi = gete_ms(channel,thsi,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate,RelativeMS);
eeg_cr = gete_ms(channel,cr,durationMS,offsetMS,bufferMS,filtfreq,filttype,filtorder,resampledRate,RelativeMS);

plot(linspace(offsetMS,durationMS,size(eeg_shsc,2)),mean(eeg_shsc,1),'b');
hold on
plot(linspace(offsetMS,durationMS,size(eeg_thsc,2)),mean(eeg_thsc,1),'r');
hold on
plot(linspace(offsetMS,durationMS,size(eeg_cr,2)),mean(eeg_cr,1),'k');
hold on
plot(linspace(offsetMS,durationMS,size(eeg_shsi,2)),mean(eeg_shsi,1),'g');
hold on
plot(linspace(offsetMS,durationMS,size(eeg_thsi,2)),mean(eeg_thsi,1),'m');
hold off
xlim([offsetMS durationMS]);
