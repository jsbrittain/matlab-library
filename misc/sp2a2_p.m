function [f,t,cl,sp] = sp2a2_p(dat1,dat2,i1,i2,i3,i4,i5,i6,i7)
% function [f,t,cl,sp] = sp2a2_p(dat1,dat2,samp_rate,seg_pwr,opt_str,[dat3])
%                        sp2a2_p(dat1,dat2,start_times,duration,samp_rate,seg_pwr,opt_str,[dat3])
%                        sp2a2_p(dat1,dat2,start_times,offset,points,samp_rate,seg_pwr,opt_str,[dat3])
%
% For type 0,1,2 analysis respectively
%

if (nargin==5)          % Type 0 Analysis without partial coherence
    samp_rate=i1;
    seg_pwr=i2;
    opt_str=i3;
    spec_type=0;
    partial=0;
elseif (nargin==6)      % Type 0 Analysis with partial coherence
    samp_rate=i1;
    seg_pwr=i2;
    opt_str=i3;
    spec_type=0;
    dat3=i4;
    partial=1;
elseif (nargin==7)      % Type 1 Analysis
    start_times=i1;
    duration=i2;
    samp_rate=i3;
    seg_pwr=i4;
    opt_str=i5;
    spec_type=1;
    partial=0;
elseif (nargin==8)
    if (size(i6,1)>1)
        start_times=i1;
        duration=i2;
        samp_rate=i3;
        seg_pwr=i4;
        opt_str=i5;
        spec_type=1;
        dat3=i6;
        partial=1;
    elseif (size(i6,1)==1)
        start_times=i1;
        samp_rate=i4;
        points_ms=i3;
        offset_ms=i2;
        points=round((i3/1000)*samp_rate);
        offset=round((i2/1000)*samp_rate);
        seg_pwr=i5;
        opt_str=i6;
        spec_type=2;
        partial=0;
    end
elseif (nargin==9)
    start_times=i1;
    samp_rate=i4;
    points_ms=i3;
    offset_ms=i2;
    points=round((i3/1000)*samp_rate);
    offset=round((i2/1000)*samp_rate);
    seg_pwr=i5;
    opt_str=i6;
    spec_type=2;
    dat3=i7;
    partial=1;
else
    error(' Incorrect number of input arguments');
end;

if (nargout<3)
    error(' Not enough output arguments');
end
if (nargout>4)
    error(' Too many output arguments');
end;
if (nargout==4)
    spect_out_flag=1;
else
    spect_out_flag=0;
end;

if (partial)
    disp(['Performing a Partial Coherence, Type: ' num2str(spec_type)]);
else
    disp(['Performing a type ' num2str(spec_type),' analysis...']);
end;

pts_tot=length(dat1);
if (length(dat2)~=pts_tot)
    error(['Unequal length data arrays. Lengths: ' num2str(pts_tot) ' and ' num2str(length(dat2))]);
end;
if (partial & length(dat3)~=pts_tot)
    error(['Predictor data array is not the correct length.   Should be: ' num2str(pts_tot) ' but is: ' num2str(length(dat3))]);
end;

seg_size=2^seg_pwr;     % DFT segment length (T)
seg_samp_min=0;         % Can specify an absolute value - overrides seg_min_fac
seg_min_fac=0.05;       % Define what fraction of samples MUST be present for us to calculate a segment (Type 1 analysis)
burst_limit_min=0;      % Specify a minimum data length (points) for type 2 analysis

% Decide on minimum number of samples allowed in a segment
if (spec_type==1)
    if (seg_samp_min==0)
        seg_samp_min=round(seg_size*seg_min_fac);
    end;
end;

if (spec_type==2)
    if (points>seg_size)        % If we are trying to look at more than one DFT length
        disp(['Data samples per segment: ' num2str(points) ', DFT segment size: ' num2str(seg_size)]);
        error('sp2a2_p - Type 2 analysis only uses one segment - points variable has been set too long');
    end;
    if (points<burst_limit_min) % If 'points' is less than a specified minimum length
        disp(['Data samples per segment: ' num2str(points) ', DFT segment size: ' num2str(seg_size)]);
        error('sp2a2_p - Points variable is below minimum allowed');
    end;
end;

% Arrange data into separate columns, one column per segment (T rows)

% For type 0 analysis, arrange into L columns, T rows
if (spec_type==0)       % For a disjoint sections analysis
    seg_tot=fix(pts_tot/seg_size);      % Number of complete segments
    samp_tot=seg_tot*seg_size;          % Number of samples to analyse R=LT
    
    % Arrange data into L columns each with T rows
    rd1=reshape(dat1(1:samp_tot),seg_size,seg_tot);
    rd2=reshape(dat2(1:samp_tot),seg_size,seg_tot);
    no_of_samples(1:seg_tot)=seg_size;  % Record how many samples are used in each segment
    md1=mean(rd1);
    md2=mean(rd2);                      % Calculate mean for each segment
    if (partial)
        rd3=reshape(dat3(1:samp_tot),seg_size,seg_tot);
        md3=mean(rd3);
    end;
    
elseif (spec_type==1)   % For type 1 analysis
    seg_num=1;          % Use as many segments as required for each burst.  This counts them.
    samp_tot=0;
    for n=1:length(start_times)     % Do this for every burst start time
        burst_offset=0;             % A counter to move through the burst
        while ((duration(n)-burst_offset)>=seg_samp_min)    % While there are enough samples left in the burst for another segment
            st = start_times(n) + burst_offset;             % Start time for this section
            et = st + seg_size - 1;                         % End time for this section
            if (et > start_times(n)+duration(n))            % If we've gone too far
                et = start_times(n)+duration(n)-1;          % then set et to be the end of the duration
            end;
            burst_offset=burst_offset+(et-st)+1;            % Keep record of the burst offset. Move on to the end of the segment, then on one more for start of next segment
            no_of_samples(seg_num) = et-st+1;               % From st to et inclusive, so add one.
            rd1(:,seg_num)=[dat1(st:et)', zeros(1,(seg_size-no_of_samples(seg_num)) )]';    % Store segment in its own column, length st:et, pad with zeros if shorter than T.
            rd2(:,seg_num)=[dat2(st:et)', zeros(1,(seg_size-no_of_samples(seg_num)) )]';
            samp_tot = samp_tot + no_of_samples(seg_num);   % Keep a count of the number of samples involved
            md1(seg_num)=mean(rd1(1:no_of_samples(seg_num),seg_num));
            md2(seg_num)=mean(rd2(1:no_of_samples(seg_num),seg_num));
            
            if (partial)
                rd3(:,seg_num)=[dat3(st:et)', zeros(1,(seg_size-no_of_samples(seg_num)) )]';
                md3(seg_num)=mean(rd3(1:no_of_samples(seg_num),seg_num));
            end;
            seg_num=seg_num+1;      % Count on to next segment
        end;
    end;
    
elseif (spec_type==2)   % For type 2 analysis. Just take a single segment
    rd1=zeros(seg_size,length(start_times));                % Preallocate memory for speed
    rd2=zeros(seg_size,length(start_times));
    if (partial)
        rd3=zeros(seg_size,length(start_times));
    end;
    samp_tot=0;
    seg_num=0;                      % No segments filled in yet
    for n=1:length(start_times)
        st = start_times(n) + offset;                       % Stat time for this section
        if ((st>0) & (st+points-1 <= length(dat1)))         % Check section is within range - ignore those that aren't
            seg_num=seg_num+1;                              % Count how many segments we've filled in...
            et = st+points-1;                               % st+points is the next start time, so st+points-1 is the end of this one
            no_of_samples(seg_num)=points;                  % Count the number of samples used in each segment (ignoring zero padding)
            samp_tot=samp_tot+no_of_samples(seg_num);
            
            rd1(1:no_of_samples(seg_num),seg_num)=dat1(st:et);
            rd2(1:no_of_samples(seg_num),seg_num)=dat2(st:et);  % Store the segment in its own column, length st:et
                                                                % The rest are already filled up with zeros
            md1(seg_num)=mean(rd1(1:no_of_samples(seg_num),seg_num));   % Find the mean, excluding zero padding, of each segment
            md2(seg_num)=mean(rd2(1:no_of_samples(seg_num),seg_num));
            
            if (partial)
                rd3(1:no_of_samples(seg_num),seg_num)=dat3(st:et);
                md3(seg_num)=mean(rd3(1:no_of_samples(seg_num),seg_num));
            end;
        end;
    end;
    
    if (size(rd1,2)>seg_num)        % If we preallocated too much memory (i.e. a start time was not used)
        rd1(:,(seg_num+1):size(rd1,2))=[];      % Delete any excess columns
        rd2(:,(seg_num+1):size(rd2,2))=[];
        if (partial)
            rd3(:,(seg_num+1):size(rd3,2))=[];
        end;
    end;
end;

% Process options
trend_chan_1=0;     % Set defaults - options off
trend_chan_2=0;
trend_chan_3=0;
mains_flag=0;
rect_chan_1=0;
rect_chan_2=0;
rect_chan_3=0;
if (nargin<5)       % No options supplied
    opt_str='';
end;
options=deblank(opt_str);
while (any(options))    % Parse individual options from string
    [opt,options]=strtok(options);
    optarg=opt(2:length(opt));  % Determine option argument
    switch (opt(1))
        case 'r'                % Rectification option
            i=str2num(optarg);
            if (i<0 | i>6)
                error(['error in option argument -- r' optarg]);
            end;
            if (i>=3)
                rect_chan_3=1;  % Rectify channel 3
                i=i-4;
            end;
            if (i==0 | i==2)
                rect_chan_1=1;  % Rectify channel 1
            end;
            if (i>=1)
                rect_chan_2=1;  % Rectify channel 2
            end;
        case 't'        % Linear de-trend option
            i=str2num(optarg);
            if (i<0 | i>6)
                error(['error in option argument -- t' optarg]);
            end;
            if (i>=3)
                trend_chan_3=1;  % De-trend channel 3
                i=i-4;
            end;
            if (i==0 | i==2)
                trend_chan_1=1;  % De-trend channel 1
            end;
            if (i>=1)
                trend_chan_2=1;  % De-trend channel 2
            end;
        case 'm'        % Mains suppression option
            mains_flag=1;
        otherwise
            error(['Illegal option -- ' opt]);  % Illegal option
    end;
end;

if (trend_chan_1 | trend_chan_2)        % Index for fitting data with polynomial
    trend_x=(1:seg_size)';
end;

for ind=1:size(rd1,2)                   % Loop across columns/segments
    if rect_chan_1
        rd1(1:no_of_samples(ind),ind)=rd1(1:no_of_samples(ind),ind)-md1(ind);   % Subtract mean from ch 1. NOT from zero padding!!!
        rd1(:,ind)=abs(rd1(:,ind));     % Rectification of ch 1 (Full wave)
    end;
    if rect_chan_2
        rd2(1:no_of_samples(ind),ind)=rd2(1:no_of_samples(ind),ind)-md2(ind);   % Subtract mean from ch 2.
        rd2(:,ind)=abs(rd2(:,ind));     % Rectification of ch 2 (Full wave)
    end;
    if (rect_chan_3 & partial)
        rd1(1:no_of_samples(ind),ind)=rd1(1:no_of_samples(ind),ind)-md1(ind);   % Subtract mean from ch 3.
        rd1(:,ind)=abs(rd1(:,ind));     % Rectification of ch 3 (Full wave)
    end;
    
    % Must remove trend from all data except zero padding... or make zero mean again
	if trend_chan_1                     % Linear trend removal
        p=polyfit(trend_x(1:no_of_samples(ind)),rd1(1:no_of_samples(ind),ind),1);   % Fit 1st order polynomial
        rd1(1:no_of_samples(ind),ind)=rd1(1:no_of_samples(ind),ind)-p(1)*trend_x(1:no_of_samples(ind))-p(2);    % Subtract from ch 1.
    else                                % Otherwise make it zero mean
        mrd1=mean(rd1(1:no_of_samples(ind),ind));   % Find the new mean (excluding zero padding) of each segment
        rd1(1:no_of_samples(ind),ind)=rd1(1:no_of_samples(ind),ind)-mrd1;   % Remove the new mean - but NOT on the zero padding!
    end;
    if trend_chan_2                     % Linear trend removal
        p=polyfit(trend_x(1:no_of_samples(ind)),rd2(1:no_of_samples(ind),ind),1);   % Fit 1st order polynomial
        rd2(1:no_of_samples(ind),ind)=rd2(1:no_of_samples(ind),ind)-p(1)*trend_x(1:no_of_samples(ind))-p(2);    % Subtract from ch 2.
    else                                % Otherwise make it zero mean
        mrd2=mean(rd2(1:no_of_samples(ind),ind));   % Find the new mean (excluding zero padding) of each segment
        rd2(1:no_of_samples(ind),ind)=rd2(1:no_of_samples(ind),ind)-mrd2;   % Remove the new mean - but NOT on the zero padding!
    end;
    if (partial)
        if trend_chan_3                     % Linear trend removal
            p=polyfit(trend_x(1:no_of_samples(ind)),rd3(1:no_of_samples(ind),ind),1);   % Fit 1st order polynomial
            rd3(1:no_of_samples(ind),ind)=rd3(1:no_of_samples(ind),ind)-p(1)*trend_x(1:no_of_samples(ind))-p(2);    % Subtract from ch 3.
        else                                % Otherwise make it zero mean
            mrd3=mean(rd3(1:no_of_samples(ind),ind));   % Find the new mean (excluding zero padding) of each segment
            rd3(1:no_of_samples(ind),ind)=rd3(1:no_of_samples(ind),ind)-mrd3;   % Remove the new mean - but NOT on the zero padding!
        end;
    end;
end;

% Before here, we need to:
% 0. Put data into segments
% 1. Find mean of each segment. Take away --> zero mean
% 2. Rectify each segment
% 3. Find mean and remove it again if we're not doing trend removal
% 4. Optionally remove linear trend (which removes another zero mean)
% 5. Pad with zeros
% We've actually padded with zeros at the beginning, but only remove means
% on selections of data, NOT zero padding

% Call sp2_p_fn()   New periodogram based spectral estimation routine.
if (partial)        % Send the predictor for partial coherence
    [f,t,cl,sp_struct]=sp2_p_fn(rd1,rd2,samp_rate,seg_samp_min,seg_size,samp_tot,mains_flag,rd3);
else                % Or don't
    [f,t,cl,sp_struct]=sp2_p_fn(rd1,rd2,samp_rate,seg_samp_min,seg_size,samp_tot,mains_flag);
end;

% sp_struct will be returned from sp2_p_fn, but sp will only be given to
% the caller if they have specified the sp output parameter
if (spect_out_flag)     % Setup sp if necessary (4 columns)
    sp(:,1)=sp_struct.f11;
    sp(:,2)=sp_struct.f22;
    sp(:,3)=sp_struct.f21r;
    sp(:,4)=sp_struct.f21i;
end;

% Shift values to remove values below minimum frequency resolution.
deltaf=f(1,1);
if (spec_type==2)
    min_freq=1000/points_ms;        % 1/sec
    min_index=fix(min_freq/deltaf); % Work out lowest index to use
    if (min_index>1)
        f(1:(min_index-1),:)=[];
    end;
end;

% Work out burst_tot_var (basically 1/L) for the coherence confidence limit:
burst_tot_var=0;
for n=1:size(rd1,2)
    burst_weight(n)=no_of_samples(n)/samp_tot;
    burst_tot_var=burst_tot_var + burst_weight(n)*burst_weight(n);
end;

cl.type=spec_type;
cl.partial=partial;
cl.seg_tot_var=1/burst_tot_var;         % Effective no. of segments: for multivariate pooled spectra
                                        % i.e. this is the effective L
if (partial)        % For partial coherence limits, replace L with L-1 (1 predictor)
    cl.seg_tot_var=cl.seg_tot_var-1;
end;

ch_pwr=1/(cl.seg_tot_var-1);
cl.f_c95=0.8512*sqrt(1/cl.seg_tot_var); % Confidence limit for spectral estimates, PBMB (6.2)
cl.ch_c95=1-0.05^ch_pwr;                % Confidence limit for coherence PBMB (6.6)
% N.B. Confidence interval for log spectra is TWICE this value

% Set additional elements in cl structure
cl.N1=0;                % N1, No of events in ch 1. (zero for TS data)
cl.N2=0;                % N2, No of events in ch 2. (zero for TS data)
cl.P1=0;                % P1, mean intensity ch 1.
cl.P2=0;                % P2, mean intensity ch 2.
cl.opt_str=opt_str;     % Copy of options string
cl.what='';             % Field for plot label

% Display No of segments & resolution
disp(['Segments: ' num2str(size(rd1,2)) ', Segment length: ' num2str(seg_size/samp_rate) ' sec, Resolution: ' num2str(cl.df) 'Hz']);
switch spec_type
    case {0,1}
        disp(['Type: ' num2str(spec_type) ', Samples: ' num2str(samp_tot)]);
    case 2
        disp(['Type: ' num2str(spec_type) ', Samples: ' num2str(samp_tot) ', Offset: ' num2str(offset_ms) 'ms, Duration: ' num2str(points_ms) ' ms']);
end;
