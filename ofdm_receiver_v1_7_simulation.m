
%%


%ofdm_receiver
clear all;
clc;
close all;
addpath('');% path of the transmitter code  
files_path = ''; % path of experiment data

%% temporary parameter
max_pathnum = 25;
power_threshold = 40;
err_locate_max_variance = 20;
max_path_delay = 200;
freq_com_type = 'NON';% frequency offset compensation method, 'NON': no compensation; 'Current Frame': using frequency offset detected in this frame; 'Freq Linear Interp': using frequency offset detected in current and next frame with linear interpolace
freq_com_2 = 'OFF';% OFF ON
freq_com_2_adjust = 'OFF';% OFF ON
rm_in_cs = 'FALSE';
rm_in4sync = 'FALSE';%'TURE' 'FALSE'



para.frame_type = 'type5';
para.channel_type = 'type1';
para.coding_type = 'LTE_turbo_0.33';
para.data_m_type = 'QPSK';
%para.channel_est_type = 'linear_LS';
para = simulation_setting(para);
para = simulation_source(para);

% passband filter
N_readframe=2;
frame_type = para.frame_type;
N_cp_ofdm = para.N_cp_ofdm;
cell_type_flag = para.cell_type_flag;
frame_structure = para.frame_structure;
sync_pilot_length = para.sync_pilot_length;
data_m_type = para.data_m_type;
N_fft = para.N_fft;
N_pilot_ofdm = para.N_pilot_ofdm; %512/4 for channel estimation
N_null_ofdm = 64; %512/8 for frequency offset estimation
N_data_ofdm = para.N_data_ofdm;
coding_type = para.coding_type;
N_coded_bit_frame = para.N_coded_bit_frame;
N_uncoded_bit_frame = para.N_uncoded_bit_frame;
carrier_signal =para.carrier_signal;
N_sample_sym = para.N_sample_sym;
N_sample_frame = para.N_sample_frame;
symbol_interval = para.symbol_interval;
sample_interval = para.sample_interval;
fc = para.fc;

ce_equ_type = 'scce&rmin_lsequ';% method of channel estimation and equalization 

%
% parameter calculating
N_ofdm_frame = length( find(frame_structure == 1) ); 
freq_pass_gate = 0.5; % min frequency offset

itera=1;
ori_err_bit_num_coded = -1*ones(N_readframe,itera);

% data generate



% frame tone mapping



freq_frame = zeros(512,5);
pilot_posi = [1:4:para.N_fft]+1;
freq_frame( pilot_posi, : ) = 101;
null_posi = [1:8:para.N_fft]+3;
freq_frame( null_posi, : ) = 102;
pilot_index = find(freq_frame == cell_type_flag(2));
null_index = find(freq_frame == cell_type_flag(3));
data_index = find(freq_frame == 0);
freq_frame_ofdm=freq_frame;

% data generate end




ori_total_raw_ber=zeros(16,itera);
for snr=1:16  
for i_read = 1:N_readframe
frame_ofdm_freq=zeros(512,5);   
pilot_tone_bits = 1 - 2*(rand(N_pilot_ofdm*2*N_ofdm_frame,1)>0.5);
pilot_frame = (pilot_tone_bits(1:2:end) + 1j*pilot_tone_bits(2:2:end))/sqrt(2);
null_frame = zeros(N_null_ofdm*N_ofdm_frame,1); 
freq_frame(pilot_index) = pilot_frame;
freq_frame(null_index) = null_frame; 
coded_bits_frame = rand(1,3200)>0.5;
if strcmpi(para.data_m_type, 'QPSK')
    Mrank = 2;
elseif strcmpi(para.data_m_type, 'QAM16')
    Mrank = 4;
end
data_frame = f_modulation(coded_bits_frame,Mrank);
data_frame_original = data_frame;
freq_frame(data_index) = data_frame;
para.freq_frame = freq_frame;   
%yt = rayleigh(reshape(freq_frame,1,2560));
%freq_frame=reshape(yt,512,5);
    snr_bd=snr-1-10;
    thegma_sqare=1/power(10,(snr_bd/10));
    thegmaw=sqrt(thegma_sqare/8.98);
    thegmai=sqrt((400*thegma_sqare)/8.98);
    noise_w_r=randn(2560,1)*(thegmaw/sqrt(2));
    noise_w_i=randn(2560,1)*(thegmaw/sqrt(2));
    noise_i_r=randn(2560,1)*(thegmai/sqrt(2));
    noise_i_i=randn(2560,1)*(thegmai/sqrt(2));
    noise_w=noise_w_r+1j*noise_w_i;
    noise_i=noise_i_r+1j*noise_i_i;
    noise=0.98*noise_w+0.02*noise_i;
    noise_frame=reshape(noise,512,5);
    frame_ofdm_block=freq_frame+noise_frame;
    %frame_ofdm_block=freq_frame;
    frame_ofdm_freq = fftshift( fft( frame_ofdm_block, N_fft ), 1 ) / sqrt(N_fft);
    pilot_frame_freq=pilot_frame;
    
    %Channel est and equ
    %[frame_equ snr_est channel_est] = ce_equ_snr(frame_ofdm_freq, para, 40);
N_ofdm_frame = para.N_ofdm_frame;
N_fft = para.N_fft;

channel_est_pilot = zeros(N_fft,1);
snr_est = zeros( N_fft, N_ofdm_frame );
tap_num=40;

if strcmpi(ce_equ_type, 'lsce_lsequ')
    

    pilot_frame_block = reshape(pilot_frame, [], N_ofdm_frame );

    channel_est_freq = zeros( N_fft, N_ofdm_frame );  

    for i_ofdm = 1:N_ofdm_frame
        ofdm_freq_i = frame_ofdm_freq( :,i_ofdm );
        pilot_tone_i = pilot_frame_block( :, i_ofdm );

        pilot_tone_recv = ofdm_freq_i(pilot_index);

        channel_est_LS = pilot_tone_recv ./ pilot_tone_i;
        channel_est_zero = [channel_est_LS; zeros( N_fft-length(channel_est_LS),1 )];
        channel_est_pilot(pilot_index) = channel_est_LS;

        channel_est_time = ifft( channel_est_zero, N_fft ) * sqrt(N_fft);
        noise_delta = channel_est_time(201:400)' * channel_est_time(201:400) / 200 * 512/length(channel_est_LS);

        channel_est_i = interp1(pilot_index,channel_est_LS, [1:N_fft],'PCHIP', 'extrap' );
        channel_est_freq(:,i_ofdm) = channel_est_i.';
        snr_est(:,i_ofdm) = max( 0.01, (channel_est_i.' .* channel_est_i'-noise_delta) / noise_delta );

    end
    
    
    frame_equ = frame_ofdm_freq .* conj( channel_est_freq ) ./ abs(channel_est_freq).^2;
    channel_est = channel_est_freq;

elseif strcmpi(ce_equ_type, 'scce&rmin_lsequ') % est channel with in and then rmin
        
    polit_index_ofdm=find((freq_frame_ofdm(:,1) == 101));
    pilot_data_index=find((freq_frame_ofdm(:,1) == 101)|(freq_frame_ofdm(:,1) == 0));
    change_frame_ofdm_freq=frame_ofdm_freq;
    P=zeros(128,512);
    for a=1:128
        P(a,polit_index_ofdm(a))=1;
    end
    pilot_frame_block = reshape(pilot_frame_freq, [], N_ofdm_frame );
    
    N_fft_p = length( polit_index_ofdm );
    fft_matrix_pilot = fft( eye(N_fft_p) ) / sqrt( N_fft_p );
    
    fft_matrix = fft( eye(N_fft) ) / sqrt( N_fft );
    
    channel_est_freq = zeros( N_fft, N_ofdm_frame );
    itera1_h_est = zeros( N_fft, N_ofdm_frame );
    frame_ofdm_freq_rmin = zeros( N_fft, N_ofdm_frame );
    ofdm_freq_rmin= zeros( N_fft, N_ofdm_frame );
    noise_delta=zeros(5,1);
    pilot_noise=zeros(128,5);
    itera1_pilot=zeros(128,5);
            itera1_rt=zeros( N_fft, N_ofdm_frame );
        itera1_rp=zeros( N_fft, N_ofdm_frame );
        
for u=1:itera
    for i_ofdm = 1:N_ofdm_frame
        
        ofdm_freq_i = frame_ofdm_freq( :,i_ofdm );
        change_ofdm_freq_i=change_frame_ofdm_freq( :,i_ofdm );
        pilot_tone_i = pilot_frame_block( :, i_ofdm );

        pilot_tone_recv = change_ofdm_freq_i(polit_index_ofdm);
        
        %measure_matrix = [ diag( pilot_tone_i ) * fft_matrix_pilot fft_matrix_pilot ];
        %[tap_index x_vector] = omp_algorithm(pilot_tone_recv,measure_matrix,tap_num,1e-5);
        %channel_est_cs = fft_matrix_pilot * x_vector(1:end/2);
        %in_est_cs = fft_matrix_pilot * x_vector(end/2+1:end);%频域
        
        %         in_freq( pilot_index ) = in_est_cs;
%         in_time = ifft( fftshift(in_freq,1), N_fft ) * sqrt(N_fft) * N_fft / length(pilot_index); % BUG
%         ofdm_time = ifft( fftshift(ofdm_freq_i,1), N_fft ) * sqrt(N_fft);
        ofdm_time = ifft( change_ofdm_freq_i, N_fft ) * sqrt(N_fft);
        ofdm_time_time_power = ofdm_time .* conj(ofdm_time);
        ofdm_time_time_power_mean = mean( ofdm_time_time_power );
        in_index = (find(ofdm_time_time_power > ofdm_time_time_power_mean*5 ));
        Pi=zeros(512,length(in_index));
        for b=1:length(in_index)
            Pi(in_index(b),b)=1;
        end
        %end estimate the IN position using received signal
        Fi=P*fft_matrix*Pi;
        LS_in=LS_est(pilot_tone_recv,Fi);
        %[123;182;183;185;197;199;200;201;203;204;392];
        %est_matrix = fft_matrix(pilot_index,in_index);
        %in_time_select = pinv(est_matrix) * in_est_cs;
        %if(LS_in==in_time_select)
        %    correct=1;
        %end
        %ofdm_time_rmin = ofdm_time;
%         ofdm_time_rmin(in_index) = ofdm_time(in_index) - in_time(in_index);
%          ofdm_time_rmin(in_index) = ofdm_time(in_index) - in_time_select;
      
        ofdm_freq_rmin(:,i_ofdm) = ofdm_freq_i - fft_matrix*Pi*LS_in;
       
        %frame_ofdm_freq_rmin( :,i_ofdm ) = fft( ofdm_time_rmin, N_fft ) / sqrt(N_fft);
        itera1_pilot(:,i_ofdm)=ofdm_freq_rmin(polit_index_ofdm,i_ofdm);
        PD=zeros(128,512);
        for h=1:128
            PD(h,polit_index_ofdm(h))=pilot_tone_i(h);
        end
        pdf= PD * fft_matrix;
        %[~,itera1_h_est(:,i_ofdm)]=omp_algorithm(itera1_pilot(:,i_ofdm),pdf,tap_num,1e-5);

        fre_channel=fft_matrix *itera1_h_est(:,i_ofdm);
        fre_channel=ones(512,1);
%         frame_ofdm_freq_rmin( :,i_ofdm ) = fftshift( fft( ofdm_time_rmin, N_fft ), 1 ) / sqrt(N_fft);
        
        fre_pilot_channel=fre_channel(polit_index_ofdm);
        %pilot_tone_recv = frame_ofdm_freq_rmin( pilot_index,i_ofdm );
        pilot_noise(:,i_ofdm) = itera1_pilot(:,i_ofdm) - fre_pilot_channel .* pilot_tone_i;
        %pilot_noise = pilot_tone_recv - channel_est_cs .* pilot_tone_i;
        noise_delta(i_ofdm) = mean( pilot_noise(:,i_ofdm) .* conj(pilot_noise(:,i_ofdm)) );


        %channel_est_i = interp1(pilot_index,channel_est_cs, [1:N_fft],'PCHIP', 'extrap' );
        channel_est_freq(:,i_ofdm) = fre_channel.';
        snr_est(:,i_ofdm) = max( 0.01, (fre_channel.' .* fre_channel'-noise_delta(i_ofdm)) / noise_delta(i_ofdm) );
        
    end
    
    frame_equ =ofdm_freq_rmin .* conj( channel_est_freq ) ./ abs(channel_est_freq).^2;
    channel_est = channel_est_freq;
    

    

    % data detect
    data_index = find( freq_frame_ofdm == 0 );

    data_equ = frame_equ(data_index);

    snr_est_data = snr_est( data_index );

    if strcmpi(data_m_type, 'QPSK')
        N_bit_symbol = 2;
    elseif strcmpi(data_m_type, 'QAM16')
        N_bit_symbol = 4;
    end
    demodu_xm_out = f_demodulation_symbol (data_equ,length(data_equ),snr_est_data,N_bit_symbol);

    bits_rcv_hard = demodu_xm_out>0;
    
     %decoding    
  
    % find error position
    ori_err_index_coded = find(bits_rcv_hard ~= coded_bits_frame);
    ori_err_bit_num_coded(i_read,u) = length( ori_err_index_coded );
    ori_total_raw_ber(snr,u)=ori_total_raw_ber(snr,u)+ori_err_bit_num_coded(i_read,u);
    %for u=1:3  
%    data_frame = f_modulation(bits_rcv_hard,N_bit_symbol);
%    frame_fb = freq_frame;
%    frame_fb(data_index) = data_frame;
%    for i_ofdm = 1:N_ofdm_frame
%        change_frame_ofdm_freq(:,i_ofdm)=frame_ofdm_freq(:,i_ofdm)-diag(frame_fb(:,i_ofdm))*fft_matrix*itera1_h_est(:,i_ofdm);
%    end
end
    % re-generate the frame
end
    if strcmpi(rm_in_cs, 'TURE')
        
        frame_fb = freq_frame;
        data_frame = f_modulation(bits_rcv_hard,N_bit_symbol);
        data_frame(rm_index) = 0;

        frame_fb(data_index) = data_frame;

        [frame_rm_in channel_est snr_est] = channel_in_est( frame_ofdm_freq, frame_fb, para );
        frame_equ = frame_rm_in ./ channel_est;

        snr_est_data = snr_est( data_index );
        data_equ_rm_in = frame_equ(data_index);
        demodu_xm_out = f_demodulation_symbol (data_equ_rm_in,length(data_equ_rm_in),snr_est_data,N_bit_symbol);
        bits_rcv_hard_rm_in = demodu_xm_out>0;
        err_index_coded_rm_in = find(bits_rcv_hard_rm_in ~= coded_bits_frame);
        bits_rcv_hard = bits_rcv_hard_rm_in;

    end
    %end   
    %decoding    
end
end

color=['r','g','b'];
for u=1:1
t=-10:5;
y=ori_total_raw_ber(:,u)/(N_readframe*3200);
%1/2码率是3250；1/3码率是3250
%ori_total_coded_ber(u)/(ori_actual_frame(u)*1088)%1/2码率是1632；1/3码率是1088
%ori_frame_ber(u)/ori_actual_frame(u)
plot(t,y,color(u))
hold on
end