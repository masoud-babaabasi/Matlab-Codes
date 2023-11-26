function  set_freq(s,f_cpu,freq)
% s = serialport(serial_com,115200);
TXbuf = zeros(1,2);
TXbuf(1) = double('F');
TXbuf(2) = round( f_cpu*1e6 / ( freq * 1e3 ) / 2 ) ;
write(s , TXbuf , "uint8");
% clear s
end

