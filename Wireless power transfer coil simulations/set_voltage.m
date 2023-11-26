function set_voltage(s,voltage)
% s = serialport(serial_com,115200);
TXbuf = zeros(1,2);
TXbuf(1) = double('V');
TXbuf(2) = floor( (voltage - 1.9) * 255 / 10 ) ;
write(s , TXbuf , "uint8");
% clear s
end

