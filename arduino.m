clear all
a = arduino('COM3','Mega2560','Libraries',{'RotaryEncoder','ShiftRegister'})
configurePin(a, 'D22', 'pullup');
channelA = 'D18';
channelB = 'D19';
encoder = rotaryEncoder(a,channelA,channelB)
dataPin = 'D29';
clockPin = 'D30';
latchPin = 'D28';
register = shiftRegister(a,'74HC595',dataPin,clockPin,latchPin)
digitTable = {...     
    '11000000', ...  % 0
    '11001111', ...  % 1
  	'10100100', ...  % 2
  	'10110000', ...  % 3
  	'10011001', ...  % 4
  	'10010010', ...  % 5
  	'10000010', ...  % 6
  	'11111000', ...  % 7
  	'10000000', ...  % 8
 	'10010000'  ...  % 9
};
for k=1:inf
   
      pinstatus = readDigitalPin(a, 'D22');
      if pinstatus == 0
       %    writeDigitalPin(a, 'D13', 0);
      %      writeDigitalPin(a, 'D12', 0)
     %       writeDigitalPin(a, 'D7', 0)
     %       writeDigitalPin(a, 'D10', 0)
     %       writeDigitalPin(a, 'D11', 0)
            writeDigitalPin(a, 'D8', 0)

      else
     %    writeDigitalPin(a, 'D13', 1);
     %     writeDigitalPin(a, 'D12', 1)
          writeDigitalPin(a, 'D8', 1)
     %     writeDigitalPin(a, 'D7', 1)
     %     writeDigitalPin(a, 'D10', 1)
     %     writeDigitalPin(a, 'D11', 1)
      end
      
  
   count = readCount(encoder);
   pos = mod(count,48);
   fprintf('position: %d\n',pos);
   
   
   
   voltage = readVoltage(a, 'A0');
   writePWMVoltage(a, 'D12', voltage);
   
   voltage = readVoltage(a, 'A1');
   writePWMVoltage(a, 'D13', voltage);
   
   digit = digitTable{4};
   write(register,digit);
   
end
    
