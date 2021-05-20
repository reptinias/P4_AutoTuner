clear all
a = arduino('COM3','Mega2560','Libraries','ExampleLCD/LCDAddon','ForceBuildOn',true);
lcd = addon(a,'ExampleLCD/LCDAddon','RegisterSelectPin','D36','EnablePin','D38','DataPins',{'D40','D42','D44','D46'});
configurePin(a, 'D22', 'pullup');
configurePin(a, 'D24', 'pullup');
initializeLCD(lcd);
Fs = 44000;                    %# sampling frequency in Hz
T = 1;                        %# length of one interval signal in sec


%# prepare audio recording
recObj = audiorecorder(Fs,16,1);



printLCD(lcd,'Start singing'); 

for i=1:30
    
    recordblocking(recObj, T);

    %# get data and compute FFT
    sig = getaudiodata(recObj);
    
    win = kbdwin(1024);
    overlapLength = 0.75*numel(win);

    S = stft(sig, ...
    "Window",win, ...
    "OverlapLength",overlapLength, ...
    "Centered",false);
    
    pinstatus1 = readDigitalPin(a, 'D22');
    pinstatus2 = readDigitalPin(a, 'D24');
    
    if pinstatus1 == 0 && pinstatus2 == 1
        
        nsemitones = 2;
        lockPhase = true;
        audioOut = shiftPitch(S,nsemitones, ...
                     "Window",win, ...
                     "OverlapLength",overlapLength, ...
                     "LockPhase",lockPhase);
        printLCD(lcd,'pitch up'); 
        sound(audioOut,Fs)
    elseif pinstatus1 == 1 && pinstatus2 == 0
        
        nsemitones = -2;
        lockPhase = true;
        audioOut = shiftPitch(S,nsemitones, ...
                     "Window",win, ...
                     "OverlapLength",overlapLength, ...
                     "LockPhase",lockPhase);
        printLCD(lcd,'Pitch decrease');
        sound(audioOut,Fs)
        
    elseif pinstatus1 == 0 && pinstatus2 == 0
       printLCD(lcd,'original'); 
       sound(sig,Fs)
    end
end
