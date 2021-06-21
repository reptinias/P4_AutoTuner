%Variables
clear all
a = arduino('COM3','Mega2560','Libraries','ExampleLCD/LCDAddon'); %start arduino object  with LCD lilbrary included
lcd = addon(a,'ExampleLCD/LCDAddon','RegisterSelectPin','D36','EnablePin','D38','DataPins',{'D40','D42','D44','D46'}); %start lcd object and declare pins
configurePin(a, 'D22', 'pullup'); %button 1 declare the pin
configurePin(a, 'D24', 'pullup');% button 2 declare the pin
initializeLCD(lcd); %initializeLCD object 
Fs = 44000;                    %#Fs in HZ
T = 1;                        %# signal length

%begin the recording
recObj = audiorecorder(Fs,16,1);
printLCD(lcd,'Start singing'); %print on LCD screen


%begin for loop for 30 seconds only (this is for testing, can just be set to inf for final product)
for i=1:30
    
    recordblocking(recObj, T); %stop the recording at interval T
    audio = getaudiodata(recObj); %retrieve the audio data and save
    
    win = kbdwin(1024); %declare the window size of the kbdwin
    overlapLength = 0.75*numel(win); %declare the overlaplength in relation to win size

    %Perform short-time fourirer transform on the signal and save the
    %processed signal
    S = stft(audio, ...           
    "Window",win, ...
    "OverlapLength",overlapLength, ...
    "Centered",false);
    
    pinstatus1 = readDigitalPin(a, 'D22'); %Read the status of button 1
    pinstatus2 = readDigitalPin(a, 'D24'); %read the status of button 2
    
    %if button 2 is pressed and button 1 is not pressed then increase pitch
    if pinstatus1 == 0 && pinstatus2 == 1 
        nsemitones = 2; %set the change in semitones
        lockPhase = true; %lockPhase to true for higher fidelity in signal
        %shiftPitch command and save the new signal as audioOut
        audioOut = shiftPitch(S,nsemitones, ...  
                     "Window",win, ...
                     "OverlapLength",overlapLength, ...
                     "LockPhase",lockPhase);
        printLCD(lcd,'pitch up'); %display on LCD screen
        sound(audioOut,Fs) %play the modified audio signal
      
      %if button 1 is pressed and button 2 is not then decrease pitch  
    elseif pinstatus1 == 1 && pinstatus2 == 0
        nsemitones = -2;
        lockPhase = true;
        audioOut = shiftPitch(S,nsemitones, ...
                     "Window",win, ...
                     "OverlapLength",overlapLength, ...
                     "LockPhase",lockPhase);
        printLCD(lcd,'Pitch decrease');
        sound(audioOut,Fs)
      %if no buttons are pressed the play original sound  
    elseif pinstatus1 == 0 && pinstatus2 == 0
       printLCD(lcd,'original'); 
       sound(audio,Fs)
    end
end
