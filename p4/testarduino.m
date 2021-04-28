clear all;

a = arduino('com3','Mega2560','libraries','ExampleLCD/LCDAddon','ForceBuildOn',true);

lcd = addon(a,'ExampleLCD/LCDAddon','RegisterSelectPin','D36','EnablePin','D38','DataPins',{'D40','D42','D44','D46'});

initializeLCD(lcd);

printLCD(lcd,'test');

