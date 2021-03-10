import matplotlib
import pyaudio as pya
import struct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
import time

matplotlib.use('TkAgg')

# Autotuner baseret på vesten musik teori.
# Derfor arbejdes der med "The 12 temperment"
# Frekvensen for det 12 forskellige noder i oktav 0
c = 16.351
cS = 17.324
d = 18.354
dS = 19.445
e = 20.6
f = 21.827
fS = 23.124
g = 24.499
gS = 25.956
a = 27.5
aS = 29.135
b = 30.868

keyNumber = 3

# laver en liste med "grund" frekvenserne for hver node
noteListNumber = [c, cS, d, dS, e, f, fS, g, gS, a, aS, b]
noteListLetter = ["C", "CS", "D", "DS", "E", "F", "FS", "G", "GS", "A", "AS", "B"]
# laver en tom liste der kommer til at indholde alle frekvenser i 9 oktaver.
allNoteList = []

# Lister over de forskellige tone arter en sang kan være i
cMaM = [c, d, e, f, g, a, b]
gMeM = [c, d, e, fS, g, a, b]
dMbM = [cS, d, e, fS, g, a, b]
aMfSM = [cS, d, e, fS, gS, a, b]
eMcSM = [cS, dS, e, fS, gS, a, b]
bMgSM = [cS, dS, e, fS, gS, aS, b]
gbMdSM = [cS, dS, f, fS, gS, aS, b]
dbMaSM = [cS, dS, f, fS, gS, aS, c]
cbMabM = [b, cS, dS, e, fS, gS, aS]
fSMebM = [b, cS, dS, f, fS, gS, aS]
cSMbbM = [c, cS, dS, f, fS, gS, aS]
abMfM = [c, cS, dS, f, g, gS, aS]
ebMcM = [c, d, dS, f, g, gS, aS]
bbMgM = [c, d, dS, f, g, a, aS]
fMdM = [c, d, e, f, g, a, aS]

# En lister der indeholder alle tone arter
keyList = [cMaM, gMeM, dMbM, aMfSM, eMcSM, bMgSM, gbMdSM, dbMaSM, cbMabM, fSMebM, cSMbbM, abMfM, ebMcM, bbMgM, fMdM]
# Print grund frekvensen for tone arten
print(keyList[keyNumber])

# Test variable. Skal bruges som en "dummy" frekvens
testFre = 2341

# Vi fylder nu allNoteList med værdier, som er de frekvenser der er ind key.
# Vi starter med at gå igennem 9 oktaver så kan finde frekvenser mellem c0 - b8
for x in range(9):
    # Her går vi igennem alle frekvenserne for de tolv grundnoder og ganger dem med 2^n, hvor n  er oktaven
    # man vil finde noden i. Forklaringen kan også læses i "Indtroduction to audioprocessing" side 69(nice)
    # og side 70
    for y in keyList[keyNumber]:
        allNoteList.append(y * 2 ** x)

# value for at finde node betegnelsen
fixedNote = min(allNoteList, key=lambda z: abs(z - testFre))
noteIndex = allNoteList.index(fixedNote)

octaveValue = 0
while noteIndex > 7:
    octaveValue += 1
    noteIndex -= 7

if fixedNote - testFre < 0:
    print("you are sharp by " + str(abs(fixedNote - testFre)))
else:
    print("you are flat by " + str(fixedNote - testFre))

print("autotuned note is " + noteListLetter[noteIndex - 1] + str(octaveValue))


class AudioHandler(object):
    # RECORD_SECONDS = 5

    def __init__(self):
        self.FORMAT = pya.paInt16
        self.CHANNELS = 1
        self.RATE = 44100
        self.CHUNK = 1024 * 2
        self.stream = None
        self.p = None

    def start(self):
        print("start")
        # Instantiate PyAudio
        self.p = pya.PyAudio()

        # Open stream using callback
        self.stream = self.p.open(format=self.FORMAT,
                                  channels=self.CHANNELS,
                                  rate=self.RATE,
                                  frames_per_buffer=self.CHUNK,
                                  input=True,
                                  output=True,
                                  #stream_callback=self.callback
                                  )

    def stop(self):
        print("stop")
        self.stream.close()
        self.p.terminate()

    def callback(self, in_data, frame_count, time_info, flag):
        print("in data")
        print(in_data)
        numpy_array = np.frombuffer(in_data, dtype=np.float32)
        print("numpy array")
        print(numpy_array)
        return None, pya.paContinue

    def playAudio(self):
        root = tk.Tk()

        fig, ax = plt.subplots()

        canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.

        x = np.arange(0, 2 * self.CHUNK, 2)
        line, = ax.plot(x, np.random.rand(self.CHUNK))

        ax.set_ylim(0, 255)
        ax.set_xlim(0, self.CHUNK)

        # for loopet hvis man vil have programmet til at køre i et bestemt stykke tid
        # for i in range(0, int(self.RATE / self.CHUNK * self.RECORD_SECONDS)):

        # while loopet hvis man vil køre programmet i et ubestemt tidsinterval
        while True:
            data = self.stream.read(self.CHUNK)  # læser audio streamen og gemmer den som en stringbuffer af bits
            # print(len(data))

            # omdanner buffer fra mikrofonen til et float tuple og derefter laver det til et numpy array
            numpy_array = np.array(struct.unpack(str(2 * self.CHUNK) + 'B', data), dtype='b')[::2] + 127
            # print(len(numpy_array))

            fourierTransform = np.fft.fft(numpy_array)/len(numpy_array)
            print(np.abs(fourierTransform))

            # plt.plot(fourierTransform)
            # plt.show()

            line.set_ydata(numpy_array)
            canvas.draw()
            canvas.flush_events()

            self.stream.write(data, self.CHUNK)  # afspiller audioen fra numpy arrayet
        # self.stop()


audio = AudioHandler()
audio.start()       # open the the stream
audio.playAudio()

# fourierTransform = np.fft.fft(numpy_array)/len(numpy_array)
# print(fourierTransform)