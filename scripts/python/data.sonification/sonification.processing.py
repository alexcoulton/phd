f = open("/mnt/e/phd.project.main/rotation1scripts_v4/original_data/datasonification/chr1a.fa", "r")
random_seq = open("/mnt/e/phd.project.main/rotation1scripts_v4/original_data/datasonification/random.sequence.fa", "r")

# f = open("/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/datasonification/chr1a.fa", "r")



# notes = f.read(100000)
notes2 = notes[7:-1]



notes2 = random_seq.read()


import re
notes3 = re.sub("\n", "", notes2)
notes4 = list(notes3)

notes_trans = []



for i in notes4:
	if(i == "T"):
		notes_trans.append(re.sub("T", "60", i))
	if(i == "A"):
		notes_trans.append(re.sub("A", "62", i))
	if(i == "G"):
		notes_trans.append(re.sub("G", "64", i))
	if(i == "C"):
		notes_trans.append(re.sub("C", "65", i))

notes_trans2 = [ int(x) for x in notes_trans ]





##### PARALLEL PROCESSING ATTEMPT ######

# def f(nucleotide):
# 	if(nucleotide == "T"):
# 		g = re.sub("T", "60", nucleotide)
# 	if(nucleotide == "A"):
# 		g = re.sub("A", "62", nucleotide)
# 	if(nucleotide == "G"):
# 		g = re.sub("G", "64", nucleotide)
# 	if(nucleotide == "C"):
# 		g = re.sub("C", "65", nucleotide)

# 	return g

# import multiprocessing
# pool = multiprocessing.Pool()

pool.map(f, notes4)





### MIDI CREATION ####

from midiutil import MIDIFile

degrees  = notes_trans2 # MIDI note number
# degrees  = [60, 62, 64, 65, 67, 69, 71, 72]  # MIDI note number
track    = 0
channel  = 0
time     = 0    # In beats
duration = 1    # In beats
tempo    = 60   # In BPM
volume   = 100  # 0-127, as per the MIDI standard

MyMIDI = MIDIFile(1)  # One track, defaults to format 1 (tempo track is created
                      
MyMIDI.addTempo(track, time, tempo)

for i, pitch in enumerate(degrees):
    MyMIDI.addNote(track, channel, pitch, time + i, duration, volume)

with open("wheat1test.random.mid", "wb") as output_file:
    MyMIDI.writeFile(output_file)


##### PROCESS GENE SEQUENCES #### 

aaseq = open("/mnt/e/phd.project.main/rotation1scripts_v4/original_data/datasonification/genes.fa", "r")
aaseq2 = aaseq.read()
aaseq3 = aaseq2.split(">")

import re

def removenewlines(seq):
	return re.sub("\n", "", seq)

aaseq4 = map(removenewlines, aaseq3)



def aamidi(aa):
	g = ""
	if(aa == "A"):
		g = re.sub("A", "36", aa)
	if(aa == "R"):
		g = re.sub("R", "38", aa)
	if(aa == "N"):
		g = re.sub("N", "40", aa)
	if(aa == "D"):
		g = re.sub("D", "41", aa)
	if(aa == "B"):
		g = re.sub("B", "43", aa)
	if(aa == "C"):
		g = re.sub("C", "45", aa)
	if(aa == "E"):
		g = re.sub("E", "47", aa)
	if(aa == "Q"):
		g = re.sub("Q", "48", aa)
	if(aa == "Z"):
		g = re.sub("Z", "50", aa)
	if(aa == "G"):
		g = re.sub("G", "52", aa)
	if(aa == "H"):
		g = re.sub("H", "53", aa)
	if(aa == "I"):
		g = re.sub("I", "55", aa)
	if(aa == "L"):
		g = re.sub("L", "57", aa)
	if(aa == "K"):
		g = re.sub("K", "59", aa)
	if(aa == "M"):
		g = re.sub("M", "60", aa)
	if(aa == "F"):
		g = re.sub("F", "62", aa)
	if(aa == "P"):
		g = re.sub("P", "64", aa)
	if(aa == "S"):
		g = re.sub("S", "65", aa)
	if(aa == "T"):
		g = re.sub("T", "67", aa)
	if(aa == "W"):
		g = re.sub("W", "69", aa)
	if(aa == "Y"):
		g = re.sub("Y", "71", aa)
	if(aa == "V"):
		g = re.sub("V", "72", aa)
	if(aa == "*"):
		g = re.sub("\\*", "74", aa)
	return g


aaseq5 = map(lambda x: list(x), aaseq4)

aaseq6 = map(lambda x: 
	map(lambda y: int(aamidi(y)), x), 
	aaseq5)


from midiutil import MIDIFile

for i in range(0, 100):
	number1 = i
	degrees  = aaseq6[i] # MIDI note number
	# degrees  = [60, 62, 64, 65, 67, 69, 71, 72]  # MIDI note number
	track    = 0
	channel  = 0
	time     = 0    # In beats
	duration = 1    # In beats
	tempo    = 60   # In BPM
	volume   = 100  # 0-127, as per the MIDI standard
	MyMIDI = MIDIFile(4)  # One track, defaults to format 1 (tempo track is created	                      
	MyMIDI.addTempo(track, time, tempo)
	MyMIDI.addTempo(1, time, tempo)
	
	for i, pitch in enumerate(degrees):	
		if pitch == 52 or pitch == 36 or pitch == 72 or pitch == 57 or pitch == 55 or pitch == 62 or pitch == 69 or pitch == 60 or pitch == 64:
			MyMIDI.addNote(0, channel, pitch, time + (i * 3), duration, volume)
		elif pitch == 65 or pitch == 67 or pitch == 45 or pitch == 71 or pitch == 40 or pitch == 48:
			MyMIDI.addNote(1, channel, pitch, time + (i * 3), duration, volume)		
		elif pitch == 41 or pitch == 47:
			MyMIDI.addNote(2, channel, pitch, time + (i * 3), duration, volume)
		elif pitch == 53 or pitch == 59 or pitch == 38:
			MyMIDI.addNote(3, channel, pitch, time + (i * 3), duration, volume)
	
	with open("aaseqs/wheat1test.aa" + str(number1) + ".mid", "wb") as output_file:
	    MyMIDI.writeFile(output_file)