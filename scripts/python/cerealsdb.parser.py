import requests
from BeautifulSoup import BeautifulSoup
import pdb
url = 'https://www.cerealsdb.uk.net/cerealgenomics/CerealsDB/quicksearch.php'
payload = {'SNP_id' : 'BS00021827'}


g = requests.post(url, data = payload)

coord = g.text.find("Contig")

g.text[coord:coord+1000]


#soup = BeautifulSoup(g.text, "html.parser")

file1 = open("/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/luzie_imputation/markers.not.in.cereals.db.txt")
markers1 = file1.read()
markers1 = markers1.split("\r\n")

contigs = []
for i in markers1:
	payload = {'SNP_id' : i}
	g = requests.post(url, data = payload)
	coord = g.text.find("Contig")	
	contigs.append(g.text[coord:coord+1000])

contigs_trunc = []
for i in contigs:
	coord1 = i.find('<div class="contig_sequence">')
	coord2 = i.find('</div>')
	contigs_trunc.append(i[coord1:coord2])

coord1 = contigs[1].find('<div class="contig_sequence">')
coord2 = contigs[1].find('</div>')

def cleaner1(x):
	x = x.replace('<span class="red" "bold">', '')
	x = x.replace('</span>', '')
	x = x.replace('<div class="contig_sequence">', '')
	x = x.replace('[', '')
	x = x.replace(']', '')
	return x

contigs_trunc2 = map(cleaner1, contigs_trunc)

