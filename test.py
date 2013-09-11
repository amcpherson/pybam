
import pybam

fa = pybam.fasta()

fa.open('/Users/amcphers/Scratch/Homo_sapiens.GRCh37.71.dna.chromosome.fa')

print fa.get('1', 100000)

pup = pybam.pileup()
pup.open('/Users/amcphers/Scratch/patient7/bam.bylibrary-normal_blood.hg19.sorted.bam')

print pup.refnames

pup.jump('20', 100170)

p = pup.next()

print pup.refnames[p[12]], p

pup.jump('13', 100170)

p = pup.next()

print pup.refnames[p[12]], p

