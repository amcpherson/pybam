
import pybam

pup = pybam.pileup()
pup.open('/Users/amcphers/Scratch/destruct_simulation/bwa.sorted.bam')

for a in range(10):
    info = pup.next()
    print info.position
    print info.bases

