
import pybam

pup = pybam.pileup()
pup.open('/Users/amcphers/Scratch/destruct_simulation/bwa.sorted.bam')

pup.jump('20', 100000)

for a in range(10):
    info = pup.next()
    print info

