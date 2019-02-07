
from Generator import Generator
from Sequence import Sequence

gn = Generator()

print(gn.generateSequence(10))

sequence = Sequence(gn.getSequence())

print(sequence.sequence)
print(sequence.randomModify())








