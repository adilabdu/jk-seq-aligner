class DataProcess:

    def __init__(self, referenceFile, queryFile, queryType='single'):
        self.referenceSequence = ''
        self.referenceHeader = ''
        self.querySequences = ''
        self.queryHeader = ''

        self.referenceFile = open(referenceFile, 'r')
        self.readSequences(self.referenceFile, seq='R')

        self.queryFile = open(queryFile, 'r')
        self.readSequences(self.queryFile, seq='Q')

    def readSequences(self, file, seq='R'):
        header = ''
        sequence = ''

        for r, line in enumerate(file):
            if r != 0:
                sequence += line.rstrip()
            else:
                header = line.rstrip()

        if seq == 'R':
            self.referenceSequence = sequence
            self.referenceHeader = header
        else:
            self.querySequences = sequence
            self.queryHeader = header
