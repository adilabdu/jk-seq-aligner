from typing import List
import sys
import time


class Node(object):

    def __init__(self, key, data=None, sequence='R'):
        self.key = key
        self.children = {}
        self.notVisited = True

    def __repr__(self):
        # return "[" + str(self.data) + ", " + str(self.key) + ", " + str(self.sequence) + ", " + str(
        #     self.notVisited) + "]"
        return self.key


class Overlap(object):

    def __init__(self, length):
        self.length = length
        self.queryIndex = [None] * 2
        self.kMer = ''
        self.referenceIndex = [None] * 2
        self.anchorChild = False

    def setReferenceIndex(self, referenceIndex, position):
        self.referenceIndex[position] = referenceIndex

    def setQueryIndex(self, queryIndex, position):
        self.queryIndex[position] = queryIndex

    def setKMer(self, kMer):
        self.kMer = kMer

    def joinWith(self, overlap):
        self.queryIndex[1] = overlap.queryIndex[1]
        self.referenceIndex[1] = overlap.referenceIndex[1]
        self.kMer += overlap.kMer[2:]

    def setAsChild(self):
        self.anchorChild = True

    def isChild(self):
        return self.anchorChild

    def __repr__(self):
        # return 'KMer: ' + self.KMer + \
        #        ' Query: ' + str(self.queryIndex[0]) + ', ' + str(self.queryIndex[1]) + \
        #        ' Reference: ' + str(self.referenceIndex[0]) + ', ' + str(self.referenceIndex[1])
        return self.kMer


class Anchor:

    def __init__(self):
        self.overlaps = []
        self.counter = 0

    def append(self, o):
        self.overlaps.append(o)

    def __repr__(self):
        return str(self.overlaps)


def kMerConstruct(sequence, kMerSize):
    for k in range(0, len(sequence) - (kMerSize - 1)):
        yield sequence[k:k + kMerSize]


def anchor(listOfOverlaps: List[Overlap], currentAnchor: Anchor, anc, lastFlag=False):

    if listOfOverlaps:
        for over in listOfOverlaps:

            if not currentAnchor.overlaps:      # Are we starting a new Anchor?
                if over.isChild():              # If so, is the currently selected Overlap a child?
                    continue                    # Then, immediately iterate to the new Overlap

            currentAnchor.append(over)                           # Append this Overlap to the current Anchor
            potentialAnchors = appropriateOverlaps(over)         # Generate the Potential Anchoring Overlaps

            if over == listOfOverlaps[-1]:
                anchor(potentialAnchors, currentAnchor, anc, True)
            else:
                anchor(potentialAnchors, currentAnchor, anc)     # Recall the function with updated parameters

    else:   # If there aren't Potential Anchoring Overlaps, go here
        fullAnchor = Anchor()                               # Initiate a fresh Anchor Object
        fullAnchor.overlaps = list(currentAnchor.overlaps)  # Copy currentAnchor's overlaps to the new Object
        anc.append(fullAnchor)                              # Dump the new object into the final Anchors list

        currentAnchor.overlaps.pop()
        if lastFlag and currentAnchor.overlaps:
            currentAnchor.overlaps.pop()

        return anc


def appropriateOverlaps(appOverlap):

    listOfOvers = []

    for over in overlaps:
        if over.referenceIndex[0] > appOverlap.referenceIndex[1] and over.queryIndex[0] > appOverlap.queryIndex[1]:
            listOfOvers.append(over)
            over.setAsChild()

    return listOfOvers


if __name__ == '__main__':

    REFERENCE_FILE = str(sys.argv[1])
    QUERY_FILE = str(sys.argv[2])
    KMER_SIZE = int(sys.argv[3])
    OUTPUT_FILE = str(sys.argv[4])

    TIME_START = time.time()

    # output file
    file = open('Final/' + OUTPUT_FILE + '.jres', 'w')
    fileIndex = 0

    # read Reference
    referenceSequence = ''
    referenceHeader = ''
    referenceFile = open(REFERENCE_FILE, 'r')
    for r, line in enumerate(referenceFile):
        if r != 0:
            referenceSequence += line.rstrip()
        else:
            referenceHeader = line.rstrip()

    # read Query
    sequences = ''
    queryHeaders = []
    querySequences = []
    queryFile = open(QUERY_FILE, 'r')
    for q, line in enumerate(queryFile):
        if line[0] != '>':
            sequences += line.rstrip()
        else:
            queryHeaders.append(line.rstrip())

            if sequences != '':
                querySequences.append(sequences)
                sequences = ''

    if sequences != '':
        querySequences.append(sequences)

    # kMer Reference
    kMerRef = []
    referenceIndex = 0
    for i, kMer in enumerate(kMerConstruct(referenceSequence, KMER_SIZE)):
        kMerRef.append(Node(kMer, i))
        referenceIndex = i

    # iterate through each query sequence
    for q, query in enumerate(querySequences):

        # kMer Query
        kMerQuery = []
        queryIndex = 0
        for i, kMer in enumerate(kMerConstruct(query, KMER_SIZE)):
            kMerQuery.append(Node(kMer, i))
            queryIndex = i

        # overlap
        overlaps = []
        currentOverlap = Overlap(KMER_SIZE)
        startOfOverlap = [0, 0]
        rIndex = 0

        while rIndex < len(kMerRef):
            qIndex = 0
            while qIndex < len(kMerQuery):

                if kMerQuery[qIndex].key == kMerRef[rIndex].key:

                    if currentOverlap.kMer == '':  # if new overlap

                        startOfOverlap[0] = rIndex
                        startOfOverlap[1] = qIndex

                        currentOverlap = Overlap(KMER_SIZE)
                        currentOverlap.kMer = kMerRef[rIndex].key
                        currentOverlap.setQueryIndex(qIndex, 0)
                        currentOverlap.setReferenceIndex(rIndex, 0)

                        currentOverlap.setQueryIndex(qIndex, 1)
                        currentOverlap.setReferenceIndex(rIndex, 1)

                        kMerQuery[qIndex].notVisited = False
                        kMerRef[rIndex].notVisited = False

                    else:  # if appending to overlap

                        currentOverlap.kMer += kMerRef[rIndex].key[-1]
                        currentOverlap.length += 1
                        currentOverlap.setQueryIndex(qIndex, 1)
                        currentOverlap.setReferenceIndex(rIndex, 1)

                        kMerQuery[qIndex].notVisited = False
                        kMerRef[rIndex].notVisited = False

                    if rIndex + 1 < referenceIndex and qIndex + 1 < queryIndex:
                        rIndex = rIndex + 1
                        qIndex = qIndex + 1

                    else:
                        if currentOverlap.length >= KMER_SIZE:
                            overlaps.append(currentOverlap)
                        currentOverlap = Overlap(KMER_SIZE)
                        rIndex = startOfOverlap[0]
                        qIndex = startOfOverlap[1] + 1

                else:  # if overlap NOT found

                    if currentOverlap.kMer == '':  # if overlap hasn't occurred
                        qIndex += 1

                    else:  # if overlap had occurred
                        if currentOverlap.length >= KMER_SIZE:
                            overlaps.append(currentOverlap)
                        currentOverlap = Overlap(KMER_SIZE)
                        rIndex = startOfOverlap[0]
                        qIndex = startOfOverlap[1] + 1

            if currentOverlap.kMer != '':
                if currentOverlap.length >= KMER_SIZE:
                    overlaps.append(currentOverlap)
                currentOverlap = Overlap(KMER_SIZE)
            rIndex += 1

        overlaps.sort(key=lambda x: x.referenceIndex[1], reverse=False)
        overlaps.sort(key=lambda x: x.queryIndex[1], reverse=False)

        index = 0
        while index < len(overlaps) - 1:
            if overlaps[index].referenceIndex[1] == overlaps[index + 1].referenceIndex[1]:
                if overlaps[index].queryIndex[1] == overlaps[index + 1].queryIndex[1]:
                    del overlaps[index + 1]
                    continue
            index += 1

        overlaps.sort(key=lambda x: x.referenceIndex[0], reverse=False)

        # anchor
        anchors = []
        bestAnchor = None
        if overlaps:
            anchor(overlaps, Anchor(), anchors)

            bestAnchor = Anchor()
            longestAnchor = 0
            for a in anchors:
                anchorLength = KMER_SIZE
                for overlap in a.overlaps:
                    anchorLength += (len(overlap.kMer) - KMER_SIZE) * KMER_SIZE
                if longestAnchor == 0:
                    longestAnchor = anchorLength
                    bestAnchor = a
                elif anchorLength > longestAnchor:
                    longestAnchor = anchorLength
                    bestAnchor = a

            anchors = []

        else:
            kMerQuery = []

        # output result
        file.write('Query: ' + repr(fileIndex) + '\n')
        file.write('Reference: ' + repr(referenceHeader) + '\n')
        file.write('Query:     ' + repr(queryHeaders[q]) + '\n')
        file.write('\n')
        file.write(referenceSequence)
        file.write('\n')

        anchorWithMismatch = ['-'] * (len(referenceSequence))
        symbols = [' '] * (len(referenceSequence))

        matchedFirst = -1
        matchedLast = 0
        if bestAnchor is not None:
            for overlap in bestAnchor.overlaps:

                matchedLast = overlap.referenceIndex[0] + len(overlap.kMer)
                if matchedFirst == -1:
                    matchedFirst = overlap.referenceIndex[0]

                symbols[overlap.referenceIndex[0]:overlap.referenceIndex[0] + len(overlap.kMer) + 1] \
                    = ['|'] * (overlap.referenceIndex[1] - overlap.referenceIndex[0] + KMER_SIZE)
                anchorWithMismatch[overlap.referenceIndex[0]:overlap.referenceIndex[0] + len(overlap.kMer) + 1] \
                    = overlap.kMer[:]

                if matchedFirst != -1:
                    anchorWithMismatch[:matchedFirst] = [' '] * matchedFirst

                anchorWithMismatch[matchedLast:] = [' '] * (len(query) - (len(query) - matchedLast))

        for symbol in symbols:
            file.write(symbol)
        file.write('\n')

        match = 0
        for element in anchorWithMismatch:
            if element == '-' or element == ' ':
                continue
            match += 1

        mismatch = len(query) - match
        aligned = round((match / (mismatch + match)) * 100, 2)

        for char in anchorWithMismatch:
            file.write(char)

        file.write('\n')
        file.write('\n')

        if bestAnchor is not None:
            for ove in bestAnchor.overlaps:
                file.write('Reference: [' + repr(ove.referenceIndex[0]) + ', ' + repr(
                    ove.referenceIndex[1] + KMER_SIZE - 1) + '] and Query: [' + repr(
                    ove.queryIndex[0]) + ', ' + repr(ove.queryIndex[1] + KMER_SIZE - 1) + ']\n')
        file.write('\n')
        file.write('Match: ' + repr(match) + ' | ' + 'Mismatch: ' + repr(mismatch) + '\n')
        file.write('Aligned: ' + repr(aligned) + '(%)' + '\n')
        file.write('\n')
        file.write('\n')
        fileIndex += 1

    print('Total Execution Time: ' + str(time.time() - TIME_START) + ' sec.')
