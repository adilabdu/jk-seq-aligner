from typing import List

from Core.Node import Node
from Core.Overlap import Overlap
from Core.Anchor import Anchor


class SeqAligner:

    def __init__(self, referenceSequence, querySequences, kMerSize, genMers=True):
        self.referenceSequence = referenceSequence
        self.querySequences = querySequences
        self.kMerSize = kMerSize
        self.referenceIndex = 0
        self.queryIndex = 0
        self.references = []
        self.queries = []
        self.overlaps = []
        self.anchors = []

        # to generate k-mers on initiate
        if genMers:
            self.kMer(self.references, self.referenceSequence, self.referenceIndex, seq='R')
            self.kMer(self.queries, self.querySequences, self.queryIndex, seq='Q')

    def kMer(self, kList, sequence, index, seq='R'):
        for i in self.kMerConstruct(sequence, self.kMerSize):
            kList.append(Node(i, index))
            index += 1

        if seq == 'R':
            self.referenceIndex = index
        else:
            self.queryIndex = index

    def overlap(self):

        currentOverlap = Overlap(self.kMerSize)
        startOfOverlap = [0, 0]
        rIndex = 0

        while rIndex < len(self.references):
            qIndex = 0
            while qIndex < len(self.queries):

                if self.queries[qIndex].key == self.references[rIndex].key:

                    if currentOverlap.kMer == '':  # if new overlap

                        startOfOverlap[0] = rIndex
                        startOfOverlap[1] = qIndex

                        currentOverlap = Overlap(self.kMerSize)
                        currentOverlap.kMer = self.references[rIndex].key
                        currentOverlap.setQueryIndex(qIndex, 0)
                        currentOverlap.setReferenceIndex(rIndex, 0)

                        currentOverlap.setQueryIndex(qIndex, 1)
                        currentOverlap.setReferenceIndex(rIndex, 1)

                        self.queries[qIndex].notVisited = False
                        self.references[rIndex].notVisited = False

                    else:  # if appending to overlap

                        currentOverlap.kMer += self.references[rIndex].key[-1]
                        currentOverlap.length += 1
                        currentOverlap.setQueryIndex(qIndex, 1)
                        currentOverlap.setReferenceIndex(rIndex, 1)

                        self.queries[qIndex].notVisited = False
                        self.references[rIndex].notVisited = False

                    if rIndex + 1 < self.referenceIndex and qIndex + 1 < self.queryIndex:
                        rIndex = rIndex + 1
                        qIndex = qIndex + 1

                    else:
                        if currentOverlap.length >= self.kMerSize:
                            self.overlaps.append(currentOverlap)
                        currentOverlap = Overlap(self.kMerSize)
                        rIndex = startOfOverlap[0]
                        qIndex = startOfOverlap[1] + 1

                else:  # if overlap NOT found

                    if currentOverlap.kMer == '':  # if overlap hasn't occurred
                        qIndex += 1

                    else:  # if overlap had occurred
                        if currentOverlap.length >= self.kMerSize:
                            self.overlaps.append(currentOverlap)
                        currentOverlap = Overlap(self.kMerSize)
                        rIndex = startOfOverlap[0]
                        qIndex = startOfOverlap[1] + 1

            if currentOverlap.kMer != '':
                if currentOverlap.length >= self.kMerSize:
                    self.overlaps.append(currentOverlap)
                currentOverlap = Overlap(self.kMerSize)
            rIndex += 1

        self.overlaps.sort(key=lambda x: x.referenceIndex[1], reverse=False)
        self.overlaps.sort(key=lambda x: x.queryIndex[1], reverse=False)

        index = 0
        while index < len(self.overlaps) - 1:
            if self.overlaps[index].referenceIndex[1] == self.overlaps[index + 1].referenceIndex[1]:
                if self.overlaps[index].queryIndex[1] == self.overlaps[index + 1].queryIndex[1]:
                    del self.overlaps[index + 1]
                    continue
            index += 1

        self.overlaps.sort(key=lambda x: x.referenceIndex[0], reverse=False)

        return self.overlaps

    def appropriateOverlaps(self, overlap):

        listOfOverlaps = []

        for over in self.overlaps:
            if over.referenceIndex[0] > overlap.referenceIndex[1] and over.queryIndex[0] > overlap.queryIndex[1]:
                listOfOverlaps.append(over)
                over.setAsChild()

        return listOfOverlaps

    def anchor(self, listOfOverlaps: List[Overlap], currentAnchor: Anchor, lastFlag=False):

        if listOfOverlaps:
            for ove in listOfOverlaps:

                if not currentAnchor.overlaps:      # Are we starting a new Anchor?
                    if ove.isChild():               # If so, is the currently selected Overlap a child?
                        continue                    # Then, immediately iterate to the new Overlap

                currentAnchor.append(ove)                           # Append this Overlap to the current Anchor
                potentialAnchors = self.appropriateOverlaps(ove)    # Generate the Potential Anchoring Overlaps

                if ove == listOfOverlaps[-1]:
                    self.anchor(potentialAnchors, currentAnchor, True)
                else:
                    self.anchor(potentialAnchors, currentAnchor)     # Recall the function with updated parameters

        else:   # If there aren't Potential Anchoring Overlaps, go here
            fullAnchor = Anchor()                               # Initiate a fresh Anchor Object
            fullAnchor.overlaps = list(currentAnchor.overlaps)  # Copy currentAnchor's overlaps to the new Object
            self.anchors.append(fullAnchor)                     # Dump the new object into the final Anchors list

            currentAnchor.overlaps.pop()
            if lastFlag and currentAnchor.overlaps:
                currentAnchor.overlaps.pop()

            return

    def anchorScore(self, listOfOverlaps, currentAnchor):

        self.anchor(listOfOverlaps, currentAnchor)
        bestAnchor = Anchor()

        longestAnchor = 0
        for anchor in self.anchors:
            anchorLength = self.kMerSize
            for overlap in anchor.overlaps:
                anchorLength += (len(overlap.kMer) - self.kMerSize) * self.kMerSize
            if longestAnchor == 0:
                longestAnchor = anchorLength
                bestAnchor = anchor
            elif anchorLength > longestAnchor:
                longestAnchor = anchorLength
                bestAnchor = anchor

        return bestAnchor

    def bestAnchor(self, overlap):
        if overlap:
            return self.anchorScore(overlap, Anchor())
        else:
            self.queries.clear()
            self.anchors.clear()
            # index += 1
        return

    @staticmethod
    def kMerConstruct(sequence, kMerSize):
        for k in range(0, len(sequence) - (kMerSize - 1)):
            yield sequence[k:k + kMerSize]
