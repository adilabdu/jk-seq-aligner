from Core.DataProcess import DataProcess
from Core.SeqAligner import SeqAligner
from Core.OutputFile import OutputFile

if __name__ == "__main__":

    data = DataProcess('Dataset/reference', 'Dataset/query')

    aligner = SeqAligner(data.referenceSequence, data.querySequences, kMerSize=5, genMers=True)
    overlap = aligner.overlap()
    anchor = aligner.bestAnchor(overlap)

    output = OutputFile('output-test' + '.jres',
                        data.referenceSequence,
                        data.referenceHeader, 5).writeFile(data.queryHeader,
                                                           len(data.querySequences),
                                                           anchor)
