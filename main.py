from Core.DataProcess import DataProcess
from Core.SeqAligner import SeqAligner
from Core.OutputFile import OutputFile

if __name__ == "__main__":

    data = DataProcess('Dataset/reference', 'Dataset/query')

    output = OutputFile('Final/output-test-multi' + '.jres',
                        data.referenceSequence,
                        data.referenceHeader, kMerSize=5)

    aligner = SeqAligner(kMerSize=5)
    kMerRef = aligner.kMer([], data.referenceSequence, 0, seq='R')

    for q, query in enumerate(data.querySequences):
        kMerQuery = aligner.kMer([], query, 0, seq='Q')
        overlap = aligner.overlap(kMerRef['kMers'], kMerRef['index'], kMerQuery['kMers'], kMerQuery['index'])
        anchor = aligner.bestAnchor(overlap)

        output.writeFile(data.queryHeaders[q], len(data.querySequences[q]), anchor)
        aligner.clear()
