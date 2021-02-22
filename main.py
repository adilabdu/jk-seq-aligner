from Core.DataProcess import DataProcess
from Core.SeqAligner import SeqAligner
from Core.OutputFile import OutputFile

if __name__ == "__main__":

    data = DataProcess('Dataset/reference', 'Dataset/query')

    output = OutputFile('Final/output-test-multi' + '.jres',
                        data.referenceSequence,
                        data.referenceHeader, 5)

    aligner = SeqAligner(kMerSize=5)
    kMerRef = aligner.kMer([], data.referenceSequence, 0, seq='R')

    queryList = [
        'ACGTAGGTCCTACGGTCCAAAGAT',
        'TTTTTTTTTTTT',
        'ACGTACGTCCGGTCCACGTAAAAGAT',
        'ACGTACGTCCGGTCCATTAGATATAAATATA',
    ]

    for query in queryList:
        kMerQuery = aligner.kMer([], query, 0, seq='Q')
        overlap = aligner.overlap(kMerRef['kMers'], kMerRef['index'], kMerQuery['kMers'], kMerQuery['index'])
        anchor = aligner.bestAnchor(overlap)

        output.writeFile(data.queryHeader, len(data.querySequences), anchor)
        aligner.clear()
        print(overlap)
