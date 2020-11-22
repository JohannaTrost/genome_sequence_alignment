import numpy as np
import argparse
import matplotlib.pyplot as plt


def align(seq1, seq2, gap_score, extcn_score, sub_matrix, local):
    """
    Aligning seq2 to seq1 using the specified or default score values and the recurrence function according to the
    alignment strategy which is local or global. Filling the score and traceback matrix.

    :param string seq1:
        A string containing a nucleotide or protein sequence
    :param string seq2:
        Like seq1
    :param int gap_score:
        Score for a gap 'opening' in case of deletion or insertion or a nucleotide
    :param int extcn_score:
        Score for a gap following an other gap
    :param dictionary sub_matrix:
        2D dictionary containing scores (int) for nucleotide matches and mismatches
    :param boolean local:
        True if the alignment is local False if it is global
    :return: tuple score_matrix, trcback_matrix:
        Tuple with 2 2D numpy arrays containing all scores and all directions (up, left, diag)
    """

    def recurrence():
        """
        this recurrence function computes
                         | D(i-1, j-1) + s(xi,xj)
            D(i, j) = max| D(i-1, j) + g
                         | D(i, j-1) + g
                         | 0   (only for local alignment)
            where D is the score matrix,
            s is the score for xi and xj and
            g is the score for gap extensions or openings
        :return: tuple candidates[max_score_pos], directions[max_score_pos]:
            Tuple with optimal score for cell (i, j) and a string containing the direction according to type of
            alignment (gap, match/mismatch) at (i, j)
        """

        # score in case of match
        if nt1 in sub_matrix.keys() and nt2 in sub_matrix[nt1].keys():  # matrix might only be half filled
            match = score_matrix[i, j] + int(sub_matrix[nt1][nt2])
        elif nt2 in sub_matrix.keys() and nt1 in sub_matrix[nt2].keys():  # given pair is e.g. not in matrix as 'AG' but as 'GA'
            match = score_matrix[i, j] + int(sub_matrix[nt2][nt1])
        else:
            raise KeyError('{} or {} is missing in the substitution matrix'.format(nt1, nt2))

        # score in case of gap
        if trcback_matrix[i, j + 1] == 'up' or trcback_matrix[i, j + 1] == 'left':  # checking for gap in upper cell
            gap_seq1 = score_matrix[i, j + 1] + extcn_score
        else:  # it is a gap opening
            gap_seq1 = score_matrix[i, j + 1] + gap_score
        if trcback_matrix[i + 1, j] == 'up' or trcback_matrix[i + 1, j] == 'left':  # checking for gap in left cell
            gap_seq2 = score_matrix[i + 1, j] + extcn_score
        else:  # it a gap opening
            gap_seq2 = score_matrix[i + 1, j] + gap_score

        # find the the highest score
        if local:  # smith and waterman
            candidates = np.array([match, gap_seq1, gap_seq2, 0])
            directions = ['diag', 'up', 'left', 'done']
        else:  # global alignment needleman wunsch
            candidates = np.array([match, gap_seq1, gap_seq2])
            directions = ['diag', 'up', 'left']

        max_score_pos = np.argmax(candidates)  # optimal score

        return candidates[max_score_pos], directions[max_score_pos]

    lg_seq1 = len(seq1)
    lg_seq2 = len(seq2)

    # Init score matrix such that first row (seq1) and column (seq2) count only gaps (gap extension)
    score_matrix = init_score_matrix(gap_score, extcn_score, lg_seq1, lg_seq2, local)

    # Init traceback matrix
    trcback_matrix = init_trcback_matrix(lg_seq1, lg_seq2, local)

    # Filling score matrix along with traceback
    for i, nt1 in enumerate(seq2):
        for j, nt2 in enumerate(seq1):
            score, direction = recurrence()
            score_matrix[i + 1, j + 1] = score
            trcback_matrix[i + 1, j + 1] = direction

    return score_matrix, trcback_matrix


def init_score_matrix(gap_score, extcn_score, lg_seq1, lg_seq2, local):
    """
    Filling of the first row and column of the score matrix such that they only contain gap or gap extension scores.
    The rest will be initialized with zeros. The shape of the matrix is lg_seq2+1 x lg_seq1+1. In case of local
    alignment the matrix will only contain zeros.

    :param int gap_score:
        Score for a gap 'opening' in case of deletion or insertion or a nucleotide
    :param int extcn_score:
        Score for a gap following an other gap
    :param int lg_seq1:
        Length of the first sequence
    :param int lg_seq2:
        Length of the second sequence
    :param boolean local:
        True if the alignment is local False if it is global
    :return: array score_matrix:
        The initialized array
    """
    # create score matrix
    score_matrix = np.zeros(shape=(lg_seq2 + 1, lg_seq1 + 1), dtype=int)

    if not local:  # global alignment needleman wunsch
        if extcn_score != 0:
            start = abs(gap_score)
            stop_seq1 = abs(gap_score) + lg_seq1 * abs(extcn_score)
            stop_seq2 = abs(gap_score) + lg_seq2 * abs(extcn_score)
            step = abs(extcn_score)

            score_matrix[0, 1:] = range(start, stop_seq1, step)
            score_matrix[1:, 0] = range(start, stop_seq2, step)
        else:
            score_matrix[0, 1:] = gap_score
            score_matrix[1:, 0] = gap_score

        if gap_score < 0:
            score_matrix[0, 1:] *= -1
            score_matrix[1:, 0] *= -1
    else:  # smith and waterman: first row and column contains 0s
        pass
    return score_matrix


def init_trcback_matrix(lg_seq1, lg_seq2, local):
    """
    Initializing the traceback matrix such that the value at (0,0) is 'done' and the rest of the first row contains
    'left' and the rest of the first column 'up'. In case of local alignment first row and column will contain 'done'.
    The matrix will be of shape lg_seq2+1 x lg_seq1+1. The other posititons will contain null values

    :param boolean local:
        True if the alignment is local False if it is global
    :return: array trcback_matrix:
        The initialized array
    """
    # create traceback matrix
    trcback_matrix = np.full(shape=(lg_seq2 + 1, lg_seq1 + 1), dtype=object, fill_value=None)

    if not local:  # needleman wunsch
        trcback_matrix[0, 0] = 'done'  # at this point we are done
        trcback_matrix[0, 1:] = 'left'  # corresponds to gap (extension)
        trcback_matrix[1:, 0] = 'up'  # corresponds to gap (extension)
    else:  # smith and waterman
        trcback_matrix[0, 0:] = 'done'
        trcback_matrix[0:, 0] = 'done'
    return trcback_matrix


def get_substitution_matrix(filename=None):
    """
    creation of substitution matrix from csv file or default matrix
    :param string filename:
        csv filename (string) or None (default)
    :return: dictionary sub_matrix:
        Dictionary [nuc1][nuc2] : score (int)
    """

    sub_matrix = {}
    if filename is None:  # default substitution matrix as specified in 'UE Bioinformatique BISM
        # Projet Needleman-Wunsch - L3 BISM 2019-2020'
        sub_matrix = {'A': {'A': '2'},
                      'C': {'A': '-1', 'C': '2'},
                      'G': {'A': '1', 'C': '-1', 'G': '2'},
                      'T': {'A': '-1', 'C': '1', 'G': '-1', 'T': '2'}}

    else:  # take matrix from input
        with open(filename) as sub_matrix_csv:
            # read csv file
            lines = sub_matrix_csv.read().splitlines()

            csv_format_check(lines)  # format check

            # init matrix with nucleotids from first line of csv
            nucleotids = lines[0].split(',')[1:]

            for nt in nucleotids:
                sub_matrix[nt] = {}
            try:
                # fill matrix of substitution
                for line in lines[1:]:
                    line = line.strip()
                    line = clean(line.split(','))
                    key = line[0]
                    for i, score in enumerate(line[1:]):
                        sub_matrix[key][nucleotids[i]] = int(score)
            except KeyError:
                print('The substitution matrix file is not correctly formatted, check the README for correctly '
                      'formatted examples')
                exit()
            except ValueError as err:
                print('The substitution matrix file is not correctly formatted, check the README for correctly '
                      'formatted examples.\n', err)
                exit()

    return sub_matrix


def csv_format_check(lines):
    if lines[0][0] != ',':
        raise ValueError('The first character of the csv-formatted substitution matrix file needs to be a '
                         'comma!')
    nb_comma = lines[0].count(',')
    for line in lines:
        if line.count(',') != nb_comma:
            raise ValueError('Error in input csv-file: the amount of commas needs to be the same for every line')



def traceback(i, j, trcback_matrix, seq1, seq2):
    """
    tracing back the 'path' from the first aligned sequences pair untill the last. 'left' is being interpreted as
    moving one cell to the left as well as a gap on the 2. sequence, 'up' means moving one cell up or a gap on the first
    sequence and 'diag' means moving diagonally which is a match if the nucleotides at (i, j) are equal and a mismatch
    otherwise. 'done' indicates the end of the alignment.

    :param int i:
        Index of row in matrix
    :param int j:
        Index of column in matrix
    :param array trcback_matrix:
        Traceback matrix of strings 'up', 'left', 'diag', 'done'
    :param string seq1:
        A string containing a nucleotide or protein sequence
    :param string seq2:
        Like seq1
    :return:
        A list containing the start position of the aligned sequence (i, j) followed by the
        type of alignment for each aligned pair (match, mismatch, insertion or deletion)
    """

    if trcback_matrix[i, j] == 'done':
        return [i, j]  # start positions for sequences
    elif trcback_matrix[i, j] == 'left':
        return traceback(i, j - 1, trcback_matrix, seq1, seq2) + ['gap seq2']
    elif trcback_matrix[i, j] == 'up':
        return traceback(i - 1, j, trcback_matrix, seq1, seq2) + ['gap seq1']
    elif trcback_matrix[i, j] == 'diag' and seq1[j - 1] == seq2[i - 1]:
        return traceback(i - 1, j - 1, trcback_matrix, seq1, seq2) + ['match']
    elif trcback_matrix[i, j] == 'diag' and seq1[j - 1] != seq2[i - 1]:
        return traceback(i - 1, j - 1, trcback_matrix, seq1, seq2) + ['mismatch']
    else:
        pass


def fasta2nuc(filename):
    """
    opens and reads a fasta file to extract is identifier and sequence

    :param string filename:
        String containing the filename (with relative path) for the fasta sequence
    :return:
        Id and nucleotide sequence from file
    """

    if filename.split('.')[1] != 'fasta':
        raise ValueError('wrong file format .' + filename.split('.')[1] + ' => .fasta expected')

    with open(filename) as file:
        id = file.readline()[1:].strip()  # identifier of sequence
        nuc = ''  # sequence of nucleotides
        for line in file:
            if line[0] == '>':  # take only first seq in case file contains more than 1
                break
            nuc += line.strip()
    return id, nuc


def clean(raw_list):
    """
    cleaning up a list such that it has no None, empty string or blank space as element
    :param list raw_list:
        A list of strings
    :return:
        Clean_list which is the input list without None, empty string or blank space
    """

    clean_list = []
    for strng in raw_list:
        if strng is not None and strng != '' and strng != ' ':
            clean_list.append(strng)
    return clean_list


def max_pos(matrix):
    """
    Finds all positions of maximum value of the matrix

    :param array matrix:
        2D matrix of integers
    :return: array positions:
        Matrix containing all indices of positions of max value
    """
    max = np.max(matrix)
    positions = []
    for i, arr in enumerate(matrix):
        for j, val in enumerate(arr):
            if val == max:
                positions.append((i,j))
    return np.asarray(positions)


def show_results(id1, id2, seq1, seq2, final_score, align_types, start_pos, end_pos, local):
    """
    printing results of alignment
    e.g.:
        The Alignment of Seq1 and Seq2
        ______________________________________________________

                CGTAC*
                :||:|
                AGTCCT
        ______________________________________________________

        Length: 6
        Number of matches: 3
        Number of mismatches: 2
        Number of gaps: 1
        ______________________________________________________

        Total Score: -6

    :param string id1:
        Identifier string of sequence 1
    :param string id2:
        Identifier string of sequence 2
    :param string seq1:
        A string containing a nucleotide or protein sequence
    :param string seq2:
        Like seq1
    :param int final_score:
        Optimal score for alignment
    :param list align_types:
        List of strings indicating a match, mismatch or gap
    :param list start_pos:
        List of indices of start position of alignment
    :param list end_pos:
        List of indices of end position of alignment
    :param boolean local:
        True if the alignment is local False if it is global
    """

    # matrix of symbols for alignment
    symbols_align = {'gap seq1': ' ',
                     'gap seq2': ' ',
                     'mismatch': ':',
                     'match': '|'}

    align_symbs = ''.join(map(str, [symbols_align[elem] for elem in align_types]))
    if local:
        seq2 = seq2[start_pos[0]: end_pos[0] + 1]
        seq1 = seq1[start_pos[1]: end_pos[1] + 1]
    for i, type in enumerate(align_types):
        if type == 'gap seq1':
            seq1 = seq1[:i] + '*' + seq1[i:]
        if type == 'gap seq2':
            seq2 = seq2[:i] + '*' + seq2[i:]
    line = '_' * max(int((max(len(seq1), len(seq2)))) + 8, len('The Alignment of {} and {}'.format(id1, id2))) + '\n'
    print('\nThe Alignment of {} and {}'.format(id1, id2))
    print(line)
    print('\t{}\n\t{}\n\t{}'.format(seq1, align_symbs, seq2))
    print(line)
    print('Length: {}'.format(len(seq1)))
    print('Number of matches: {}'.format(align_types.count('match')))
    print('Number of mismatches: {}'.format(align_types.count('mismatch')))
    print('Number of gaps: {}'.format((align_types.count('gap seq1') + align_types.count('gap seq2'))))
    print(line)
    print('Total Score: {}\n'.format(final_score))


def show_matrix(matrix, title, seq1, seq2):
    seq1 = ' ' + seq1
    seq2 = ' ' + seq2
    plt.figure()
    table = plt.table(cellText=matrix,
                          colWidths=[0.1] * matrix.shape[1],
                          colLabels=[' ' + nt + ' ' for nt in seq1],
                          rowLabels=[' ' + nt + ' ' for nt in seq2],
                          rowColours=['lightgrey'] * matrix.shape[0],
                          colColours=['lightgrey'] * matrix.shape[1],
                          loc='center')
    table.auto_set_font_size(False)
    table.scale(1, 2)

    # Removing ticks and spines to get the figure only with table
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)

    plt.title(title, size='20')
    plt.draw()


def main():
    """
    Runs the alignment algorithm (Needlemann Wunsch approach) with possibility to do a
    local alignment (Smith and Waterman). Provides all necessary parameters like sequences to be aligned,
    the substitution matrix, gap and gap extension score according to the users input. Finally, prints the results

    Run in commandline as ``python </path/to/> alignment.py [args]``

    Type ``python </path/to/> alignment.py -h`` to see all available options.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('fasta1', type=str,
                        help='The </path/to/> name of the first fasta file of which the first sequence will be aligned')
    parser.add_argument('fasta2', type=str,
                        help='The </path/to/> name of the first fasta file of which the first sequence will be aligned '\
                             'to the sequence of fasta1')
    parser.add_argument('-s', '--substitution_matrix', type=str,
                        help='csv file <path/to/> name containing the substitution matrix. Note that both a half '\
                             'complete or complete matrix are possible. For protein sequences default is Blosum62 '\
                             'matrix. For nucleotide sequences find default matrix in README under Arguments explained')
    parser.add_argument('-g', '--gap_score', default=-10, help='Score for a gap. Default score is -10', type=int)
    parser.add_argument('-e', '--extension_score', default=-1,
                        help='Score for a gap extension, for a gap directly following after another gap. Default '
                             'score is -1', type=int)
    parser.add_argument('-p', '--protein_alignment', action='store_true',
                        help='2 protein sequences will be aligned. Take into account option -s' )
    parser.add_argument('-l', '--local_alignment', action='store_true',
                        help='local alignment technique using Smith and Waterman approach')
    parser.add_argument('-sm', '--show_matrices', action='store_true',
                        help='flag to plot traceback and score matrix')

    args = parser.parse_args()

    # get the sequences from fasta file inputs
    id1, seq1 = fasta2nuc(args.fasta1)
    id2, seq2 = fasta2nuc(args.fasta2)

    # Get substitution matrix from input
    if args.substitution_matrix:
        csv_sub_matrix = args.substitution_matrix
    else:
        if args.protein_alignment:  # Blosum62 is default for protein alignment
            csv_sub_matrix = 'BLOSUM62.csv'
        else:  # default will be set later
            csv_sub_matrix = None

    sub_matrix = get_substitution_matrix(csv_sub_matrix)  # create dict from csv

    # transform sequences so that they are compatible with substitution matrix entries
    if 'U' not in sub_matrix.keys():
        seq1 = seq1.upper().replace('U', 'T')
        seq2 = seq2.upper().replace('U', 'T')
    else:
        seq1 = seq1.upper()
        seq2 = seq2.upper()

    score_matrix, trcback_matrix = align(seq1, seq2, args.gap_score, args.extension_score, sub_matrix,
                                         args.local_alignment)

    # determine position of the final score in score matrix
    if not args.local_alignment:
        pos_final_score = np.asarray([(trcback_matrix.shape[0] - 1, trcback_matrix.shape[1] - 1)])
    else:
        pos_final_score = max_pos(score_matrix)

    for i, pos in enumerate(pos_final_score):
        # saves for each aligned nucleotide couple if it is a gap for seq1 or seq2, mismatch or match
        trace = traceback(pos[0], pos[1], trcback_matrix, seq1, seq2)
        align_types, start_pos = trace[2:], trace[:2]
        final_score = score_matrix[pos[0], pos[1]]
        # generating output
        if len(pos_final_score) > 0:
            factor = min(int((max(len(seq1), len(seq2))/2)) * 15, 40)
            print('#' * factor + ' ' + str(i+1) + ' ' + '#' * factor + '\n')
        show_results(id1, id2, seq1, seq2, final_score, align_types, start_pos,
                     np.asarray(pos) - 1, args.local_alignment)

    if args.show_matrices:
        show_matrix(score_matrix, 'Score Matrix', seq1, seq2)
        show_matrix(trcback_matrix, 'Traceback Matrix', seq1, seq2)
        plt.show()


if __name__ == '__main__':
    main()
