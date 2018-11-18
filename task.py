#task 1.1
def Pattern_Count_Problem(Genome, Pattern):

    count = 0
    for i in range(len(Genome) - len(Pattern)):
        if Pattern == Genome[i : i+len(Pattern)]:
            count = count + 1

    return(count)

#task 1.2
def Frequent_Words_Problem(Text, k):

    Pattern = []
    Count = []
    for i in range(len(Text) - k):
        pattern = Text[i : i+k]
        count = Pattern_Count_Problem(Text, pattern)

        Pattern.append(pattern)
        Count.append(count)

    answer = []
    MAX = max(Count)
    for i in range(len(Count)):
        if Count[i] == MAX:
            answer.append(Pattern[i])

    return(set(answer))

#task 1.3
def Reverse_Complement_Problem(DNA):
    dictionary = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    RNA = ''
    for i in range(len(DNA)):
        RNA += dictionary[DNA[i]]

    return(RNA[::-1])

#task 2.1
def Protein_Translation_Problem(RNA):
    dictionary = {'AAA':'K','AAC':'N','AAG':'K','AAU':'N','ACA':'T','ACC':'T','ACG':'T','ACU':'T','AGA':'R','AGC':'S','AGG':'R','AGU':'S','AUA':'I','AUC':'I','AUG':'M','AUU':'I','CAA':'Q','CAC':'H','CAG':'Q','CAU':'H','CCA':'P','CCC':'P','CCG':'P','CCU':'P','CGA':'R','CGC':'R','CGG':'R','CGU':'R','CUA':'L','CUC':'L','CUG':'L','CUU':'L','GAA':'E','GAC':'D','GAG':'E','GAU':'D','GCA':'A','GCC':'A','GCG':'A','GCU':'A','GGA':'G','GGC':'G','GGG':'G','GGU':'G','GUA':'V','GUC':'V','GUG':'V','GUU':'V','UAA':'','UAC':'Y','UAG':'','UAU':'Y','UCA':'S','UCC':'S','UCG':'S','UCU':'S','UGA':'','UGC':'C','UGG':'W','UGU':'C','UUA':'L','UUC':'F','UUG':'L','UUU':'F'}

    amino = ''
    for i in range(0, len(RNA) - 2, 3):
        amino += dictionary[RNA[i : i+3]]

    return(amino)

#task 2.2
def Peptide_Encoding_Problem(DNA, amino):

    answer = []

    for i in range(0, len(DNA) - len(amino)*3 + 1):
        dna1 = DNA[i : i + len(amino)*3]
        dna2 = Reverse_Complement_Problem(dna1)

        str1 = dna1.replace('T', 'U')
        str2 = dna2.replace('T', 'U')

        if amino == Protein_Translation_Problem(str1) or amino == Protein_Translation_Problem(str2):
            answer.append(DNA[i: i+len(amino)*3])

    return(answer)

#task 2.3
def Subpeptides_Count_Problem(n):
    return(n * (n - 1))

#task 2.4
def Generating_Theoretical_Spectrum_Problem(Peptide, cyclic=True):
    dictionary = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                  'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

    if cyclic == True:
        elem = [Peptide]
        Cyclic_peptide = Peptide * 2
        for i in range(len(Peptide)): #step
            for j in range(1,len(Peptide)): #shift
                elem.append(Cyclic_peptide[i : i + j])
    else:
        elem = []
        for i in range(len(Peptide)):
            for j in range(1,len(Peptide) - i + 1):
                elem.append(Peptide[i : i + j])

    answer = [0]
    for i in range(len(elem)):
        sum = 0
        for j in list(elem[i]):
            sum += dictionary[j]
        answer.append(sum)
    answer.sort()

    return(answer)

#task 3.1
def Cyclopeptide_Sequencing_Problem(Spectrum):
    Theory = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
                  'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


    Peptides = []
    for Peptide in Theory.keys():
        if Theory[Peptide] in Spectrum:
            Peptides.append(Peptide)

    Peptides_result = []
    while Peptides != []:
        Peptide = Peptides.pop(0)

        for next_peptide in Theory.keys():
            Expand_peptide = Peptide + next_peptide
            Cyclospectrum = Generating_Theoretical_Spectrum_Problem(Expand_peptide)
            Nonsyclospectrum = Generating_Theoretical_Spectrum_Problem(Expand_peptide, False)[:-1]
            if Cyclospectrum[-1] == Spectrum[-1]:
                if Cyclospectrum == Spectrum:
                    Peptides_result.append(Expand_peptide)
            elif [mass for mass in Cyclospectrum[:-1] if mass not in Spectrum] == []:
                for mass in Spectrum:
                    if mass in Nonsyclospectrum:
                        Nonsyclospectrum.remove(mass)
                if Nonsyclospectrum == []:
                    Peptides.append(Expand_peptide)

    result = []
    for Peptide in Peptides_result:
        mass = []
        for c in Peptide:
            mass.append(Theory[c])
        if mass not in result:
            result.append(mass)

#task 3.2
def Cyclopeptide_Scoring_Problem(Peptide, Spectrum, cyclic=True):
    Theory_spectrum = Generating_Theoretical_Spectrum_Problem(Peptide, cyclic=cyclic)

    score = 0
    Same_mass = list(set(Theory_spectrum) & set(Spectrum))
    for mass in Same_mass:
        score += min(Theory_spectrum.count(mass), Spectrum.count(mass))

    return(score)

#task 3.3
def Trim(Leaderboard, Spectrum, N):
    dict_pepide_score = dict()
    for Peptide in Leaderboard:
        dict_pepide_score[Peptide] = Cyclopeptide_Scoring_Problem(Peptide, Spectrum, False)

    sorted_dict_pepide_score = sorted(dict_pepide_score.items(), key=lambda x: x[1], reverse=True)

    num = N
    if num < len(sorted_dict_pepide_score):
        while sorted_dict_pepide_score[num - 1][1] == sorted_dict_pepide_score[num][1]:
            num += 1
            if num == len( sorted_dict_pepide_score):
                break
    else:
        num = len(sorted_dict_pepide_score)

    Leaderboard.clear()
    trim_dict_pepide_score = sorted_dict_pepide_score[:num][:]
    for elem in trim_dict_pepide_score:
        Leaderboard.add(elem[0])

    return Leaderboard

def Expand(Leaderboard):
    dict = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

    result = set()
    for peptide in dict:
        for elem in Leaderboard:
            new_peptide = elem + peptide
            result.add(new_peptide)
    return result

def Leaderboard_Cyclopeptide_Sequencing(Spectrum, N):
    dict = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
    Leaderboard = {''}
    LeaderPeptide = ''
    while len(Leaderboard) !=0:
        Leaderboard = Expand(Leaderboard)
        leaderboard_temp = Leaderboard.copy()
        for peptide in Leaderboard:
            pep_spect = Generating_Theoretical_Spectrum_Problem(peptide, False)
            if pep_spect[-1] == Spectrum[-1]:
                if Cyclopeptide_Scoring_Problem(peptide, Spectrum, False) > Cyclopeptide_Scoring_Problem(LeaderPeptide, Spectrum, False):
                    LeaderPeptide = peptide
            elif int(pep_spect[-1]) > int(Spectrum[-1]):
                leaderboard_temp.remove(peptide)
        leaderboard_temp = Trim(leaderboard_temp, Spectrum, N)
        Leaderboard = leaderboard_temp

    result = []
    for char in LeaderPeptide:
        result.append(dict[char])
    return result

#task 4.1
def Mismatches_Count(Pattern1, Pattern2):
    count = 0
    for i in range(len(Pattern1)):
        if Pattern1[i] != Pattern2[i]:
            count += 1

    return count

def Peptidew_Expand(Peptide):
    Elements = ['A', 'T', 'G', 'C']

    result = set()
    for peptide in Peptide:
        for elem in Elements:
            new_peptide = peptide + elem
            result.add(new_peptide)

    return result

def Patterns_Difference(Peptide1, d):
    Peptides = {''}
    for i in range(len(Peptide1)):
        Peptides = Peptidew_Expand(Peptides)

    result = set()
    for Peptide2 in Peptides:
        if Mismatches_Count(Peptide1, Peptide2) <= d:
            result.add(Peptide2)

    return result

def Motif_Enumeration_Problem(DNA, k, d):
    Patterns = set()
    for str in DNA:
        for i in range(len(str) - k + 1):
            Pattern = str[i:i+k]
            Patterns_new = Patterns_Difference(Pattern, d)
            add = False
            for Pattern_new in Patterns_new:
                for str_check in DNA:
                    add = False
                    for j in range(len(str_check) - k + 1):
                        if Mismatches_Count(Pattern_new, str_check[j:j+k]) <= d:
                            add = True
                    if add == False: break
                if add: Patterns.add(Pattern_new)

#task 4.2
def Hamming_Distance(Pattern, DNA):
    HEMMING_DIST = 0
    # Patterns = []
    for str in range(len(DNA)):
        MIN_SCORE = 100000
        for i in range(len(DNA[0]) - len(Pattern) + 1):
            score = Mismatches_Count(Pattern, DNA[str][i:i+len(Pattern)])
            if score < MIN_SCORE:
                MIN_SCORE = score
                # Patterns.append(DNA[str][i:i+len(Pattern)])
        HEMMING_DIST += MIN_SCORE

    return HEMMING_DIST

def Median_String_Problem(k, DNA):
    distance = 100000
    Median = ''

    Patterns = {''}
    for i in range(k):
        Patterns = Peptidew_Expand(Patterns)

    for Pattern in Patterns:
        distance_new = Hamming_Distance(Pattern, DNA)
        if distance > distance_new:
            distance = distance_new
            Median = Pattern

    return Median