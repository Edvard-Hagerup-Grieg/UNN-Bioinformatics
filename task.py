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
def Generating_Theoretical_Spectrum_Problem(Peptide):
    dictionary = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

    elem = [Peptide]
    Cyclic_peptide = Peptide * 2
    for i in range(len(Peptide)): #step
        for j in range(1,len(Peptide)): #shift
            elem.append(Cyclic_peptide[i : i + j])

    answer = [0]
    for i in range(len(elem)):
        sum = 0
        for j in list(elem[i]):
            sum += dictionary[j]
        answer.append(sum)
    answer.sort()

    return(answer)