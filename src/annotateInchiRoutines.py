# Porthmeus
# 16.06.23

from src.matchInchi import matchInchi

def findOptimalInchi(inchis:list[str]) -> str:
    ''' 
    Find the "best" inchi string of equivalently annotated inchi strings.

    Sometimes it happens that the same molecule gets annotated with more than one inchi string. We need to handle this in some way. Usually it is just some form of isomer, which is not a big deal to take care of, but in case the Inchis are fundamentally different here we set following rules. 
    1. Use the inchi which has most trues after src.matchInchi.matchInchi()
    2. If there are equals - use the most complex InChI string, that is the one with the most "/" in the string
    3. If this number is also equal, simply use the longest inchi string
    4. If this is also equal, use the first entry by chance
    '''
    
    # rule 1
    matches = []
    for i in range(len(inchis)):
        k = 0
        for j in range(len(inchis)):
            if j != i:
                k = k+ int(matchInchi(inchis[i],inchis[j])[0])
        matches.append(k)

    max_match = max(matches)
    inchis = [x for x,y  in zip(inchis, matches) if y == max_match]

    # rule 2
    counts = []
    for inchi in inchis:
        counts.append(inchi.count("/"))
    max_count = max(counts)
    inchis = [x for x,y  in zip(inchis, counts) if y == max_count]

    # rule 3
    counts = []
    for inchi in inchis:
        counts.append(len(inchi))
    max_count = max(counts)
    inchis = [x for x,y  in zip(inchis, counts) if y == max_count]

    # rule 4
    
    return(inchis[0])
