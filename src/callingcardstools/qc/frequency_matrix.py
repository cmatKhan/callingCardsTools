def frequency_matrix(dna_list):
    """Tally the frequency of each base by position

    :param dna_list: a list of dna sequences of the same length

    this is hand coded because there isn't an easy way to do this in BioPython?
    In R, it is like so --

    library(Biostrings)
    start <- Sys.time()
    x = DNAStringSet(df$srt)
    freq_mat = consensusMatrix(x, baseOnly = TRUE)
    print(freq_mat)
    end <- Sys.time()
    print(paste0("Runtime: ", end-start))

    and takes .004 seconds. Could not find an easy function to do the same
    with biopython -- if one exists consider replace the function below with
    that
    """
    # CITE: https://hplgit.github.io/bioinf-py/doc/pub/html/main_bioinf.html
    n = len(dna_list[0])
    A = [0]*n
    T = [0]*n
    G = [0]*n
    C = [0]*n
    other = [0]*n
    for dna in dna_list:
        try:
            for index, base in enumerate(dna):
                if base == 'A':
                    A[index] += 1
                elif base == 'C':
                    C[index] += 1
                elif base == 'G':
                    G[index] += 1
                elif base == 'T':
                    T[index] += 1
                else:
                    other[index] +=1
        except IndexError:
            print("dna: %s; index: %s; n: %s" %(base, index, n))
        except TypeError as e:
            print("error in DNA list. Value of what was \
                expected to be a sequence/char is:  %s. Error: %s" %(dna, e))
            raise


    # make position frequency matrix
    df = pd.DataFrame({'A':A,'C':C,'G':G,'T':T,'other':other}).transpose()
    # transform to position probability matrix
    df = (df/df.sum(axis=0)*100)

    return(df)
