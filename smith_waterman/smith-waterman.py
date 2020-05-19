def alignmentPlot(seq1, seq2, fill_matrix, trace):
	""" This function is optional if you would like a plot of score and trace matrices """
	seq1, seq2 = "3" + seq1[::-1], "5" + seq2
	import numpy as np
	import matplotlib.pyplot as plt
	fig, ax = plt.subplots(figsize=(8,8))
	im = ax.imshow(fill_matrix)
	ax.xaxis.tick_top()
	ax.set_xticks(np.arange(len(seq2)))
	ax.set_yticks(np.arange(len(seq1)))
	ax.set_xticklabels(seq2)
	ax.set_yticklabels(seq1)
	for i in range(0,len(seq1)):
	    for j in range(0,len(seq2)):
	        if [i,j] in trace:
	            text = ax.text(j, i, np.int_(fill_matrix)[i,j], ha="center", va="center", color="r")
	        else:
	            text = ax.text(j, i, np.int_(fill_matrix)[i,j], ha="center", va="center", color="w")
	plt.xlabel("Smith-Waterman Local Alignment")
	fig.tight_layout()
	plt.show()

def alignScore(base1, base2):
    if (base1 == 'G' and base2 == 'C') or (base1 == 'C' and base2 == 'G'):
        return 3
    elif (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A'):
        return 2
    elif (base1 == 'G' and base2 == 'U') or (base1 == 'U' and base2 == 'G'):
        return 1
    else:
        return -1
    
def matrices(seq1, seq2):
    import numpy as np
    seq1, seq2 = "3"+seq1, "5"+seq2
    n, m = len(seq1), len(seq2)
    scores, traces = (np.zeros([n,m]), np.zeros([n,m]))
    for i in range(1, n):
        for j in range(1, m):
            diagonal = scores[i-1,j-1] + alignScore(seq1[i], seq2[j])
            left = scores[i,j-1] - 1
            up = scores[i-1,j] - 1
            scores[i,j] = max(0, left, up, diagonal)
            whereFrom = {1:diagonal, 3:left, 5:up, 0:0}
            traces[i,j] = sum([i for i in whereFrom.keys() if whereFrom[i] == max(whereFrom.values())])
    return scores, traces
    
def trace(seq1,seq2,i,j,store,ends,tM,sM):
    while i>0 and j>0 and sM[i,j] != 0:
        value=tM[i,j]
        if value == 1:         
            i, j = i-1,j-1
            store = [(seq1[i]+s1, seq2[j]+s2) for s1, s2 in store]
        elif value == 3:
            j = j-1
            store = [("-"+s1, seq2[j]+s2) for s1, s2 in store]
        elif value == 5:
            i = i-1
            store = [(seq1[i]+s1, "-"+s2) for s1, s2 in store]
        elif value == 6:
            store_diag = [(seq1[i-1]+s1, seq2[j-1]+s2) for s1, s2 in store]
            store_up = [(seq1[i-1]+s1, "-"+s2) for s1, s2 in store]
            return trace(seq1,seq2,i-1,j-1,store_diag,ends,tM,sM),trace(seq1,seq2,i-1,j,store_up,ends,tM,sM)
        elif value == 4:
            store_diag = [(seq1[i-1]+s1, seq2[j-1]+s2) for s1, s2 in store]
            store_left = [("-"+s1, seq2[j-1]+s2) for s1, s2 in store]
            return trace(seq1,seq2,i-1,j-1,store_diag,ends,tM,sM),trace(seq1,seq2,i,j-1,store_left,ends,tM,sM)
        elif value == 8:
            store_left = [("-"+s1, seq2[j-1]+s2) for s1, s2 in store]
            store_up = [(seq1[i-1]+s1, "-"+s2) for s1, s2 in store]
            return trace(seq1,seq2,i,j-1,store_left,ends,tM,sM),trace(seq1,seq2,i,j-1,store_up,ends,tM,sM)
        elif value == 9:
            store_diag = [(seq1[i-1]+s1, seq2[j-1]+s2) for s1, s2 in store]
            store_left = [("-"+s1, seq2[j-1]+s2) for s1, s2 in store]
            store_up = [(seq1[i-1]+s1, "-"+s2) for s1, s2 in store]
            return trace(seq1,seq2,i-1,j-1,store_diag,ends,tM,sM),trace(seq1,seq2,i,j-1,store_left,ends,tM,sM),trace(seq1,seq2,i,j-1,store_up,ends,tM,sM)
    starts = [len(seq1)-i,j+1]
    return store,starts,ends


def flatten(alist):
    import collections.abc
    for l in alist:
        if isinstance(l, collections.abc.Iterable) and not isinstance(l, (str, bytes)):
            yield from flatten(l)
        else:
            yield l
            
def symbolGenerator(seq1, seq2):
    symbols = ""
    for i in range(len(seq1)):
        base1 = seq1[i]
        base2 = seq2[i]
        if (base1 == 'G' and base2 == 'C') or (base1 == 'C' and base2 == 'G'):
            symbols = symbols + "|"
        elif (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A'):
            symbols = symbols + "|"
        elif (base1 == 'G' and base2 == 'U') or (base1 == 'U' and base2 == 'G'):
            symbols = symbols + ":"
        else:
            symbols = symbols + " "
    return symbols
    
def formatOutput(store, max_score, file):
    for i in range(int(len(store)/6)):
        seq3 = store[6*i + 0]
        seq5 = store[6*i + 1]
        symbols = symbolGenerator(seq3, seq5)
        start3 = store[6*i + 2]
        start5 = store[6*i + 3]
        end3 = store[6*i + 4]
        end5 = store[6*i + 5]
        print(f"Alignment {i+1}:", file=file)
        print("%-7.f %s %7.f" % (start3, "3\'  "+seq3+"  5\'", end3), file=file)
        print("%-7s %s %7s" % ("","    "+symbols+"    ",""),file=file)
        print("%-7.f %s %7.f" % (start5, "5\'  "+seq5+"  3\'", end5),file=file)
        print("",file=file)

def readFasta(file):
    """ Removes all comments (#) and sequence headers (>) and any white space """
    fasta = [line.rstrip() for line in open(file, "r") if ">" not in line and "#" not in line]
    fasta = list(filter(None, fasta))
    return fasta

if __name__ == '__main__':
	# Read in fasta and produce list of sequences where each test case can be found at
	# indices of 2n,2n+1 for n = 0,1,2,...... e.g [seqT1_1,seqT1_2,seqT2_1, seqT2_2.....]
	import os
	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
	fasta = [fasta.rstrip() for fasta in os.listdir() if ".fasta" in fasta]
	print("()@#$&)@#$)(@#$)@#$@#)$@#*$&*@#$@#_$&@#*$&@$(@#)$*@#)$_&(*@#)*$&)_$&#_@#")
	print("()@#$&)@#$)(@#$)@#$@#)$@#*$&*@#$@#_$&@#*$&@$(@#)$*@#)$_&(*@#)*$&)_$&#_@#")
	print("()@#               @#)$@                   $(@#)                       #")
	print("()@#               @#)$@                   $(@#)                       #")
	print("()@#$&)@      $)@#$@#)$@     @#$@#_$&@#    $(@#)                       #")
	print("()@#$&)@      $)@#$@#)$@     @#$@#_$&@#    $(@#)     $_&     *$&)_     #")
	print("()@#$&)@      $)@#$@#)$@     @#$@#_$&@#    $(@#)     $_&     *$&)_     #")
	print("()@#$&)@      $)@#$@#)$@     @#$@#_$&@#    $(@#)     $_&     *$&)_     #")
	print("()@#$&)@      $)@#$@#)$@     @#$@#_$&@#    $(@#)     $_&     *$&)_     #")
	print("()@#$&)@      $)@#$@#)$@                   $(@#)     $_&     *$&)_     #")
	print("()@#$&)@      $)@#$@#)$@                   $(@#)     $_&     *$&)_     #")
	print("()@#$&)@#$)(@#$)@#$@#)$@#*$&*@#$@#_$&@#*$&@$(@#)$*@#)$_&(*@#)*$&)_$&#_@#")
	print("()@#$&)@#$)(@#$)@#$@#)$@#*$&*@#$@#_$&@#*$&@$(@#)$*@#)$_&(*@#)*$&)_$&#_@#")
	print("========================================================================")
	print("")
	print("          Smith-Waterman Implementation with recursive traceback        ")
	print("                     modified for RNA hybridization                     ")
	print("")
	print("========================================================================")
	import timeit
	start = timeit.default_timer()
	# Compute alignments for all fasta files in directory
	for f in fasta:
	    test_cases = readFasta(f)
	    n_cases = int(len(test_cases)/2)
	    out = open("".join([f.split(".fasta")[0],"_output.txt"]), "w")
	    for n in range(n_cases):
	        import numpy as np
	        seq1 = test_cases[2*n][::-1]
	        seq2 = test_cases[2*n+1]
	        scoreMatrix, traceMatrix = matrices(seq1, seq2)
	        max_positions = np.where(scoreMatrix == np.amax(scoreMatrix))
	        alignments = []
	        print("====================================", file=out)
	        print(f"            Test case {n+1}", file=out)
	        print(f"            Score = {int(np.amax(scoreMatrix))}", file=out)
	        print("====================================", file=out)
	        print("", file=out)
	        for pos in range(len(max_positions[0])):
	            i = max_positions[0][pos]
	            j = max_positions[1][pos]
	            store = [("","")]
	            starts = [(i,j)]
	            ends=[len(seq1)-i+1,j]
	            traced = trace(seq1, seq2, i, j, store, ends,traceMatrix, scoreMatrix)
	            traced = list(flatten(traced))
	            alignments = alignments + traced
	        formatOutput(alignments, np.amax(scoreMatrix),out)
	    out.close()
	stop = timeit.default_timer()
	
	print("")
	print("              Wow, your script finished in", stop - start, "seconds!    ")
	print("")
	print("========================================================================")
