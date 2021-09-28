import numpy as np

def compute(a, b, subDict):
	score = 0
	for i in range(len(a)):
		score += subDict[a[i]][b[i]]
	return score

def matToDict(chars, scoring_matrix):
	subDict = {}
	cList = []
	for i in range(len(chars)):
		cList += chars[i]
	cList += '-'

	for i in range(0, len(cList)):
		for j in range(0, len(cList)):
			if cList[i] not in subDict:
				subDict[cList[i]] = {}
			subDict[cList[i]][cList[j]] = scoring_matrix[i][j]
	return subDict

def findPathExhaustive(sm, sc, a, b, s, t):
	alignments = [["","",[],[]]]
	alignment = ["","",[],[]]
	if sc[a-1][b-1]+sm[s[a]][t[b]] == sc[a][b]:
		path = findPath(sm, sc, a-1, b-1, s, t)
		sAlgn = path[0] + s[a]
		tAlgn = path[1] + t[b]
		sIndices = path[2] + [a]
		tIndices = path[3] + [b]
		alignments += [[sAlgn,tAlgn,sIndices,tIndices]]
	if sc[a-1][b]+sm[s[a]]["-"] == sc[a][b]:
		path = findPath(sm, sc, a-1, b, s, t)
		tAlgn = path[1] + "-"
		sAlgn = path[0] + s[a]
		alignments += [[sAlgn,tAlgn,path[2],path[3]]]
	if sc[a][b-1]+sm["-"][t[b]] == sc[a][b]:
		path = findPath(sm, sc, a, b-1, s, t)
		tAlgn = path[1] + t[b]
		sAlgn = path[0] + "-"
		alignments += [[sAlgn,tAlgn,path[2],path[3]]]
	for a in alignments:
		if len(a[0]) > len(alignment[0]):
			alignment = a
	return alignment

def backtrack(sm, sc, a, b, s, t):
	sAlgn,tAlgn,sIndices,tIndices = "","",[],[]
	while a > 0 or b > 0:
		if a > 0 and b > 0 and sc[a-1][b-1]+sm[s[a]][t[b]] == sc[a][b]:
			sAlgn = s[a] + sAlgn
			tAlgn = t[b] + tAlgn
			sIndices = [a-1] + sIndices
			tIndices = [b-1] + tIndices
			a -= 1
			b -= 1
		elif a > 0 and sc[a-1][b]+sm[s[a]]["-"] == sc[a][b]:
			tAlgn = "-" + tAlgn
			sAlgn = s[a] + sAlgn
			a -= 1
		elif b > 1 and sc[a][b-1]+sm["-"][t[b]] == sc[a][b]:
			tAlgn = t[b] + tAlgn
			sAlgn = "-" + sAlgn
			b -= 1
		else:
			a = 0
			b = 0
	return [sAlgn,tAlgn,sIndices,tIndices]

def SWalign(s,t,subDict):
	s = "-" + s
	t = "-" + t
	maxScore = 0
	position = [0,0]
	score = np.zeros((len(s),len(t)))
	for i in range(1, len(s)):
		for j in range(1, len(t)):
			ms = max(score[i-1][j-1]+subDict[s[i]][t[j]], score[i-1][j]+subDict[s[i]]["-"], score[i][j-1]+subDict["-"][t[j]], 0)
			score[i][j] = ms
			if ms >= maxScore:
				position = [i,j]
				maxScore = ms
	alignment = backtrack(subDict,score,position[0],position[1],s,t)

	return [maxScore, alignment[2], alignment[3]]

def dynprog(alphabet, scoring_matrix, sequence1, sequence2):
	subDict = matToDict(alphabet, scoring_matrix)
	alignment = SWalign(sequence1, sequence2, subDict)
	return alignment

def localScore(s, t, subDict):
	s = "-" + s
	t = "-" + t
	score = np.zeros((2,len(t)))
	startPoints = np.zeros((2,len(t),2))
	start = [0,0]
	maxScore = -1000
	end = [0,0]
	for i in range(len(s)):
		for j in range(len(t)):
			if i == 0 or j == 0:
				score[i%2][j] = 0
				startPoints[i%2][j][0] = i
				startPoints[i%2][j][1] = j
			else:
				diag = score[i%2-1][j-1]+subDict[s[i]][t[j]]
				hor = score[i%2-1][j]+subDict[s[i]]["-"]
				vert = score[i%2][j-1]+subDict["-"][t[j]]
				ms = max(diag, hor, vert, 0)
				
				if ms == diag:
					startPoints[i%2][j][0] = startPoints[i%2-1][j-1][0]
					startPoints[i%2][j][1] = startPoints[i%2-1][j-1][1]
				elif ms == hor:
					startPoints[i%2][j][0] = startPoints[i%2-1][j][0]
					startPoints[i%2][j][1] = startPoints[i%2-1][j][1]
				elif ms == vert:
					startPoints[i%2][j][0] = startPoints[i%2][j-1][0]
					startPoints[i%2][j][1] = startPoints[i%2][j-1][1]
				elif ms == 0:
					startPoints[i%2][j][0] = i
					startPoints[i%2][j][1] = j
				score[i%2][j] = ms
				if ms >= maxScore:
					end = [i,j]
					start = list(startPoints[i%2][j])
					maxScore = ms
	return [maxScore,end,start]

def globalScore(s, t, lastRow, subDict):
	s = "-" + s
	t = "-" + t
	score = np.zeros((2,len(t)))

	for i in range(len(s)):
		for j in range(len(t)):
			if i == 0 and j == 0:
				score[i][j] = 0
			elif i == 0 and j != 0:
				score[i][j] = score[i][j-1] + subDict["-"][t[j]]
			elif i != 0 and j == 0:
				score[i%2][j] = score[i%2-1][j] + subDict[s[i]]["-"]
			elif i > 0 and j > 0:
				score[i%2][j] = max(score[i%2-1][j-1]+subDict[s[i]][t[j]], score[i%2-1][j]+subDict[s[i]]["-"], score[i%2][j-1]+subDict["-"][t[j]])
	if lastRow:
		return score[i%2]
	else:
		return score[i%2][-1]

def NWalign(s, t, subDict):
	s = "-" + s
	t = "-" + t
	score = np.zeros((len(s),len(t)))
	for i in range(len(s)):
		for j in range(len(t)):
			if i == 0 and j == 0:
				score[i][j] = 0
			elif i == 0:
				score[i][j] = score[i][j-1] + subDict["-"][t[j]]
			elif j == 0:
				score[i][j] = score[i-1][j] + subDict[s[i]]["-"]
			else:
				score[i][j] = max(score[i-1][j-1]+subDict[s[i]][t[j]], score[i][j-1]+subDict["-"][t[j]], score[i-1][j]+subDict[s[i]]["-"])
	i,j = len(s)-1,len(t)-1
	S = ""
	T = ""
	cont = True
	while i > 0 or j > 0:
		if i > 0 and j > 0 and score[i-1][j-1]+subDict[s[i]][t[j]] == score[i][j]:
			S = s[i] + S
			T = t[j] + T
			i -= 1
			j -= 1
		elif i > 0 and score[i-1][j]+subDict[s[i]]["-"] == score[i][j]:
			T = "-" + T
			S = s[i] + S
			i -= 1
		elif j > 0 and score[i][j-1]+subDict["-"][t[j]] == score[i][j]:
			T = t[j] + T
			S = "-" + S
			j -= 1
	return [S,T]

def Hirschberg(s, t, subDict):
	if len(s) < 2 or len(t) < 2:
		return NWalign(s,t,subDict)
	else:
		mid = len(s)//2
		h1 = globalScore(s[:mid],t,True,subDict)
		h2 = globalScore(s[mid:][::-1],t[::-1],True,subDict)
		maxScore = globalScore(s,t,False,subDict)
		sumLists = [h1[i] + h2[::-1][i] for i in range (len(h1))]
		split = np.where(sumLists == np.amax(sumLists))[0]
		left = Hirschberg(s[:mid],t[:split[0]],subDict)
		right = Hirschberg(s[mid:],t[split[0]:],subDict)
		return [left[i] + right[i] for i in range(len(left))]

def linear(s, t, subDict):
	scoreData = localScore(s, t, subDict)
	targetScore = scoreData[0]
	start = scoreData[2]
	end = scoreData[1]
	alignment = Hirschberg(s[int(start[0]):end[0]], t[int(start[1]):end[1]], subDict)
	aIndex = len(alignment[0]) - 1
	sPos = scoreData[1][0] - 1
	tPos = scoreData[1][1] - 1
	sIndices = []
	tIndices = []

	while aIndex >= 0:
		if alignment[0][aIndex] != "-" and alignment[1][aIndex] != "-":
			sIndices = [sPos] + sIndices
			tIndices = [tPos] + tIndices
			sPos = sPos - 1
			tPos = tPos - 1
		elif alignment[0][aIndex] == "-":
			tPos = tPos - 1
		elif alignment[1][aIndex] == "-":
			sPos = sPos - 1
		aIndex = aIndex - 1
	
	return [targetScore,sIndices,tIndices]

def dynproglin(alphabet, scoring_matrix, sequence1, sequence2):
	subDict = matToDict(alphabet, scoring_matrix)
	alignment = linear(sequence1, sequence2, subDict)
	return alignment

def FASTA(s,t,subDict,k,w):
	# index all k-length substrings in s
	swapped = False
	if len(s) > len(t):
		s, t = t, s
		swapped = True
	kList = []
	kDict = {}
	seeds = []
	for i in range(len(s)-k):
		if s[i:i+k] not in kDict:
			kDict[s[i:i+k]] = []
		kDict[s[i:i+k]] += [i]

	for i in range(len(t)-k):
		if t[i:i+k] in kDict:
			sIndices = kDict[t[i:i+k]]
			for seed in sIndices:
				seeds += [[seed,i,k]]

	diags = {}
	for seed in seeds:
		diag = seed[0]-seed[1]
		if diag not in diags:
			diags[diag] = []
		diags[diag] += [seed]

	for d in diags:
		diags[d].sort()
		i = 0
		changed = False
		while i < len(diags[d])-1:
			u, v = diags[d][i], diags[d][i+1]
			#print(u,v)
			#print(s[u[0]:u[0] + u[2]], t[u[1]:u[1]+u[2]])
			#print(s[v[0]:v[0] + v[2]], t[v[1]:v[1]+v[2]])
			length = v[0] - u[0] + k
			if compute(s[u[0]:u[0]+length],t[u[1]:u[1]+length],subDict) > compute(s[u[0]:u[0]+u[2]],t[u[1]:u[1]+u[2]],subDict):
				del(diags[d][i+1])
				diags[d][i][2] = length
			else:
				i += 1

	highScore = 0
	highScorers = []
	for d in diags:
		for seed in diags[d]:
			newScore = compute(s[seed[0]:seed[0]+seed[2]],t[seed[1]:seed[1]+seed[2]],subDict)
			if newScore > highScore:
				highScore = newScore
				highScorers = [seed]
			elif newScore == highScore:
				highScorers += [seed]

	if len(highScorers) > 0:
		bestAlignment = []
		for h in highScorers:
			diff = h[0] - h[1]
			out = bandedDP(s,t,diff,subDict,w)
			if len(bestAlignment) == 0 or out[0] > bestAlignment[0]:
				bestAlignment = out
	else:
		bestAlignment = bandedDP(s,t,0,subDict,w)
	if swapped:
		bestAlignment[1], bestAlignment[2] = bestAlignment[2], bestAlignment[1]
	return bestAlignment

def bandedDP(s,t,bandCenter,subDict,w):
	#region: diff finds center point at start
	s = "-" + s
	t = "-" + t
	score = np.zeros((len(s),len(t)))
	maxScore = 0
	position = [0,0]
	bandLeft = int(bandCenter - w/2)
	for j in range(len(t)):
		for i in range(w+1):
			if bandLeft + i == 0 or j == 0:
				score[bandLeft+i][j] = 0
			elif bandLeft + i > 0 and bandLeft + i < len(s):
				if i == 0:
					ms = max(score[bandLeft+i-1][j-1] + subDict[s[bandLeft+i]][t[j]], score[bandLeft+i][j-1] + subDict["-"][t[j]],0)
				elif i == w:
					ms = max(score[bandLeft+i-1][j-1] + subDict[s[bandLeft+i]][t[j]], score[bandLeft+i-1][j] + subDict[s[bandLeft+i]]["-"],0)
				else:
					ms = max(score[bandLeft+i-1][j-1] + subDict[s[bandLeft+i]][t[j]], score[bandLeft+i][j-1] + subDict["-"][t[j]], score[bandLeft+i-1][j] + subDict[s[bandLeft+i]]["-"],0)
				if ms > maxScore:
					maxScore = ms
					position = [bandLeft+i,j]
				score[bandLeft+i][j] = ms
		bandLeft += 1
	alignment = backtrack(subDict, score, position[0], position[1], s, t)
	score = float(compute(alignment[0],alignment[1],subDict))
	return [score, alignment[2], alignment[3]]

def heuralign(alphabet,scoring_matrix,sequence1,sequence2):
	subDict = matToDict(alphabet,scoring_matrix)
	k = 3		# key length
	w = 8		# band width
	alignment = FASTA(sequence1,sequence2,subDict,k,w)
	return alignment
