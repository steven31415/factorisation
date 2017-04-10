import sys
import itertools
import pysais
import numpy as np

# To install PySAIS module needed for Suffix Arrays and LCP
# follow README.md instructions on https://github.com/AlexeyG/PySAIS

# class RMQ as taken from:
# https://gist.github.com/m00nlight/1f226777a49cfc40ed8f
#
# This is most likely not an efficient implementation of RMQ
# Use as a placeholder
#
# When performing queries on a 0-indexed string
# the range considered is [a, b)
#
# class Vertex, shortestPathBFS, traverseShortestPath taken from:
# http://pythonfiddle.com/shortest-path-bfs/



class Vertex:
	"""
	Adjacency List implementation of a graph vertex
	We create a very simple class to represent or Graph nodes so we can use it in our Graph Traversal Algorithms
	Just the bare essentials were included here
	"""
	def __init__(self, vert_id):
		"""
		Constructor
		:param vert_id: The id that uniquely identifies the vertex.
		"""
		self.vert_id = vert_id          # simple type
		self.neighbors = []             # type List[Vertex]
		self.status = 'undiscovered'    # undiscovered | discovered | explored

		self.distance = -1              # shortest distance from source node in shortest path search
		self.previous = None            # previous vertex in shortest path search

	def addVertex(self, vertex):
		"""
		Adds a new vertex as an adjacent neighbor of this vertex
		:param vertex: new Vertex() to add to self.neighbors
		"""
		self.neighbors.append(vertex)

	def getNeighbors(self):
		"""
		Returns a list of all neighboring vertices
		:return: list of vertexes
		"""
		return self.neighbors



def shortestPathBFS(vertex):
	"""
	Shortest Path - Breadth First Search
	:param vertex: the starting graph node
	:return: does not return, changes in place
	"""
	if vertex is None:
		return

	queue = []                                  # our queue is a list with insert(0) as enqueue() and pop() as dequeue()
	queue.insert(0, vertex)

	while len(queue) > 0:
		current_vertex = queue.pop()                    # remove the next node in the queue
		next_distance = current_vertex.distance + 1     # the hypothetical distance of the neighboring node

		for neighbor in current_vertex.getNeighbors():

			if neighbor.distance == -1 or neighbor.distance > next_distance:    # try to minimize node distance
				neighbor.distance = next_distance       # distance is changed only if its shorter than the current
				neighbor.previous = current_vertex      # keep a record of previous vertexes so we can traverse our path
				queue.insert(0, neighbor)



def traverseShortestPath(target):
	"""
	Traverses backward from target vertex to source vertex, storing all encountered vertex id's
	:param target: Vertex() Our target node
	:return: A list of all vertexes in the shortest path
	"""
	vertexes_in_path = []

	while target.previous:
		vertexes_in_path.append(target.vert_id)
		target = target.previous

	return vertexes_in_path



class RMQ:
	def __init__(self, n):
		self.sz = 1
		self.inf = (1 << 31) - 1
		while self.sz <= n: self.sz = self.sz << 1
		self.dat = [self.inf] * (2 * self.sz - 1)

	def update(self, idx, x):
		idx += self.sz - 1
		self.dat[idx] = x
		while idx > 0:
			idx = (idx - 1) >> 1
			self.dat[idx] = min(self.dat[idx * 2 + 1], self.dat[idx * 2 + 2])

	def query(self, a, b):
		return self.query_help(a, b, 0, 0, self.sz)

	def query_help(self, a, b, k, l, r):
		if r <= a or b <= l:
			return sys.maxsize
		elif a <= l and r <= b:
			return self.dat[k]
		else:
			return min(self.query_help(a, b, 2 * k + 1, l, (l + r)>>1),
					self.query_help(a, b, 2 * k + 2, (l + r) >> 1, r))



# converts string into degenerate sequence
def convert_to_degenerate(sequence):

	degenerate_sequence = []
	i = 0

	while i < len(sequence):

		if sequence[i] != "[":
			degenerate_sequence.append([sequence[i]])
		else:
			non_solid = []

			i += 1
			while sequence[i] != "]":
				non_solid.append(sequence[i])
				i += 1

			degenerate_sequence.append(non_solid)

		i += 1

	return degenerate_sequence



################################################################################

# factorising function
def get_maximal_palindromic_factorisation(degenerate_sequence, k, alphabet):

	lambda_sequence = ""
	lambda_locations = []
	alphabet_value = {}

	# assign numerical values to alphabet characters
	for i in range(len(alphabet)):
		alphabet_value[alphabet[i]] = i

	# replace non solid symbols with lambda characters
	ascii_offset = 128
	lambda_value = ascii_offset

	for i in range(len(degenerate_sequence)):
		char = degenerate_sequence[i]

		if len(char) == 1:
			lambda_sequence = lambda_sequence + char[0]
		else:
			lambda_locations.append(i)
			lambda_sequence = lambda_sequence + chr(lambda_value)
			lambda_value += 1


	# MATCH TABLE FORMAT:
	#
	# Usage: match_table[i][j]
	# i = l1...lk + s1...  (lambdas + alphabet)
	# j = l1...lk          (lambdas)

	###########################
	#   l1 l2 ... lk s1 ...   #
	# l1                      #
	# l2                      #
	# ...                     #
	# lk                      #
	# s1                      #
	# ...                     #
	###########################

	# store match table
	match_table = []

	for i in range(len(lambda_locations) + len(alphabet)):
		match_table.append([False] * len(lambda_locations))

	for i in range(len(lambda_locations)):
		for j in range(len(lambda_locations)):
			p1 = 0
			p2 = 0

			char1 = degenerate_sequence[lambda_locations[i]]
			char2 = degenerate_sequence[lambda_locations[j]]

			while p1 < len(char1) and p2 < len(char2):
				if char1[p1] == char2[p2]:
					match_table[i][j] = True
					break
				elif char1[p1] < char2[p2]:
					p1 += 1
				else:
					p2 += 1

	for i in range(len(alphabet)):
		for j in range(len(lambda_locations)):
			if alphabet[i] in degenerate_sequence[lambda_locations[j]]:
				match_table[len(lambda_locations) + i][j] = True


	# create match checking funtion
	def real_match(char1, char2):           #CHANGE THIS SO THAT LAMBDAS ARE CONSTANT TIME IDENTIFIED
		if char1 == '#' or char2 == '#' or char1 == '$' or char2 == '$':
			return char1 == char2

		if char1 in alphabet:
			if char2 in alphabet:
				return char1 == char2
			else:
				return match_table[len(lambda_locations) + alphabet_value[char1]][ord(char2) - ascii_offset]
		else:
			if char2 in alphabet:
				return match_table[len(lambda_locations) + alphabet_value[char2]][ord(char1) - ascii_offset]
			else:
				return match_table[ord(char1) - ascii_offset][ord(char2) - ascii_offset]


	# create string for suffix tree
	suffix_tree_string = lambda_sequence + '#' + lambda_sequence[::-1] + '$'

	# store length of suffix tree string
	n_st = len(suffix_tree_string)

	# store length of sequence
	n = len(degenerate_sequence)

	# create suffix array and LCP array
	sa = pysais.sais(suffix_tree_string)
	lcp, lcp_lm, lcp_mr = pysais.lcp(suffix_tree_string, sa)

	# create inverse suffix array
	isa = [0] * n_st
	for i in range(0, n_st):
		isa[sa[i]] = i

	# preprocess LCP array for RMQ
	rmq = RMQ(n_st)
	for i in range(n_st):
		rmq.update(i, lcp[i])


	# define generic LCP function
	def LCP(i, j):
		if i == j:
			return n_st - i

		a = isa[i]
		b = isa[j]

		if a > b:
			a, b = b, a

		return rmq.query(a, b)


	# define real LCP function taking into account degenerate positions
	def real_LCP(i, j, k):
		real_lcp = 0
		mismatch_count = 0

		# perform up to k + 1 LCP queries for each suffix
		while mismatch_count < k + 1:

			# perform an LCP and update next location to check
			real_lcp = real_lcp + LCP(i + real_lcp, j + real_lcp)

			# stop if end of string is reached
			if i + real_lcp >= n_st or j + real_lcp >= n_st:
				break

			# stop if real mismatch encountered
			char1 = suffix_tree_string[i + real_lcp]
			char2 = suffix_tree_string[j + real_lcp]

			if not real_match(char1, char2):
				break
			else:
				real_lcp += 1

		# increment mismatch count
		mismatch_count += 1

		# remove 1 to ensure final mismatch is not included
		return real_lcp


	# determine locations of even palindromes
	even_palindromes = [(0, 0)] * (n - 1)
	
	for i in range(1, n):
		j = 2 * n - i + 1

		e = real_LCP(i, j, k)

		left = i - e
		right = i + e - 1

		if left <= right:
			even_palindromes[i - 1] = (left, right)
		else:
			even_palindromes[i - 1] = None


	# determine locations of odd palindromes
	odd_palindromes = [(0, 0)] * n

	for i in range(1, n + 1):
		j = 2 * n - i + 2

		e = real_LCP(i, j, k)

		left = i - e - 1
		right = i + e - 1

		if left <= right:
			odd_palindromes[i - 1] = (left, right)
		else:
			odd_palindromes[i - 1] = None


	# build graph of maximal palindromes
	vertices = [0] * (n + 1)
	for i in range(n + 1):
		vertices[i] = Vertex(i)

	for palindrome in even_palindromes:
		if palindrome != None:
			start = palindrome[0]
			end = palindrome[1]
			vertices[start].addVertex(vertices[end + 1])

	for palindrome in odd_palindromes:
		if palindrome != None:
			start = palindrome[0]
			end = palindrome[1]
			vertices[start].addVertex(vertices[end + 1])

	source = vertices[0]
	target = vertices[n]


	# determine shortest path corresponding to factorisation
	shortestPathBFS(source)
	vertices_in_path = traverseShortestPath(target)
	if len(vertices_in_path) > 0:
		return [0] + vertices_in_path[::-1][:-1]
	else:
		return None

################################################################################



# perform tests
line_count = 0
correct_result_count = 0

with open('tests.txt') as f:
	next(f)

	for line in f:

		values = str.split(line)

		k = int(values[0])

		if k > 128:
			print("This basic implementation does not accept k greater than 128")
			break

		alphabet = values[1]
		sequence = values[2]
		correct_result = eval(''.join(values[3:]))

		degenerate_sequence = convert_to_degenerate(sequence)
		factorisation = get_maximal_palindromic_factorisation(degenerate_sequence, k, alphabet)

		# print result of each test
		print("TEST: " + str(degenerate_sequence))
		print("FACTORISATION: " + str(factorisation))

		# increment test result counts
		if factorisation == correct_result:
			correct_result_count += 1
			print("Correct")
		else:
			print("Incorrect")
		print("\n")

		line_count += 1

# print summary of test results
print(str(correct_result_count) + " out of " + str(line_count) + " tests passed")
	
