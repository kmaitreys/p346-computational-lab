#Summation of first N natural numbers

def sumN(n):
	N = n
	s_n = 0
	while n > 0:
		s_n = s_n + n
		n = n - 1;
	print(str(s_n) + " is the sum of first " + str(N) + " natural numbers.")


#Summation of first N odd numbers

def sumOdd(n):
	N = n
	s_n = 0
	while n > 0:
		s_n = s_n + (2*n - 1)
		n = n - 1;
	print(str(s_n) + " is the sum of first " + str(N) + " odd numbers.")



#Arithmetic Progression

def sumAP(a, d, n):
	s_ap = 0
	A = a
	i = 0
	while i < n:
		s_ap = s_ap + a
		a = a + d
		i = i + 1;
	print(str(float(s_ap)) + " is the sum of first " + str(n) + " terms of the arithmetic progression with first term " + str(A) + " and common difference " + str(d) + ".")



#Geometric Progression

def sumGP(a, r, n):
	s_gp = 0
	A = a
	i = 0
	while i < n:
		s_gp = s_gp + a
		a = a*r
		i = i + 1;
	print(str(float(s_gp)) + " is the sum of first " + str(n) + " terms of the geometric progression with first term " + str(A) + " and common ratio " + str(r) + ".")



#Harmonic Progression

def sumHP(a, d, n):
	if a == 0:
		print("A harmonic progression cannot have zero as one of its term!") #taking care of exceptional case in HP
	else:
		s_hp = 0
		A = a
		i = 0
		while i < n:
			s_hp = s_hp + 1/a
			a = a + d
			i = i + 1;
		print(str(float(s_hp)) + " is the sum of first " + str(n) + " terms of the harmonic progression with first term " + str(A) + " and common difference " + str(d) + ".")

#Factorial

def fact(n):
	f = 1
	if n == 0:
		return 1 #taking care of factorial of zero
	else:
		for i in range(1, n + 1):
			f = f*i
	return f
