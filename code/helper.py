

def splitInterval (a: float, b: float, n: int) -> list[float] :
	'''Returns a list of n ordered equally-spaced real numbers in [a, b]'''
	delta = (b - a) / n
	return  [a + i * delta for i in range(n)]