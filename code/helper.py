'''Contains functions that can be useful in many situations and fit in no specific file'''



def splitInterval (a: float, b: float, n: int) -> list[float] :
	'''Returns a list of n ordered equally-spaced real numbers in [a, b]'''
	delta = (b - a) / n
	return  [a + i * delta for i in range(n)]