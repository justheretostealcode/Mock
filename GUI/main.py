""" Main Modul zum GUI starten"""

def fib(start: int) -> float:
    """
    Fibonacci recursiv berechnen

    Args:
        start (int): der Wert der Fibonacci Folge

    Returns:
        float: Das Ergebnis der Fib reihe
    """
    if start == 0:
        return 0
    if start == 1:
        return 1
    return fib(start-1) + fib(start-2)

if __name__ == "__main__":
    print("Hello world")
    print(f"Fib von 10: {fib(10)}")
    input("..")
