import gmpy2

def compute_sum(filename):
    total_sum = gmpy2.mpz(0)
    with open(filename, 'r') as file:
        for line in file:
            _, _, third_value = line.strip().split()
            total_sum += gmpy2.mpz(third_value)
    return total_sum

if __name__ == "__main__":
    filename = "PartialSums9.txt"
    result = compute_sum(filename)
    print("Total Sum:", result)
