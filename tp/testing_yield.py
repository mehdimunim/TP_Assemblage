# generator to print even numbers
def print_even(test_list):
    for i in test_list:
        if i % 2 == 0:
            yield i


x = print_even([1, 2, 2, 2, 19, 38, 29])
for i in x:
    print(i)

