def test_sum_two_positive_numbers():
    assert sum_two_numbers(5, 3) == 8

def test_sum_positive_and_negative_numbers():
    assert sum_two_numbers(10, -7) == 3

def test_sum_two_negative_numbers():
    assert sum_two_numbers(-4, -6) == -10
