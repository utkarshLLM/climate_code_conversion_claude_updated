from config import client
from utils import completion_with_backoff
import re


def get_code_analysis(code):
    prompt = f"""
    Perform a static analysis on the following code. Identify independent units of code that function as a whole and are not dependent on other parts of the code. For each unit, provide:
    1. The line numbers it spans
    2. A brief description of its functionality
    3. An explanation of why it's considered independent

    Here's the code:

    {code}

    Human: Analyze the code above and provide the requested information.

    Assistant: Certainly! I'll analyze the code and identify independent units as requested.
    """

    response = completion_with_backoff(
        prompt=prompt,
        model="claude-2.1",
        max_tokens_to_sample=1000,
        stop_sequences=["\nHuman:"],
    )

    return response.completion

def parse_analysis(analysis):
    units = re.findall(r'Unit \d+:(.*?)(?=Unit \d+:|$)', analysis, re.DOTALL)
    parsed_units = []

    for unit in units:
        lines = re.search(r'Lines: (\d+)-(\d+)', unit)
        description = re.search(r'Description: (.*)', unit)
        independence = re.search(r'Independence: (.*)', unit)

        if lines and description and independence:
            parsed_units.append({
                'lines': (int(lines.group(1)), int(lines.group(2))),
                'description': description.group(1).strip(),
                'independence': independence.group(1).strip()
            })

    return parsed_units

# Example usage
code_to_analyze = """
def factorial(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n * factorial(n-1)

def fibonacci(n):
    if n <= 1:
        return n
    else:
        return fibonacci(n-1) + fibonacci(n-2)

def main():
    print(factorial(5))
    print(fibonacci(7))

if __name__ == "__main__":
    main()
"""

analysis = get_code_analysis(code_to_analyze)
parsed_units = parse_analysis(analysis)

for i, unit in enumerate(parsed_units, 1):
    print(f"Unit {i}:")
    print(f"Lines: {unit['lines'][0]}-{unit['lines'][1]}")
    print(f"Description: {unit['description']}")
    print(f"Independence: {unit['independence']}")
    print()