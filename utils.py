from typing import List, Dict
import csv
import difflib
import re
from loguru import logger
import requests
from tenacity import (
    retry,
    stop_after_attempt,
    wait_random_exponential,
    wait_exponential, 
    wait_fixed, 
    before_sleep_log,
    retry_if_exception_type, 
)
from config import client
from anthropic import RateLimitError, APIError, APITimeoutError
import logging
import sys

logger.remove()
logger.add(sys.stderr, level="TRACE")


def save_to_csv(dict: List[Dict], outfile: str):
    with open(outfile, "w") as f:
        w = csv.DictWriter(f, fieldnames=dict[0].keys())
        w.writeheader()
        w.writerows(dict)


def find_diffs(source_code: str, previous_source_code: str) -> str:
    """
    Finds the diffs between the newly generated source_code and unit_tests and the corresponding python_function and python_unit_tests.
    """

    # Find the diff between the unit tests
    diff = difflib.unified_diff(
        source_code.splitlines(), previous_source_code.splitlines()
    )

    # Convert the diff to a string
    diff_str = "\n".join(diff)

    return diff_str


def remove_ansi_escape_codes(s):
    ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
    return ansi_escape.sub("", s)


def extract_code_block(completion):
    # Extract the code block from the completion
    code = completion.content[0].text.split("```")[1]
    code = code.replace("python\n", "")

    return code

def extract_code_block_new(full_reponse):
    # Extract the code block from the completion
    code = full_reponse.content[0].text.split("<translated_code>")[1]
    code = code.replace("</translated_code>", "")
    
    return code


logging_logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def get_last_fortran_point_info(completion):
    text = completion.content[0].text.strip()
    match = re.search(r'<stopping_information>(.*?)</stopping_information>', text, re.DOTALL)
    return match.group(1) if match else ''

# @retry(
    # wait=wait_fixed(60), #Since the rate limit is per minute
    # stop=stop_after_attempt(20),
    # before_sleep=before_sleep_log(logging_logger, log_level=logging.INFO),)

@retry(
    retry=retry_if_exception_type((RateLimitError, APIError, APITimeoutError)),
    stop=stop_after_attempt(20),
    wait=wait_fixed(60),
    before_sleep=before_sleep_log(logging_logger, log_level=logging.INFO),
    reraise=True, )
def completion_with_backoff(**kwargs):
    model_ = kwargs.get('model')
    temperature_ = kwargs.get('temperature')
    messages_ = kwargs.get('messages')
    max_tokens = kwargs.get('max_tokens')
    system_ = kwargs.get('system')

    try:
        return client.messages.create(
            max_tokens=max_tokens, 
            system = system_, 
            model=model_, 
            messages=messages_,
            temperature=temperature_, 
            extra_headers={
            "anthropic-beta": "max-tokens-3-5-sonnet-2024-07-15"
            }, 
        )
    except (RateLimitError, APIError, APITimeoutError) as e:
        print(f"Rate limit exceeded. This error was raised by 'completion_with_backoff' Retrying... ({str(e)})")
        raise

def find_nth(haystack, needle, n):
    """
    Finds index of the *n*'th occurrence of *needle* within *haystack*.
    Returns -1 when the *n*'th occurrence is not found.
    """
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start + len(needle))
        n -= 1
    return start


def extract_unit_test_code(message):
    start_marker = "UNIT TESTS:"
    end_marker = "```"

    if message.find(start_marker) < 0:
        return None

    # find the position of start and end markers
    start = message.find(start_marker) + len(start_marker)
    end = message.rfind(end_marker)

    # extract the unit test code
    unit_test_code = message[start:end].strip()

    # Remove ``` and python
    unit_test_code = unit_test_code.replace("```", "")
    unit_test_code = unit_test_code.replace("python\n", "")

    return unit_test_code


def extract_source_code(message: str):
    start_marker = "SOURCE CODE:"
    end_marker = "```"

    if message.find(start_marker) < 0 or find_nth(message, end_marker, 2) < 0:
        return None

    # find the position of start and end markers
    start = message.find(start_marker) + len(start_marker)
    end = find_nth(message, end_marker, 2)

    # extract the source code
    source_code = message[start:end].strip()

    # Remove ``` and python
    source_code = source_code.replace("```", "")
    source_code = source_code.replace("python\n", "")

    return source_code

def extract_source_code_new(message: str):
    match = re.search(r'<translated_code>(.*?)</translated_code>', message, re.DOTALL)
    return match.group(1) if match else ''


    return source_code

def options_menu(options: list[str]):
    response = ""

    print("Please select one of the following options:")
    for i, option in enumerate(options):
        print(f"{i+1}. {option}")

    while True:
        response = input(f"Select an option [1-{len(options)}]: ")
        try:
            choice = int(response)
            if choice >= 1 and choice <= len(options):
                return options[choice - 1]
        except:
            pass


def write_to_file(source_code: str, outfile: str):
    with open(outfile, "a") as f:
        f.write("\n")
        f.write(source_code)
        f.write("\n")


