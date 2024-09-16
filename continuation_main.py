from testing import run_tests, TestResult, check_if_docker_is_running
from utils import logger
from debug import fix_problem
from generate import generate_unit_tests, generate_python
from continuation import get_continuation, get_continuation_from_file
from config import model_name
import messages


def continue_(incomplete_python_function, 
              fortran_code: str):
    if not check_if_docker_is_running():
        logger.error("Docker is not running. Please start Docker and try again.")
        return ""

    original_prompt = messages.translate_to_python_messages(fortran_code, )
    
    # instead of the incomplete python code, it just starts where it left off and
    # then we can manually assemble it
    # Write a code to assmeble it rather
    python_code = get_continuation_from_file(original_prompt, #in its default form 
                                   incomplete_python_function, #completion will be replaced by incomplete_python_function 
                                   max_tokens = 8192,
                                   system = messages.system_translate,
                                   model=model_name,
                                   temperature=0,)
    
    return python_code

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Continue Incomplete Fortran to Python translation code."
    )
    parser.add_argument(
        "--fortranfile", 
        type=str, 
        help="the original fortran code", #no default
    )
    parser.add_argument(
        "--infile", type=str, help="input incomplete Python file path" #no default
    )
    parser.add_argument(
        "--outfile",
        type=str,
        help="output Python file path", #no default
    )

    args = parser.parse_args()

    ## the command for this would be such that infile and outfile are the same
    ## In other words its gonna overwrite the inputfile.
    fortranfile = args.fortranfile
    infile = args.infile
    outfile = args.outfile

    with open(fortranfile, "r") as f:
        logger.info(f"Reading from {fortranfile}")
        fortran_code = f.read()

    with open(infile, "r") as f:
        logger.info(f"Reading from {infile}")
        incomplete_python_function = f.read()
        
        python_code = continue_(incomplete_python_function, fortran_code) # won't really be the correct python code, 
                                                                          # would have LLM replies and some XML tags etc.

    with open(outfile, "w") as f:
        f.write(python_code)
        logger.info(f"Continued Python translation written to {outfile}")

