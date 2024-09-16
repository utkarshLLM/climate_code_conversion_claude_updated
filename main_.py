from testing import check_if_docker_is_running
from utils import logger
from generate import generate_python_from_chunks
import openai


def translate(fortran_code: str):
    if not check_if_docker_is_running():
        logger.error("Docker is not running. Please start Docker and try again.")
        return ""

    python_code = generate_python_from_chunks(fortran_code)

    return python_code

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Translate Fortran code to Python code."
    )
    parser.add_argument(
        "--infile", type=str, default="source/input.f90", help="input Fortran file path"
    )
    parser.add_argument(
        "--outfile",
        type=str,
        default="source/output.py",
        help="output Python file path",
    )
    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile

    with open(infile, "r") as f:
        logger.info(f"Reading from {infile}")
        fortran_function = f.read()
        python_code = translate(fortran_function)
        
    with open(outfile, "w") as f:
        f.write(python_code)
        logger.info(f"Python translation written to {outfile}")
