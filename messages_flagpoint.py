from utils import logger


system_generate = """You're a programmer proficient in Python and unit testing.\ 
You can write and execute Python code by enclosing it in triple backticks, e.g. ```code goes here```"""
system_translate = """You're a programmer proficient in both Fortran and Python.\ 
More importantly you know how to connect these languages very well. Thus you are extremely good and thorough at exactly translating\ 
fortran code to python code."""

def generate_python_test_messages(python_function: str):
    prompt = f"""
Generate 3 unit tests for the following Python function using pytest. No need to import the module under test. ```python\n{python_function}\n```
    """

    logger.debug(f"PROMPT: {prompt}")

    return [
        {
            "role": "user",
            "content": prompt,
        },
    ]

example_unidentified_processes = """ In this fortran code, its calling a process called phys_getopts(). Something that is\
not present in this script and would be present in the larger atmospheric model instead. Even then, in the translation\
you should have a process by the same name with the same parameters; you don't need to define it separately, since it probably\
exists somewhere in the atmospheric model.
<fortran_code>
! Return error message
errstring = ' '

call phys_getopts(do_clubb_sgs_out = do_clubb_sgs)

! initialize  output fields for number conc qand ice nucleation
ncai(1:ncol,1:pver)=0._r8 
ncal(1:ncol,1:pver)=0._r8  
</fortran_code>
"""

def translate_to_python_messages(fortran_code: str, feedback: str = ""):
    prompt = f"""
    The following Fortran code is a part of the Earth Atmospheric model.\
    In this Fortran code, there is a subroutine by the name "micro_mg_tend". Your task is to translate subroutine "micro_mg_tend"\
    from Fortran code into its equivalent in the Python programming language.\
    This translation should be thorough and complete, ensuring that the resulting
    Python code maintains the original functionality and logic of the Fortran code.\n
    
    Here is the Fortran code you need to translate:
    <fortran>
    ```\n{fortran_code}```\n
    </fortran>
    \n\n

    To begin the translation process, carefully analyze the Fortran code and identify the following key
    elements:\n
    - Variables and their data types \n
    - Control structures (e.g., loops, conditionals) \n
    - Functions and subroutines \n
    - Input/output operations \n
    - Mathematical operations and expressions \n\n

    Once you have identified these elements, start translating each one into its Python equivalent. Keep
    in mind the following guidelines: \n
    - Fortran uses implicit variable typing, while Python uses explicit typing. Ensure that you declare\
    variables with the appropriate data types in Python.\n
    - Fortran uses different syntax for control structures compared to Python. For example, Fortran uses \
    "DO" for loops, while Python uses "for" and "while" loops. Make sure to use the correct Python \
    syntax for each control structure.\n
    - Fortran functions and subroutines should be translated into Python functions. Pay attention to the \
    arguments and return values of each function.\n
    - Fortran input/output operations, such as "READ" and "WRITE," should be replaced with their Python \
    equivalents, such as "input()" and "print()."\n
    - If there are functions or variables that you don't recognize, \
    include them in your translations anyway, take this for example: {example_unidentified_processes}.\n
    - Mathematical operations and expressions in Fortran should be translated to their Python \
    counterparts. Be aware of any differences in operator precedence and function names.\n\n
    
    Under no condtion, should you use placeholders like the following:\
    <placeholders>
    "# ... (initialize other constants and parameters)",  "# ... (perform calculations for each column)", \
    "# ... (perform calculations for each vertical level)", "# ... (more calculations)", "# ... (prepare output variables)", \
    "# ... (more post-processing)".
    </placeholders>
    
    Don't omit or simplify any parts of the code during translation, instead be thorough with your translations.\n\n

    As you translate each element, organize the Python code in a clear and readable format. Follow \
    Python best practices and conventions, such as using meaningful variable names, proper indentation,\
    and comments where necessary.\n\n

    When you are done, provide your Python code translation in <translated_code> tags.\n\n

    Whenever the translation stops, you have to describe the place in the fortran code where it stopped.\
    Enclose this information in <stopping_information> tags.\n\n

    Your only output should only be the Python code and this information about the stopping point. You should provide no\
    explanations or clarifications apart from that.\n\n

    Remember, the goal is to create a Python version of the Fortran code that is functionally equivalent \
    and follows Python's syntax and best practices. Take your time and carefully consider each element \
    of the Fortran code to ensure an accurate and totally complete translation.\n\n

    
    """
    logger.debug(f"PROMPT: {prompt}")

    return [
        {
            "role": "user",
            "content": prompt,
        },
    ]


## placeholder prompt:
# Instead of using placeholders described above, translate each and every initialization, calculation and post-processing thoroughly.\n\n