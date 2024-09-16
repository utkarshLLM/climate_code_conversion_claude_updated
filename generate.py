import testing as testing
import utils as utils
from utils import logger, completion_with_backoff
import messages
from messages import system_generate, system_translate
from config import model_name
from continuation import get_continuation
from tree_sitter_script import chunk_code


###########################

def generate_unit_tests(python_function: str):
    logger.debug("Generating unit tests based on python code...")

    completion = completion_with_backoff(
        model=model_name,
        system = system_generate,
        max_tokens=8192,  
        temperature=0.0,
        messages=messages.generate_python_test_messages(python_function)
    )

    logger.trace(f'COMPLETION: {completion.content[0].text}')  
    # type: ignore

    unit_tests = utils.extract_code_block(completion)

    return unit_tests

###########################

def generate_python(fortran_code: str):
    logger.debug("Translating function to Python...")

    # The original prompt
    prompt_json = messages.translate_to_python_messages(fortran_code)

    completion = completion_with_backoff(
        max_tokens = 8192,
        system = system_translate,
        model=model_name,
        messages=prompt_json,
        temperature=0,
    )


    # Extract the code block from the completion
    python_function = utils.extract_code_block(completion)

    # trace
    logger.trace(f'PYTHON_FUNCTION_EXTRACTED: {python_function}')

    return python_function

###########################

def generate_python_with_continuation(fortran_code: str):
    
    logger.debug("Translating function to Python...")

    # The original prompt
    prompt_json = messages.translate_to_python_messages(fortran_code)

    completion = completion_with_backoff(
        max_tokens = 8192,
        system = system_translate,
        model=model_name,
        messages=prompt_json,
        temperature=0,
    )
    
    continued_completion = get_continuation(prompt_json, 
                                           completion, 
                                           max_tokens = 8192,
                                           system = system_translate,
                                           model=model_name, 
                                           temperature=0,)
    
    
    # Extract the code block from the completion
    python_function = utils.extract_code_block_new(continued_completion) ## WHENEVER ITS SWITCHED ON, THE OUTPUT CODE BEHAVES BETTER

    # trace
    logger.trace(f'PYTHON_FUNCTION_EXTRACTED: {python_function}')

    return python_function

###########################
   
## Chunk the code and then translate it.
def generate_python_from_chunks(fortran_code: str, algorithm: str = "greedy"):
    logger.debug("Translating function to Python...")

    # chunked code
    code_chunks = chunk_code(fortran_code)
      
    if algorithm=="greedy":
        full_response = ""
        for code_chunk in code_chunks:

            prompt_chunk_json = messages.translate_to_python_messages(code_chunk)

            completion = completion_with_backoff(
                max_tokens = 8192,
                system = system_translate,  ## import this from messages.
                model=model_name,
                messages=prompt_chunk_json, ## this needs to be changed
                temperature=0, ## can play around with this for fine-tuning     
            )

            logger.trace(f'COMPLETION: {completion.content[0].text}')

            translated_chunk = completion.content[0].text.strip()
            full_response = full_response + "\n" + translated_chunk  # Starts in a new line

    else:
        for code_chunk in code_chunks:
                
                print("")
            
                # prompt_chunk_json = messages.translate_to_python_from_chunk_messages(code_chunk)
            
                # completion = completion_with_backoff(
                #     max_tokens = 8192,
                #     system = system_translate_chunks,  ## import this from messages.
                #     model=model_name,
                #     messages=prompt_chunk_json,
                #     temperature=0, ## can play around with this for fine-tuning
                # )
  
    
    # Extract the code block from the completion
    python_function = full_response

    # trace
    logger.trace(f'PYTHON_FUNCTION_EXTRACTED: {python_function}')

    return python_function


############################################################################################################



