from config import client
from utils import completion_with_backoff, get_last_fortran_point_info
from loguru import logger


def get_continuation_prompt(original_prompt, last_fortran_point):
    
    continuation_prompt = f"""
    You are tasked with continuing the translation of Fortran code to Python code. The fortran code provided to you\
    was in process of translation to python before it stopped. You are provided with details about at what exact point\
    in the Fortran code did the translation stop, and your job is to use this information and continue translation\
    from Fortan to Python from this point on.\
    Here's what you need to know:

    1. Original Translation Prompt:
    <original_prompt>
    {original_prompt}
    </original_prompt>

    2. Information about the last point where the python code stopped:
    <info_last_fortran_code_point>
    {last_fortran_point}
    </info_last_fortran_code_point>

    3. Output Format:
    Provide your continued translation in the following format:
    <continued_translation>
    [Your continued Python code here]
    </continued_translation>

    Remember, your task is to seamlessly continue the translation as if it had never been interrupted.
    Maintain consistency with the original translation prompt and the previously generated Python code.
    """

    return continuation_prompt

def get_continuation(prompt_json, stopping_point, **kwargs):
    max_attempts = 3
    full_response = ""
    attempts = 0

    original_prompt = prompt_json
    
    model_ = kwargs.get('model')
    temperature_ = kwargs.get('temperature')
    max_tokens = kwargs.get('max_tokens')
    system_ = kwargs.get('system')

    while attempts < max_attempts:
        try:
            
            ###### Make a breaking condition
            # if completion.stop_reason == "stop":
            #     # Task completed
            #     return full_response
            
            
            prompt = get_continuation_prompt(original_prompt, stopping_point)

            new_completion = client.messages.create(max_tokens=max_tokens, 
                                                     system = system_, 
                                                     model=model_, 
                                                     messages=[
                                                         {
                                                             "role":"user", 
                                                             "content":prompt
                                                             }
                                                             ],
                                                    temperature=temperature_, 
                                                    extra_headers={
                                                        "anthropic-beta": "max-tokens-3-5-sonnet-2024-07-15"
                                                        },)

            logger.trace(f'COMPLETION: {new_completion.content[0].text}') # I'll get the trace anyway

            stopping_point = get_last_fortran_point_info(new_completion)
            
            attempts += 1
        
        except Exception as e:
            print(f"An error occurred, Continuation didn't work: {e}")
            attempts += 1
    
    return full_response

def get_continuation_from_file(prompt_json, completion, **kwargs):
    max_attempts = 1
    full_response = ""
    attempts = 0

    original_prompt = prompt_json
    
    model_ = kwargs.get('model')
    temperature_ = kwargs.get('temperature')
    max_tokens = kwargs.get('max_tokens')
    system_ = kwargs.get('system')

    full_response += completion # completion is already a string

    while attempts < max_attempts:
        try:
            
            ###### Make a breaking condition
            # if completion.stop_reason == "stop":
            #     # Task completed
            #     return full_response
            
            
            prompt = get_continuation_prompt(original_prompt, full_response)

            new_completion = client.messages.create(
            max_tokens=max_tokens, 
            system = system_, 
            model=model_, 
            messages=[
                {
                    "role":"user", 
                    "content":prompt
                 }
                 ],
            temperature=temperature_, 
            extra_headers={
            "anthropic-beta": "max-tokens-3-5-sonnet-2024-07-15"
            }, )

            chunk = new_completion.content[0].text.strip() #current response
            ## add something that removes the tags or the LLM response text.
            full_response += chunk # current response added to the 'total' response
            
            attempts += 1
        
        except Exception as e:
            print(f"An error occurred, Continuation didn't work: {e}")
            attempts += 1
    
    return full_response



