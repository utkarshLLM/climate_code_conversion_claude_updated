from config import client
from utils import completion_with_backoff
from loguru import logger

def get_continuation_prompt(original_prompt, full_response):
    
    continuation_prompt = f"""
    You are an AI assistant tasked with continuing the translation of Fortran code to Python code. The
    original translation was incomplete, and your job is to pick up where it left off. Here's what you
    need to know:

    1. Original Translation Prompt:
    <original_prompt>
    {original_prompt}
    </original_prompt>

    2. Last Generated Python Code:
    <last_python_code>
    {full_response}
    </last_python_code>

    3. Instructions for Continuing the Translation:
    - Carefully review the last generated Python code to understand where the translation stopped.
    - Continue the translation from the exact point where it left off.
    - Maintain the same translation style and conventions used in the previous part of the code.

    4. Guidelines for Maintaining Consistency:
    - Use the same indentation style as in the previously translated code.
    - Keep variable and function names consistent with those already used.
    - If specific Python libraries or modules were imported in the earlier part, continue to use them as
    appropriate.

    5. Output Format:
    Provide your continued translation in the following format:
    <continued_translation>
    [Your continued Python code here]
    </continued_translation>

    Remember, your task is to seamlessly continue the translation as if it had never been interrupted.
    Maintain consistency with the original translation prompt and the previously generated Python code.
    """

    return continuation_prompt

def get_continuation(prompt_json, completion, **kwargs):
    max_attempts = 3
    full_response = ""
    attempts = 0

    original_prompt = prompt_json
    
    model_ = kwargs.get('model')
    temperature_ = kwargs.get('temperature')
    max_tokens = kwargs.get('max_tokens')
    system_ = kwargs.get('system')

    chunk = completion.content[0].text.strip()
    full_response += chunk

    while attempts < max_attempts:
        try:
            
            ###### Make a breaking condition
            # if completion.stop_reason == "stop":
            #     # Task completed
            #     return full_response
            
            
            prompt = get_continuation_prompt(original_prompt, full_response)

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

            logger.trace(f'COMPLETION: {new_completion.content[0].text}')

            chunk = new_completion.content[0].text.strip() #current response
            ## add something that removes the tags or the LLM response text.
            full_response += chunk # current response added to the 'total' response
            
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
            
            new_completion = completion_with_backoff(max_tokens=max_tokens, 
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
        

            chunk = new_completion.content[0].text.strip() #current response
            ## add something that removes the tags or the LLM response text.
            full_response += chunk # current response added to the 'total' response
            
            attempts += 1
        
        except Exception as e:
            print(f"An error occurred, Continuation didn't work: {e}")
            attempts += 1
    
    return full_response







##### The prompt used before:
# prompt = f"""
#                 The task given to you previously was this: 
#                 <Instructions> 
#                 {original_prompt} 
#                 </Instructions>
                
#                 Here's what has been generated so far from the instructions provided above:
#                 <response_so_far>
#                 {full_response}
#                 </response_so_far>
                
#                 Please continue from where you left off, keeping in mind the previous instructions \
#                 and while also using this as feedback from what has been generated so far: 
#                 <additional_instructions>
#                 {add_inst}
#                 </additional_instructions>
#                 """


####### Additional instructions try and except block
# try:
#     text_from_response = response.content[0].text.strip()
#     add_inst = text_from_response.split('</translated_code>')[1]

#     prompt = get_continuation_prompt(original_prompt, full_response, add_inst)
    
# except:
#    prompt = get_continuation_prompt(original_prompt, full_response)