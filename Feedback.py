# updated to have a feedback loop, which works for 3 iterations.
# def generate_python(fortran_code: str):
#     logger.debug("Translating function to Python...")

#     python_function = ""
#     feedback = ""

#     for iteration in range(3):
#         logger.debug(f"Starting iteration {iteration + 1}")

        
#         # Include feedback in the messages for subsequent iterations
#         if iteration == 0:
#             current_messages = messages.translate_to_python_messages(fortran_code)
#             # Splitting the string by the closing ''' delimiter
#             parts = current_messages.split("'''")
#             ### Add something for if the 'parts' if empty
#             # The later characters will be in the second part after the split
#             current_feedback = parts[-1].strip()
#         else:
#             current_messages = messages.translate_to_python_messages(fortran_code)
#             # Splitting the string by the closing ''' delimiter
#             parts = current_messages.split("'''")
#             # The later characters will be in the second part after the split
#             current_feedback = parts[-1].strip()

#         completion = completion_with_backoff(
#             max_tokens = 4096, 
#             system = system_translate, 
#             model=model_name,
#             messages=current_messages,
#             temperature=0,
#         )

#         logger.trace(f'COMPLETION: {completion.content[0].text}')

#         # Extract the code block from the completion
#         python_function = utils.extract_code_block(completion)

#         # Generate feedback for the next iteration
#         if iteration < 2:  # No need for feedback on the last iteration
#             feedback = generate_feedback(fortran_code, python_function)

#         logger.debug(f"Completed iteration {iteration + 1}")

#     return python_function

# def generate_feedback(fortran_code: str, python_code: str):
#     logger.debug("Generating feedback for the translation...")
#     
#     message = """this is a message from the LLM, extract anything that sounds like feedback, or sounds like \
#     things that need to be done.\
#     """
#     completion = completion_with_backoff(
#         max_tokens = 1024, 
#         system = system_feedback,
#         model=model_name,
#         messages=messages.generate_feedback_messages(fortran_code, python_code),
#         temperature=0,
#     )

#     logger.trace(f'FEEDBACK: {completion.content[0].text}')

#     return completion.content[0].text