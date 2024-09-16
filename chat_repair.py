import os
import argparse
import time
import requests
from pathlib import Path, PurePath
import re
import openai

from raptor.utils import RichLog, QueueClient
from raptor.modules.prompt_composing import clean_up
from raptor.data_helpers import CodeNetVerifier

from transformers import AutoTokenizer
import torch

OPENAI_KEY          = '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'
lang_to_extension   = {"c": "c", "c++": "cpp", "rust": "rs", "java": "java", "python": "py"}

def dir_iter(dirname, model_str, overwrite, pattern=".*"):

    load_dir = Path(os.path.dirname(__file__)).joinpath(dirname)
    assert load_dir.exists()

    # technique = PurePath(load_dir).name.split('_')[1]

    for problem_dir in load_dir.iterdir():
        if not problem_dir.is_dir():
            continue

        for prompt_file in problem_dir.iterdir():

            if '_prompt' not in prompt_file.stem:
                continue
            
            if not re.match(pattern, prompt_file.stem):
                continue

            prompt      = prompt_file.read_text()
            sol_id      = prompt_file.stem.split('_')[0]
            technique   = prompt_file.stem.split('_')[2]
            source_lang = prompt_file.stem.split('_')[3]
            target_lang = prompt_file.stem.split('_')[4]
            problem_id  = PurePath(problem_dir).name
            target_ext  = lang_to_extension[target_lang.lower()]

            save_dir = load_dir
            save_dir.mkdir(exist_ok=True)
            save_dir.joinpath(problem_id).mkdir(exist_ok=True)
            save_path = save_dir.joinpath(f"{problem_id}/{sol_id}_{technique}_{source_lang}_{target_lang}_{model_str}-chatrepair.{target_ext}")
            baseline_path = save_dir.joinpath(f"{problem_id}/{sol_id}_{technique}_{source_lang}_{target_lang}_{model_str}.{target_ext}")

            if save_path.is_file() and not overwrite:
                continue

            yield {'prompt': prompt,
                   'save_path': save_path,
                   'baseline_path': baseline_path,
                   'target_language': target_lang,
                   'problem_id': problem_id,
                   'filename': prompt_file.name}

class OpenAIException(Exception):
    pass

def openai_gen(messages):
    count = 0
    while True:
        try:
            chat = openai.ChatCompletion.create(
                model="gpt-3.5-turbo-0125", messages=messages, temperature=0, seed=42
            )
            break
        except openai.error.InvalidRequestError as e:
            raise OpenAIException(f"Encountered an error with OpenAI API {e}")
        except openai.error.RateLimitError as e:
            count += 1
            if count >= 5:
                raise OpenAIException("Too many retries")
            print("OpenAI Rate Limit Error. Waiting 10 seconds and retrying")
            time.sleep(10)
        except:
            raise OpenAIException("Unknown error")
        
    return chat.choices[0].message.content

def wrap_format(code: str, lang: str):
    return f'```{lang.lower()}\n{code}\n```'

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Chat Repair')
    parser.add_argument('-overwrite',       dest="overwrite",       action="store_true",                help="Whether to re-translate if a translation alreay exists")
    parser.add_argument('-dir',             dest="dir",             default="",             type=str,   help="Path to directory containing one sub-folder per problem")
    parser.add_argument('-pattern',         dest="pattern",         default=".*",                        type=str,   help="Pattern to match the files")
    args = parser.parse_args()

    possible_languages = ["c", "rust", "java", "python", "go"]
    openai.api_key = OPENAI_KEY
    cnv = CodeNetVerifier()

    model_str = 'gpt'

    for data in dir_iter(args.dir, model_str, args.overwrite, args.pattern):

        messages = [{"role": "system", "content": "You are a intelligent code assistant. You generate only code, without any accompanying messages or text."},
                    {"role": "user",   "content": data['prompt']}]
        
        replies = []

        for attempt in range(1, 5+1):

            try:
                raw_reply   = openai_gen(messages)
            except OpenAIException as e:
                RichLog.error(f"Encountered an error with OpenAI API {e}")
                break

            reply       = raw_reply.strip()

            reply_head, reply_body, reply_tail = "", reply, ""

            # If there is a ```rust anywhere in the code, remove everything before it
            if "```rust" in reply_body:
                reply_head = reply_body[:reply_body.index("```rust") + len(f"```rust")]
                reply_body = reply_body[reply_body.index("```rust") + len(f"```rust"):]

            # If there is a ``` anywhere in the code, remove everything after it
            if "```" in reply_body:
                reply_tail  = reply_body[reply_body.index("```"):]
                reply_body  = reply_body[:reply_body.index("```")]
            
            if attempt == 1:
                reply_to_write = f"/*\n{reply_head}\n*/\n{reply_body}\n/*\n{reply_tail}\n*/"

                if args.overwrite or not data['baseline_path'].is_file():
                    data['baseline_path'].write_text(reply_to_write)
                    RichLog.info(f"Wrote translated file to {data['baseline_path'].absolute()}")

            result = cnv.verify(reply_body, data['target_language'], data['problem_id'])

            if (result['status'] == "Success") or (result['message'] == 'Timeout') or (result['status'] == "No tests"):
                break
            else:
                if attempt == 5:
                    break

                RichLog.info(f"Translating this file didn't work due to {result['status']}. Giving feedback and trying again. Attempt {attempt+1}")
                repair_message = f"The above code produced the error:\n\n{result['message']}\n\n"
                repair_message += "Can you fix that and re-generate the code? Please re-generate the *entire* program, not just the part that caused the error."

                messages = [{"role": "system", "content": "You are a intelligent code assistant. You generate only code, without any accompanying messages or text."},
                            {"role": "user",  "content": wrap_format(reply_body, "rust") + "\n\n" + repair_message}]
            
                replies += [{"role": "Assistant", "content": raw_reply}, {"role": "User", "content": repair_message}]
            
        try:
            intermediate_messages = ["/*{}:\n\n{}\n*/".format(message['role'], message['content']) for message in replies]
            intermediate_messages = "\n".join(intermediate_messages)
            reply_to_write = f"{intermediate_messages}/*\n{reply_head}\n*/\n{reply_body}\n/*\n{reply_tail}\n*/"

            data['save_path'].write_text(reply_to_write)
            RichLog.info(f"Wrote translated file to {data['save_path'].absolute()}")
        except NameError:
            RichLog.info(f"Failed to translate {data['filename']}")

        try:
            del reply_body
        except NameError:
            pass