import os
import dotenv
import promptlayer

dotenv.load_dotenv()

# Only use promptlayer if the API key is provided
if os.environ.get("PROMPTLAYER_API_KEY"):
    promptlayer.api_key = os.environ.get("PROMPTLAYER_API_KEY")
    anthropic = promptlayer.anthropic
    claude_api_key = os.environ.get("CLAUDE_API_KEY")
    client = anthropic.Anthropic(api_key=claude_api_key)
    
else:
    import anthropic

    claude_api_key = os.environ.get("CLAUDE_API_KEY")

    client = anthropic.Anthropic(
    # defaults to os.environ.get("ANTHROPIC_API_KEY")
    api_key=claude_api_key,)


model_name = "claude-3-5-sonnet-20240620"